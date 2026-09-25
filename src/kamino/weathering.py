import numpy as np
import numpy.typing as npt

from kamino.chemistry import get_b_eq, get_k, stoichiometry, K_FUNCTIONS, elements, ION_CHARGE, alk_idx, ca_idx, mg_idx, si_idx, na_idx, al_idx, fe_idx, c_idx
from kamino.precipitation import get_precipitation
from kamino.mineral_info import (basalt_49, clay_minerals, carbonate_minerals,
                                 lt_equilibrium_buffer_minerals, MINERAL_MOLAR_MASS)
from kamino.constants import (
    YR,
    EARTH_ATM,
    EARTH_CRUST_PRODUCTION_RATE_PER_AREA,
    EARTH_HYDROTHERMAL_FLUX_PER_AREA,
    EARTH_CONTINENTAL_WEATHERING_REF,
    ALPHA_REF,
    A_SEAFLOOR_EARTH,
    EARTH_DUST_FLUX_TO_OCEAN,
    COSMIC_DUST_FLUX_PER_AREA,
    SEDIMENT_GRAIN_DENSITY,
    SEDIMENT_POROSITY,
)

def seafloor_reactive_area(T: float, pH: float, rate: float, alpha: float, clog: bool=True, cover: bool=True, sedimentation_rate: float | None = None) -> float:

    T_ref = 280 # Pore space temperature
    T_c = 7 # Activation energy 92 kJ
    pH_ref = 8.1
    t_clog_ref = 20e6 * YR # Coogan & Gillis (2018), Fig 6
    beta = 0 # no pH dependence
    h_cover = 100 # m
    S_ref = EARTH_DUST_FLUX_TO_OCEAN / A_SEAFLOOR_EARTH / SEDIMENT_GRAIN_DENSITY / (1 - SEDIMENT_POROSITY)  # Earth's dust alone, bulk; only if no rate is passed
    S_min = COSMIC_DUST_FLUX_PER_AREA / SEDIMENT_GRAIN_DENSITY / (1 - SEDIMENT_POROSITY)  # floor: cosmic dust, bulk

    t_clog = t_clog_ref * np.exp(- (T - T_ref) / T_c) * (pH / pH_ref) ** beta
    S = max(sedimentation_rate, S_min) if sedimentation_rate is not None else S_ref
    t_cover = h_cover / S

    if clog and cover:
        t_reduce = (1 / t_clog + 1 / t_cover) ** -1
    elif clog:
        t_reduce = t_clog
    elif cover:
        t_reduce = t_cover
    else:
        return alpha

    return alpha * rate * t_reduce * (1 - np.exp(- 1 / (rate * t_reduce)))

def get_weathering_flux(
        P: float,
        T: float,
        P_CO2: float,
        b_input: npt.NDArray[np.float64],
        alpha: float | None=None, rate: float | None=None,
        J: float | None=None,
        sedimentation_rate: float | None=None,
        high_temperature: bool=False,
        crust_composition: dict[str, float]=basalt_49,
        precipitating_minerals: list[str] = clay_minerals + carbonate_minerals,
        cover: bool=True,
        clog: bool=False,
        fO2: float=0,
        water_rock_ratio: float | None=None,
        pe: float | None=None,
        ) -> tuple[npt.NDArray[np.float64], dict[str, float]]:

    if alpha is None:
        alpha = ALPHA_REF

    if rate is None:
        rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA

    if J is None:
        # hydrothermal flux proportional to crust production rate
        J = EARTH_HYDROTHERMAL_FLUX_PER_AREA * (rate / EARTH_CRUST_PRODUCTION_RATE_PER_AREA)

    if len(b_input) == 0:
        b_input = np.zeros(elements.shape)

    if precipitating_minerals is None:
        precipitating_minerals = []

    primary_equilibrium_buffer = ([] if (high_temperature or water_rock_ratio is None) else list(lt_equilibrium_buffer_minerals)) 

    b_eq_primary, pH = get_b_eq(P, T, P_CO2, crust_composition, b_input=b_input, precipitating_minerals=primary_equilibrium_buffer, high_temperature=high_temperature, fO2=fO2, water_rock_ratio=water_rock_ratio, dissolve_only=True, pe=pe)
    k_primary = get_k(P, T, pH, crust_composition)
    k_nonzero = k_primary != 0  # save before replacing zeros with inf
    k_primary = np.where(k_nonzero, k_primary, np.inf)

    A_reactive = seafloor_reactive_area(T, pH, rate, alpha, clog, cover, sedimentation_rate=sedimentation_rate)

    with np.errstate(invalid='ignore', divide='ignore'):

        # Flux formula: F = A*(b_eq - b_in) / (b_eq/k + A/J)
        F_primary = A_reactive * (b_eq_primary - b_input) / (b_eq_primary / k_primary + A_reactive / J)
        F_primary = np.where(k_nonzero, F_primary, 0.0)  # zero out elements with no basalt stoichiometry (e.g. C)

        # Da = k*A / (J*b_eq): dimensionless ratio of transport to kinetic resistance.
        b_eq_safe = np.where(b_eq_primary > 0, b_eq_primary, np.inf)
        Da_primary = (k_primary * A_reactive) / (J * b_eq_safe)

        b_pore = b_input + F_primary / J

    if 'Calcite' not in crust_composition:
        b_pore[c_idx] = b_input[c_idx]

    flux = F_primary

    supply_efficiency = 1 - np.exp(-Da_primary[0])

    weathering_diagnostics = {
        'A_reactive': A_reactive,
        'supply_efficiency': supply_efficiency,
        'b_pore': b_pore,
        'secondary_SI': {},
    }

    if precipitating_minerals and not high_temperature:
        d_b_secondary, _, SI_dict = get_precipitation(P, T, b_pore, precipitating_minerals, [], pe=pe)
        flux += J * d_b_secondary
        b_pore = b_pore + d_b_secondary
        weathering_diagnostics['secondary_SI'] = SI_dict
        weathering_diagnostics['b_pore'] = b_pore

    flux[alk_idx] = float(np.dot(ION_CHARGE, flux))

    # Charge-consistent alkalinity Damkohler
    q_alk = np.where(ION_CHARGE > 0, ION_CHARGE, 0.0)
    q_alk[al_idx] = 0.0
    k_finite = np.where(k_nonzero, k_primary, 0.0)
    k_alk = float(np.dot(q_alk, k_finite))
    b_alk = float(np.dot(q_alk, b_eq_primary))

    Da_alk = (k_alk * A_reactive) / (J * b_alk) if b_alk > 0 else np.nan

    weathering_diagnostics['Da'] = Da_alk

    return flux, weathering_diagnostics


_MG_FRACTION = 0.28  # Mg/(Ca+Mg) in continental silicate weathering — Gaillardet et al. (1999)
_NA_CA_FRACTION = 0.67  # Na/Ca from silicate weathering — Gaillardet et al. (1999), global rivers

def get_continental_weathering_flux(
    T: float,
    P_CO2: float,
    F_alk_ref: float = EARTH_CONTINENTAL_WEATHERING_REF,
    T_ref: float = 288.0,
    P_CO2_ref: float = EARTH_ATM * 280e-6,
    beta: float = 0.3,
    T_e: float = 17.0,
) -> npt.NDArray[np.float64]:
    """Walker-Hays-Kasting continental silicate weathering parameterization.

    Returns flux per unit land area [mol/m²/s] using mixed CaSiO3+MgSiO3 stoichiometry.
    Total Alk flux calibrated to Earth; Ca/Mg split from Gaillardet et al. (1999) river data.

    Parameters
    ----------
    T        : surface temperature [K]
    P_CO2    : atmospheric CO2 partial pressure [Pa]
    F_alk_ref: reference alkalinity flux at (T_ref, P_CO2_ref) [mol_eq/m²/s]
    beta     : CO2 sensitivity exponent (default 0.3)
    T_e      : temperature e-folding scale [K] (default 17 K ~ 70 kJ/mol)
    """
    f = (P_CO2 / P_CO2_ref) ** beta * np.exp((T - T_ref) / T_e)
    F_alk = F_alk_ref * f   # mol_eq / m² / s

    flux = np.zeros(len(elements))

    # Mixed CaSiO3 + MgSiO3 (WHK thermostat sets the Ca+Mg silicate weathering rate),
    # cation split Ca/Mg; F_alk is the Ca+Mg alkalinity (2 per cation -> Ca+Mg = F_alk/2).
    flux[ca_idx]  = F_alk / 2 * (1 - _MG_FRACTION)
    flux[mg_idx]  = F_alk / 2 * _MG_FRACTION
    flux[si_idx]  = F_alk / 2
    # Na from continental feldspar (albite) weathering.
    flux[na_idx]  = flux[ca_idx] * _NA_CA_FRACTION

    # Alkalinity = charge the weathering delivers (2 Ca + 2 Mg + Na). Na is now INCLUDED
    # (was excluded under the balanced-Na assumption), so continental weathering is
    # charge-consistent with the seafloor weathering and ocean alkalinity stays equal to
    # the ion charge balance. Na-silicate weathering is thus an alkalinity source here, to
    # be balanced by the alkalinity removed when Na is sunk (F_na_cont, reverse weathering).
    flux[alk_idx] = float(np.dot(ION_CHARGE, flux))

    return flux


# Graham & Pierrehumbert (2020), ApJ 896, 115, Table 1 and eqs 41-45.
RUNOFF_FRACTION = 0.2             # Gamma, runoff / precipitation (Oki et al. 2001)
PRECIP_REF = 0.99 / YR            # m/s, modern global-mean precipitation (Xie & Arkin 1997)
PRECIP_SENSITIVITY = 0.03         # epsilon, fractional precipitation change per K
T_PRECIP_REF = 288.0              # K
RHO_WATER = 1000.0                # kg/m^3

# Reactive-area scaling for continents; placeholder until calibrated against Earth.
ALPHA_CONT_REF = 1.0

ROCK_DENSITY = 2900.0             # kg/m^3, basalt; ~2700 for granitic crust


def latent_heat_volumetric(T: float) -> float:
    """Latent heat of vaporisation [J/m^3], Henderson-Sellers (1984) via G&P (2020) eq 44."""
    return 1.918e9 * (T / (T - 33.91)) ** 2


def continental_runoff(T: float, S_avg: float | None = None, land_fraction: float | None = None) -> float:
    """Runoff water flux per unit land area [kg/m^2/s], G&P (2020) eqs 41-45.

    S_avg         : globally averaged absorbed instellation (1 - a) S / 4 [W/m^2]; None disables the energetic limit.
    land_fraction : if given, the energetic limit is scaled by (1 - land_fraction) (G&P eq 45, their speculative variant).
    """
    p = max(PRECIP_REF * (1 + PRECIP_SENSITIVITY * (T - T_PRECIP_REF)), 0.0)
    if S_avg is not None:
        p_lim = S_avg / latent_heat_volumetric(T)
        if land_fraction is not None:
            p_lim *= (1 - land_fraction)
        p = min(p, p_lim)
    return RUNOFF_FRACTION * p * RHO_WATER


def continental_reactive_area(rate: float, alpha: float) -> float:
    """Reactive area per unit land area, proportional to crust production rate."""
    return alpha * rate / EARTH_CRUST_PRODUCTION_RATE_PER_AREA


def supply_limited_k(T: float, pH: float, composition: dict[str, float], A: float, erosion_rate: float, rho: float = ROCK_DENSITY) -> npt.NDArray[np.float64]:
    """get_k with each mineral capped by erosional supply: 1/r_i = 1/(A k_i x_i) + 1/(E rho x_i / M_i), returned per unit A."""
    k = np.zeros(elements.shape)
    for mineral, x in composition.items():
        r_kin = A * float(K_FUNCTIONS[mineral](T, pH)) * x
        r_sup = erosion_rate * rho * x / (MINERAL_MOLAR_MASS.get(mineral, 150.0) / 1000.0)
        if r_kin > 0 and r_sup > 0:
            k += stoichiometry[mineral] / (1 / r_kin + 1 / r_sup)
    return k / A if A > 0 else k


def get_continental_weathering_flux_mac(
        P: float,
        T: float,
        P_CO2: float,
        alpha: float | None = None,
        rate: float | None = None,
        J: float | None = None,
        S_avg: float | None = None,
        land_fraction: float | None = None,
        crust_composition: dict[str, float] = basalt_49,
        precipitating_minerals: list[str] = clay_minerals,
        kinetic: bool = False,
        retain_al_fe: bool = True,
        erosion_rate: float | None = None,
        rock_density: float = ROCK_DENSITY,
        fO2: float = 0,
        water_rock_ratio: float | None = None,
        pe: float | None = None,
        ) -> tuple[npt.NDArray[np.float64], dict[str, float]]:
    """Lithology-aware continental weathering, per unit land area [mol/m^2/s].

    Same law as get_weathering_flux, with the seafloor A_r replaced by a crust-production
    proportionality and the hydrothermal J replaced by G&P (2020) runoff. Input water is rain
    and k is evaluated at its pH in both modes. kinetic=True drops the runoff /
    thermodynamic limit, F = A k, i.e. the J -> inf limit of the default.
    erosion_rate [m/s], if given, caps each mineral at its erosional supply (West et al. 2005).
    """
    if alpha is None:
        alpha = ALPHA_CONT_REF

    if rate is None:
        rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA

    if J is None:
        J = continental_runoff(T, S_avg, land_fraction)

    if precipitating_minerals is None:
        precipitating_minerals = []

    A_reactive = continental_reactive_area(rate, alpha)

    weathering_diagnostics = {'A_reactive': A_reactive, 'runoff': J, 'secondary_SI': {}}

    if J <= 0:
        # No runoff (T below p = 0): nothing is carried off the continents.
        weathering_diagnostics.update({'b_pore': np.zeros(elements.shape), 'Da': np.nan, 'pH': np.nan})
        return np.zeros(elements.shape), weathering_diagnostics

    # Input is rainwater (pure water + P_CO2). k is taken at its pH: at Da << 1 the fluid on the rock
    # stays near input chemistry, and this gives k the positive pCO2 dependence of G&P (2020) eq 40.
    b_input, pH = get_b_eq(P, T, P_CO2, {}, b_input=np.zeros(elements.shape), pe=pe)
    if erosion_rate is None:
        k_primary = get_k(P, T, pH, crust_composition)
    else:
        k_primary = supply_limited_k(T, pH, crust_composition, A_reactive, erosion_rate, rock_density)
    k_nonzero = k_primary != 0

    b_eq_primary, pH_eq = (None, np.nan) if kinetic else get_b_eq(P, T, P_CO2, crust_composition, b_input=b_input, fO2=fO2, water_rock_ratio=water_rock_ratio, dissolve_only=True, pe=pe)

    with np.errstate(invalid='ignore', divide='ignore'):
        if kinetic:
            F_primary = A_reactive * k_primary
        else:
            # F = A*b_eq / (b_eq/k + A/J): harmonic mean of kinetic (A k) and runoff (J b_eq) limits.
            k_safe = np.where(k_nonzero, k_primary, np.inf)
            F_primary = A_reactive * (b_eq_primary - b_input) / (b_eq_primary / k_safe + A_reactive / J)
        F_primary = np.where(k_nonzero, F_primary, 0.0)

    b_pore = b_input + F_primary / J
    flux = F_primary.copy()

    # Secondary clays in the regolith (Kaolinite, Goethite) remove Al and Fe, and with them their charge.
    if precipitating_minerals:
        d_b_secondary, _, SI_dict = get_precipitation(P, T, b_pore, precipitating_minerals, [], pe=pe)
        flux += J * d_b_secondary
        b_pore = b_pore + d_b_secondary
        weathering_diagnostics['secondary_SI'] = SI_dict

    # Al and Fe stay in the regolith (Nesbitt & Young 1982); the clay step alone misses them in dilute runoff.
    if retain_al_fe:
        flux[[al_idx, fe_idx]] = 0.0

    flux[alk_idx] = float(np.dot(ION_CHARGE, flux))

    # Charge-consistent alkalinity Damkohler, as in get_weathering_flux; undefined in kinetic mode.
    q_alk = np.where(ION_CHARGE > 0, ION_CHARGE, 0.0)
    q_alk[al_idx] = 0.0
    if retain_al_fe:
        q_alk[fe_idx] = 0.0
    k_alk = float(np.dot(q_alk, np.where(k_nonzero, k_primary, 0.0)))
    if kinetic:
        Da_alk = np.nan
    else:
        b_alk = float(np.dot(q_alk, b_eq_primary))
        Da_alk = (k_alk * A_reactive) / (J * b_alk) if b_alk > 0 else np.nan

    weathering_diagnostics.update({'b_pore': b_pore, 'Da': Da_alk, 'pH': pH, 'pH_eq': pH_eq, 'F_kinetic_alk': A_reactive * k_alk})

    return flux, weathering_diagnostics
