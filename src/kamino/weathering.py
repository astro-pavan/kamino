import numpy as np
import numpy.typing as npt
from scipy.optimize import least_squares 

from kamino.chemistry import get_b_eq, get_k, elements, alk_idx, c_idx, ca_idx, mg_idx, si_idx, na_idx, cl_idx, so4_idx
from kamino.precipitation import get_precipitation
from kamino.mineral_info import *
from kamino.constants import YR, R_EARTH, EARTH_ATM

rate_ref = 1 / (50e6 * YR)
J_ref = 1.4e15 / YR # kg / yr
A_seafloor = 0.7 * 4 * np.pi * R_EARTH ** 2
J_ref_normalised = J_ref / A_seafloor

def seafloor_reactive_area(T: float, pH: float, rate: float, alpha: float, clog: bool=True, cover: bool=True, sedimentation_rate: float | None = None) -> float:

    T_ref = 280 # Pore space temperature
    T_c = 7 # Activation energy 92 kJ
    pH_ref = 8.1
    t_clog_ref = 20e6 * YR # Coogan & Gillis (2018), Fig 6
    beta = 0 # no pH dependence
    h_cover = 100 # m
    S_ref = 5 / (1e6 * YR) # 5 m / Myr
    S_min = 0.3 / (1e6 * YR) # 0.3 m / Myr background (dust)

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
        ) -> tuple[npt.NDArray[np.float64], dict[str, float]]:

    if alpha is None:
        alpha = ALPHA_REF

    if rate is None:
        rate = rate_ref

    if J is None:
        J = J_ref_normalised * (rate / rate_ref) # hydrothermal flux proportional to crust production rate

    if len(b_input) == 0:
        b_input = np.zeros(elements.shape)

    if precipitating_minerals is None:
        precipitating_minerals = []

    b_eq_primary, pH = get_b_eq(P, T, P_CO2, crust_composition, b_input=b_input, high_temperature=high_temperature, fO2=fO2, water_rock_ratio=water_rock_ratio)
    k_primary = get_k(P, T, pH, crust_composition)
    k_nonzero = k_primary != 0  # save before replacing zeros with inf
    k_primary = np.where(k_nonzero, k_primary, np.inf)

    A_reactive = seafloor_reactive_area(T, pH, rate, alpha, clog, cover, sedimentation_rate=sedimentation_rate)

    # Flux formula: F = A*(b_eq - b_in) / (b_eq/k + A/J)
    F_primary = A_reactive * (b_eq_primary - b_input) / (b_eq_primary / k_primary + A_reactive / J)
    F_primary = np.where(k_nonzero, F_primary, 0.0)  # zero out elements with no basalt stoichiometry (e.g. C)

    # Da = k*A / (J*b_eq): dimensionless ratio of transport to kinetic resistance.
    # Da>>1 → thermodynamically limited (pore fluid near equilibrium with basalt);
    # Da<<1 → kinetically limited.
    b_eq_safe = np.where(b_eq_primary > 0, b_eq_primary, np.inf)
    Da_primary = (k_primary * A_reactive) / (J * b_eq_safe)
    b_pore = b_input + (b_eq_primary - b_input) * (1 - np.exp(-Da_primary))

    if 'Calcite' not in crust_composition:
        b_pore[c_idx] = b_input[c_idx]

    flux = F_primary

    supply_efficiency = 1 - np.exp(-Da_primary[0])

    weathering_diagnostics = {
        'Da': Da_primary[0],
        'A_reactive': A_reactive,
        'supply_efficiency': supply_efficiency,
        'b_pore': b_pore,
        'secondary_SI': {},
    }

    if precipitating_minerals and not high_temperature:
        d_b_secondary, _, SI_dict = get_precipitation(P, T, b_pore, precipitating_minerals, [])
        flux += J * d_b_secondary
        b_pore = b_pore + d_b_secondary
        weathering_diagnostics['secondary_SI'] = SI_dict
        weathering_diagnostics['b_pore'] = b_pore
    
    return flux, weathering_diagnostics

# Calibration at actual Earth pore conditions: T_surface=288K gives T_sf=277K, T_pore=286K.
# Using modern seawater concentrations as b_input ensures ALPHA_REF produces the observed
# ~1 Tmol/yr seafloor Alk flux at Earth steady state, not at an empty-ocean reference.
# Secondary precipitation excluded (precipitating_minerals=[]) to isolate primary dissolution.
_T_ref_calib = 286
_P_ref_calib = 1000 * 10 * 3000
_P_CO2_ref_calib = EARTH_ATM * 280e-6
_seafloor_flux = 1e12 / YR
_seafloor_flux_normalised = _seafloor_flux / A_seafloor

_calib_b_input = np.zeros(elements.shape)
_calib_b_input[alk_idx]  = 2.3e-3
_calib_b_input[ca_idx]   = 10.3e-3
_calib_b_input[mg_idx]   = 52.8e-3
_calib_b_input[na_idx]   = 480e-3
_calib_b_input[cl_idx]   = 550e-3
_calib_b_input[so4_idx]  = 28e-3
_calib_b_input[si_idx]   = 0.1e-3
_calib_b_input[c_idx]    = 2.0e-3

def _calibration_residual(a_array):
    a_val = a_array[0]
    flux, _ = get_weathering_flux(
        _P_ref_calib, _T_ref_calib, _P_CO2_ref_calib,
        _calib_b_input, alpha=a_val, rate=rate_ref, precipitating_minerals=[],
    )
    alkalinity_flux = flux[0]
    return (alkalinity_flux - _seafloor_flux_normalised) / _seafloor_flux_normalised

# Run the optimization
_root = least_squares(_calibration_residual, [100.0])
ALPHA_REF = float(_root.x[0])

# Reference alkalinity flux per unit land area at modern Earth conditions.
# Calibrated so that 0.3 land fraction gives ~8 Tmol eq/yr total, which
# balances modern volcanic outgassing after accounting for seafloor weathering.
_CONT_LAND_AREA_EARTH = 0.3 * 4 * np.pi * R_EARTH ** 2   # m²
EARTH_CONTINENTAL_WEATHERING_REF = (8e12 / YR) / _CONT_LAND_AREA_EARTH  # mol_eq / m² / s

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

    # Mixed CaSiO3 + MgSiO3 — same total Alk, cation split Ca/Mg.
    flux[alk_idx] = F_alk
    flux[ca_idx]  = F_alk / 2 * (1 - _MG_FRACTION)
    flux[mg_idx]  = F_alk / 2 * _MG_FRACTION
    flux[si_idx]  = F_alk / 2
    # Na from continental feldspar (albite) weathering. No Alk added: balanced Na cycle assumption
    # (Coogan 2022) — the HCO3- produced when Na-silicate weathers is consumed when Na is
    # eventually removed from the ocean via reverse weathering or subduction.
    flux[na_idx]  = flux[ca_idx] * _NA_CA_FRACTION

    return flux