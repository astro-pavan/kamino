import numpy as np

import kamino.H21.chili.equilibrium as eq
import kamino.H21.chili.kinetics as ki
import kamino.H21.chili.parameters as pr
import kamino.H21.chili.climate as cl
from kamino.utils import *

from kamino.constants import YR, POROSITY, R_EARTH

KeqFuncs   = eq.import_thermo_data(eq.DATABASE_DIR / 'species.csv')
DICeqFuncs = eq.get_DICeq(pr.xCO2, pr.T, pr.Pfull, KeqFuncs)

logkDict   = ki.import_kinetics_data()
kFuncs     = ki.get_keff(pr.T, pr.pHfull, logkDict)

basalt_composition = ['woll','enst','ferr','anoh','albh']
granite_composition = ['albi', 'kfel', 'phlo', 'anni', 'quar']

# Mineral fractions by mass, mirroring ocean_chem's basalt/granite compositions.
# Used when weighted_sum=True in get_k_eff.
basalt_composition_weighted = {
    'woll': 0.1,  # Wollastonite
    'enst': 0.2,  # Enstatite
    'anoh': 0.3,  # Anorthite (high-T)
    'albh': 0.2,  # Albite (high-T)
    'fors': 0.1,  # Forsterite
    'faya': 0.1,  # Fayalite
}

granite_composition_weighted = {
    'albi': 0.25,  # Albite
    'kfel': 0.25,  # K-feldspar
    'phlo': 0.20,  # Phlogopite
    'anni': 0.20,  # Annite
    'quar': 0.10,  # Quartz
}

# mol alkalinity produced per mol mineral dissolved.
# Derived from H+ stoichiometry of dissolution reactions (Phreeqc / ocean_chem database).
alk_stoich = {
    'woll': 2.0,  # CaSiO3 + 2H+ -> Ca2+ + SiO2 + H2O
    'enst': 2.0,  # MgSiO3 + 2H+ -> Mg2+ + SiO2 + H2O
    'ferr': 2.0,  # FeSiO3 + 2H+ -> Fe2+ + SiO2 + H2O
    'anoh': 8.0,  # CaAl2Si2O8 + 8H+ -> Ca2+ + 2Al3+ + 2SiO2 + 4H2O
    'albh': 4.0,  # NaAlSi3O8 + 4H+ -> Na+ + Al3+ + 3SiO2 + 2H2O
    'albi': 4.0,
    'fors': 4.0,  # Mg2SiO4 + 4H+ -> 2Mg2+ + SiO2 + 2H2O
    'faya': 4.0,  # Fe2SiO4 + 4H+ -> 2Fe2+ + SiO2 + 2H2O
    'kfel': 4.0,  # KAlSi3O8 + 4H+ -> K+ + Al3+ + 3SiO2 + 2H2O
    'phlo': 8.0,  # KMg3AlSi3O10(OH)2 + 8H+ -> K+ + 3Mg2+ + Al3+ + 3SiO2 + 5H2O
    'anni': 8.0,  # KFe3AlSi3O10(OH)2 + 8H+ -> K+ + 3Fe2+ + Al3+ + 3SiO2 + 5H2O
    'quar': 0.0,  # SiO2: no alkalinity production
}

def get_C_eq(P: float, T: float, x_CO2: float, granite=False):
    
    lithology = 'grah' if granite else 'bash'
    
    arg = np.array((x_CO2, T, P))
    C_eq = DICeqFuncs[lithology]['ALK'](arg) * 1000 # convert to mol/m^3
    return C_eq

def get_k_eff(P: float, T: float, x_CO2: float, granite=False,
              weighted_sum=False, alk_correction=False):
    """
    Returns the effective kinetic rate coefficient k_eff [mol m⁻² yr⁻¹].

    Parameters
    ----------
    weighted_sum : bool
        If False (default), use the original CHILI approach: k_eff = min over
        all minerals (slowest mineral is rate-limiting).
        If True, use a composition-weighted sum over minerals (parallel reactions),
        matching the approach in ocean_chem's get_k.
    alk_correction : bool
        If False (default), k_eff is in mol mineral m⁻² yr⁻¹ (no stoichiometry).
        If True, multiply by mol_alk / mol_mineral so that k_eff is in
        mol alkalinity m⁻² yr⁻¹, making it dimensionally consistent with C_eq.
        For weighted_sum=True: each mineral's rate is multiplied by its own stoich.
        For weighted_sum=False: the minimum k is multiplied by the
        composition-weighted mean stoichiometry.
    """

    lithology = 'grah' if granite else 'bash'

    arg = np.array((x_CO2, T, P))
    pH = DICeqFuncs[lithology]['pH'](arg)

    if weighted_sum:
        composition = granite_composition_weighted if granite else basalt_composition_weighted
        k_eff = 0.0
        for mineral, fraction in composition.items():
            k_mineral = kFuncs[mineral](T, pH)
            stoich = alk_stoich[mineral] if alk_correction else 1.0
            k_eff += fraction * k_mineral * stoich
    else:
        composition = granite_composition if granite else basalt_composition
        k_eff = -1
        for mineral in composition:
            if k_eff == -1:
                k_eff = kFuncs[mineral](T, pH)
            else:
                k_eff = smooth_min(kFuncs[mineral](T, pH), k_eff)
        if alk_correction:
            # Weighted-average stoichiometry over the weighted composition
            ref_composition = granite_composition_weighted if granite else basalt_composition_weighted
            mean_stoich = sum(f * alk_stoich[m] for m, f in ref_composition.items())
            k_eff = k_eff * mean_stoich

    return k_eff

def get_Dw(P: float, T: float, x_CO2: float, flow_path_length: float, rock_age: float, granite=False):

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 100 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 1 # all minerals are considered reactive
    porosity = POROSITY

    C_eq = get_C_eq(P, T, x_CO2)
    k_eff = get_k_eff(P, T, x_CO2)
    
    psi = flow_path_length * (1 - porosity) * rock_density * specific_surface_area * fresh_mineral_fraction
    
    Dw = psi / (C_eq * (k_eff ** -1 + mean_molar_mass * specific_surface_area * rock_age))

    return Dw

def w_kinetic(P: float, T: float, x_CO2: float, flow_path_length: float, rock_age: float, granite=False):

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 100 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 1
    porosity = POROSITY

    k_eff = get_k_eff(P, T, x_CO2)
    psi = flow_path_length * (1 - porosity) * rock_density * specific_surface_area * fresh_mineral_fraction

    return (k_eff * psi) / (1 + mean_molar_mass * specific_surface_area * k_eff * rock_age)

def w_thermodynamic(P: float, T: float, x_CO2: float, runoff: float, granite=False):

    C_eq = get_C_eq(P, T, x_CO2)

    return runoff * C_eq[0]

def w_supply(flow_path_length: float, rock_age: float):

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 100 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 1
    porosity = POROSITY

    psi = flow_path_length * (1 - porosity) * rock_density * specific_surface_area * fresh_mineral_fraction

    return psi / (mean_molar_mass * specific_surface_area * rock_age)


def get_weathering_rate(P: float, T: float, x_CO2: float, runoff: float, flow_path_length: float, rock_age: float, granite=False,
                        weighted_sum=False, alk_correction=False) -> float:
    """
    Calculates the basalt seafloor weathering rate, giving an alkalinity production rate.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    x_CO2 : float
        Atmospheric CO2 fraction.
    runoff : float
        Fluid flow rate in seafloor pore space in m/yr.
    flow_path_length : float
        Length of flow path through seafloor pore space.
    rock_age : float
        Age of rocks in the pore space in yr.
    weighted_sum : bool
        Passed to get_k_eff. If True, use a composition-weighted sum of mineral
        rates instead of the minimum (see get_k_eff for details).
    alk_correction : bool
        Passed to get_k_eff. If True, scale k_eff by mol_alk / mol_mineral so
        it is dimensionally consistent with C_eq (see get_k_eff for details).

    Returns
    -------
    float
        Alkalinity production rate in mol / m^2 / yr
    """

    P = float(P / 1e5) # convert to bar
    P = np.clip(P, pr.P.min(), pr.P.max())
    x_CO2 = np.clip(x_CO2, pr.xCO2.min(), pr.xCO2.max())
    T = np.clip(T, pr.T.min(), pr.T.max())

    C_eq = get_C_eq(P, T, x_CO2, granite)
    k_eff = get_k_eff(P, T, x_CO2, granite, weighted_sum=weighted_sum, alk_correction=alk_correction)

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 0.005 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 0.1
    porosity = POROSITY
    
    psi = flow_path_length * (1 - porosity) * rock_density * specific_surface_area * fresh_mineral_fraction
    
    Dw = psi / (C_eq * (k_eff ** -1 + mean_molar_mass * specific_surface_area * rock_age))
    
    C = C_eq / (1 + (runoff / Dw))

    w = runoff * C

    return float(w[0])

def get_weathering_rate_KT18(P: float, T: float, x_CO2: float) -> float:
    """
    Calculates the basalt seafloor weathering rate, giving an alkalinity production rate with the older KT18 method.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    x_CO2 : float
        Atmospheric CO2 fraction.

    Returns
    -------
    float
        Alkalinity production rate in mol / m^2 / yr
    """

    P_CO2 = (P / 1e5) * x_CO2 # P_CO2 in bar

    return cl.seaf_brad1997(T, P_CO2) * 1e12 / (4 * np.pi * R_EARTH ** 2)

def get_weathering_rate_WHAK(P: float, T: float, x_CO2: float) -> float:
    """
    Calculates the basalt seafloor weathering rate, giving an alkalinity production rate with the older KT18 method.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    x_CO2 : float
        Atmospheric CO2 fraction.

    Returns
    -------
    float
        Alkalinity production rate in mol / m^2 / yr
    """

    P_CO2 = (P / 1e5) * x_CO2 # P_CO2 in bar

    return cl.cont_walk1981(T, P_CO2) * 1e12 / (4 * np.pi * R_EARTH ** 2)

# MAC weathering

mu = np.exp(2)
# A: specific surface area (m^2.kg^-1)
A = 100
# L: reactive flow path length (m)
L = 1
# phi: porosity, of soil
phi = 0.1
# t_s: soil age (year)
t_s = 1e5
# rho_sf: mineral mass to fluid volume ratio (kg.m^-3)
rho_sf = 12728
# keff_ref: reference rate constant for mineral dissolution (mol.m^-2.year^-1)
keff_ref = 8.7e-6
# T_e: Kinetic weathering temperature dependence (K)
Te = 11.1
# T_ref: modern global average temperature (K)
T_ref = 288
# pco2_ref: pre-industrial co2 (bar)
pco2_ref = 280e-6
# beta: kinetic weathering pco2 dependence
beta = 0.2
# m: mineral molar mass (kg.mol^-1)
m = 0.27
# X_r: reactive mineral concentration in fresh rock
X_r = 0.36
# lmbda: thermodynamic coefficient for C_eq, variable by orders of magnitude
lmbda = 1.4e0
# 
# n: thermodynamic pco2 dependence of C_eq
n = 0.316

q_ref = 0.2
# p_ref: global average precipitation per unit area (m.yr^-1)
p_ref = 0.99
# Gamma: fraction of precipitation that becomes runoff = q_ref/p_ref = 0.2
Gamma = q_ref / p_ref
# eps: fractional change in precipitation per K change in temperature (1/K)
eps = 0.03

def get_weathering_rate_MAC(T, pco2):

    pco2 = pco2 * 1e-5

    def Ceq(pco2):
        return lmbda * pco2**n
    
    q = Gamma * p_ref * (1 + eps * (T - T_ref))

    alpha = L * phi * rho_sf * A * X_r * mu
    top = alpha
    keff = keff_ref * np.exp((T-T_ref)/Te) * (pco2/pco2_ref)**beta
    bottom = (keff**(-1) + m*A*t_s + alpha/(q*Ceq(pco2)))

    return top/bottom
