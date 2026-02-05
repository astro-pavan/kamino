import numpy as np

import kamino.seafloor_weathering.chili.equilibrium as eq
import kamino.seafloor_weathering.chili.kinetics as ki
import kamino.seafloor_weathering.chili.parameters as pr
import kamino.seafloor_weathering.chili.climate as cl
from kamino.utils import *

from kamino.constants import YR, POROSITY, R_EARTH

KeqFuncs   = eq.import_thermo_data(eq.DATABASE_DIR / 'species.csv')
DICeqFuncs = eq.get_DICeq(pr.xCO2, pr.T, pr.Pfull, KeqFuncs)

logkDict   = ki.import_kinetics_data()
kFuncs     = ki.get_keff(pr.T, pr.pHfull, logkDict)

basalt_composition = ['woll','enst','ferr','anoh','albh']
granite_composition = ['albi', 'kfel', 'phlo', 'anni', 'quar']

def get_C_eq(P: float, T: float, x_CO2: float, granite=False):
    
    lithology = 'grah' if granite else 'bash'
    
    arg = np.array((x_CO2, T, P))
    C_eq = DICeqFuncs[lithology]['ALK'](arg) * 1000 # convert to mol/m^3
    return C_eq

def get_k_eff(P: float, T: float, x_CO2: float, granite=False):

    lithology = 'grah' if granite else 'bash'

    arg = np.array((x_CO2, T, P))
    pH = DICeqFuncs[lithology]['pH'](arg)
    k_eff = -1

    composition = granite_composition if granite else basalt_composition

    for mineral in composition:
        if k_eff == -1:
            k_eff = kFuncs[mineral](T, pH)
        else:
            k_eff = smooth_min(kFuncs[mineral](T, pH), k_eff)

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


def get_weathering_rate(P: float, T: float, x_CO2: float, runoff: float, flow_path_length: float, rock_age: float, granite=False) -> float:
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
    k_eff = get_k_eff(P, T, x_CO2, granite)

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 100 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 1
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