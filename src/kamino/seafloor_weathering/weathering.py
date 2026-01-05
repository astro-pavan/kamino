import numpy as np

import kamino.seafloor_weathering.chili.equilibrium as eq
import kamino.seafloor_weathering.chili.kinetics as ki
import kamino.seafloor_weathering.chili.parameters as pr

from kamino.constants import YR

KeqFuncs   = eq.import_thermo_data(eq.DATABASE_DIR / 'species.csv')
DICeqFuncs = eq.get_DICeq(pr.xCO2, pr.T, pr.Pfull, KeqFuncs)

logkDict   = ki.import_kinetics_data()
kFuncs     = ki.get_keff(pr.T, pr.pHfull, logkDict)

basalt_composition = ['woll','enst','ferr','anoh','albh']

def weathering_rate(P: float, T: float, x_CO2: float, runoff: float, flow_path_length: float, rock_age: float) -> float:
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
        Fluid flow rate in seafloor pore space in m/s.
    flow_path_length : float
        Length of flow path through seafloor pore space.
    rock_age : float
        Age of rocks in the pore space in s.

    Returns
    -------
    float
        Alkalinity production rate in mol / m^2 / s
    """

    P = P / 1e5 # convert to bar

    arg = np.array((x_CO2, T, P))
    pH = DICeqFuncs['bash']['pH'](arg)
    C_eq = DICeqFuncs['bash']['ALK'](arg) * 1000 # convert to mol/m^3

    k_eff = 1e99
    for mineral in basalt_composition:
        k_eff = np.minimum(kFuncs[mineral](T, pH), k_eff)

    k_eff = k_eff / YR # convert to per s rather than per year

    mean_molar_mass = 0.216 # kg / mol
    specific_surface_area = 100 # m^2 / kg
    rock_density = 2700 # kg / m^3
    fresh_mineral_fraction = 1
    porosity = 0.1
    
    psi = flow_path_length * (1 - porosity) * rock_density * specific_surface_area * fresh_mineral_fraction
    
    Dw = psi / (C_eq * (k_eff ** -1 + mean_molar_mass * specific_surface_area * rock_age))
    
    C = C_eq / (1 + (runoff / Dw))

    w = runoff * C

    return w