from kamino.ocean_chemistry.aqueous_geochemistry import *
from kamino.utils import *

import PyCO2SYS as pyco2

def get_P_CO2(P: float, T: float, alkalinity: float, DIC: float, Ca: float=0, Mg: float=0, Fe: float=0) -> float:
    """_summary_

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    alkalinity : float
        Alkalinity in mol/kgw.
    DIC : float
        Dissolved inorganic carbon in mol/kgw.
    Ca : float
        Calcium concentration in mol/kgw, by default 0.
    Mg : float, optional
        Magnesium concentration in mol/kgw, by default 0.
    Fe : float, optional
        Iron concentration in mol/kgw, by default 0.

    Returns
    -------
    float
        Partial pressure of CO2 in Pa.
    """

    composition = {
        'Ca' : Ca,
        'Mg' : Mg,
        'Fe' : Fe,
        'C'  : DIC,
        'Alkalinity' : alkalinity
    }

    # if T < -ABSOLUTE_ZERO + 0.1:
    #     return 0 # if the temperature is 0 C, the planet has frozen and no CO2 is in the atmosphere? 

    T = smooth_max(T, 273.5)
    T = smooth_min(T, 500)

    input_lines = solution_block(P, T, composition, None) + output_block(saturation_indexes=['CO2(g)'])
    output = run_PHREEQC(input_lines, single_output=True)

    SI = get_output_saturation_indexes(output, get_initial=True)
    CO2_SI = SI['CO2(g)']

    return (10 ** CO2_SI) * P

def get_P_CO2_v2(P, T, alk, DIC):

    known_alkalinity = alk * 1e6  # umol/kg (example ocean value)
    known_dic = DIC * 1e6         # umol/kg (example ocean value)

    temperature_surface = T + ABSOLUTE_ZERO # degrees Celsius
    pressure_surface = 0     # dbar (0 dbar is roughly surface pressure)
    salinity = 35            # PSU (Standard ocean salinity)

    temperature_surface = smooth_max(0.5, temperature_surface)
    temperature_surface = smooth_min(99.5, temperature_surface)

    # --- The Calculation ---
    results = pyco2.sys(
        par1=known_alkalinity,  # Value of first parameter
        par2=known_dic,         # Value of second parameter
        par1_type=1,            # 1 tells code par1 is Alkalinity
        par2_type=2,            # 2 tells code par2 is DIC
        salinity=salinity,      # Required for equilibrium constants (K values)
        temperature=temperature_surface,
        pressure=pressure_surface
    )

    # --- Extracting the Answer ---
    # The results dictionary contains calculated values for everything else.
    # 'pCO2' is the partial pressure in micro-atmospheres (uatm).

    pCO2_ocean = results['pCO2'] * EARTH_ATM

    return pCO2_ocean