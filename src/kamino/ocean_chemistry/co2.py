from kamino.ocean_chemistry.aqueous_geochemistry import *

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

    if T < -ABSOLUTE_ZERO + 0.1:
        return 0 # if the temperature is 0 C, the planet has frozen and no CO2 is in the atmosphere? 

    input_lines = solution_block(P, T, composition, None) + output_block(saturation_indexes=['CO2(g)'])
    output = run_PHREEQC(input_lines, single_output=True)

    SI = get_output_saturation_indexes(output, get_initial=True)
    CO2_SI = SI['CO2(g)']

    return (10 ** CO2_SI) * P