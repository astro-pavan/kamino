from kamino.ocean_chemistry.aqueous_geochemistry import *

def calcite_precipitation_rate(P: float, T: float, alkalinity: float, DIC: float, Ca: float, Mg: float=0, Fe: float=0) -> tuple[float, float]:
    """
    Calculates calcite preciptation rate and calcite stauration index.

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
        Calcium concentration in mol/kgw.
    Mg : float, optional
        Magnesium concentration in mol/kgw, by default 0.
    Fe : float, optional
        Iron concentration in mol/kgw, by default 0.

    Returns
    -------
    tuple[float, float]
        Calcite precipitation rate in mol/s/kgw, Calcite saturation index
    """

    composition = {
        'Ca' : Ca,
        'Mg' : Mg,
        'Fe' : Fe,
        'C'  : DIC,
        'Alkalinity' : alkalinity
    }

    output = kinetics(P, T, composition, None, {'Calcite': 1e-10}, 1, 100)
    k = get_output_kinetics_phases(output)['Calcite']
    SI = get_output_saturation_indexes(output)['Calcite']

    return np.maximum(k, 0), SI

def magnesite_precipitation_rate(P: float, T: float, alkalinity: float, DIC: float, Ca: float, Mg: float, Fe: float) -> tuple[float, float]:
    """
    Calculates magnestite preciptation rate and calcite stauration index.

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
        Calcium concentration in mol/kgw.
    Mg : float, optional
        Magnesium concentration in mol/kgw.
    Fe : float, optional
        Iron concentration in mol/kgw.

    Returns
    -------
    tuple[float, float]
        Magnesite precipitation rate in mol/s/kgw, Magnesite saturation index
    """

    composition = {
        'Ca' : Ca,
        'Mg' : Mg,
        'Fe' : Fe,
        'C'  : DIC,
        'Alkalinity' : alkalinity
    }

    output = kinetics(P, T, composition, None, {'Magnesite': 1e-10}, 1, 100)
    k = get_output_kinetics_phases(output)['Magnesite']
    SI = get_output_saturation_indexes(output)['Magnesite']

    return np.maximum(k, 0), SI

def siderite_precipitation_rate(P: float, T: float, alkalinity: float, DIC: float, Ca: float, Mg: float, Fe: float) -> tuple[float, float]:
    """
    Calculates siderite preciptation rate and calcite stauration index.

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
        Calcium concentration in mol/kgw.
    Mg : float, optional
        Magnesium concentration in mol/kgw.
    Fe : float, optional
        Iron concentration in mol/kgw.

    Returns
    -------
    tuple[float, float]
        Siderite precipitation rate in mol/s/kgw, Siderite saturation index
    """

    composition = {
        'Ca' : Ca,
        'Mg' : Mg,
        'Fe' : Fe,
        'C'  : DIC,
        'Alkalinity' : alkalinity
    }

    output = kinetics(P, T, composition, None, {'Siderite': 1e-10}, 1, 100)
    k = get_output_kinetics_phases(output)['Siderite']
    SI = get_output_saturation_indexes(output)['Siderite']

    return np.maximum(k, 0), SI