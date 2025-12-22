from kamino.ocean_chemistry.aqueous_geochemistry import *

def calcite_precipitation_rate(P: float, T: float, alkalinity: float, DIC: float, Ca: float, Mg: float=0, Fe: float=0) -> tuple[float, float]:

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

def siderite_precipitation_rate(P: float, T: float, alkalinity: float, DIC: float, Ca: float, Mg: float=0, Fe: float=0) -> tuple[float, float]:

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