peridotite_42 = {
    'Forsterite': 0.7,
    'Fayalite': 0.3
} # 100% Olivine

peridotite_44 = {
    'Wollastonite': 0.333 * 0.3,
    'Enstatite': 0.666 * 0.3,
    'Forsterite': 0.7 * 0.7,
    'Fayalite': 0.3 * 0.7
} # 70% Olivine 30% Pyroxene

peridotite_47 = {
    'Wollastonite': 0.333 * 0.7,
    'Enstatite': 0.666 * 0.7,
    'Forsterite': 0.7 * 0.3,
    'Fayalite': 0.3 * 0.3
} # 30% Olivine 70% Pyroxene

basalt_49 = {
    'Anorthite': 0.1,
    'Wollastonite': 0.333 * 0.8,
    'Enstatite': 0.666 * 0.8,
    'Forsterite': 0.7 * 0.1,
    'Fayalite': 0.3 * 0.1
} # 10% Olivine 80% Pyroxene 10% Plagioclase

basalt_51 = {
    'Anorthite': 0.3,
    'Wollastonite': 0.333 * 0.7,
    'Enstatite': 0.666 * 0.7,
} # 70% Pyroxene 30% Plagioclase

basalt_composition = {
    'Wollastonite': 0.1,
    'Enstatite': 0.2,
    'Anorthite': 0.5,
    # 'Albite': 0.2,
    'Forsterite': 0.18,
    'Fayalite': 0.02
}

hydrothermal_minerals = [
    'Albite', 'Anorthite', 'Enstatite', 'Forsterite',
    'Fayalite', 'Wollastonite', 'K-Feldspar', 'Diopside',
] # All minerals in hydrothermal database

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)

carbonate_minerals = ['Calcite', 'Siderite'] # Ca and Fe sinks
clay_minerals = ['Kaolinite', 'Goethite'] # Al, Fe sinks
silica_minerals = ['SiO2(am)'] # Si sink
reverse_weathering_minerals = ['Saponite-Mg'] # Mg sink (via reverse weathering)