basalt_composition = {
    'Wollastonite': 0.1,
    'Enstatite': 0.2,
    'Anorthite': 0.5,
    # 'Albite': 0.2,
    'Forsterite': 0.18,
    'Fayalite': 0.02
}

secondary_composition = {
    'Smectite-high-Fe-Mg': 0.35,
    'Kaolinite': 0.25,
    'Calcite': 0.15,
    'Goethite': 0.10,
    'Sepiolite': 0.10,
    'Magnesite': 0.05,
}

secondary_sinks = [
    'Calcite', 'Magnesite', 'Siderite', 'Dolomite', 'Gypsum', 
    'Kaolinite', 'Goethite', 'Smectite-high-Fe-Mg', 'Sepiolite'
]

hydrothermal_minerals = [
    'Albite', 'Anorthite', 'Enstatite', 'Forsterite',
    'Fayalite', 'Wollastonite', 'K-Feldspar', 'Diopside',
]

# Hydrothermal composition: basalt minerals present in hydrothermal.dat (bl-1kb.dat).
# Kaolinite is not in hydrothermal.dat, so it is omitted as a precipitating Al-sink.
hydrothermal_composition = {m: v for m, v in basalt_composition.items()
                             if m in hydrothermal_minerals}

carbonate_minerals = ['Calcite', 'Magnesite', 'Siderite', 'Dolomite']
secondary_sink_minerals = ['Kaolinite', 'Quartz', 'Goethite']
reverse_weathering_minerals = ['Smectite-high-Fe-Mg', 'Sepiolite']

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)