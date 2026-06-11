komatiite_42 = {
    'Forsterite': 0.7,
    'Fayalite': 0.3
} # 100% Olivine

komatiite_44 = {
    'Wollastonite': 0.333 * 0.3,
    'Enstatite': 0.666 * 0.3,
    'Forsterite': 0.7 * 0.7,
    'Fayalite': 0.3 * 0.7
} # 70% Olivine 30% Pyroxene

basalt_47 = {
    'Wollastonite': 0.333 * 0.7,
    'Enstatite': 0.666 * 0.7,
    'Forsterite': 0.7 * 0.3,
    'Fayalite': 0.3 * 0.3
} # 30% Olivine 70% Pyroxene

basalt_49 = {
    'Anorthite': 0.7 * 0.1,
    'Albite': 0.3 * 0.1,
    'Wollastonite': 0.333 * 0.8,
    'Enstatite': 0.666 * 0.8,
    'Forsterite': 0.7 * 0.1,
    'Fayalite': 0.3 * 0.1
} # 10% Olivine 80% Pyroxene 10% Plagioclase

basalt_51 = {
    'Anorthite': 0.7 * 0.3,
    'Albite': 0.1 * 0.3,
    'Wollastonite': 0.333 * 0.7,
    'Enstatite': 0.666 * 0.7,
} # 70% Pyroxene 30% Plagioclase

# HT hydrothermal composition: replaces Wollastonite+Enstatite with Diopside
# (CaMgSi2O6), which is the actual dominant pyroxene in MORB. Wollastonite
# (CaSiO3) is a metamorphic mineral absent from fresh oceanic crust and gives
# unrealistically high equilibrium Ca (~62 mM) at HT conditions.
basalt_49_HT = {
    'Anorthite': 0.7 * 0.1,
    'Albite':    0.3 * 0.1,
    'Diopside':  0.8,           # replaces Wollastonite (0.333*0.8) + Enstatite (0.666*0.8)
    'Forsterite': 0.7 * 0.1,
    'Fayalite':  0.3 * 0.1,
} # 10% Olivine 80% Diopside-pyroxene 10% Plagioclase (HT-appropriate)

hydrothermal_minerals = [
    'Albite', 'Anorthite', 'Enstatite', 'Forsterite',
    'Fayalite', 'Wollastonite', 'K-Feldspar', 'Diopside',
    'Quartz', 'Talc', 'Chrysotile',
] # All minerals in hydrothermal database

# Molar masses (g/mol) for primary crust minerals — used to convert
# water-to-rock mass ratios into per-mineral mole amounts for PHREEQC.
MINERAL_MOLAR_MASS = {
    'Forsterite':  140.69,   # Mg2SiO4
    'Fayalite':    203.77,   # Fe2SiO4
    'Enstatite':   100.39,   # MgSiO3
    'Wollastonite':116.16,   # CaSiO3
    'Diopside':    216.55,   # CaMgSi2O6
    'Anorthite':   278.21,   # CaAl2Si2O8
    'Albite':      262.22,   # NaAlSi3O8
    'K-Feldspar':  278.33,   # KAlSi3O8
}

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)

carbonate_minerals = ['Calcite', 'Siderite'] # Ca and Fe sinks
clay_minerals = ['Kaolinite', 'Goethite'] # Al, Fe sinks
silica_minerals = ['SiO2(am)'] # Si sink
reverse_weathering_minerals = ['Sepiolite(d)', 'Smectite-Na', 'Greenalite'] # Mg, Na, Fe sinks (via reverse weathering)
evaporite_minerals = ['Halite'] # Cl (and Na) sink; only active when land_fraction > 0