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

hydrothermal_minerals = [
    'Albite', 'Anorthite', 'Enstatite', 'Forsterite',
    'Fayalite', 'Wollastonite', 'K-Feldspar', 'Diopside',
]

# Hydrothermal composition: basalt minerals present in hydrothermal.dat (bl-1kb.dat).
# Kaolinite is not in hydrothermal.dat, so it is omitted as a precipitating Al-sink.
hydrothermal_composition = {m: v for m, v in basalt_composition.items()
                             if m in hydrothermal_minerals}

carbonate_minerals = ['Calcite', 'Siderite']
secondary_sink_minerals = ['Kaolinite', 'SiO2(am)', 'Goethite']

# Authigenic clay minerals that form by reverse weathering (Isson & Planavsky 2018).
# Fe-rich clays dominate in anoxic, high-Si (Precambrian) oceans.
# Mg-rich clays dominate in oxygenated, low-Si (Phanerozoic) oceans.
reverse_weathering_minerals_fe = [
    'Greenalite',    # Fe₃Si₂O₅(OH)₄       Fe-serpentine;       Alk:Si = 3.00
    'Minnesotaite',  # Fe₃Si₄O₁₀(OH)₂      Fe-talc;             Alk:Si = 1.50
    'Chamosite-7A',  # Fe₂Al₂SiO₅(OH)₄     Berthierine;         Alk:Si = 4.00
    'Daphnite-14A',  # Fe₅Al₂Si₃O₁₀(OH)₈   Fe-chlorite;         Alk:Si = 3.33
    'Nontronite-Mg', # Mg₀.₁₆₅Fe₂Al₀.₃₃Si₃.₆₇·nH₂O  Fe³⁺-smectite
]
reverse_weathering_minerals_mg = [
    'Sepiolite',     # Mg₄Si₆O₁₅(OH)₂      Mg-chain silicate;   Alk:Si = 1.33
    'Clinochlore-14A', # Mg₅Al₂Si₃O₁₀(OH)₈ Mg-chlorite;         Alk:Si = 3.33
    'Saponite-Mg',   # Mg₃.₁₆₅Al₀.₃₃Si₃.₆₇O₁₀(OH)₂  Mg-smectite
]
reverse_weathering_minerals = reverse_weathering_minerals_fe + reverse_weathering_minerals_mg

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)