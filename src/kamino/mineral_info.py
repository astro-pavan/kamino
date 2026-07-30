komatiite_42 = {
    'Forsterite': 0.7,
    'Fayalite': 0.3
} # 100% Olivine

# Pyroxene is represented by Diopside (CaMgSi2O6), the real dominant clinopyroxene
# in MORB. It replaces the earlier Wollastonite(0.333·X)+Enstatite(0.666·X) proxy:
# Wollastonite (CaSiO3) is a metamorphic mineral absent from fresh oceanic crust and
# floods the equilibrium fluid with Ca (~20 mM LT / ~62 mM HT). See project memory
# "HT weathering fixes" and scratchpad/verify_beq.py.
komatiite_44 = {
    'Diopside': 0.3,
    'Forsterite': 0.7 * 0.7,
    'Fayalite': 0.3 * 0.7
} # 70% Olivine 30% Pyroxene

basalt_47 = {
    'Diopside': 0.7,
    'Forsterite': 0.7 * 0.3,
    'Fayalite': 0.3 * 0.3
} # 30% Olivine 70% Pyroxene

basalt_49 = {
    'Anorthite': 0.7 * 0.1,
    'Albite': 0.3 * 0.1,
    'Diopside': 0.8,
    'Forsterite': 0.7 * 0.1,
    'Fayalite': 0.3 * 0.1
} # 10% Olivine 80% Pyroxene 10% Plagioclase

basalt_51 = {
    'Anorthite': 0.7 * 0.3,
    'Albite': 0.1 * 0.3,
    'Diopside': 0.7,
} # 70% Pyroxene 30% Plagioclase

hydrothermal_minerals = [
    'Albite', 'Anorthite', 'Enstatite', 'Forsterite',
    'Fayalite', 'Wollastonite', 'K-Feldspar', 'Diopside',
    'Quartz', 'Talc', 'Chrysotile',
    'Clinochlore', 'Epidote', 'Clinozoisite',
] # All minerals in hydrothermal database

# Secondary phases allowed to PRECIPITATE in the high-temperature (greenschist)
# equilibrium. Quartz buffers Si; Clinochlore (Mg-chlorite) is the Mg sink;
# Epidote / Clinozoisite (Ca-Al) buffer Ca. These replace Forsterite as the
# unphysical Mg-sink proxy and add the Ca sinks the assemblage was missing.
# Sulfate phases (anhydrite) are deliberately excluded (sulfur stays decoupled).
ht_secondary_minerals = ['Quartz', 'Clinochlore', 'Epidote', 'Clinozoisite']

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
    'Quartz':       60.08,   # SiO2 (only emitted for silica-oversaturated rocks)
}

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)

carbonate_minerals = ['Calcite', 'Siderite'] # Ca and Fe sinks
clay_minerals = ['Kaolinite', 'Goethite'] # Al, Fe sinks
silica_minerals = ['SiO2(am)'] # Si sink
reverse_weathering_minerals = ['Sepiolite(d)', 'Saponite-Na', 'Greenalite'] # Mg, (Mg+Na), Fe sinks (via reverse weathering). Saponite-Na (trioctahedral Mg-smectite, basalt alteration clay) is the Na sink; least Al-limited Na-smectite.
evaporite_minerals = ['Halite'] # Cl (and Na) sink; only active when land_fraction > 0

# Primary (igneous) rock-forming minerals — the minerals that make up the crust
# compositions above. In the model these are a SOURCE ONLY: they dissolve, but they do
# not crystallise out of cold seawater. Several of them (Enstatite, Forsterite, Albite)
# are supersaturated at seawater pH, so if PHREEQC is allowed to equilibrate them freely
# it back-precipitates them and wrecks the cation budget (Mg -> ~0, Ca -> ~200 mM).
# chemistry.py therefore tags exactly these with PHREEQC's `dissolve_only` modifier.
#
# Secondary/authigenic phases (carbonates, clays, silica, evaporites, reverse-weathering
# clays) are deliberately NOT listed here: they must be free to both precipitate and
# dissolve. Quartz is likewise excluded — it is used as a precipitating Si buffer at HT.
primary_minerals = {
    'Anorthite', 'Albite', 'K-Feldspar',
    'Wollastonite', 'Enstatite', 'Diopside',
    'Forsterite', 'Fayalite',
}

# The CIPW norm and its oxide molar masses now live in crust_composition.py, which
# builds crust_composition dicts from bulk-rock analyses / planetary parameters.