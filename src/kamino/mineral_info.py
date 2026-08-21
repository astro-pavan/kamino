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
# Epidote/Clinozoisite (Ca-Al) REMOVED: they consumed ~all the Ca released by
# Diopside/Anorthite dissolution at HT (dCa/-dMg -> 0), starving calcite burial of its
# Ca feedstock and breaking the intended Mg->Ca exchange. Without them the HT releases
# Ca (flux ~0.02 -> ~3.7 Tmol/yr). Clinochlore keeps the Mg sink; Quartz buffers Si.
ht_secondary_minerals = ['Quartz', 'Clinochlore']

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
    'Nepheline':   142.05,   # NaAlSiO4 (feldspathoid; CIPW desilication product of Albite)
    'K-Feldspar':  278.33,   # KAlSi3O8
    'Quartz':       60.08,   # SiO2 (only emitted for silica-oversaturated rocks)
}

hydrothermal_mineral_string = ' '.join(hydrothermal_minerals)

carbonate_minerals = ['Calcite', 'Siderite', 'Nahcolite'] # Ca, Fe, Na sinks
clay_minerals = ['Kaolinite', 'Goethite'] # Al, Fe sinks
silica_minerals = ['SiO2(am)'] # Si sink
reverse_weathering_minerals = ['Sepiolite(d)', 'Saponite-Na', 'Greenalite'] # Mg, (Mg+Na), Fe sinks (via reverse weathering). Saponite-Na (trioctahedral Mg-smectite, basalt alteration clay) is the Na sink; least Al-limited Na-smectite.
evaporite_minerals = ['Halite'] # Cl (and Na) sink; only active when land_fraction > 0

# Secondary phases allowed to precipitate DURING the low-temperature primary equilibrium
# (get_b_eq) rather than only afterwards. The HT path has always done this (get_b_eq appends
# ht_secondary_minerals); the LT path never has. The wiring exists so the two paths CAN be made
# consistent, but the list is EMPTY -- an empty list is behaviourally identical to not passing
# it at all, so this is off, not merely defaulted.
#
# History (fast_13). ['Kaolinite'] was tried and reverted. Al is the trace species in Albite's
# solubility product (~1e-7 mol/kgw), so removing it un-saturates the feldspars, and sweep-wide
# it worked: median ocean Na went 6.6e-7 -> 0.911 mM, runs with Na > 1 mM went 1 -> 434, and
# clamp-pinned oceans halved (23% -> 9.5%).
#
# It was reverted because it does not help WHERE IT MATTERS. At the states near Da = 1 -- where
# the kinetic/thermodynamic transition and the pCO2 reversal have to happen -- it moves the
# delivered F[Na] essentially not at all (measured -0.013 to +0.009 Tmol/yr across w/r = 3-600,
# with and without it). The sweep-wide gain came from the high water/rock ratio it required,
# not from the buffer. And it cost convergence: it is the only configuration that fails
# outright at water_rock_ratio = None (see the guard in weathering.py, and fast_13's first run,
# which failed 100% of its chemistry and still reported 95% clean timeouts).
#
# If it is reinstated: ONLY the Al sink belongs here. Adding the Mg sinks
# (reverse_weathering_minerals) forces the equilibrium to strip Mg from the through-flowing
# seawater, and each mole of Mg2+ removed takes 2 eq of alkalinity with it -- measured
# F[Alk] +2.5 -> -207 Tmol/yr. Those are ocean-facing sinks and stay in the post-equilibrium
# precipitation step.
lt_equilibrium_buffer_minerals = []  # OFF -- see above

# Primary (igneous) rock-forming minerals - they generally only dissolve
primary_minerals = {
    'Anorthite', 'Albite', 'Nepheline', 'K-Feldspar',
    'Wollastonite', 'Enstatite', 'Diopside',
    'Forsterite', 'Fayalite',
}