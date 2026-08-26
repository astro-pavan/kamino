"""make_database.py — assemble the low-temperature seafloor-weathering PHREEQC database.

Design (see project notes / paper methods):

  Base ("basic model"): the ThermoChimie SIT database (``sit.dat``). This supplies
    the SIT activity model — valid to the ionic strengths of the modelled oceans
    (I up to ~3-4 mol/kg) — together with the aqueous speciation (including the
    Mg/Ca carbonate and bicarbonate ion pairs and the high-pH silica species
    H3SiO4-/H2SiO4-2) and the authigenic precipitating assemblage that is already
    present in sit.dat (Calcite, Siderite, Kaolinite, Goethite, Greenalite, Halite,
    Smectite, and the Na-carbonate Na sink Nahcolite [+ Trona/Natron]).

  Added PHASES: the primary basalt silicates absent from sit.dat, taken from a
    SUPCRT-lineage source (llnl.dat, same lineage as the SUPCRTBL high-temperature
    database). Their dissolution is kinetically- and transport-limited and stays
    far from equilibrium, so SUPCRT-lineage log K used inside the SIT reference
    frame is acceptable for the dissolution step. A few gap-fill precipitating
    phases (e.g. SiO2(am), Sepiolite) are added the same way.

  RATES: mineral dissolution kinetics from Hermanska et al. (2022, 2023), via
    Kinec_v3.dat.

Citations for the paper:
  * Activity model / aqueous speciation: ThermoChimie v12a (Giffaut et al. 2014,
    Applied Geochemistry 49, 225; Grivé et al. 2015, Applied Geochemistry 55, 85);
    SIT model (Ciavatta 1980); PHREEQC v3 (Parkhurst & Appelo 2013).
  * Primary silicate thermodynamics: SUPCRTBL (Zimmer et al. 2016) / llnl.dat.
  * Kinetics: Hermanska et al. (2022, 2023).
"""

import os
import re


# ── Locate the databases shipped with the installed `phreeqc` package ──────────
def _phreeqc_database_dir():
    """Return the databases directory of the installed ``phreeqc`` package, or None.

    This resolves the path from the active Python environment rather than a
    hard-coded location, so it follows whichever venv the code runs in.
    """
    try:
        import phreeqc
    except ImportError:
        return None
    d = os.path.join(os.path.dirname(phreeqc.__file__), "databases")
    return d if os.path.isdir(d) else None


phreeqc_path = _phreeqc_database_dir()
SUPCRTBL_path = "https://github.com/HydrogeoIU/SupPHREEQC"

# Output directory (this file's directory: src/kamino/data/)
DATA_DIR = os.path.dirname(os.path.abspath(__file__))


# ── Source configuration ───────────────────────────────────────────────────────
# The "basic model": everything downstream is layered onto this complete database.
# Two activity models are supported. The choice is a real modelling decision, so it is a
# parameter rather than a hard-coded file -- see make_database(base=...).
#
#   "pitzer"  pitzer.dat. The Pitzer specific-ion-interaction model, valid to high ionic
#             strength and the model this project's results were produced on. DEFAULT.
#   "sit"     sit.dat (ThermoChimie v12a). Specific ion Interaction Theory, quoted to
#             I ~ 3-4 mol/kgw, which covers every ocean state the model reaches (I = 0.7-2.8).
#             Far fuller Al/Fe/Ca complexation -- 87 aqueous species involving the model's
#             elements against Pitzer's 28.
#
# MEASURED TRADE-OFF (identical PHREEQC solve, this machine):
#   speed          Pitzer 2.2 ms/call, SIT 20.2 ms/call -- 9.1x. The cost is not overhead: it is
#                  the larger complexation network for the model's OWN elements (87/28 = 3.1,
#                  squared = 9.6, matching the observed ratio for a Newton solve). It cannot be
#                  trimmed away -- see the note above _restrict_base's removal.
#   weathering     alkalinity flux agrees to within 10-20% (median SIT/Pitzer 1.09x) and the
#                  feedback strength d ln(F)/dT to within ~15%. The model sits in the
#                  kinetically-limited regime, where F -> A_r*k and the equilibrium
#                  concentrations the activity model sets divide out of the flux entirely.
#
# So Pitzer is the default: 9x faster for a ~10% difference in the quantity that matters. SIT
# remains available for a sensitivity test, and is the better choice if the model is ever pushed
# into the thermodynamically-limited (high Damkohler) regime, where b_eq stops dividing out.
BASE_DATABASES = {"pitzer": "pitzer.dat", "sit": "sit.dat"}
DEFAULT_BASE = "pitzer"
BASE_DATABASE = BASE_DATABASES[DEFAULT_BASE]   # retained for backwards compatibility

# ── Aluminium graft (pitzer base only) ─────────────────────────────────────────
# pitzer.dat defines 24 elements and ALUMINIUM IS NOT ONE OF THEM. That is not an oversight to
# work around quietly: it means the Pitzer framework has no Al interaction parameters, so grafted
# Al gets the long-range electrostatic term only. Al is trace here (b_eq[Al] ~ 2.5e-4 mM) but
# thermodynamically load-bearing -- albite's solubility product, the nepheline argument, and
# Kaolinite as the Al sink all run through it. Document this limitation wherever Al speciation is
# quoted from a Pitzer-based run.
#
# The graft reproduces exactly what the previous hand-built runtime database (hybrid_ocean.dat,
# generated by a make_hybrid.py that was never committed) contained: one master species and five
# hydrolysis species. Source is Kinec_v3.dat, which is bundled with the `phreeqc` package and is
# the same source ocean_chem.dat cites for its thermodynamics -- so the values are traceable and
# the build needs no file outside the installed package. The master-species line is byte-identical
# in Kinec_v3.dat, llnl.dat and the old hybrid_ocean.dat.
#
# Deliberately NOT grafted: Kinec_v3.dat's organic Al complexes (Al(CH3COO)2+ and similar).
# Acetate is not a tracked species, and importing them would add ligands the model cannot balance.
GRAFT_SOURCE = "Kinec_v3.dat"

# Master species and aqueous species to graft when the base lacks them, reproducing exactly what
# the previous hand-built hybrid_ocean.dat carried. Deliberately NOT grafted: Kinec_v3.dat's
# organic complexes (Al(CH3COO)2+ and similar) -- acetate is not a tracked species and importing
# them would add ligands the model cannot balance.
GRAFTS = {
    # pitzer.dat defines 24 elements and aluminium is not among them.
    "Al": dict(
        masters=["Al"],
        species=[
            "Al+3 = Al+3",
            "H2O + Al+3 = Al(OH)+2 + H+",
            "2 H2O + Al+3 = Al(OH)2+ + 2 H+",
            "3 H2O + Al+3 = Al(OH)3 + 3 H+",
            "4 H2O + Al+3 = Al(OH)4- + 4 H+",
        ]),
    # pitzer.dat defines Fe only as Fe+2. Goethite (FeOOH + 3 H+ = Fe+3 + 2 H2O) and the ferric
    # side of the iron system need the redox pair, so graft the valence master species and the
    # Fe+3 aqueous species. The Fe(+2)/Fe(+3) lines are byte-identical in Kinec_v3.dat and in the
    # old hybrid_ocean.dat.
    "Fe(+3)": dict(
        # Kinec_v3.dat writes the ferrous/ferric couple against O2
        # (H+ + Fe+2 + 0.25 O2 = Fe+3 + 0.5 H2O); the model's redox bookkeeping and the old
        # hybrid_ocean.dat both use the electron form, which phreeqc.dat and wateq4f.dat carry.
        source="phreeqc.dat",
        masters=["Fe(+2)", "Fe(+3)"],
        species=[
            "Fe+2 = Fe+3 + e-",
            "Fe+3 + H2O = FeOH+2 + H+",
            "Fe+3 + 2 H2O = Fe(OH)2+ + 2 H+",
            "Fe+3 + 3 H2O = Fe(OH)3 + 3 H+",
            "Fe+3 + 4 H2O = Fe(OH)4- + 4 H+",
            "Fe+3 + Cl- = FeCl+2",
        ]),
}

# Phases the model needs that pitzer.dat does not define. The sit.dat base already carries all of
# these, so they are added only for the pitzer base -- adding them unconditionally would silently
# override ThermoChimie's versions in the SIT build.
PITZER_EXTRA_PHASES = {
    "Siderite":    [("llnl.dat", "Siderite")],
    "Kaolinite":   [("llnl.dat", "Kaolinite")],
    "Goethite":    [("llnl.dat", "Goethite")],
    "Saponite-Na": [("llnl.dat", "Saponite-Na")],
    "Greenalite":  [("llnl.dat", "Greenalite")],
}

# Primary basalt silicates + gap-fill precipitating phases that are NOT already
# in sit.dat under the model's canonical name. Each entry lists candidate
# (database, source_name) pairs, tried in order; the first that resolves is used
# and renamed to the canonical (model) name. Names differ between databases, so
# several fallbacks are given.
ADDED_PHASES = {
    "Forsterite":   [("llnl.dat", "Forsterite")],
    "Fayalite":     [("llnl.dat", "Fayalite")],
    "Enstatite":    [("llnl.dat", "Enstatite")],
    "Diopside":     [("llnl.dat", "Diopside")],
    "Wollastonite": [("llnl.dat", "Wollastonite"),
                     ("phreeqc.dat", "Wollastonite"),
                     ("wateq4f.dat", "Wollastonite"),
                     ("PHREEQC_ThermoddemV1.10_15Dec2020.dat", "Wollastonite")],
    "Anorthite":    [("llnl.dat", "Anorthite")],
    # Melilite endmember produced by the CIPW desilication cascade when albite->nepheline
    # cannot clear the silica deficit (silica-poor, Ca-rich melts at high mantle Mg/Si).
    # NOTE: this phase was already present in hybrid_ocean.dat but was NOT recorded here,
    # so a rebuild would silently have dropped it. Thermodynamics only -- Kinec_v3 has no
    # RATES block for it, so mineral_rates.akermanite_k is a documented PROXY.
    "Akermanite":   [("llnl.dat", "Akermanite")],
    # Fe endmembers of the pyroxenes. Adding them lets the norm keep iron in pyroxene
    # instead of forcing all of it into fayalitic olivine.
    "Hedenbergite": [("llnl.dat", "Hedenbergite")],
    "Ferrosilite":  [("llnl.dat", "Ferrosilite")],
    "Albite":       [("llnl.dat", "Albite"), ("sit.dat", "Albite-low")],
    # Feldspathoid produced by the CIPW desilication cascade in silica-undersaturated (high
    # mantle Mg/Si) crusts. Its solubility product carries silica to the FIRST power, not the
    # third as albite's does, so unlike albite it is not shut down by a silica-flooded pore
    # fluid — see crust_composition.cipw_norm step 5b.
    "Nepheline":    [("llnl.dat", "Nepheline")],
    "K-Feldspar":   [("llnl.dat", "K-Feldspar"), ("llnl.dat", "K-feldspar")],
    "SiO2(am)":     [("llnl.dat", "SiO2(am)"),
                     ("PHREEQC_ThermoddemV1.10_15Dec2020.dat", "SiO2(am)"),
                     ("phreeqc.dat", "SiO2(a)")],
    "Sepiolite(d)": [("llnl.dat", "Sepiolite(d)"), ("llnl.dat", "Sepiolite")],
    # NB: the Na/Mg reverse-weathering smectite is Saponite-Na, which is already
    # present in the sit.dat base (SIT-consistent, trioctahedral Mg-smectite — the
    # basalt alteration clay), so it is NOT added here.
    #
    # NB: the Na sink. Na-aluminosilicate clays (Saponite-Na, analcime, Na-zeolites)
    # are all Na/Al ~ 1 (Na charge-balanced by Al-for-Si substitution), so their Na
    # uptake is capped by the ~1 nM dissolved Al and is negligible — verified in
    # PHREEQC (scratchpad/rw_na_uptake.py): no clay escapes the Al bottleneck. The
    # only Al-FREE Na sink is a Na-carbonate, whose saturation is set by Na + carbonate
    # (concentration-driven, not Al-limited). These are already in the sit.dat base
    # (SIT-consistent), so like Saponite-Na they are NOT added here — they are listed
    # for provenance and to mark them as the intended Na sink:
    #   Nahcolite  NaHCO3            — the appropriate phase for high-pCO2 / bicarbonate
    #                                   ("carbonated") oceans; forms at moderate pH and
    #                                   the high DIC of the CO2-rich worlds Kamino targets.
    #   Trona      Na3H(CO3)2:2H2O   — intermediate-alkalinity soda oceans (higher pH).
    #   Natron     Na2CO3:10H2O      — carbonate-dominated soda oceans (highest pH).
    # Nahcolite is the one that actually saturates first at high-CO2 conditions
    # (scratchpad/na_carbonate_sink.py); Trona/Natron cover the higher-pH tail if a
    # world is alkalinity- rather than CO2-dominated. All act as a saturation CEILING
    # (they clamp Na once Na*carbonate hits saturation), so on Earth-like low-DIC
    # oceans they stay inert and Na is only weakly regulated; on high-alkalinity
    # oceans they self-limit Na. To activate the sink, add these to the model's
    # ocean precipitating-mineral set (mineral_info.py) — the phases are DB-ready here.
}

# Na-carbonate Na sink (Al-free), already present in the sit.dat base — see the NB in
# ADDED_PHASES. Listed here so the intended sink assemblage is programmatically
# discoverable; Nahcolite is the primary (high-pCO2) phase, Trona/Natron the higher-pH
# soda tail. Not added to the DB (would duplicate the base); wire into the model's
# ocean precipitating set (mineral_info.py) to activate.
NA_CARBONATE_SINK = ["Nahcolite", "Trona", "Natron"]

# Cross-database species compatibility shim.
# The base (ThermoChimie) names dissolved silica H4(SiO4). Phases sourced from
# llnl / SUPCRTBL / Thermoddem write it as SiO2 or H4SiO4. Dissolved monomeric
# silica is the same species in either convention (the two hydration waters are
# formal, a_H2O ~ 1), so we alias the foreign names onto the base master species
# with log_k 0 rather than rewriting every foreign reaction. This lets the added
# phase reactions resolve while keeping the base's silica thermodynamics.
# Dissolved monomeric silica is spelled three ways across these databases: `H4(SiO4)` in sit.dat,
# `H4SiO4` in pitzer.dat (and the old hybrid_ocean.dat), and `SiO2` in llnl/Kinec/Thermoddem
# phases. They are the same species -- the two hydration waters are formal, a_H2O ~ 1 -- so the
# foreign spellings are aliased onto whichever one the chosen base uses as its master species,
# with log_k 0. Section 11 of the development history calls this the naming trap; getting the
# direction wrong yields "Reaction for species has not been defined".
_SILICA_MASTER = {"sit": "H4(SiO4)", "pitzer": "H4SiO4"}


def _compat_solution_species(base):
    master = _SILICA_MASTER[base]
    others = [sp for sp in ("H4(SiO4)", "H4SiO4") if sp != master]
    # PHREEQC defines a species by writing it on the RIGHT of a reaction whose left-hand side is
    # already-known species. So the base's master goes on the LEFT and each foreign spelling is
    # defined from it -- writing these the other way round yields
    # "Reaction for species has not been defined".
    lines = ["SOLUTION_SPECIES",
             "# --- make_database.py compatibility aliases (silica naming) ---",
             f"# Base master species is {master}; define the other spellings from it."]
    for other in others:
        lines += [f"    {master} = {other}", "        log_k 0.0"]
    lines += [f"    {master} = SiO2 + 2 H2O", "        log_k 0.0"]
    return "\n".join(lines) + "\n"

# Dissolution kinetics (RATES). Hermanska et al. (2022/2023) via Kinec_v3.dat.
RATES_MINERALS = {
    "Forsterite":   [("Kinec_v3.dat", "Forsterite")],
    "Fayalite":     [("Kinec_v3.dat", "Fayalite")],
    "Enstatite":    [("Kinec_v3.dat", "Enstatite")],
    "Diopside":     [("Kinec_v3.dat", "Diopside")],
    "Wollastonite": [("Kinec_v3.dat", "Wollastonite")],
    "Anorthite":    [("Kinec_v3.dat", "Anorthite")],
    "Albite":       [("Kinec_v3.dat", "Albite")],
    "Nepheline":    [("Kinec_v3.dat", "Nepheline")],
    "K-Feldspar":   [("Kinec_v3.dat", "K-Feldspar")],
}


# ── Low-level PHREEQC section / block parsers (proven in create_db.py) ──────────
def _read(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def _extract_section_text(text, section_kw):
    """Return the text of a PHREEQC section, from its keyword to the next keyword."""
    m = re.search(rf"^{section_kw}\s*$", text, re.MULTILINE)
    if not m:
        return ""
    rest = text[m.end():]
    nxt = re.search(r"^[A-Z][A-Z_]+\s*$", rest, re.MULTILINE)
    end = m.end() + nxt.start() if nxt else len(text)
    return text[m.start():end]


# Column-0 keywords that are phase DATA, not a new phase name (defensive: in most
# databases these are indented, but some place them flush-left).
_PHASE_DATA_KEYWORDS = {
    "log_k", "logk", "delta_h", "delta_g", "vm", "t_c", "p_c", "omega",
    "gflag", "add_logk", "add_log_k", "-analytic", "-vm", "-log_k", "-delta_h",
}


def _looks_like_phase_name(bare):
    """True if a line begins a new PHASES entry (its phase-name line).

    A name line is flush-left, non-blank, not a comment, carries no reaction
    ('=' -> reaction line), is not a '-qualifier' (-analytic, -Vm, ...), and its
    first token is not a phase-data keyword (log_k, delta_h, ...).
    """
    if not bare or bare[0] in (" ", "\t"):      # indented -> data line
        return False
    s = bare.strip()
    if not s or s.startswith("#"):               # blank / comment
        return False
    if "=" in s:                                  # reaction line
        return False
    if s.startswith("-"):                        # -analytic, -Vm, -log_k, ...
        return False
    if s.split()[0].lower() in _PHASE_DATA_KEYWORDS:
        return False
    return True


def _extract_named_block(section_text, name, is_rates):
    """Extract a single named block (mineral) from a PHASES or RATES section.

    PHASES: name line through everything up to (but not including) the next
            phase-name line (see ``_looks_like_phase_name``). Handles databases
            that put the reaction equation flush-left (e.g. sit.dat).
    RATES:  name line through the terminating ``-end``.
    """
    lines = section_text.splitlines()
    result, inside = [], False
    for line in lines:
        bare = line.rstrip()
        is_name = (
            bare == name
            or bare.startswith(name + " ")
            or bare.startswith(name + "\t")
            or bare.startswith(name + "#")
        )
        if is_name and not inside:
            inside = True
            result = [line]
        elif inside:
            if is_rates:
                result.append(line)
                if bare.strip() == "-end":
                    break
            else:
                if _looks_like_phase_name(bare):
                    break
                result.append(line)
    while result and not result[-1].strip():
        result.pop()
    return "\n".join(result) if result else None


def _rename_block(block, source_name, model_name):
    """Rename a block's leading identifier (and any ``name$ = "..."`` reference)
    from ``source_name`` to ``model_name``, so blocks pulled from databases that
    use a different label resolve to the canonical names in kamino.mineral_info.
    """
    if block is None or source_name == model_name:
        return block
    lines = block.splitlines()
    lines[0] = re.sub(rf"^{re.escape(source_name)}\b", model_name, lines[0], count=1)
    return "\n".join(ln.replace(f'"{source_name}"', f'"{model_name}"') for ln in lines)


# ── Public retrieval functions ─────────────────────────────────────────────────
def retrieve_mineral_thermodynamic_data(mineral_name, database_path):
    """Return the PHASES block (dissolution reaction + log_K / analytic expression)
    for ``mineral_name`` from the PHREEQC database at ``database_path``, or None if
    the phase is not defined there."""
    phases = _extract_section_text(_read(database_path), "PHASES")
    return _extract_named_block(phases, mineral_name, is_rates=False)


def retrieve_kinetic_data(mineral_name, database_path):
    """Return the RATES block (BASIC rate program, name line … ``-end``) for
    ``mineral_name`` from the PHREEQC database at ``database_path``, or None if no
    rate program is defined there."""
    rates = _extract_section_text(_read(database_path), "RATES")
    return _extract_named_block(rates, mineral_name, is_rates=True)


def _resolve(model_name, candidates, retriever):
    """Try each (database, source_name) candidate in order; return (block, note).

    ``block`` is renamed to ``model_name``; ``note`` records the resolved source
    or a warning if none matched.
    """
    for db_file, source_name in candidates:
        db_path = os.path.join(phreeqc_path, db_file) # type: ignore
        if not os.path.isfile(db_path):
            continue
        block = retriever(source_name, db_path)
        if block is not None:
            block = _rename_block(block, source_name, model_name)
            note = f"{model_name}: {db_file} ({source_name})"
            return block, note
    tried = ", ".join(f"{d}:{n}" for d, n in candidates)
    return None, f"# WARNING: {model_name} not found in any candidate [{tried}]"


# ── On trimming the base database (measured, rejected) ─────────────────────────
# sit.dat carries ~1373 aqueous species and ~1769 phases spanning actinides, ore metals and
# organics, against the ten elements kamino tracks, so restricting it looks like free speed.
# It is not. Measured on an identical PHREEQC solve:
#
#     pure water          Pitzer 0.05 ms   SIT  1.17 ms   (21.7x -- fixed overhead)
#     full kamino ocean   Pitzer 2.22 ms   SIT 20.22 ms   ( 9.1x)
#
# The fixed, database-size-dependent overhead is only ~1.1 ms of the 20.2 ms; the rest is
# chemistry-dependent, because PHREEQC speciates only the elements actually present. The cost is
# the complexation network for the model's OWN elements -- 87 species against the Pitzer
# database's 28 -- and 87/28 = 3.1, squared = 9.6, which matches the observed 8.8x for a Newton
# solve on a larger species system.
#
# So trimming other elements recovers ~5% and a filter that drops master species also breaks the
# database's dependency graph ("Elements in species have not been tabulated"). The 9x is the
# price of the fuller Al/Fe/Ca complexation SIT was chosen for. Do not retry this.
#
def _base_elements(base_text):
    """Elements defined in a database's SOLUTION_MASTER_SPECIES block."""
    els, in_block = set(), False
    for raw in base_text.split("\n"):
        st = raw.split("#")[0].strip()
        if st == "SOLUTION_MASTER_SPECIES":
            in_block = True
            continue
        if st in ("SOLUTION_SPECIES", "PHASES", "EXCHANGE_MASTER_SPECIES", "RATES", "END"):
            in_block = False
            continue
        if in_block and st:
            name = st.split()[0]
            els.add(name)                 # full name, e.g. "Fe(+3)"
            els.add(name.split("(")[0])   # and the bare element, e.g. "Fe"
    return els


# ── Extracting the aluminium graft ─────────────────────────────────────────────
def _extract_master_species(db_text, element):
    """The SOLUTION_MASTER_SPECIES line defining `element`, or None."""
    in_block = False
    for raw in db_text.split("\n"):
        line = raw.split("#")[0]
        st = line.strip()
        if st == "SOLUTION_MASTER_SPECIES":
            in_block = True
            continue
        if st in ("SOLUTION_SPECIES", "PHASES", "EXCHANGE_MASTER_SPECIES", "RATES", "END"):
            in_block = False
            continue
        if in_block and st and st.split()[0] == element:
            return raw.rstrip()
    return None


def _extract_species_block(db_text, reaction):
    """A SOLUTION_SPECIES entry (its reaction line plus indented parameters), matched on the
    reaction with whitespace normalised, or None."""
    want = " ".join(reaction.split())
    lines = db_text.split("\n")
    in_block = False
    for i, raw in enumerate(lines):
        st = raw.split("#")[0].strip()
        if st == "SOLUTION_SPECIES":
            in_block = True
            continue
        if st in ("PHASES", "EXCHANGE_MASTER_SPECIES", "RATES", "END",
                  "SOLUTION_MASTER_SPECIES"):
            in_block = False
            continue
        if in_block and st and " ".join(st.split()) == want:
            out = [raw.rstrip()]
            for nxt in lines[i + 1:]:
                if nxt.strip() and not nxt[:1].isspace():
                    break
                out.append(nxt.rstrip())
            return "\n".join(out).rstrip()
    return None


def _graft(element):
    """Master species + aqueous species for `element`, for bases that lack it.

    Returns (text, provenance_notes, warnings). See GRAFTS for what is taken and why.
    """
    spec = GRAFTS[element]
    source = spec.get("source", GRAFT_SOURCE)
    src = os.path.join(phreeqc_path, source)
    if not os.path.exists(src):
        return "", [], [f"# WARNING: {source} not found; {element} NOT grafted"]
    text = _read(src)
    notes, warns, masters, blocks = [], [], [], []

    for m in spec["masters"]:
        line = _extract_master_species(text, m)
        if line is None:
            warns.append(f"# WARNING: no {m} master species in {source}")
        else:
            masters.append(line)
            notes.append(f"MASTER  {m}: {source}")
    if not masters:
        return "", notes, warns

    for reaction in spec["species"]:
        block = _extract_species_block(text, reaction)
        if block is None:
            warns.append(f"# WARNING: {element} species not found in {source}: {reaction}")
        else:
            blocks.append(block)
            notes.append(f"SPECIES {reaction}: {source}")

    body = (
        f"# ═══ {element} graft ═══\n"
        f"# The chosen base does not define {element}. Master and aqueous species are taken from\n"
        f"# {source}, bundled with the `phreeqc` package and the source ocean_chem.dat cites\n"
        f"# for its thermodynamics, so the build needs no file outside the installed package.\n"
        "# NOTE for the Pitzer base: this supplies THERMODYNAMICS only. The Pitzer framework has\n"
        "# no interaction parameters for these ions, so they receive the long-range term alone.\n"
        "SOLUTION_MASTER_SPECIES\n"
        + "\n".join(masters) + "\n\n"
        "SOLUTION_SPECIES\n"
        + "\n\n".join(blocks) + "\n"
    )
    return body, notes, warns


# ── Assemble the database ──────────────────────────────────────────────────────
def make_database(name=None, database_path=None, base=DEFAULT_BASE):
    """Assemble the low-temperature weathering database and write it to disk.

    Imports the basic model (ThermoChimie SIT), appends the primary-silicate and
    gap-fill PHASES and the Hermanska RATES, and writes the result.

    Parameters
    ----------
    name : str
        Database name; used in the header and as the default output filename.
    database_path : str or None
        Output path. Defaults to ``<this dir>/<name>.dat``.

    Returns
    -------
    (output_path, provenance, warnings)
    """
    if phreeqc_path is None:
        raise FileNotFoundError(
            "Could not locate the installed `phreeqc` package's databases directory. "
            "Install `phreeqc` in the active environment (e.g. big-venv)."
        )
    if base not in BASE_DATABASES:
        raise ValueError(f"base must be one of {sorted(BASE_DATABASES)}, not {base!r}")
    if name is None:
        name = f"lt_weathering_{base}"
    if database_path is None:
        database_path = os.path.join(DATA_DIR, f"{name}.dat")

    base_path = os.path.join(phreeqc_path, BASE_DATABASES[base])
    base_text = _read(base_path).rstrip()

    # PHREEQC STOPS READING A DATABASE AT `END`. pitzer.dat closes its PITZER block with one, so
    # anything appended after it -- every added PHASE, every RATES block, the aluminium graft --
    # is silently invisible: the file loads with rc=0 and no error, and the first symptom is
    # "Phase not found in database" at run time. This is the same trap section 23.9 records for
    # hand-edited databases, and it is why the previous hybrid_ocean.dat interleaved its phases
    # into the base (PHASES at line 259, END at 1145) rather than appending them.
    #
    # Databases do not need a terminating END, so strip any standalone one. sit.dat has none,
    # which is why the SIT build worked while the Pitzer build did not.
    _stripped = [ln for ln in base_text.split("\n") if ln.strip().upper() != "END"]
    if len(_stripped) != len(base_text.split("\n")):
        base_text = "\n".join(_stripped).rstrip()

    provenance, warnings = [f"BASE    {BASE_DATABASES[base]} ({base} activity model)"], []

    # pitzer.dat has no aluminium; sit.dat does. Graft only where it is missing.
    graft_text = ""
    _have = _base_elements(base_text)
    for element in GRAFTS:
        if element in _have:
            continue
        body, notes, warns = _graft(element)
        graft_text += body + "\n"
        provenance.extend(notes)
        warnings.extend(warns)

    # Phases the pitzer base lacks. Adding these unconditionally would override ThermoChimie's
    # own versions in the SIT build, so they are base-conditional.
    extra_phases = PITZER_EXTRA_PHASES if base == "pitzer" else {}

    # Added primary-silicate + gap-fill PHASES
    phase_blocks = []
    for model_name, candidates in {**ADDED_PHASES, **extra_phases}.items():
        block, note = _resolve(model_name, candidates,
                               retrieve_mineral_thermodynamic_data)
        if block is not None:
            phase_blocks.append(block)
            provenance.append("PHASES  " + note)
        else:
            phase_blocks.append(note)
            warnings.append(note)

    # Dissolution kinetics (RATES)
    rate_blocks = []
    for model_name, candidates in RATES_MINERALS.items():
        block, note = _resolve(model_name, candidates, retrieve_kinetic_data)
        if block is not None:
            rate_blocks.append(block)
            provenance.append("RATES   " + note)
        else:
            rate_blocks.append(note)
            warnings.append(note)

    header = _build_header(name, provenance, warnings)

    with open(database_path, "w", encoding="utf-8") as f:
        f.write(header)
        f.write(base_text + "\n\n")
        if graft_text:
            f.write(graft_text + "\n")
        f.write("# ═══ Cross-database species compatibility (silica naming) ═══\n")
        f.write(_compat_solution_species(base) + "\n")
        f.write("# ═══ Added primary-silicate & gap-fill PHASES (SUPCRT-lineage) ═══\n")
        f.write("PHASES\n\n")
        f.write("\n\n".join(phase_blocks) + "\n\n")
        f.write("# ═══ Dissolution kinetics — Hermanska et al. (2022/2023) ═══\n")
        f.write("RATES\n\n")
        f.write("\n\n".join(rate_blocks) + "\n")

    print(f"Written: {database_path}")
    if warnings:
        print(f"  {len(warnings)} unresolved mineral(s):")
        for w in warnings:
            print("   " + w)
    else:
        print("  All minerals resolved.")
    return database_path, provenance, warnings


def _build_header(name, provenance, warnings):
    lines = [
        f"# {name}.dat",
        "# Low-temperature seafloor-weathering PHREEQC database.",
        "#",
        "# Base (activity model + aqueous speciation + authigenic precipitating",
        "#   assemblage): ThermoChimie SIT database, sit.dat (v12a).",
        "#   Giffaut et al. (2014); Grivé et al. (2015); SIT: Ciavatta (1980);",
        "#   PHREEQC: Parkhurst & Appelo (2013).",
        "# Added primary-silicate PHASES: SUPCRT-lineage (llnl.dat / SUPCRTBL,",
        "#   Zimmer et al. 2016) — dissolution step only (far from equilibrium).",
        "# RATES: Hermanska et al. (2022, 2023) via Kinec_v3.dat.",
        "#",
        f"# Generated by make_database.py from databases in: {phreeqc_path}",
        "#",
        "# Provenance (phase/rate -> source):",
    ]
    lines += [f"#   {p}" for p in provenance]
    if warnings:
        lines.append("#")
        lines.append("# UNRESOLVED (need a source; database incomplete until fixed):")
        lines += [f"#   {w.lstrip('# ')}" for w in warnings]
    return "\n".join(lines) + "\n\n"


if __name__ == "__main__":
    make_database()
