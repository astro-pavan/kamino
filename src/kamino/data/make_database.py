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
BASE_DATABASE = "sit.dat"

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
COMPAT_SOLUTION_SPECIES = """\
SOLUTION_SPECIES
# --- make_database.py compatibility aliases (silica naming) ---
# Foreign phases (llnl/SUPCRTBL/Thermoddem) reference SiO2 / H4SiO4; the SIT base
# master species is H4(SiO4). Treat them as identical dissolved silica.
    H4(SiO4) = H4SiO4
        log_k 0.0
    H4(SiO4) = SiO2 + 2 H2O
        log_k 0.0
"""

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


# ── Assemble the database ──────────────────────────────────────────────────────
def make_database(name="lt_weathering_sit", database_path=None):
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
    if database_path is None:
        database_path = os.path.join(DATA_DIR, f"{name}.dat")

    base_path = os.path.join(phreeqc_path, BASE_DATABASE)
    base_text = _read(base_path).rstrip()

    provenance, warnings = [], []

    # Added primary-silicate + gap-fill PHASES
    phase_blocks = []
    for model_name, candidates in ADDED_PHASES.items():
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
        f.write("# ═══ Cross-database species compatibility (silica naming) ═══\n")
        f.write(COMPAT_SOLUTION_SPECIES + "\n")
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
