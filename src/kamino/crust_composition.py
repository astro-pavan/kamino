"""Generate a seafloor-crust mineral composition from astronomical composition parameters.

Pipeline (two observationally-motivated knobs -> a `crust_composition` dict):

    mantle molar Mg/Si, core-formation dIW
      --feo_from_delta_iw()-----------> mantle FeO wt% (and core mass fraction)
      --oxide_composition()-----------> primary-melt oxide vector (wt%), from the MAGEMin table
      --mineral_composition()---------> CIPW norm -> DB-valid mineral fractions

The oxide table `data/crust_compositions.csv` is generated offline by
`data/make_crust_compositions.jl`: MAGEMin (Riel et al. 2022) with the Holland, Green & Powell
(2018) igneous dataset, isentropic decompression melting of a McDonough & Sun (1995) pyrolite
whose MgO/SiO2 and FeO have been re-set to the target axes, stopped at a fixed melt fraction
F = 0.20 (the value Guimond et al. 2024 adopt, being where clinopyroxene leaves the melting
assemblage for Earth's mantle). This module only interpolates that table and runs the norm.

The two axes
------------
`mantle_mg_si`  molar Mg/Si of the mantle. Earth = 1.25. Controls olivine vs orthopyroxene, and
                through the melt, the normative feldspar/feldspathoid balance. Mantles are
                olivine-free below ~0.8 and orthopyroxene-free above ~1.6.

`delta_iw`      CORE-FORMATION oxygen fugacity, log10 units relative to iron-wustite. Sets how
                much iron stayed in the mantle as FeO rather than going to the core as metal,
                via the metal-silicate equilibrium dIW = 2 log10(a_FeO_sil / a_Fe_met).
                Earth = -2. THE AXIS IS LOGARITHMIC IN FeO: -5 -> 0.26 wt%, -2 -> 8.05 wt%,
                -1 -> 24 wt%, which is already at the ~25 wt% ceiling beyond which the
                thermodynamic models are unreliable.

                This is NOT the fO2 of a modern melt -- see `melt_delta_iw` and the
                DELTA_IW_SELF_OXIDATION note in constants.py before handing a value to the
                outgassing model.

What is deliberately NOT an axis
--------------------------------
* T_p -- not observable, and not free: a mantle that cannot melt cannot transport heat
  magmatically. Solved per composition by the generator and read back from the table.
* C/O -- a regime switch, not a continuum. Carbides/nitrides need C/O > 0.9 (~1% of FGK
  dwarfs), and those planets are not silicate at all, so HGP18 cannot model them. Below 0.9
  the major elements do not care. Kamino's carbon input is the `outgassing` parameter.
* Ca/Al -- real, but both elements are refractory so stellar Ca/Al varies least of the
  candidate ratios. Held at pyrolite.
* Na2O -- the strongest UN-SWEPT lever, and it should be reported as such: alkalis move the
  solidus more than MgO/FeO does (Guimond et al. 2024 section 4.1), and albite vs nepheline is
  what sets this model's ocean Na. Excluded because Na is moderately volatile and
  devolatilisation is stochastic, so it does not belong on a grid indexed by observables.

Citations
---------
  * Guimond, Wang, Seidler, Sossi, Mahajan & Shorttle (2024), Rev. Mineral. Geochem. 90, 259.
  * Riel, Kaus, Green & Berlie (2022), G3 23, e2022GC010427 (MAGEMin).
  * Holland, Green & Powell (2018), J. Petrol. 59, 881 (thermodynamic dataset).
  * McDonough & Sun (1995), Chem. Geol. 120, 223 (pyrolite).
  * Frost & McCammon (2008), Annu. Rev. Earth Planet. Sci. 36, 389 (IW buffer).
"""

import functools
import os
import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from kamino.constants import DELTA_IW_SELF_OXIDATION, EARTH_DELTA_IW, EARTH_MANTLE_MG_SI
from kamino.mineral_info import MINERAL_MOLAR_MASS

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
CRUST_TABLE = os.path.join(_DATA_DIR, 'crust_compositions.csv')

# Major-element oxides carried through the pipeline (the ones the CIPW norm uses). Cr2O3 is in
# the generator's output but not here: the norm has no chromite and Cr is not a tracked ion.
PIPELINE_OXIDES = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O']

# Column in the CSV that supplies each pipeline oxide (the table reports total iron as FeO).
_TABLE_COLUMN = {ox: ('FeO' if ox == 'FeOt' else ox) for ox in PIPELINE_OXIDES}

# Oxide molar masses (g/mol) — convert a weight-percent oxide analysis to cation moles.
OXIDE_MOLAR_MASS = {
    'SiO2':   60.084,
    'TiO2':   79.866,
    'Al2O3': 101.961,
    'Fe2O3': 159.688,
    'FeO':    71.844,
    'FeOt':   71.844,   # total iron reported as FeO
    'MnO':    70.937,
    'MgO':    40.304,
    'CaO':    56.077,
    'Na2O':   61.979,
    'K2O':    94.196,
    'P2O5':  141.945,
    'Cr2O3': 151.990,
}

# Fe2O3 -> FeO mass-conversion factor (2 * M_FeO / M_Fe2O3): combine split iron into
# total FeO, since the CIPW norm routes all iron to fayalitic olivine.
_FE2O3_TO_FEO = 2 * OXIDE_MOLAR_MASS['FeO'] / OXIDE_MOLAR_MASS['Fe2O3']

# ---------------------------------------------------------------------------------------------
# Redox: dIW <-> mantle FeO, and dIW <-> absolute fO2
# ---------------------------------------------------------------------------------------------

# McDonough & Sun (1995) pyrolite, NON-Fe oxides only (wt%). Iron is set by the dIW axis and the
# rest of the budget is renormalised against it, so these are proportions, not absolutes. Must
# stay identical to PYROLITE_NON_FE in make_crust_compositions.jl.
PYROLITE_NON_FE = {
    'SiO2': 45.00, 'Al2O3': 4.45, 'CaO': 3.55, 'MgO': 37.80,
    'K2O': 0.029, 'Na2O': 0.36, 'TiO2': 0.201, 'Cr2O3': 0.384,
}
_CATIONS_PER_OXIDE = {'SiO2': 1, 'Al2O3': 2, 'CaO': 1, 'MgO': 1,
                      'K2O': 2, 'Na2O': 2, 'TiO2': 1, 'Cr2O3': 2}

# Cation moles contributed by 1 wt% of the (renormalised) non-Fe budget.
_K_CATIONS = sum(
    PYROLITE_NON_FE[ox] / sum(PYROLITE_NON_FE.values()) * _CATIONS_PER_OXIDE[ox] / OXIDE_MOLAR_MASS[ox]
    for ox in PYROLITE_NON_FE
)

# Earth's BSE FeO, the anchor for the whole redox axis (McDonough & Sun 1995).
EARTH_MANTLE_FEO = 8.05  # wt%

# Average mid-ocean ridge basalt, for reference ONLY -- nothing in the model consumes it. It is
# the observational anchor the figures compare the computed crusts against: real oceanic crust,
# put through the same norm.
#
# Gale, Dalton, Langmuir, Su & Schilling (2013), G3 14, 489 -- "The mean composition of ocean
# ridge basalts", the global "All MORB" average.
#
# !! THESE NUMBERS ARE UNVERIFIED AGAINST THE PAPER. It is paywalled and no accessible secondary
# source quotes the table, so they were not confirmed at the primary source. They are internally
# consistent with what this repository already asserts about MORB (CaO/Al2O3 = 0.775 against the
# ~0.78 quoted in docs/crust_composition.md section 7; Mg# 0.564, in the usual 0.55-0.60 range;
# oxides summing to 99.56) and the norm returns a textbook basalt (53 wt% plagioclase, 25 wt%
# clinopyroxene) -- but CHECK THEM AGAINST GALE ET AL. BEFORE PUBLISHING A FIGURE THAT USES THEM.
MORB_OXIDES = {
    'SiO2': 50.47, 'TiO2': 1.68, 'Al2O3': 14.70, 'FeOt': 10.43, 'MnO': 0.18,
    'MgO': 7.58, 'CaO': 11.39, 'Na2O': 2.79, 'K2O': 0.16, 'P2O5': 0.18,
}

# Metal-silicate activity constant, X_Fe(metal) / gamma_FeO, in
#     X_FeO(silicate, cation mole fraction) = _FEO_ACTIVITY_CONST * 10 ** (dIW / 2)
#
# CALIBRATED, not measured: fixed so that dIW = EARTH_DELTA_IW reproduces EARTH_MANTLE_FEO
# exactly. The first-principles value carries ~0.4 log units of slop -- gamma_FeO = 1.0 gives
# Earth dIW = -2.35, gamma_FeO = 1.5 gives -1.99 -- so the absolute number is a label and only
# the Earth anchor is meaningful. Same lesson as the T_p offset in section 24.2 of the
# development history: anchor on the composition, not on the number.
#
# The mapping is also OURS, not Putirka & Rarick's: they define alpha_Fe = Fe_BSP/Fe_BP on a
# cation weight basis with explicit core Ni and Si. Adopt their formulation before citing them.
def _x_feo_from_feo_wt(feo_wt: float) -> float:
    """Cation mole fraction of FeO in a pyrolite mantle carrying `feo_wt` wt% FeO."""
    n_fe = feo_wt / OXIDE_MOLAR_MASS['FeO']
    return n_fe / (n_fe + (100.0 - feo_wt) * _K_CATIONS)


_FEO_ACTIVITY_CONST = _x_feo_from_feo_wt(EARTH_MANTLE_FEO) / 10 ** (EARTH_DELTA_IW / 2)


def feo_from_delta_iw(delta_iw: float) -> float:
    """Mantle FeO (wt%) implied by a core-formation oxygen fugacity `delta_iw`.

    Inverts the metal-silicate equilibrium Fe + 1/2 O2 = FeO, for which
    dIW = 2 log10(a_FeO_silicate / a_Fe_metal), with the activity constant calibrated on Earth.

    Because the relation is logarithmic, a one-log-unit step in dIW is roughly a factor of
    three in FeO: -5 -> 0.26, -4 -> 0.82, -3 -> 2.59, -2 -> 8.05, -1 -> 24.1 wt%.
    """
    x = _FEO_ACTIVITY_CONST * 10 ** (delta_iw / 2)
    if not 0.0 < x < 1.0:
        raise ValueError(f'delta_iw={delta_iw} implies X_FeO={x:.3g}, outside (0, 1)')
    feo = 100.0 * x * _K_CATIONS / ((1.0 - x) / OXIDE_MOLAR_MASS['FeO'] + x * _K_CATIONS)
    if feo > 25.0:
        warnings.warn(
            f'delta_iw={delta_iw} gives mantle FeO = {feo:.1f} wt%; thermodynamic models are '
            f'unreliable above ~25 wt% (Guimond et al. 2024, section 3.1.1)', stacklevel=2)
    return feo


def delta_iw_from_feo(feo_wt: float) -> float:
    """Inverse of `feo_from_delta_iw`: core-formation dIW implied by a mantle FeO in wt%."""
    if not 0.0 < feo_wt < 100.0:
        raise ValueError('feo_wt must be in (0, 100)')
    return 2 * np.log10(_x_feo_from_feo_wt(feo_wt) / _FEO_ACTIVITY_CONST)


def mantle_composition(mantle_mg_si: float = EARTH_MANTLE_MG_SI,
                       delta_iw: float = EARTH_DELTA_IW) -> dict[str, float]:
    """Bulk mantle oxides (wt%) on the two axes -- the INPUT to the melting calculation.

    Mirrors `mantle_composition` in make_crust_compositions.jl and must stay identical to it: iron
    first from dIW, renormalising the non-Fe oxides to (100 - FeO) at their pyrolite proportions,
    then MgO/SiO2 re-split at fixed (MgO + SiO2) mass to hit the target molar ratio. Applied the
    other way round the two axes would not be orthogonal.

    Iron is reported as 'FeOt' to match `oxide_composition`, though the mantle is all-ferrous by
    construction (the MAGEMin `O` component is 0.0), so FeOt == FeO here.
    """
    if mantle_mg_si <= 0:
        raise ValueError('mantle_mg_si must be positive')
    feo = feo_from_delta_iw(delta_iw)
    scale = (100.0 - feo) / sum(PYROLITE_NON_FE.values())
    ox = {o: v * scale for o, v in PYROLITE_NON_FE.items()}
    total_mg_si = ox['MgO'] + ox['SiO2']
    si = total_mg_si / (1 + mantle_mg_si * OXIDE_MOLAR_MASS['MgO'] / OXIDE_MOLAR_MASS['SiO2'])
    ox['SiO2'], ox['MgO'] = si, total_mg_si - si
    ox['FeOt'] = feo
    return ox


# Iron-wustite buffer, log10 fO2_IW = A/T + B + C*(P-1)/T with T in K and P in bar
# (Frost & McCammon 2008 / O'Neill 1987). These coefficients are duplicated verbatim from
# m-class/outgassing/outgassing_model.fo2_IW so the two modules agree by construction --
# change them in both places or not at all.
_IW_A, _IW_B, _IW_C = -27215.0, 6.57, 0.055


def log10_fo2_iw(T: float, P_bar: float = 1.0) -> float:
    """log10 of the absolute oxygen fugacity (bar) of the iron-wustite buffer at T (K), P (bar)."""
    return _IW_A / T + _IW_B + _IW_C * (P_bar - 1.0) / T


def delta_iw_to_fo2(delta_iw: float, T: float, P_bar: float = 1.0) -> float:
    """Absolute oxygen fugacity (bar) for a buffer offset `delta_iw` at T (K), P (bar).

    The buffer position moves ~6 orders of magnitude between 1600 K / 1 bar and 1400 K / 10 GPa
    while the redox state of the rock is unchanged, which is the whole reason fO2 is quoted as an
    offset. A dIW is meaningless without the P and T it was evaluated at.
    """
    return 10.0 ** (delta_iw + log10_fo2_iw(T, P_bar))


def fo2_to_delta_iw(fo2_bar: float, T: float, P_bar: float = 1.0) -> float:
    """Buffer offset for an absolute oxygen fugacity `fo2_bar` at T (K), P (bar)."""
    return np.log10(fo2_bar) - log10_fo2_iw(T, P_bar)


def melt_delta_iw(delta_iw: float = EARTH_DELTA_IW,
                  self_oxidation: float = DELTA_IW_SELF_OXIDATION) -> float:
    """fO2 of a MODERN melt, given the CORE-FORMATION fO2 -- the handoff to the outgassing model.

    These are different quantities and conflating them is the easy mistake here. Core formation
    fixes how much iron the mantle keeps (Earth: IW-2, FeO 8.05 wt%). The mantle then
    self-oxidises by Fe disproportionation, raising Fe3+/SigmaFe at constant total Fe until the
    modern upper mantle sits near FMQ ~ IW+3.5. Volatile speciation in a degassing melt responds
    to the latter, mantle FeO to the former.

    `self_oxidation` is an Earth anchor (+5.5), not a law -- the true offset grows with mantle
    pressure. Pass your own if you have a better estimate. The return value is in the same units
    and sign convention as `outgassing_model.IW_offset`, so it can be handed across unchanged.
    """
    return delta_iw + self_oxidation


# ---------------------------------------------------------------------------------------------
# The melt-composition table
# ---------------------------------------------------------------------------------------------

# Module-level cache built by load_crust_table().
_INTERPOLATOR: RegularGridInterpolator | None = None
_TABLE_AXES: dict | None = None

# FeO spans nearly two orders of magnitude across the dIW axis, so it is interpolated in log
# space; everything else is near-linear in the two axes and is not.
_LOG_INTERPOLATED = {'FeOt'}


def load_crust_table(path: str = CRUST_TABLE) -> tuple[RegularGridInterpolator, dict]:
    """Read `crust_compositions.csv` and build the (mg_si, delta_iw) -> oxides interpolator.

    The generator writes a full regular grid, so this is a `RegularGridInterpolator` rather than
    a scattered one. Also returns the axes and the auxiliary columns (T_p, melt fraction,
    residual assemblage, CaO/Al2O3, warnings) for diagnostics.

    Cached module-level; call again to force a re-read after regenerating the table.
    """
    import pandas as pd  # local import: only needed to (re)build the interpolator

    df = pd.read_csv(path, comment='#')
    mg_si = np.array(sorted(df['mg_si'].unique()))
    d_iw = np.array(sorted(df['delta_iw'].unique()))
    if len(df) != len(mg_si) * len(d_iw):
        raise ValueError(
            f'{path} is not a complete regular grid: {len(df)} rows for a '
            f'{len(mg_si)} x {len(d_iw)} grid. Rerun make_crust_compositions.jl -- a partial '
            f'table means grid points failed, and interpolating across the hole would hide it.')

    df = df.sort_values(['mg_si', 'delta_iw'])
    shape = (len(mg_si), len(d_iw))
    values = np.stack([
        (np.log(df[_TABLE_COLUMN[ox]].to_numpy()) if ox in _LOG_INTERPOLATED
         else df[_TABLE_COLUMN[ox]].to_numpy()).reshape(shape)
        for ox in PIPELINE_OXIDES
    ], axis=-1)

    # bounds_error=True on purpose. Extrapolating this table is what section 23.2 of the
    # development history punished: silently clipped compositions look exactly like an
    # insensitive model.
    interpolator = RegularGridInterpolator((mg_si, d_iw), values, method='linear',
                                           bounds_error=True)
    axes = {
        'mg_si': mg_si, 'delta_iw': d_iw,
        # T_melt in isobaric tables, T_p in the older isentropic ones.
        'T_melt': (df['T_melt'] if 'T_melt' in df else df['T_p']).to_numpy().reshape(shape),
        'melt_fraction': df['melt_fraction'].to_numpy().reshape(shape),
        'mantle_feo': df['mantle_feo'].to_numpy().reshape(shape),
        'CaO_Al2O3': df['CaO_Al2O3'].to_numpy().reshape(shape),
        'residual_phases': df['residual_phases'].to_numpy().reshape(shape),
        'table': df,
    }

    global _INTERPOLATOR, _TABLE_AXES
    _INTERPOLATOR, _TABLE_AXES = interpolator, axes
    return interpolator, axes


def oxide_composition(mantle_mg_si: float = EARTH_MANTLE_MG_SI,
                      delta_iw: float = EARTH_DELTA_IW) -> dict[str, float]:
    """Primary-melt oxides (wt%) for a mantle of molar Mg/Si `mantle_mg_si` at `delta_iw`.

    Bilinear interpolation of the MAGEMin table. Iron is returned as total FeO ('FeOt'), which is
    what the norm consumes. Raises if either axis is off the table -- see load_crust_table.
    """
    if _INTERPOLATOR is None:
        load_crust_table()

    values = _INTERPOLATOR([[mantle_mg_si, delta_iw]])[0]  # type: ignore[misc]
    oxides = {}
    for ox, v in zip(PIPELINE_OXIDES, values):
        oxides[ox] = float(np.exp(v)) if ox in _LOG_INTERPOLATED else max(float(v), 0.0)

    total = sum(oxides.values())
    return {ox: v / total * 100.0 for ox, v in oxides.items()}


@functools.lru_cache(maxsize=512)
def _mineral_composition_cached(mantle_mg_si: float, delta_iw: float,
                                cipw_items: tuple) -> tuple:
    """Cached inner call. pyrolite's CIPW costs ~1.6 s against 0.03 ms for a hand-rolled norm, and
    a sweep holds the composition fixed while varying instellation/outgassing/crust production --
    so without this every run in the grid would pay the full cost again for an identical answer.

    Every input that changes the assemblage MUST arrive through `cipw_items` so it is part of the
    cache key. A module-level flag read inside `cipw_norm` would not be, and flipping it would
    silently return the previously-cached assemblage."""
    oxides = oxide_composition(mantle_mg_si, delta_iw)
    return tuple(cipw_norm(oxides, **dict(cipw_items)).items())


def mineral_composition(mantle_mg_si: float = EARTH_MANTLE_MG_SI,
                        delta_iw: float = EARTH_DELTA_IW, **cipw_kwargs) -> dict[str, float]:
    """Full pipeline: (mantle Mg/Si, core-formation dIW) -> weight-fraction crust mineral dict.

    Extra keyword arguments (emit_quartz, kfeldspar, verbose) pass through to cipw_norm. The
    returned dict is a valid `crust_composition` for get_weathering_flux / get_b_eq.

    `emit_quartz` defaults to True here, unlike in cipw_norm itself: below Mg/Si ~0.8 the melts
    are genuinely silica-oversaturated and discarding the excess is neither mass-conservative nor
    honest -- it relabels a rhyolite as a basalt. Quartz is in `primary_minerals`, so it dissolves
    only and cannot act as a low-temperature silica buffer.
    """
    cipw_kwargs.setdefault('emit_quartz', True)
    return dict(_mineral_composition_cached(float(mantle_mg_si), float(delta_iw),
                                            tuple(sorted(cipw_kwargs.items()))))



# ---------------------------------------------------------------------------------------------
# CIPW norm — pyrolite, with corrections for the phases the database cannot express
# ---------------------------------------------------------------------------------------------

# Oxides the database's phase set can express. Everything else (TiO2, K2O, MnO, P2O5, Cr2O3) is
# removed from the input and the remainder renormalised, BEFORE the norm runs. Removing them
# afterwards does not work: standard CIPW allocates Fe to magnetite and ilmenite and K to
# orthoclase/leucite ahead of the ferromagnesian phases, so deleting those products strands their
# cations and leaves the silica balance — which is what sets the olivine/pyroxene split — wrong.
_NORM_OXIDES = ['SiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O']

# pyrolite applies a ferric-iron correction by DEFAULT, and `Fe_correction=None` does not disable
# it (normative.py: `if Fe_correction is None: Fe_correction = "LeMaitre"`). The default assigns
# Fe2O3/FeO from a TAS rock-type classification (_MiddlemostTASRatios, ratios 0.1-0.5), which
# would invent ~3.6 wt% normative magnetite from an all-ferrous melt and silently contradict the
# dIW axis. The correction is skipped only where BOTH FeO and Fe2O3 are > 0, so Fe2O3 must be a
# tiny POSITIVE number rather than zero. This is undocumented and version-fragile: the guard in
# `cipw_norm` fails loudly if it ever stops working.
_FERRIC_SENTINEL = 1e-9

# pyrolite endmembers that map straight onto a database phase.
_PYROLITE_DIRECT = {
    'quartz': 'Quartz', 'albite': 'Albite', 'nepheline': 'Nepheline', 'anorthite': 'Anorthite',
    'clinoenstatite': 'Diopside',      # CaO MgO (SiO2)2 -- diopside proper, despite the name
    'enstatite': 'Enstatite',          # MgO SiO2
    'forsterite': 'Forsterite', 'fayalite': 'Fayalite',
    # Fe endmembers of the pyroxenes. Both are in the database (grafted from llnl.dat by
    # make_database.py) with a measured Fe-cpx proxy rate (mineral_rates.augite_k), so the norm's
    # own iron assignment is kept rather than exchanged into olivine.
    'clinoferrosilite': 'Hedenbergite',   # CaO FeO (SiO2)2
    'ferrosilite': 'Ferrosilite',         # FeO SiO2
}

# Endmembers with no database counterpart, converted by the reactions in `cipw_norm`.
_PYROLITE_CORRECTED = {
    'dicalcium silicate': 172.239,     # Ca2SiO4          -- larnite
}
# Aggregate rows pyrolite also returns; counting them would double-count their endmembers.
_PYROLITE_GROUPS = {'olivine', 'hypersthene', 'diopside'}

_M_SIO2 = 60.084
_CIPW_IMPL = None


def _pyrolite_cipw():
    """Import pyrolite's CIPW norm, neutralising its matplotlib side effect. Cached."""
    global _CIPW_IMPL
    if _CIPW_IMPL is None:
        import matplotlib
        from pyrolite.mineral.normative import CIPW_norm
        # pyrolite writes pyrolite.mplstyle into the matplotlib config dir on import and applies
        # it globally; its `legend.bbox_to_anchor` line is not a valid rcParam, so every
        # matplotlib import in every sweep worker prints a "Bad key" warning. This cost real time
        # twice (development history sections 5 and 23.8). Strip the line: pyrolite only writes
        # the file when it is ABSENT, so editing it is durable where deleting it is not.
        style = os.path.join(matplotlib.get_configdir(), 'stylelib', 'pyrolite.mplstyle')
        try:
            if os.path.exists(style):
                with open(style) as fh:
                    lines = fh.readlines()
                kept = [ln for ln in lines if 'bbox_to_anchor' not in ln]
                if len(kept) != len(lines):
                    with open(style, 'w') as fh:
                        fh.writelines(kept)
        except OSError:
            pass  # cosmetic only -- never let this break the norm
        _CIPW_IMPL = CIPW_norm
    return _CIPW_IMPL


def cipw_norm(oxides: dict[str, float], emit_quartz: bool = True, kfeldspar: bool = False,
              verbose: bool = False) -> dict[str, float]:
    """Convert a melt major-element analysis into a crust mineral composition.

    The norm itself is pyrolite's implementation (Williams et al. 2020) of the CIPW norm; this
    function restricts the input to the oxides the model's PHREEQC database can express, and then
    converts the one normative phase that has no counterpart in that database. Both conversions
    are balanced and use only database-available phases:

        larnite   + 2 diopside -> 2 akermanite + SiO2          (releases silica)
        nepheline + 2 SiO2     -> albite                       (reabsorbs it)

    Larnite is the textbook desilication product and is rejected as a cement clinker phase
    dissolving ~10^5x faster than diopside — see mineral_rates.akermanite_k. Routing it through
    diopside rather than through enstatite matters: none of the silica-deficient compositions
    carries free enstatite, and the only other balanced route produces wollastonite, which
    floods Ca.

    Normative hedenbergite and ferrosilite are emitted directly (see `_PYROLITE_DIRECT`). Earlier
    versions exchanged their iron into fayalite, because the database had no Fe-pyroxene; both
    endmembers are now present with a measured Fe-cpx proxy rate, so the norm's own iron
    assignment stands. That exchange overstated iron delivery — most severely below Mg/Si ~0.8,
    where the melts carry no normative olivine and every ferrous atom ended up in fayalite, the
    fastest-dissolving host in the database.

    `emit_quartz` defaults True: below Mg/Si ~0.8 the melts are genuinely silica-oversaturated and
    discarding the excess is neither mass-conservative nor honest. `kfeldspar` is accepted for
    signature compatibility and must be False -- K is not a tracked ion and K2O is removed from the
    input before the norm runs.

    Returns weight fractions summing to 1, valid as a `crust_composition`.
    """
    import pandas as pd

    if kfeldspar:
        raise ValueError('kfeldspar=True is not supported: K2O is stripped before the norm runs')

    g = lambda k: float(oxides.get(k, 0.0))
    feo = g('FeO') + g('FeOt') + _FE2O3_TO_FEO * g('Fe2O3')
    kept = {'SiO2': g('SiO2'), 'Al2O3': g('Al2O3'), 'FeO': feo,
            'MgO': g('MgO'), 'CaO': g('CaO'), 'Na2O': g('Na2O')}
    total = sum(kept.values())
    if total <= 0:
        raise ValueError('no expressible oxides in the input; check the analysis')
    row = {k: v / total * 100.0 for k, v in kept.items()}
    row.update(TiO2=0.0, Fe2O3=_FERRIC_SENTINEL, MnO=0.0, K2O=0.0, P2O5=0.0, Cr2O3=0.0)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = _pyrolite_cipw()(pd.DataFrame([row]).astype('float64'), adjust_all_Fe=False)
    wt = {c: float(out[c].iloc[0]) for c in out.columns}

    # Guard the two assumptions this path rests on: that stripping the oxides suppressed the
    # phases we cannot express, and that the ferric sentinel suppressed the Fe correction.
    unexpected = {c: v for c, v in wt.items()
                  if v > 1e-3 and c not in _PYROLITE_DIRECT and c not in _PYROLITE_CORRECTED
                  and c not in _PYROLITE_GROUPS}
    if unexpected:
        raise RuntimeError(
            f'pyrolite returned phases this pipeline cannot express: '
            f'{ {k: round(v, 3) for k, v in unexpected.items()} }. If magnetite is among them the '
            f'ferric sentinel has stopped suppressing the LeMaitre/TAS Fe correction -- check '
            f'pyrolite\'s normative.py before trusting any crust composition.')

    mol = {name: wt[c] / MINERAL_MOLAR_MASS[name] for c, name in _PYROLITE_DIRECT.items()
           if wt.get(c, 0.0) > 0}
    lar = wt.get('dicalcium silicate', 0.0) / _PYROLITE_CORRECTED['dicalcium silicate']
    get = lambda m: mol.get(m, 0.0)
    warn_msgs: list[str] = []
    free_silica = 0.0

    # 1. Larnite -> akermanite, drawing Mg and Si from diopside and releasing silica.
    if lar > 1e-12:
        if 2 * lar > get('Diopside') + 1e-12:
            raise ValueError(f'{2*lar:.4g} mol diopside needed to convert normative larnite but '
                             f'only {get("Diopside"):.4g} available')
        mol['Diopside'] = get('Diopside') - 2 * lar
        mol['Akermanite'] = get('Akermanite') + 2 * lar
        free_silica += lar

    # 2. Reabsorb that silica by reversing desilication, then keep any remainder as quartz.
    if free_silica > 1e-12 and get('Nepheline') > 1e-12:
        conv = min(free_silica / 2.0, get('Nepheline'))
        mol['Nepheline'] = get('Nepheline') - conv
        mol['Albite'] = get('Albite') + conv
        free_silica -= 2 * conv
    if free_silica > 1e-12:
        if emit_quartz:
            mol['Quartz'] = get('Quartz') + free_silica
        else:
            warn_msgs.append(f'{free_silica:.3g} mol SiO2 released by the larnite conversion was '
                             f'dropped (set emit_quartz=True to keep it)')

    if verbose:
        print(f'  pyrolite: ' + ', '.join(f'{c} {v:.2f}' for c, v in sorted(wt.items())
                                          if v > 1e-3 and c not in _PYROLITE_GROUPS))
        print(f'  corrections: larnite {lar:.4g} mol')
    for m in warn_msgs:
        warnings.warn(m, stacklevel=2)

    weights = {m: n * MINERAL_MOLAR_MASS[m] for m, n in mol.items() if n > 1e-12}
    tot = sum(weights.values())
    if tot <= 0:
        raise ValueError('CIPW norm produced no minerals; check the oxide input')
    return {m: w / tot for m, w in weights.items()}


def _cipw_norm_native(oxides: dict[str, float], emit_quartz: bool = False, kfeldspar: bool = False,
              verbose: bool = False) -> dict[str, float]:
    """Convert a bulk-rock major-element analysis into a crust mineral composition.

    Takes major-element oxides in weight-percent and returns a weight-fraction
    dictionary of igneous ("normative") minerals, restricted to the phases that
    exist in the model's PHREEQC database (see hydrothermal_minerals /
    MINERAL_MOLAR_MASS):

        Forsterite, Fayalite, Enstatite, Diopside, Anorthite, Albite
        (+ Wollastonite only as a Ca-overflow fallback, + Quartz if emit_quartz,
         + K-Feldspar only if kfeldspar=True)

    The result is normalised to sum to 1 and can be passed straight in as a
    `crust_composition` (like basalt_49 etc.). Values are WEIGHT fractions, which
    is how get_b_eq consumes a composition (grams of mineral per kgw / molar mass
    -> moles; see chemistry.py `water_rock_ratio` handling).

    This is a compact CIPW-style norm (Cross-Iddings-Pirsson-Washington) that
    allocates cations to phases in the standard order feldspar -> diopside ->
    pyroxene/olivine, using the silica balance to split the ferromagnesian
    remainder between orthopyroxene (Enstatite) and olivine (Forsterite).

    Because the database has no magnetite/ilmenite/apatite/nepheline and no
    Fe-bearing pyroxene, it makes three deliberate simplifications:
      * TiO2 and P2O5 are dropped (P removes a little Ca as apatite first);
      * ALL iron is routed to fayalitic olivine (the only Fe silicate available),
        so pass total iron as FeO/FeOt (FeO and Fe2O3 are summed);
      * Diopside and Enstatite are pure Mg endmembers.
    These match how the existing hardcoded compositions are built. K is likewise not
    a tracked ion (no sink for it), so K-Feldspar is off by default (see kfeldspar).

    Parameters
    ----------
    oxides : dict[str, float]
        Major-element oxides in wt% (missing oxides default to 0). Recognised
        keys: SiO2, TiO2, Al2O3, Fe2O3, FeO, FeOt, MnO, MgO, CaO, Na2O, K2O, P2O5.
    emit_quartz : bool
        If the rock is silica-oversaturated, keep the excess SiO2 as 'Quartz'.
        Default False -> excess silica is dropped with a warning (Quartz doubles
        as an HT precipitation buffer, so it is normally kept out of the crust).
    kfeldspar : bool
        Allocate K to K-Feldspar. Default False -> K is dropped (no K sink in the
        model, and the LT database lacks a K-Feldspar phase). Only set True if a
        K-Feldspar phase has been added to hybrid_ocean.dat.
    verbose : bool
        Print the cation budget and any saturation warnings.

    Returns
    -------
    dict[str, float]  weight-fraction crust composition (sums to 1).
    """
    g = lambda k: float(oxides.get(k, 0.0))

    # oxide wt% -> cation moles (per 100 g; only ratios matter, we normalise later)
    m_Al2O3 = g('Al2O3') / OXIDE_MOLAR_MASS['Al2O3']
    m_Fe2O3 = g('Fe2O3') / OXIDE_MOLAR_MASS['Fe2O3']
    m_FeO   = (g('FeO') + g('FeOt')) / OXIDE_MOLAR_MASS['FeO']
    m_MnO   = g('MnO')   / OXIDE_MOLAR_MASS['MnO']
    m_P2O5  = g('P2O5')  / OXIDE_MOLAR_MASS['P2O5']

    Si = g('SiO2') / OXIDE_MOLAR_MASS['SiO2']
    Al = 2 * m_Al2O3
    Fe = m_FeO + m_MnO + 2 * m_Fe2O3        # all iron (+Mn) as divalent -> fayalite
    Mg = g('MgO') / OXIDE_MOLAR_MASS['MgO']
    Ca = g('CaO') / OXIDE_MOLAR_MASS['CaO'] - (10.0 / 3.0) * m_P2O5   # apatite removes Ca
    Na = 2 * (g('Na2O') / OXIDE_MOLAR_MASS['Na2O'])
    K  = 2 * (g('K2O')  / OXIDE_MOLAR_MASS['K2O'])
    Ca = max(Ca, 0.0)

    warn_msgs: list[str] = []
    moles: dict[str, float] = {}

    def take(mineral: str, amount: float) -> None:
        if amount > 1e-12:
            moles[mineral] = moles.get(mineral, 0.0) + amount

    # 1. K-Feldspar  KAlSi3O8 -- OFF by default. The LT runtime database
    #    (hybrid_ocean.dat) has no K-Feldspar phase, and the model deliberately does
    #    not carry K as a tracked ion (no sink for it), so K is dropped and its Al/Si
    #    pass to plagioclase. MORB is K-poor (~0.1% K2O -> ~1% orthoclase) so this is
    #    a small approximation; warn if handed a genuinely K-rich rock.
    if kfeldspar:
        kf = min(K, Al)
        take('K-Feldspar', kf); Al -= kf; Si -= 3 * kf
    elif g('K2O') > 0.5:
        warn_msgs.append(f'K2O={g("K2O"):.2g} wt% dropped (K not modelled); '
                         f'pass kfeldspar=True if the LT database has a K-Feldspar phase')
    # 2. Albite  NaAlSi3O8
    ab = min(Na, Al)
    take('Albite', ab); Al -= ab; Si -= 3 * ab
    # 3. Anorthite  CaAl2Si2O8  (uses remaining Al)
    an = min(Al / 2.0, Ca)
    take('Anorthite', an); Ca -= an; Al -= 2 * an; Si -= 2 * an
    if Al > 1e-9:
        warn_msgs.append(f'peraluminous: {Al:.3g} mol excess Al (corundum-normative) dropped')
    # 4. Diopside  CaMgSi2O6  (remaining Ca paired with Mg)
    di = min(Ca, Mg)
    take('Diopside', di); Ca -= di; Mg -= di; Si -= 2 * di
    # 4b. Ca overflow (Mg exhausted before Ca) -> Wollastonite fallback
    if Ca > 1e-9:
        take('Wollastonite', Ca); Si -= Ca
        warn_msgs.append(f'Ca-oversaturated: {Ca:.3g} mol Ca to Wollastonite (Mg too low for diopside)')
        Ca = 0.0
    # 5. Ferromagnesian remainder: all Fe -> Fayalite; Mg -> Enstatite/Forsterite by silica
    take('Fayalite', Fe / 2.0); Si -= 0.5 * Fe          # Fe2SiO4 : 0.5 Si per Fe
    if Mg > 1e-12:
        if Si >= Mg:                    # silica-saturated: all orthopyroxene
            take('Enstatite', Mg); Si -= Mg
        elif Si <= 0.5 * Mg:            # silica-deficient: all olivine, then desilicate (step 5b)
            take('Forsterite', Mg / 2.0); Si -= 0.5 * Mg
        else:                          # mixed opx + olivine set by silica balance
            en = 2 * Si - Mg
            take('Enstatite', en)
            take('Forsterite', (Mg - en) / 2.0)
            Si -= en + 0.5 * (Mg - en)

    # 5b. Silica deficit -> desilicate Albite to Nepheline, the standard CIPW cascade.
    # Converting olivine-normative Mg is already done above, so feldspar is the next donor:
    # NaAlSi3O8 -> NaAlSiO4 + 2 SiO2 releases 2 mol Si per mol converted. Without this the
    # deficit was silently discarded (step 6 only fires on Si > 0), which both violated mass
    # balance -- the norm assigned more SiO2 to minerals than the rock contained -- and pinned
    # every undersaturated composition onto one identical assemblage, erasing the Mg/Si signal.
    if Si < -1e-12 and moles.get('Albite', 0.0) > 1e-12:
        needed = -Si / 2.0
        converted = min(needed, moles['Albite'])
        moles['Albite'] -= converted
        if moles['Albite'] <= 1e-12:
            moles.pop('Albite')
        take('Nepheline', converted)
        Si += 2 * converted
    # 5c. Still deficient -> desilicate Diopside to Akermanite, the melilite step.
    # 2 CaMgSi2O6 -> Ca2MgSi2O7 + 1/2 Mg2SiO4 + 3/2 SiO2 releases 0.75 mol Si per mol diopside
    # converted. Textbook CIPW would use larnite (Ca2SiO4) here, and Kinec_v3 has larnite
    # complete (thermodynamics AND kinetics) -- it is rejected deliberately. Larnite dissolves
    # 606x faster than Wollastonite and ~1e5x faster than Diopside, so at the 5-8 wt% this
    # cascade produces it would supply the entire dissolution flux; and it is a cement clinker
    # phase, rare in nature, whereas these melts are melilititic (section 24.4 via Medard et al.
    # 2004) and melilitites crystallise MELILITE. Akermanite is the phase that is actually there.
    # Its rate is a documented proxy -- see mineral_rates.akermanite_k.
    if Si < -1e-12 and moles.get('Diopside', 0.0) > 1e-12:
        converted = min(-Si / 0.75, moles['Diopside'])
        moles['Diopside'] -= converted
        if moles['Diopside'] <= 1e-12:
            moles.pop('Diopside')
        take('Akermanite', converted / 2.0)
        take('Forsterite', converted / 4.0)
        Si += 0.75 * converted
    if Si < -1e-12:
        # Nothing left to desilicate. Proceeding would assign more SiO2 to minerals than the rock
        # contains -- a crust that violates mass balance by a few percent of its silica, which is
        # exactly the silent failure section 23.2 was written about. Raise rather than warn: a
        # composition this far outside the phase set is not a crust, and the caller must know.
        raise ValueError(
            f'silica-deficient: {-Si:.3g} mol SiO2 deficit remains after the albite->nepheline '
            f'and diopside->akermanite cascades. This composition is beyond what the database\'s '
            f'phase set can express; it is not a mass-conservative crust.')
    # 6. Excess silica -> Quartz (optional)
    if Si > 1e-9:
        if emit_quartz:
            take('Quartz', Si)
        else:
            warn_msgs.append(f'silica-oversaturated: {Si:.3g} mol SiO2 excess dropped '
                             f'(set emit_quartz=True to keep as Quartz)')

    if verbose:
        print(f'cations  Si={g("SiO2")/OXIDE_MOLAR_MASS["SiO2"]:.4f}  Al={2*m_Al2O3:.4f}  '
              f'Fe={Fe:.4f}  Mg={g("MgO")/OXIDE_MOLAR_MASS["MgO"]:.4f}  '
              f'Ca={g("CaO")/OXIDE_MOLAR_MASS["CaO"]:.4f}  Na={Na:.4f}  K={K:.4f}')
        for msg in warn_msgs:
            print('  warning:', msg)
    for msg in warn_msgs:
        warnings.warn(msg, stacklevel=2)

    # mineral moles -> weight fractions
    weights = {m: n * MINERAL_MOLAR_MASS[m] for m, n in moles.items() if n > 0}
    total = sum(weights.values())
    if total <= 0:
        raise ValueError('CIPW norm produced no minerals; check the oxide input')
    return {m: w / total for m, w in weights.items()}
