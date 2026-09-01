import os
import sys
import glob
import json
import argparse
import functools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import cmasher as cmr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
from kamino.weathering import get_weathering_flux
from kamino.constants import (
    EARTH_HYDROTHERMAL_FLUX_PER_AREA as J_ref_normalised,
    EARTH_CRUST_PRODUCTION_RATE_PER_AREA as rate_ref,
    A_SEAFLOOR_EARTH as A_seafloor,
)
from kamino.constants import (G, EARTH_CRUST_PRODUCTION_RATE_PER_AREA, YR, SOLAR_CONSTANT,
                              STEFAN_BOLTZMANN, EARTH_MANTLE_MG_SI, EARTH_DELTA_IW)
from kamino.chemistry import alk_idx, elements
from kamino.mineral_info import (carbonate_minerals, clay_minerals, silica_minerals,
                                 reverse_weathering_minerals, evaporite_minerals)
from kamino.crust_composition import mineral_composition
from kamino.planet import KD_MG_HT, K_NA_CONT_REMOVAL
from kamino.weathering import ALPHA_REF
from kamino.planet import _S_TERR_EARTH

# Figure style. Presentation mode uses larger type and wider figures; paper mode targets the
# MNRAS column width. Set KAMINO_PRESENTATION=1 to switch without editing the file.
presentation = os.environ.get('KAMINO_PRESENTATION', '0').lower() in ('1', 'true', 'yes')

_STYLE_DIR = os.path.dirname(os.path.abspath(__file__))
STYLE_FILE = os.path.join(
    _STYLE_DIR,
    'planetary-chem-presentation.mplstyle' if presentation else 'planetary-chem-paper.mplstyle')
if not os.path.exists(STYLE_FILE):
    raise SystemExit(f"missing style file {STYLE_FILE}")
# Resolved against THIS file, not the working directory. The previous relative path
# ('experiments/planetary-chem-*.mplstyle') only resolved when the script happened to be run
# from the repository root, and raised OSError from anywhere else.
plt.style.use(STYLE_FILE)


# ---------------------------------------------------------------------------
# Figure geometry
# ---------------------------------------------------------------------------
# Publication figures are sized to the PAGE, not to whatever happened to look right on screen.
# MNRAS is A4 two-column: \columnwidth = 240 pt and \textwidth = 504 pt, at 72.27 pt/inch.
# Ask for one of those two widths by name and give the height you want in inches; nothing else
# should set a figure size directly.
#
#     fig, axes = plt.subplots(4, 1, figsize=figure_size('single', 6.0))
#
# Presentation mode scales BOTH dimensions so the aspect ratio (and therefore the layout) is the
# one the paper will use; the larger type in the presentation style then sits correctly on it.
COLUMN_WIDTH_IN = 240 / 72.27      # 3.32 in -- one MNRAS column
TEXT_WIDTH_IN = 504 / 72.27        # 6.97 in -- both columns
PAGE_WIDTHS = {'single': COLUMN_WIDTH_IN, 'double': TEXT_WIDTH_IN}

# Default height per stacked panel row, in inches, when a figure does not name its own height.
ROW_HEIGHT_IN = 1.5

_PRES_SCALE = 2.0 if presentation else 1.0


def figure_size(width='single', height=None, n_rows=1, row_height=ROW_HEIGHT_IN):
    """(width, height) in inches for a page-ready figure.

    `width` is 'single' or 'double' -- one MNRAS column or the full text width. `height` is in
    inches; omit it to get `n_rows * row_height`, which is what a stack of panels wants.
    """
    if width not in PAGE_WIDTHS:
        raise ValueError(f"width must be one of {sorted(PAGE_WIDTHS)}, not {width!r}")
    h = height if height is not None else n_rows * row_height
    return (PAGE_WIDTHS[width] * _PRES_SCALE, h * _PRES_SCALE)


def diagnostic_size(n_rows, n_cols, col_width=6.0, row_height=3.0, pad=1.5):
    """Figure size for a DIAGNOSTIC grid, which is not page-constrained.

    The wide multi-panel diagnostics (the full outgassing x crust grid, the per-crust-rate
    mineral SI panels) exist to be read on screen at whatever size they need. Forcing them into
    a column width would make them illegible, so they keep a fixed size per panel.
    """
    return (col_width * n_cols + pad, row_height * n_rows)


# fast_18's /data/pt426 path was on a machine that is no longer available.
DEFAULT_OUTPUT_PATH = os.environ.get('KAMINO_SWEEP_OUTPUT', '/home/pavan/PhD/sweep_output')

TERM_LABELS = {
    'converged':     'Converged',
    'timeout':       'Timeout (2 Gyr)',
    'wall_timeout':  'Wall-clock cap',
    'out_of_domain': 'Outside model domain',
}

# Where the run left the validity box. Only meaningful for 'out_of_domain'; this is
# recovered from the final state in Planet.time_evolve, not from a dedicated event.
WALL_LABELS = {
    'cold':     'Frozen (T → 181 K)',
    'hot':      'Runaway (T → 389 K)',
    'co2_high': 'CO₂ ceiling (10 bar)',
    'co2_low':  'CO₂ depleted (0.1 Pa)',
}

# Terminations that mean "the model ran out of validity", not "the planet did something".
# Not a habitability verdict: a run cut off at a wall has no known fate.
# Terminations that mean the model ran out of validity, and those that do not. Older sweeps used
# a larger vocabulary; plot_legacy.upgrade() maps it onto these before any figure sees it.
OUT_OF_DOMAIN = {'out_of_domain'}
# NOTE 'wall_timeout' is deliberately in NEITHER set. A run cut off by the wall-clock cap has an
# incomplete integration, so it is neither known-habitable nor known-out-of-domain; it falls
# through to the 'x' marker. 115 of the 2045 runs in the 2026-08 sweep are in this state.
HABITABLE = {'converged', 'timeout'}
T_SNOWBALL = 260.0
T_RUNAWAY  = 360.0

# Molar masses (g/mol) for the b_ocean elements, used to turn the final state into a
# salinity. C is carried as HCO₃⁻ (61) and S as SO₄²⁻ (96.06); Alkalinity is a charge
# balance rather than a mass, so it is skipped. Indices are derived from
# kamino.chemistry.elements so the mapping follows the model's element list:
# y = [P_CO2, P_H2O, *elements, r_avg], i.e. elements[i] lives at y[2 + i].
_ELEMENT_MASSES = {'C': 61.0, 'Si': 60.1, 'Al': 27.0, 'Fe': 55.8, 'Ca': 40.1,
                   'Mg': 24.3, 'Na': 23.0, 'Cl': 35.45, 'S': 96.06}
_SAL_INDICES = [2 + i for i, e in enumerate(elements) if e in _ELEMENT_MASSES]
_SAL_MASSES  = [_ELEMENT_MASSES[e] for e in elements if e in _ELEMENT_MASSES]

# Reference crust: the model derives the mineralogy from the two composition axes -- mantle
# molar Mg/Si and the core-formation oxygen fugacity dIW -- instead of a named composition like
# 'basalt_49', or the (T_p, mg_si_ratio) pair that preceded them.
REF_MG_SI = float(EARTH_MANTLE_MG_SI)
REF_DIW   = float(EARTH_DELTA_IW)

COMP_COLORS = {
    'komatiite_42': '#7b2d8b',
    'komatiite_44': '#c44bc4',
    'basalt_47':    '#e08040',
    'basalt_49':    '#d4b000',
    'basalt_51':    '#6abf69',
}
HAB_MARKERS    = {'converged': 'o', 'timeout': 's'}
FAILED_MARKERS = {'out_of_domain': 's'}

DA_LEGEND = [
    Line2D([0], [0], color='k', linestyle='-',  linewidth=1.8, label='Da < 1 (kinetic)'),
    Line2D([0], [0], color='k', linestyle='--', linewidth=1.8, label='Da ≥ 1 (thermodynamic)'),
    Line2D([0], [0], color='k', linestyle=':',  linewidth=1.8, label='$T_\\mathrm{sf}$ at floor (274 K)'),
    Line2D([0], [0], color='k', linestyle='-.', linewidth=0.8, label='Equilibrium temperature'),
]

PANEL_COLS = ['T', 'P_CO2', 'pH', 'salinity']

# Colour scale for core-formation dIW, used by EVERY figure that puts dIW on a colourbar so the
# redox axis is recognisable across them. Mg/Si and the other facet variables keep their own maps.
#
# Truncated: cmr.emerald runs from #000000, and the swept dIW range (-5 to -1) put its lower half
# into near-black, where the reduced end members were indistinguishable from each other and from
# the axis text on a white background. Trimming to the middle 70% keeps the emerald identity and
# spans luminance 0.19-0.81, so all seven grid values separate.
DIW_CMAP = cmr.get_sub_cmap('cmr.emerald', 0.25, 0.95)


def equilbrium_temperature(instellation, albedo=0.3, greenhouse=0.5):
    return (((1-albedo) * instellation * SOLAR_CONSTANT)/(4 * STEFAN_BOLTZMANN * greenhouse)) ** 0.25

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _salinity_from_y(y_list):
    try:
        return sum(
            float(y_list[i][-1]) * mass
            for i, mass in zip(_SAL_INDICES, _SAL_MASSES)
            if len(y_list) > i and len(y_list[i]) > 0
        )
    except Exception:
        return np.nan


_DIAG_NAN = {'da': np.nan, 'calcite_si': np.nan, 'ocean_si': np.nan, 'alk_flux': np.nan, 'pH': np.nan}


@functools.lru_cache(maxsize=None)
def _crust_composition(mantle_mg_si, delta_iw):
    """CIPW-norm crust mineralogy for a mantle Mg/Si and core-formation dIW (cached)."""
    return mineral_composition(mantle_mg_si, delta_iw)


def _crust_composition_of(d):
    """Crust mineralogy for a run: the stored composition if present, else the two axes.

    Runs predating the MAGEMin table stored `mantle_potential_temperature` / `mg_si_ratio`, which
    the current pipeline cannot reproduce -- T_p is no longer an input and `mg_si_ratio` was a
    multiplier on 1.23, not a Mg/Si. Those runs are only plottable via their stored
    `crust_composition`; without one there is nothing honest to draw, so say so rather than
    silently substituting Earth.
    """
    stored = d.get('crust_composition')
    if stored:
        return stored
    if 'mantle_potential_temperature' in d and 'mantle_mg_si' not in d:
        raise KeyError(
            'legacy run (mantle_potential_temperature) with no stored crust_composition: its '
            'mineralogy cannot be reconstructed by the current pipeline')
    return _crust_composition(float(d.get('mantle_mg_si', REF_MG_SI)),
                              float(d.get('delta_iw', REF_DIW)))


def _pore_conditions(d):
    """Reconstruct the seafloor/pore-space state of a run's final step from its JSON.

    Returns (b_ocean, P_pore, T_pore, T_seafloor, P_CO2, crust_rate, J_total) matching
    what Planet.dY_dt saw on the last evaluation.
    """
    y_list = d['data']['y']
    n_elements = len(y_list) - 3  # y = [P_CO2, P_H2O, *elements, r_avg]
    b_ocean = np.maximum(np.array([float(y_list[i][-1]) for i in range(2, 2 + n_elements)]), 0.0)

    mass    = float(d.get('mass',   5.972e24))
    radius  = float(d.get('radius', 6.371e6))
    gravity = G * mass / radius**2

    P_CO2     = float(d['P_CO2']) * 1e5             # bar → Pa
    T_surface = float(d['T'])
    P_H2O     = float(y_list[1][-1]) if y_list[1] else 0.0
    P_surface = float(d['background_pressure']) + P_CO2 + P_H2O

    T_seafloor = max(1.02 * T_surface - 16.7, 274.0)
    T_pore     = T_seafloor + 9
    P_pore     = P_surface + 1000 * gravity * float(d['ocean_depth'])

    # Hydrothermal flux scales with crust production (planet.py no longer splits out f_HT).
    crust_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * float(d['crust_production_rate'])
    J_total    = J_ref_normalised * (crust_rate / rate_ref)

    return b_ocean, P_pore, T_pore, T_seafloor, P_CO2, crust_rate, J_total


def _sedimentation_rate(d, b_ocean, P_pore, T_seafloor):
    """Sediment accumulation rate (m/s) from abiotic carbonate + silica burial, as in dY_dt."""
    from kamino.precipitation import get_precipitation
    from kamino.chemistry import c_idx as _c_idx, si_idx as _si_idx

    rw = bool(d.get('reverse_weathering', False))
    fast_minerals = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals
    try:
        F_prec, _, _ = get_precipitation(P_pore, T_seafloor, b_ocean, fast_minerals,
                                         precipitation_timescale=float(d.get('tau_prec', 1e5 * YR)))
        if rw:
            F_rw, _, _ = get_precipitation(P_pore, T_seafloor, b_ocean, list(reverse_weathering_minerals),
                                           precipitation_timescale=float(d.get('tau_rw', 5e6 * YR)))
            F_prec = F_prec + F_rw
    except Exception:
        return None  # fall back to the reference sedimentation rate inside the weathering law

    F_carb = max(0.0, -float(F_prec[_c_idx]))
    F_sil  = max(0.0, -float(F_prec[_si_idx]))
    ocean_water_per_area = float(d['ocean_depth']) * 1000.0
    s_terr = _S_TERR_EARTH * (float(d.get('land_fraction', 0.0)) / 0.3)
    return (F_carb * 0.100 / 2710.0 + F_sil * 0.060 / 2650.0) * ocean_water_per_area + s_terr


def _diag_from_json(d):
    """Compute weathering diagnostics from the final state stored in a JSON file.

    Returns a dict with:
      da         — Damköhler number for alkalinity
      calcite_si — calcite SI of the pore fluid *before* secondary precipitation
                   (the driving force for pore-space calcite; always valid)
      ocean_si   — calcite SI of the ocean water at seafloor T/P
                   (meaningful only when ocean Ca > 0; NaN when Ca ≈ 0)
      alk_flux   — net seafloor alkalinity flux (Tmol eq/yr)
    """
    try:
        y_list = d.get('data', {}).get('y', [])
        n_elements = len(y_list) - 3  # y_list = [P_CO2, P_H2O, *elements, r_avg]
        if not y_list or n_elements < 7:
            return _DIAG_NAN

        b_ocean, P_pore, T_pore, T_seafloor, P_CO2, crust_rate, J_total = _pore_conditions(d)

        # Pore space precipitates clays only (planet.pore_precipitating_minerals);
        # carbonates and reverse-weathering clays form in the ocean sediments.
        pore_minerals = list(clay_minerals)

        flux, diag = get_weathering_flux(
            P_pore, T_pore, P_CO2, b_ocean,
            alpha=float(d.get('alpha', 1.43)),
            rate=crust_rate, J=J_total,
            crust_composition=_crust_composition_of(d),
            sedimentation_rate=_sedimentation_rate(d, b_ocean, P_pore, T_seafloor),
            precipitating_minerals=pore_minerals,
        )

        from kamino.precipitation import get_precipitation
        from kamino.chemistry import ChemistryError, ca_idx as _ca_idx

        # pH is recomputed here rather than read from the JSON: self._pH is planet.py's
        # side-effect value, and unlike 'T' it is not re-evaluated on the accepted final state
        # before the file is written, so it can still be a Jacobian probe's value.
        # This reuses T_seafloor/P_pore from _pore_conditions, which are now correct because
        # they are built from the recomputed T_surface -- so this pH is the equilibrium pH of
        # the ACTUAL final ocean composition, matching exactly what dY_dt computes on a real
        # (non-probe) trajectory step: get_precipitation with the fast-precipitating assemblage
        # at seafloor conditions. Costs one more PHREEQC solve, on top of the several this
        # function already does, so it is only applied where diagnostics are already paid for
        # (_add_diag_columns callers), not in the cheap load_data() pass.
        try:
            fast_minerals = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals
            _, pH_recomputed, _ = get_precipitation(
                P_pore, T_seafloor, b_ocean, fast_minerals,
                precipitation_timescale=float(d.get('tau_prec', 1e5 * YR)))
            pH_recomputed = float(pH_recomputed)
        except (ChemistryError, Exception):
            pH_recomputed = np.nan

        # Pore SI: the pore space now precipitates clays only, so Calcite no longer
        # appears in secondary_SI. Evaluate it directly on the post-weathering pore
        # fluid (b_pore) — the driving force for pore-space calcite.
        try:
            b_pore = np.asarray(diag['b_pore'], dtype=float)
            _, _, si_p = get_precipitation(P_pore, T_pore, np.maximum(b_pore, 0.0), ['Calcite'],
                                           precipitation_timescale=1e6 * YR)
            calcite_si = float(si_p.get('Calcite', np.nan))
        except (ChemistryError, Exception):
            calcite_si = np.nan

        # Ocean SI: only reliable when Ca_ocean > 0; set to NaN otherwise to avoid
        # the spurious -∞ that results when ocean Ca has been depleted to the ODE floor.
        if b_ocean[_ca_idx] > 1e-6:
            try:
                _, _, si_o = get_precipitation(P_pore, T_seafloor, b_ocean, ['Calcite'],
                                               precipitation_timescale=1e6 * YR)
                ocean_si = float(si_o.get('Calcite', np.nan))
            except (ChemistryError, Exception):
                ocean_si = np.nan
        else:
            ocean_si = np.nan

        alk_flux = float(flux[alk_idx]) * A_seafloor * YR / 1e12  # Tmol eq/yr

        return {'da': float(diag['Da']), 'calcite_si': calcite_si,
                'ocean_si': ocean_si, 'alk_flux': alk_flux, 'pH': pH_recomputed}
    except Exception:
        return _DIAG_NAN.copy()


def _diag_from_run(d):
    """Diagnostics for one run, preferring the values `Planet.time_evolve` recorded.

    The model computes `da`, `calcite_si`, `ocean_si`, `alk_flux` and the seafloor pH on every
    step anyway, and since 2026-08-27 writes them for the accepted final state into a
    "diagnostics" block. Reading them costs nothing; reconstructing them here costs ~0.9 s of
    PHREEQC per run, which is ~30 minutes for a 2000-run sweep. Runs written before that change
    have no block, so `_diag_from_json` remains as the fallback and is exercised by them.
    """
    block = d.get('diagnostics')
    if isinstance(block, dict) and 'da' in block:
        out = {}
        for key, src in (('da', 'da'), ('calcite_si', 'calcite_si'), ('ocean_si', 'ocean_si'),
                         ('alk_flux', 'alk_flux'), ('pH', 'pH_seafloor')):
            v = block.get(src)
            out[key] = np.nan if v is None else float(v)
        return out
    return _diag_from_json(d)


_DIAG_CACHE = {}

# Directory the run JSONs were loaded from. `_add_diag_columns` re-reads each run to compute its
# diagnostics, and until now it used the FIGURE output directory for that -- fine in __main__,
# where they are the same directory, but silently wrong (every diagnostic NaN, and the cache
# unwritable) for any caller that renders elsewhere. load_data records the real location here.
RUN_PATH = None

# Bump when _diag_from_json changes WHAT it computes, so stale on-disk records are discarded.
_DIAG_VERSION = 1
_DIAG_CACHE_FILE = '.plot_diag_cache.json'
_diag_cache_loaded = set()
_diag_cache_dirty = set()


def _diag_cache_key(fpath):
    """Identity of a run's diagnostics: its path plus the file's size and mtime.

    Re-running a sweep rewrites the JSON, which changes both, so a stale record can never be
    served for a run that has been recomputed.
    """
    st = os.stat(fpath)
    return f'{st.st_size}:{int(st.st_mtime)}'


def _load_diag_cache(output_path):
    """Read the sidecar cache for `output_path` once per process."""
    if output_path in _diag_cache_loaded:
        return
    _diag_cache_loaded.add(output_path)
    path = os.path.join(output_path, _DIAG_CACHE_FILE)
    try:
        with open(path) as fh:
            blob = json.load(fh)
    except Exception:
        return
    if blob.get('version') != _DIAG_VERSION:
        print(f"  diagnostics cache at {path} is version {blob.get('version')}, "
              f"expected {_DIAG_VERSION} -- ignoring it.")
        return
    n = 0
    for fpath, entry in blob.get('runs', {}).items():
        try:
            if _diag_cache_key(fpath) == entry['key']:
                _DIAG_CACHE[fpath] = entry['rec']
                n += 1
        except OSError:
            continue        # run file has gone away
    if n:
        print(f"  reusing cached diagnostics for {n} run(s) from {_DIAG_CACHE_FILE}")


def _save_diag_cache(output_path):
    """Write the sidecar cache if anything new was computed for `output_path`."""
    if output_path not in _diag_cache_dirty:
        return
    _diag_cache_dirty.discard(output_path)
    # MERGE with whatever is already on disk. `_DIAG_CACHE` holds only the runs this process
    # happened to touch, so writing it verbatim would shrink a complete cache down to the subset
    # of one partial render -- silently discarding hours of work.
    runs = {}
    path = os.path.join(output_path, _DIAG_CACHE_FILE)
    try:
        with open(path) as fh:
            existing = json.load(fh)
        if existing.get('version') == _DIAG_VERSION:
            runs.update(existing.get('runs', {}))
    except Exception:
        pass
    for fpath, rec in _DIAG_CACHE.items():
        try:
            runs[fpath] = {'key': _diag_cache_key(fpath),
                           'rec': {k: (None if v is None or (isinstance(v, float) and np.isnan(v))
                                       else float(v)) for k, v in rec.items()}}
        except OSError:
            continue
    try:
        with open(path, 'w') as fh:
            json.dump({'version': _DIAG_VERSION, 'runs': runs}, fh)
        print(f"  wrote diagnostics cache for {len(runs)} run(s) -> {path}")
    except OSError as exc:
        print(f"  could not write diagnostics cache ({exc}); results are unaffected.")


def _add_diag_columns(df, output_path=None):
    """Add da, calcite_si, alk_flux and (corrected) pH columns by re-reading each JSON file.

    Each record costs ~0.9 s of PHREEQC (a full weathering equilibration plus three saturation
    solves), so 2000 runs is ~30 minutes. Results are therefore cached twice: in `_DIAG_CACHE`
    for the several figures that request overlapping subsets within one process, and in a
    sidecar JSON beside the runs so a later invocation -- re-rendering after a styling change,
    which is the common case -- pays nothing. The sidecar is keyed on each run file's size and
    mtime, so re-running a sweep invalidates its own entries automatically.

    Overwrites the 'pH' column (the JSON's stored side-effect value, which unlike 'T' is not
    re-evaluated on the final state) wherever _diag_from_json succeeded -- callers of this
    function already pay the PHREEQC cost the recompute needs, so the correction is free here.
    Rows where the recompute itself failed keep the original stored 'pH' rather than becoming
    NaN, since a stale-but-present value is more useful than none for a plot.
    """
    run_path = RUN_PATH or output_path
    _load_diag_cache(run_path)
    records = []
    todo = sum(1 for n in df['name']
               if os.path.join(run_path, f'{n}.json') not in _DIAG_CACHE)
    if todo:
        print(f"  computing diagnostics for {todo} run(s) (~{todo * 0.9 / 60:.1f} min)...",
              flush=True)
    for name in df['name']:
        fpath = os.path.join(run_path, f'{name}.json')
        if fpath in _DIAG_CACHE:
            records.append(_DIAG_CACHE[fpath])
            continue
        try:
            with open(fpath) as fh:
                d = json.load(fh)
            rec = _diag_from_run(d)
        except Exception:
            rec = _DIAG_NAN.copy()
        _DIAG_CACHE[fpath] = rec
        _diag_cache_dirty.add(run_path)
        records.append(rec)
    _save_diag_cache(run_path)
    diag_df = pd.DataFrame(records, index=df.index)
    df = df.assign(**diag_df.drop(columns=['pH']))
    df['pH'] = diag_df['pH'].where(diag_df['pH'].notna(), df['pH'])
    return df


def load_data(output_path):
    global RUN_PATH
    RUN_PATH = output_path
    files = sorted(glob.glob(os.path.join(output_path, 'planet_*.json')))
    rows = []
    for f in files:
        with open(f) as fh:
            d = json.load(fh)
        if 'termination' not in d:
            print(f"  Skipping (no termination): {os.path.basename(f)}")
            continue

        name = d.get('name', '')
        y_list = d.get('data', {}).get('y', [])
        salinity = _salinity_from_y(y_list) if y_list else np.nan

        rows.append({
            'name':               name,
            'instellation':       float(d['instellation']),
            'outgassing':         float(d['outgassing']),
            'crust_production':   float(d['crust_production_rate']),
            'reverse_weathering': bool(d.get('reverse_weathering', False)),
            'ocean_depth':        float(d['ocean_depth']),
            'mg_si':              float(d.get('mantle_mg_si', REF_MG_SI)),
            'delta_iw':           float(d.get('delta_iw', REF_DIW)),
            'f_HT':               float(d.get('f_HT', 0.0)),
            # Chemistry constants are swept axes (parameter_sweep.py). Runs differing only in
            # these must not be pooled into one line -- see _ref_chem.
            'alpha':              float(d.get('alpha', ALPHA_REF)),
            'kd_mg':              float(d.get('kd_mg_ht', KD_MG_HT)),
            'k_na':               float(d.get('k_na_cont_removal', K_NA_CONT_REMOVAL)),
            'land_fraction':      float(d.get('land_fraction', 0.0)),
            'termination':        d['termination'],
            'domain_wall':        d.get('domain_wall'),   # None for pre-domain-event runs
            'end_time_yr':        d.get('end_time_yr', np.nan),
            # Stored 'T' is trusted: Planet.time_evolve re-evaluates it on the accepted final
            # state before writing, so it corresponds to the P_CO2 in the same file. Output
            # written before that fix needs plot_legacy.upgrade(), which recomputes it.
            'T':                  float(d.get('T', np.nan)),
            'P_CO2':              d.get('P_CO2', np.nan),
            'pH':                 d.get('pH', np.nan),  # corrected in _diag_from_json when available
            'salinity':           salinity,
        })
    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} simulations.")
    return df


def _ref_crust(df):
    """Mask for the reference crust: Earth's mantle Mg/Si and core-formation dIW."""
    return np.isclose(df['mg_si'], REF_MG_SI) & np.isclose(df['delta_iw'], REF_DIW)


# Chemistry constants that parameter_sweep.py can vary, with axis labels for the sweep plots.
CHEM_KNOBS = {
    'alpha': r'Reactive area scaling $\alpha$',
    'kd_mg': r'Mg$\rightarrow$Ca exchange $k_{Mg}$',
    'k_na':  r'Na sink $k_{Na}$',
}

CHEM_SHIPPED = {'alpha': ALPHA_REF, 'kd_mg': KD_MG_HT, 'k_na': K_NA_CONT_REMOVAL}

CHEM_OVERRIDE = {}     # set from the CLI to choose which chemistry the main plots show
_chem_pinned = set()   # (column, value) already reported, so repeated _base calls print once


def _chem_reference(df, col):
    """Which value of a swept chemistry constant the main plots should show.

    Most runs wins; ties break away from a disabled term (k_mg=0 / k_na=0 are ablations and
    should never become the headline chemistry) and then toward the shipped default.
    """
    if col in CHEM_OVERRIDE:
        return CHEM_OVERRIDE[col]
    counts = df[col].value_counts()
    best = counts.max()
    tied = [v for v, n in counts.items() if n == best]
    return min(tied, key=lambda v: (v == 0, abs(v - CHEM_SHIPPED[col])))


def _ref_chem(df):
    """Mask pinning each chemistry constant to a single value.

    Everything except plot_chemistry_constants shows ONE chemistry. Without this, runs that
    differ only in alpha/kd_mg/k_na fall into the same (instellation, outgassing, crust) group
    and get drawn as a single line through several different models.
    """
    mask = pd.Series(True, index=df.index)
    for col in CHEM_KNOBS:
        if col not in df.columns or df[col].nunique() <= 1:
            continue
        ref = _chem_reference(df, col)
        mask &= (df[col] == ref)
        if (col, ref) not in _chem_pinned:
            _chem_pinned.add((col, ref))
            others = sorted(v for v in df[col].unique() if v != ref)
            print(f"  Pinning {col} = {ref:g} for the main plots "
                  f"(also present: {', '.join(f'{v:g}' for v in others)}).")
    return mask


def _base(df):
    """Sweep 1: reference crust and chemistry, rw=True, depth=3000, outgassing>0, ocean world."""
    return df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_chem(df) &
        (df['ocean_depth'] == 3000) &
        (df['land_fraction'] == 0.0) &
        (df['outgassing'] > 0)
    ]


# ---------------------------------------------------------------------------
# Shared plot helpers
# ---------------------------------------------------------------------------

def _panel_groups(split):
    """Return [(cols, filename_suffix), ...] — one entry normally, two when split."""
    if split:
        return [(['T', 'P_CO2'], '_tp'), (['pH', 'salinity'], '_chem')]
    return [(PANEL_COLS, '')]


def _style_axes(axes, cols, x_lims=(0.25, 1.45)):
    """Style a set of axes given the column names they represent."""
    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(*x_lims)
    for ax, col in zip(axes, cols):
        if col == 'T':
            ax.set_ylabel('Temperature (K)')
            ax.axhspan(T_SNOWBALL - 25, T_SNOWBALL, color='blue', alpha=0.12)
            ax.axhspan(T_RUNAWAY - 20,  T_RUNAWAY,  color='red',  alpha=0.12)
            ax.set_ylim(235, 360)
            s_eq = np.linspace(x_lims[0], x_lims[1], 300)
            ax.plot(s_eq, equilbrium_temperature(s_eq), color='k', linestyle='-.',
                    linewidth=0.8, zorder=1, alpha=0.7)
        elif col == 'P_CO2':
            ax.set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 20)
        elif col == 'pH':
            ax.set_ylabel('Ocean pH')
            ax.set_ylim(4.5, 12)
        elif col == 'salinity':
            ax.set_ylabel('Salinity (g/kg)')
            ax.set_yscale('log')
            ax.set_ylim(1e-1, 1.1e2)
        elif col == 'calcite_si':
            ax.set_ylabel('Calcite SI')
            ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
        elif col == 'alk_flux':
            ax.set_ylabel('Alk. flux (Tmol eq/yr)')
            ax.set_yscale('symlog', linthresh=0.001)
            ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    axes[-1].set_xlabel('Instellation (S/S₀)')


def _add_colorbar(fig, ax, cmap, norm, label, ticks=None, ticklabels=None, aspect=30):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, location='right', pad=0.02, aspect=aspect)
    cbar.set_label(label)
    if ticks is not None:
        cbar.set_ticks(ticks)
        if ticklabels is not None:
            cbar.set_ticklabels(ticklabels)
    return cbar


def _make_legend_handles(show_markers=True, prefix_handles=None):
    """Build legend handle list: prefix (default DA_LEGEND) + optional marker entries."""
    handles = list(prefix_handles if prefix_handles is not None else DA_LEGEND)
    if show_markers:
        handles += [plt.scatter([], [], marker=m, s=28, color='k', label=TERM_LABELS[t])
                    for t, m in HAB_MARKERS.items()]
        handles += [plt.scatter([], [], marker=m, s=50, label=TERM_LABELS[t],
                                facecolors='none', edgecolors='k', linewidths=1.4)
                    for t, m in FAILED_MARKERS.items()]
    return handles


def _legend_ncol(handles, fallback):
    """Deprecated: prefer `_add_figure_legend`, which measures instead of guessing."""
    return len(handles) if presentation else fallback


def _add_figure_legend(fig, axes, handles, loc='outside lower center', **kw):
    """Figure legend wrapped so it is never wider than the panels it describes.

    `_legend_ncol` put every entry on one row in presentation mode, which for the marker-bearing
    legends ran past the axes and widened the whole saved figure -- defeating the point of sizing
    it to a column. This measures the rendered legend and drops a column at a time until it fits
    inside the panel block, so a page-width figure stays page width.
    """
    axs = [a for a in np.ravel(np.asarray(axes, dtype=object)) if a is not None]
    fig.canvas.draw()                       # a renderer is needed to measure anything
    boxes = [a.get_window_extent() for a in axs]
    panel_w = max(b.x1 for b in boxes) - min(b.x0 for b in boxes)

    leg = None
    for ncol in range(len(handles), 0, -1):
        if leg is not None:
            leg.remove()
        leg = fig.legend(handles=handles, loc=loc, ncol=ncol, **kw)
        fig.canvas.draw()
        if leg.get_window_extent().width <= panel_w:
            break
    return leg


def _save_fig(fig, path, tight=False):
    r"""Write a figure at EXACTLY the size it was created with.

    `bbox_inches='tight'` is deliberately off. The style sets constrained_layout, which already
    fits every artist inside the canvas, and 'tight' then re-crops to the content box -- which
    for these figures came out LARGER than the requested width (3.43 in against the 3.32 in
    column). A figure wider than the column gets scaled down by \includegraphics, shrinking the
    type below the size the style chose. Pass tight=True for diagnostics, where exact width does
    not matter.
    """
    fig.savefig(path, **({'bbox_inches': 'tight'} if tight else {}))
    plt.close(fig)
    print(f"Saved {path}")


def _style_combined_col(axes_c, ci, n_cols, title='', cols=None):
    """Style one column of a multi-column grid: axis labels, tick visibility, x-label."""
    if cols is None:
        cols = PANEL_COLS
    col_axes = axes_c[:, ci]
    _style_axes(col_axes, cols)
    if title:
        axes_c[0, ci].set_title(title)
    if ci > 0:
        for ax in col_axes:
            ax.set_ylabel('')
        for row in range(axes_c.shape[0]):
            plt.setp(axes_c[row, ci].get_yticklabels(), visible=False)
    if ci != n_cols // 2:
        col_axes[-1].set_xlabel('')


def _plot_line_da_style(ax, x, y, da, color, at_floor=None, linewidth=1.8, alpha=0.8, zorder=3):
    """Draw a line: dotted where seafloor T is floored, dashed where Da≥1, solid where Da<1.

    Places an open circle at each kinetic↔thermodynamic transition (solid↔dashed only;
    transitions into/out of the dotted floor regime are not marked).
    """
    if len(x) < 2:
        return

    if at_floor is None:
        at_floor = np.zeros(len(x), dtype=bool)

    def _ls(i):
        if at_floor[i]:
            return ':'
        return '-' if (np.isfinite(da[i]) and da[i] < 1) else '--'

    seg_x = [x[0]]
    seg_y = [y[0]]
    current_ls = _ls(0)
    trans_x, trans_y = [], []

    for i in range(1, len(x)):
        ls = _ls(i)
        if ls == current_ls:
            seg_x.append(x[i])
            seg_y.append(y[i])
        else:
            ax.plot(seg_x, seg_y, color=color, linewidth=linewidth,
                    alpha=alpha, linestyle=current_ls, zorder=zorder)
            # Mark only solid↔dashed transitions (skip dotted floor segments)
            if {current_ls, ls} == {'-', '--'}:
                trans_x.append(seg_x[-1])
                trans_y.append(seg_y[-1])
            seg_x = [seg_x[-1], x[i]]
            seg_y = [seg_y[-1], y[i]]
            current_ls = ls

    if len(seg_x) >= 2:
        ax.plot(seg_x, seg_y, color=color, linewidth=linewidth,
                alpha=alpha, linestyle=current_ls, zorder=zorder)

    if trans_x:
        ax.scatter(trans_x, trans_y, facecolors='none', edgecolors=color,
                   s=30, linewidths=1.2, zorder=5)


def _best_operating_point(pool, col, what):
    """The (outgassing, crust production) pair spanning the most distinct values of `col`.

    Every one-variable sweep fixes those two at whatever defaults it was run with, and those
    defaults drifted over the project's life (outgassing 1.0 -> 0.1). Finding the pair rather
    than hardcoding it is why these figures survive that drift; each caller had its own copy.

    Returns (subset, outgassing, crust) or None when nothing varies.
    """
    counts = pool.groupby(['outgassing', 'crust_production'])[col].nunique()
    if counts.empty or counts.max() < 2:
        print(f"No {what} variation at a fixed (outgassing, crust) -- skipping.")
        return None
    best_o, best_c = counts.idxmax()
    subset = pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)]
    print(f"{what} plot: using outgassing={best_o:g}, crust={best_c:g} "
          f"({subset[col].nunique()} values).")
    return subset, best_o, best_c


def _x_limits(subset, default=(0.25, 1.45), frac=0.05):
    """Instellation axis limits with a small margin, falling back when the range is degenerate."""
    lo, hi = subset['instellation'].min(), subset['instellation'].max()
    if pd.isna(lo):
        return default
    margin = (hi - lo) * frac if hi != lo else 0.1
    return (lo - margin, hi + margin)


def _value_norm(values, pad=0.05):
    """Colour norm for a set of facet values.

    Log where they span decades AND none is zero -- the k_mg / k_na ablations set the constant to
    exactly 0, which a log norm cannot place. `pad` widens a linear range so the end members are
    not at the extreme ends of the colour map; pass 0 to match a range exactly.
    """
    lo, hi = min(values), max(values)
    if lo > 0 and hi / lo >= 10:
        return mcolors.LogNorm(vmin=lo, vmax=hi)
    span = hi - lo
    return mcolors.Normalize(vmin=lo - pad * span if span else lo - 1,
                             vmax=hi + pad * span if span else hi + 1)


def _colorbar_ticks(values, max_ticks=10):
    """Explicit ticks only when there are few enough to label individually."""
    if len(values) > max_ticks:
        return None, None
    return list(values), [f'{v:g}' for v in values]


def _faceted_lines(subset, col, values, colours, cmap, norm, cbar_label, stem, output_path,
                   split_panels=False, show_markers=False, x_lims=None,
                   ticks=None, ticklabels=None, aspect_per_row=None,
                   width='single', height=None):
    """The standard one-variable faceted line figure, one file per panel group.

    Every figure that plots T / P_CO2 / pH / salinity against instellation with one line per value
    of a single variable -- ocean depth, a chemistry constant, mantle Mg/Si, dIW -- shares this
    body. They differed only in which rows they select and how the lines are labelled, so the
    selection stays with the caller and the assembly lives here.

    `colours` is parallel to `values` rather than derived from `norm`, so a caller can colour by
    something other than the facet value itself (plot_legacy does, by melt SiO2).

    `aspect_per_row` is multiplied by the panel count rather than passed through, because the
    number of panels differs between the split and combined panel groups.
    """
    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        fig, axes = plt.subplots(n_rows, 1, sharex=True,
                                 figsize=figure_size(width, height, n_rows))
        for value, colour in zip(values, colours):
            group = subset[subset[col] == value].sort_values('instellation')
            if not group.empty:
                _plot_group_on_axes(axes, group, colour, show_markers=show_markers, cols=cols)
        _style_axes(axes, cols, **({} if x_lims is None else {'x_lims': x_lims}))
        _add_colorbar(fig, list(axes), cmap, norm, cbar_label, ticks=ticks,
                      ticklabels=ticklabels,
                      **({} if aspect_per_row is None else {'aspect': n_rows * aspect_per_row}))
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        _save_fig(fig, os.path.join(output_path, f'{stem}{sfx}.png'))


def _plot_group_on_axes(axes, group, color, linestyle='-', show_markers=True, cols=None):
    if cols is None:
        cols = PANEL_COLS
    hab = group[group['termination'].isin(HABITABLE)]
    if hab.empty:
        return
    non_hab = group[~group['termination'].isin(HABITABLE)]
    has_da = 'da' in group.columns

    for ax, col in zip(axes, cols):
        if len(group) > 1:
            if has_da and linestyle == '-':
                T_sf = 1.02 * group['T'].values - 16.7
                at_floor = T_sf <= 274.001
                _plot_line_da_style(ax, group['instellation'].values, group[col].values,
                                    group['da'].values, color, at_floor=at_floor)
            else:
                ax.plot(group['instellation'], group[col], color=color, linewidth=1.8,
                        alpha=0.8, linestyle=linestyle, zorder=3)
        if not show_markers:
            continue
        for term, hab_marker in HAB_MARKERS.items():
            sub = hab[hab['termination'] == term]
            if not sub.empty:
                ax.scatter(sub['instellation'], sub[col], color=color,
                           marker=hab_marker, s=28, zorder=4)
        for _, row in non_hab.iterrows():
            marker = FAILED_MARKERS.get(row['termination'], 'x')
            val = row[col]
            if np.isfinite(val):
                ax.scatter(row['instellation'], val, marker=marker, s=55,
                           facecolors='none', color=color, zorder=4, linewidths=1.4)


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

def plot_faceted_lines(df, output_path, all_results=True, multiple_plots=False,
                       split_panels=False, sequence=False, width='double', height=None):
    """T, P_CO2, pH, salinity vs instellation per crust rate, coloured by outgassing."""
    base = _base(df)
    base = _add_diag_columns(base, output_path)

    if all_results:
        crust_rates     = sorted(base['crust_production'].unique())
        outgassing_vals = sorted(base['outgassing'].unique())
    else:
        crust_rates = [0.1, 1, 10]
        outgassing_vals = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)

        if multiple_plots:
            for c in crust_rates:
                subset_c = base[base['crust_production'] == c]
                fig, axes = plt.subplots(n_rows, 1, sharex=True,
                                         figsize=diagnostic_size(n_rows, 1, col_width=7.0, pad=0.0))
                for o in outgassing_vals:
                    group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
                    if not group.empty:
                        _plot_group_on_axes(axes, group, cmap(norm(o)), cols=cols)
                _style_axes(axes, cols)
                _add_colorbar(fig, list(axes), cmap, norm, 'Earth Outgassing',
                              ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                              aspect=n_rows * 7.5) # type: ignore
                _h = _make_legend_handles()
                fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
                fig.suptitle(f'Crust production = {c}× Earth')
                _save_fig(fig, os.path.join(output_path, f'lines_crust_{c}{sfx}.png'))

        # Combined plot: crust rates as columns
        n_cols = len(crust_rates)
        # all_results=True is the full outgassing x crust diagnostic grid -- not page-sized.
        full_figsize = (diagnostic_size(n_rows, n_cols) if all_results
                        else figure_size(width, height, n_rows))
        fig_c, axes_c = plt.subplots(n_rows, n_cols, figsize=full_figsize,
                                      sharex=True, sharey='row', squeeze=False)
        for ci, c in enumerate(crust_rates):
            subset_c = base[base['crust_production'] == c]
            for o in outgassing_vals:
                group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes_c[:, ci], group, cmap(norm(o)),
                                        show_markers=all_results, cols=cols)
            _style_combined_col(axes_c, ci, n_cols, title=f'{c}×', cols=cols)

        _h = _make_legend_handles(show_markers=all_results)
        fig_c.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
        _add_colorbar(fig_c, list(axes_c.ravel()), cmap, norm, 'Earth Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_rows * 10)
        fig_c.suptitle('Earth crust production rate')
        fname = f'lines_combined{"_full" if all_results else ""}{sfx}.png'
        _save_fig(fig_c, os.path.join(output_path, fname), tight=all_results)

        if sequence:
            ref_crust, ref_out = 1.0, 1.0
            seq_scenarios = [
                ('single',      f'lines_seq_1{sfx}'),
                ('out_sweep',   f'lines_seq_2{sfx}'),
                ('crust_sweep', f'lines_seq_3{sfx}'),
            ]
            for scenario, seq_fname in seq_scenarios:
                fig_s, axes_s = plt.subplots(n_rows, n_cols, figsize=full_figsize,
                                              sharex=True, sharey='row', squeeze=False)
                for ci, c in enumerate(crust_rates):
                    for o in outgassing_vals:
                        show = (
                            (scenario == 'single'      and np.isclose(c, ref_crust) and np.isclose(o, ref_out)) or
                            (scenario == 'out_sweep'   and np.isclose(c, ref_crust)) or
                            (scenario == 'crust_sweep' and np.isclose(o, ref_out))
                        )
                        if not show:
                            continue
                        group = (base[(base['crust_production'] == c) & (base['outgassing'] == o)]
                                 .sort_values('instellation'))
                        if not group.empty:
                            _plot_group_on_axes(axes_s[:, ci], group, cmap(norm(o)), cols=cols, show_markers=all_results)
                    _style_combined_col(axes_s, ci, n_cols, title=f'{c}×', cols=cols)
                _add_colorbar(fig_s, list(axes_s.ravel()), cmap, norm, 'Earth Outgassing',
                              ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                              aspect=n_rows * 10)
                _h = _make_legend_handles(show_markers=all_results)
                fig_s.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
                fig_s.suptitle('Earth crust production rate')
                _save_fig(fig_s, os.path.join(output_path, seq_fname + '.png'))


def plot_ocean_depth_effect(df, output_path, show_markers=False, split_panels=False,
                            width='single', height=None):
    """T, P_CO2, pH, salinity vs instellation for Earth-like tectonics, coloured by ocean depth."""
    # The depth sweep fixes (outgassing, crust) at their sweep defaults and varies ocean_depth.
    # That default changed over time (outgassing 1.0 -> 0.1), so instead of hardcoding a value
    # (which silently produced an empty/one-depth plot), pick whichever (outgassing, crust) pair
    # actually spans the most distinct depths.
    pool = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_chem(df) &
        (df['land_fraction'] == 0.0)
    ]
    if pool.empty:
        print("No data for ocean depth sweep — skipping.")
        return
    picked = _best_operating_point(pool, 'ocean_depth', 'Ocean-depth')
    if picked is None:
        return
    subset = _add_diag_columns(picked[0], output_path)

    depths = sorted(subset['ocean_depth'].unique())
    # pad=0: depth is a physical range that should map to the colour map exactly.
    norm = _value_norm(depths, pad=0.0)
    cmap = cmr.bubblegum_r
    ticks, ticklabels = _colorbar_ticks(depths)
    _faceted_lines(subset, 'ocean_depth', depths, [cmap(norm(d)) for d in depths],
                   cmap, norm, 'Ocean Depth (m)', 'lines_ocean_depth', output_path,
                   split_panels=split_panels, show_markers=show_markers,
                   x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                   aspect_per_row=7.5, width=width, height=height)


def plot_chemistry_constants(df, output_path, show_markers=False, split_panels=False,
                             width='single', height=None):
    """Sweeps 4/5: T, P_CO2, pH, salinity vs instellation for each chemistry-constant value.

    One figure per constant that actually varies (alpha, kd_mg, k_na); the other two are held
    at their most common value so each figure isolates a single knob.
    """
    pool_all = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        (df['ocean_depth'] == 3000) &
        (df['land_fraction'] == 0.0)
    ]
    if pool_all.empty:
        print("No data for chemistry-constant sweep — skipping.")
        return

    varying = [c for c in CHEM_KNOBS if pool_all[c].nunique() > 1]
    if not varying:
        print("No chemistry-constant variation in the data — skipping.")
        return

    for col in varying:
        # Hold the other knobs fixed so this figure varies one thing only.
        held = pd.Series(True, index=pool_all.index)
        for other in varying:
            if other != col:
                held &= (pool_all[other] == _chem_reference(pool_all, other))
        pool = pool_all[held]

        picked = _best_operating_point(pool, col, col)
        if picked is None:
            continue
        subset = _add_diag_columns(picked[0], output_path)

        values = sorted(subset[col].unique())
        norm = _value_norm(values)
        cmap = cmr.ember
        ticks, ticklabels = _colorbar_ticks(values)
        _faceted_lines(subset, col, values, [cmap(norm(v)) for v in values],
                       cmap, norm, CHEM_KNOBS[col], f'lines_{col}', output_path,
                       split_panels=split_panels, show_markers=show_markers,
                       x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                       aspect_per_row=7.5, width=width, height=height)


def _composition_pool(df, ocean_depth=3000):
    """Runs usable for a composition figure: reference chemistry, land-free, one depth."""
    return df[
        df['reverse_weathering'] &
        _ref_chem(df) &
        (df['ocean_depth'] == ocean_depth) &
        (df['f_HT'] == 0.0) &
        (df['land_fraction'] == 0.0)
    ].copy()


def _composition_slice(pool):
    """Pick the (outgassing, crust) pair that spans the most crust values, and return it.

    The composition sweep fixes those two at their sweep defaults and varies the crust knob.
    That default drifted over time (outgassing 1.0 -> 0.1), so it is found rather than hardcoded.
    """
    def _spread(g):
        return max(g['mg_si'].nunique(), g['delta_iw'].nunique())
    counts = pool.groupby(['outgassing', 'crust_production']).apply(_spread, include_groups=False)
    if counts.empty or counts.max() <= 1:
        return None
    best_o, best_c = counts.idxmax()
    return pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)], best_o, best_c


def plot_crust_composition(df, output_path, split_panels=False, show_markers=False,
                           ocean_depth=3000, width='single', height=None):
    """T, P_CO2, pH, salinity vs instellation, one figure per crust-composition axis.

    Emits a SEPARATE figure for each axis that varies -- mantle Mg/Si and core-formation dIW --
    holding the other at its Earth reference. The previous version picked whichever axis varied
    first and silently dropped the rest, so a sweep varying both produced only the Mg/Si figure
    and the dIW figure could not be made at all.

    Pre-MAGEMin sweeps that varied a NAMED composition are not handled here; see
    plot_legacy.plot_named_compositions.
    """
    pool = _composition_pool(df, ocean_depth)
    if pool.empty:
        print(f"No composition data at depth {ocean_depth:g} m -- skipping.")
        return
    sliced = _composition_slice(pool)
    if sliced is None:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping.")
        return
    subset, best_o, best_c = sliced

    # Which axes actually vary here. Each becomes its own figure; the others are held at
    # reference so a line is a cut through the grid rather than a mixture of compositions.
    axes_spec = []
    if subset['mg_si'].nunique() > 1:
        axes_spec.append(('mg_si', 'Mantle Mg/Si', lambda v: f'{v:g}',
                          sorted(subset['mg_si'].unique())))
    if subset['delta_iw'].nunique() > 1:
        axes_spec.append(('delta_iw', r'Core-formation $\Delta$IW', lambda v: f'{v:+g}',
                          sorted(subset['delta_iw'].unique())))
    if not axes_spec:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping.")
        return

    subset = _add_diag_columns(subset, output_path)
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'

    for key, label, fmt, values in axes_spec:
        # Hold every OTHER composition axis at its reference value.
        cut = subset
        held = []
        for other, ref in (('mg_si', REF_MG_SI), ('delta_iw', REF_DIW)):
            if other != key and other in subset.columns and subset[other].nunique() > 1:
                cut = cut[np.isclose(cut[other], ref)]
                held.append(f'{other}={ref:g}')
        if cut.empty or cut[key].nunique() < 2:
            print(f"  {key}: fewer than 2 values once {', '.join(held)} held -- skipping.")
            continue
        values = [v for v in values if v in set(cut[key])]
        print(f"Crust composition [{key}] depth={ocean_depth:g}: {len(values)} values, "
              f"outgassing={best_o:g}, crust={best_c:g}"
              + (f", holding {', '.join(held)}" if held else ""))

        cmap = DIW_CMAP if key == 'delta_iw' else cmr.gem
        numeric = [float(v) for v in values]
        norm = _value_norm(numeric)
        cbar_label, ticklabels = label, [fmt(v) for v in values]

        _faceted_lines(cut, key, values, [cmap(norm(n)) for n in numeric],
                       cmap, norm, cbar_label, f'lines_crust_{key}{tag}', output_path,
                       split_panels=split_panels, show_markers=show_markers,
                       ticks=numeric, ticklabels=ticklabels, width=width, height=height)


def _diw_title(dw, shown):
    """Column title naming the redox end-members and the Earth reference, not just the number."""
    if np.isclose(dw, REF_DIW):
        return f'{dw:+g} (Earth-like)'
    if len(shown) >= 2 and dw == min(shown):
        return f'{dw:+g} (reduced)'
    if len(shown) >= 2 and dw == max(shown):
        return f'{dw:+g} (oxidised)'
    return f'{dw:+g}'


def _three_columns(values, ref, n=3):
    """Pick `n` representative column values: the reference, then the extremes around it.

    A six-column grid is ~11000 px wide at presentation sizing -- fine on screen, unusable in a
    paper. Three columns carry the comparison that matters (reference vs each extreme) at a
    third of the width.
    """
    vals = sorted(values)
    if n >= len(vals):
        return vals
    centre = min(vals, key=lambda v: abs(v - ref))
    lo = [v for v in vals if v < centre]
    hi = [v for v in vals if v > centre]
    picked = [centre]
    if lo:
        picked.append(lo[0])
    if hi:
        picked.append(hi[-1])
    # If one side is empty, backfill from the other so `n` columns are still returned.
    pool = [v for v in vals if v not in picked]
    while len(picked) < n and pool:
        picked.append(pool.pop(len(pool) // 2))
    return sorted(picked)[:n]


def plot_composition_grid(df, output_path, split_panels=False, show_markers=False,
                          ocean_depth=3000, min_lines=2, n_cols=3,
                          width='double', height=None):
    """The composition factorial in the style of the basic sweep: dIW as columns, Mg/Si as colour.

    Same layout as `lines_combined_full` (crust rate as columns, outgassing as colour), with the
    two crust axes substituted. This uses the WHOLE factorial rather than the one-axis cuts
    `plot_crust_composition` draws: every panel is a fixed dIW, and every line within it a fixed
    Mg/Si, so an interaction between the two axes shows up as the family of lines changing shape
    from column to column rather than merely shifting.

    Columns with fewer than `min_lines` populated Mg/Si values are dropped -- a column carrying a
    single line says nothing and costs a fifth of the figure width.
    """
    pool = _composition_pool(df, ocean_depth)
    if pool.empty:
        print(f"No composition data at depth {ocean_depth:g} m -- skipping grid.")
        return
    sliced = _composition_slice(pool)
    if sliced is None:
        print("No crust composition sweep data -- skipping grid.")
        return
    subset, best_o, best_c = sliced
    if subset['mg_si'].nunique() < 2 or subset['delta_iw'].nunique() < 2:
        print("Composition grid needs both Mg/Si and dIW to vary -- skipping.")
        return

    subset = _add_diag_columns(subset, output_path)
    mg_vals = sorted(subset['mg_si'].unique())
    usable = [d for d in sorted(subset['delta_iw'].unique())
              if subset[subset['delta_iw'] == d]['mg_si'].nunique() >= min_lines]
    if not usable:
        print("No dIW column has enough Mg/Si values -- skipping grid.")
        return
    # Three columns: the Earth reference in the middle, the most REDUCED mantle on the left and
    # the most OXIDISED on the right (dIW increases with oxidation, so ascending order puts them
    # that way round).
    diw_vals = _three_columns(usable, REF_DIW, n_cols)
    dropped = sorted(set(subset['delta_iw'].unique()) - set(diw_vals))
    print(f"Composition grid depth={ocean_depth:g}: dIW columns {diw_vals} x "
          f"{len(mg_vals)} Mg/Si lines, outgassing={best_o:g}, crust={best_c:g}"
          + (f" (not shown: dIW {dropped})" if dropped else ""))

    # Mg/Si is linear and spans 0.5-2.0, so a linear norm -- unlike outgassing, which is log.
    norm = mcolors.Normalize(vmin=min(mg_vals), vmax=max(mg_vals))
    cmap = cmr.gem
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'

    for cols, sfx in _panel_groups(split_panels):
        n_rows, ncol = len(cols), len(diw_vals)
        figsize = figure_size(width, height, n_rows)
        fig, axes = plt.subplots(n_rows, ncol, figsize=figsize,
                                 sharex=True, sharey='row', squeeze=False)
        for ci, dw in enumerate(diw_vals):
            col_df = subset[subset['delta_iw'] == dw]
            for mg in mg_vals:
                group = col_df[col_df['mg_si'] == mg].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, cmap(norm(mg)),
                                        show_markers=show_markers, cols=cols)
            _style_combined_col(axes, ci, len(diw_vals), title=_diw_title(dw, diw_vals),
                                cols=cols)
        _add_colorbar(fig, list(axes.ravel()), cmap, norm, 'Mantle Mg/Si',
                      ticks=mg_vals, ticklabels=[f'{v:g}' for v in mg_vals],
                      aspect=n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        fig.suptitle(r'Core-formation $\Delta$IW')
        _save_fig(fig, os.path.join(output_path, f'lines_composition_grid{tag}{sfx}.png'))

    # The transpose: Mg/Si as columns, dIW as colour. Same data, and it is the better
    # arrangement when the dIW effect is the one being read off.
    diw_all = sorted(subset['delta_iw'].unique())
    mg_usable = [m for m in mg_vals
                 if subset[subset['mg_si'] == m]['delta_iw'].nunique() >= min_lines]
    if len(mg_usable) < 2:
        return
    mg_cols = _three_columns(mg_usable, REF_MG_SI, n_cols)
    norm_d = mcolors.Normalize(vmin=min(diw_all), vmax=max(diw_all))
    cmap_d = DIW_CMAP
    for cols, sfx in _panel_groups(split_panels):
        n_rows, ncol = len(cols), len(mg_cols)
        figsize = figure_size(width, height, n_rows)
        fig, axes = plt.subplots(n_rows, ncol, figsize=figsize,
                                 sharex=True, sharey='row', squeeze=False)
        for ci, mg in enumerate(mg_cols):
            col_df = subset[subset['mg_si'] == mg]
            for dw in diw_all:
                group = col_df[col_df['delta_iw'] == dw].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, cmap_d(norm_d(dw)),
                                        show_markers=show_markers, cols=cols)
            _style_combined_col(axes, ci, len(mg_cols), title=f'{mg:g}', cols=cols)
        _add_colorbar(fig, list(axes.ravel()), cmap_d, norm_d, r'Core-formation $\Delta$IW',
                      ticks=diw_all, ticklabels=[f'{v:+g}' for v in diw_all],
                      aspect=n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        fig.suptitle('Mantle Mg/Si')
        _save_fig(fig, os.path.join(output_path, f'lines_composition_grid_T{tag}{sfx}.png'))


def plot_composition_map(df, output_path, s_vals=(0.7, 0.9, 1.0, 1.1), ocean_depth=3000,
                         quantity='T', relative=True, min_cells=3,
                         width='double', height=2.6):
    """Map of outcome over the Mg/Si x dIW plane, one panel per instellation.

    The composition sweep is a near-factorial, so the two axes can be shown together rather than
    as separate one-axis cuts. That answers the question the cuts cannot: whether Mg/Si and dIW
    INTERACT, or whether their effects simply add.

    `relative` (default) plots the difference from the Earth-reference composition at the SAME
    instellation, on a diverging scale centred at zero. Absolute values do not work here: the
    Mg/Si = 0.5 column runs ~60 K hotter than everything else and swallows the whole colour
    range, leaving the other seven columns as indistinguishable pale cells. The anomaly is also
    the quantity of interest -- what the composition does, not what the instellation does.

    Rows and columns with fewer than `min_cells` in-domain runs are dropped. The sweep carries a
    few cross-design extras (Mg/Si = 1.6 and dIW = -5 exist at one point each), which otherwise
    appear as near-empty stripes through the middle of the grid.

    Cells with no run are blank; runs that left the model domain are marked, not coloured -- they
    have no steady state, so colouring them by final temperature would invent a result.
    """
    pool = _composition_pool(df, ocean_depth)
    if pool.empty:
        print(f"No composition data at depth {ocean_depth:g} m -- skipping map.")
        return
    sliced = _composition_slice(pool)
    if sliced is None:
        print("No crust composition sweep data -- skipping map.")
        return
    subset, best_o, best_c = sliced
    if subset['mg_si'].nunique() < 2 or subset['delta_iw'].nunique() < 2:
        print("Composition map needs both Mg/Si and dIW to vary -- skipping.")
        return
    if quantity not in ('T', 'P_CO2', 'pH'):
        raise ValueError(f'unsupported quantity {quantity!r}')
    if quantity == 'pH':
        subset = _add_diag_columns(subset, output_path)

    live = subset[~subset['termination'].isin(OUT_OF_DOMAIN)]
    # Count DISTINCT cells along the other axis, not runs. Counting runs keeps a cross-design
    # extra like Mg/Si = 1.6 (present at one dIW but all 19 instellations), which then draws as
    # a near-empty stripe through the middle of the grid.
    mg_vals = [m for m in sorted(subset['mg_si'].unique())
               if live[live['mg_si'] == m]['delta_iw'].nunique() >= min_cells]
    diw_vals = [d for d in sorted(subset['delta_iw'].unique())
                if live[live['delta_iw'] == d]['mg_si'].nunique() >= min_cells]
    dropped = (sorted(set(subset['mg_si'].unique()) - set(mg_vals)),
               sorted(set(subset['delta_iw'].unique()) - set(diw_vals)))
    if len(mg_vals) < 2 or len(diw_vals) < 2:
        print("Composition map: too few populated rows/columns -- skipping.")
        return

    s_present = [s for s in s_vals if s in set(subset['instellation'])]
    if not s_present:
        avail = sorted(subset['instellation'].unique())
        s_present = avail[::max(1, len(avail) // 4)][:4]

    log_q = (quantity == 'P_CO2')

    def _cell(S, mg, dw):
        r = live[(live['instellation'] == S) & (live['mg_si'] == mg) & (live['delta_iw'] == dw)]
        if r.empty:
            return np.nan
        v = float(r[quantity].iloc[0])
        return np.log10(v) if log_q and v > 0 else (np.nan if log_q else v)

    grids = {}
    for S in s_present:
        g = np.full((len(diw_vals), len(mg_vals)), np.nan)
        for i, dw in enumerate(diw_vals):
            for j, mg in enumerate(mg_vals):
                g[i, j] = _cell(S, mg, dw)
        if relative:
            ref = _cell(S, REF_MG_SI, REF_DIW)
            g = g - ref if np.isfinite(ref) else g * np.nan
        grids[S] = g

    allv = np.concatenate([g[np.isfinite(g)].ravel() for g in grids.values()]) \
        if any(np.isfinite(g).any() for g in grids.values()) else np.array([])
    if allv.size == 0:
        print("  no in-domain runs (or no reference cell) -- skipping map.")
        return

    unit = {'T': 'K', 'P_CO2': 'dex', 'pH': ''}[quantity]
    qname = {'T': 'Temperature', 'P_CO2': '$P_{\\mathrm{CO_2}}$', 'pH': 'Ocean pH'}[quantity]
    # Short label: the long form clips against the figure edge even under bbox_inches='tight'.
    short = {'T': r'$\Delta T$', 'P_CO2': r'$\Delta\log P_{\mathrm{CO_2}}$',
             'pH': r'$\Delta$pH'}[quantity]
    if relative:
        lim = float(np.nanpercentile(np.abs(allv), 98)) or 1.0
        norm = mcolors.TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
        cmap = cmr.fusion_r
        cbar_label = f'{short} vs Earth crust' + (f' ({unit})' if unit else '')
    else:
        norm = mcolors.Normalize(vmin=np.nanpercentile(allv, 2), vmax=np.nanpercentile(allv, 98))
        cmap = cmr.ember if log_q else cmr.fusion_r
        cbar_label = qname + (f' ({unit})' if unit else '')

    print(f"Composition map [{quantity}, {'relative' if relative else 'absolute'}] "
          f"depth={ocean_depth:g}: {len(mg_vals)}x{len(diw_vals)} grid at S={s_present}"
          + (f"; dropped sparse Mg/Si {dropped[0]} dIW {dropped[1]}" if any(dropped) else ""))

    n = len(s_present)
    fig, axs = plt.subplots(1, n, figsize=figure_size(width, height),
                            sharey=True, squeeze=False)
    axs = axs[0]
    for ax, S in zip(axs, s_present):
        ax.pcolormesh(np.arange(len(mg_vals) + 1), np.arange(len(diw_vals) + 1), grids[S],
                      cmap=cmap, norm=norm, edgecolors='w', linewidth=0.4)
        cell = subset[subset['instellation'] == S]
        for _, r in cell.iterrows():
            if r['termination'] in OUT_OF_DOMAIN and r['mg_si'] in mg_vals \
                    and r['delta_iw'] in diw_vals:
                ax.plot(mg_vals.index(r['mg_si']) + 0.5, diw_vals.index(r['delta_iw']) + 0.5,
                        marker='x', color='0.4', markersize=5, mew=1.2)
        if relative:  # mark the reference cell the anomaly is measured against
            if REF_MG_SI in mg_vals and REF_DIW in diw_vals:
                ax.plot(mg_vals.index(REF_MG_SI) + 0.5, diw_vals.index(REF_DIW) + 0.5,
                        marker='o', mfc='none', mec='k', markersize=9, mew=1.4)
        ax.set_xticks(np.arange(len(mg_vals)) + 0.5)
        ax.set_xticklabels([f'{v:g}' for v in mg_vals], rotation=90)
        ax.set_yticks(np.arange(len(diw_vals)) + 0.5)
        ax.set_yticklabels([f'{v:+g}' for v in diw_vals])
        ax.set_xlabel('Mantle Mg/Si')
        ax.set_title(f'S = {S:g}')
    axs[0].set_ylabel(r'Core-formation $\Delta$IW')
    _add_colorbar(fig, list(axs), cmap, norm, cbar_label, aspect=18)
    handles = [Line2D([0], [0], marker='x', color='0.4', linestyle='none', markersize=5,
                      label='Outside model domain')]
    if relative:
        handles.append(Line2D([0], [0], marker='o', mfc='none', mec='k', linestyle='none',
                              markersize=8, label='Earth reference crust'))
    fig.legend(handles=handles, loc='outside lower center', ncol=len(handles))
    _save_fig(fig, os.path.join(output_path,
                                f'map_composition_{quantity}{"" if relative else "_abs"}.png'))


def plot_ratio_scatter(df, output_path, s_vals=(0.4, 0.6, 0.8, 1.0, 1.2)):
    """T and P_CO2 vs outgassing/crust-production ratio, coloured by instellation."""
    base = _base(df)
    sub = base[base['instellation'].isin(s_vals)].copy()
    if sub.empty:
        print("No runs found — skipping ratio scatter.")
        return

    sub['ratio'] = sub['outgassing'] / sub['crust_production']
    hab     = sub[sub['termination'].isin(HABITABLE)]
    non_hab = sub[~sub['termination'].isin(HABITABLE)]

    norm = mcolors.Normalize(vmin=min(s_vals), vmax=max(s_vals))
    cmap = cmr.cosmic

    figsize = figure_size('single', height=3.0)
    fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)

    for term, marker in HAB_MARKERS.items():
        grp = hab[hab['termination'] == term]
        if grp.empty:
            continue
        kw = dict(c=grp['instellation'], cmap=cmap, norm=norm,
                  s=22, alpha=0.85, zorder=4, linewidths=0)
        axes[0].scatter(grp['ratio'], grp['T'],     marker=marker, **kw)
        axes[1].scatter(grp['ratio'], grp['P_CO2'], marker=marker, **kw)

    for term, marker in FAILED_MARKERS.items():
        grp = non_hab[non_hab['termination'] == term]
        if grp.empty:
            continue
        for ax, col in zip(axes, ['T', 'P_CO2']):
            valid = grp[np.isfinite(grp[col])]
            if valid.empty:
                continue
            ax.scatter(valid['ratio'], valid[col],
                       marker=marker, s=28, zorder=3, linewidths=0.8,
                       facecolors='none',
                       edgecolors=cmap(norm(valid['instellation'].values)))

    axes[0].set_ylabel('Temperature (K)')
    axes[0].axhspan(T_SNOWBALL - 25, T_SNOWBALL, color='blue', alpha=0.12)
    axes[0].axhspan(T_RUNAWAY - 20,  T_RUNAWAY,  color='red',  alpha=0.12)
    axes[0].set_ylim(235, 360)

    axes[1].set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-8, 20)
    axes[1].set_xlabel('Outgassing / Crust production rate')
    axes[1].set_xlim([1e-3, 1e3])

    for ax in axes:
        ax.set_xscale('log')
        ax.grid(True, linestyle='--', alpha=0.4)

    _add_colorbar(fig, list(axes), cmap, norm, 'Instellation (S/S₀)', ticks=sorted(s_vals))

    marker_handles = [
        plt.scatter([], [], marker=m, s=22, color='k', label=TERM_LABELS[t])
        for t, m in HAB_MARKERS.items()
    ] + [
        plt.scatter([], [], marker=m, s=28, facecolors='none', edgecolors='k',
                    linewidths=0.8, label=TERM_LABELS[t])
        for t, m in FAILED_MARKERS.items()
    ]
    fig.legend(handles=marker_handles, loc='outside lower center', ncol=_legend_ncol(marker_handles, 2))
    _save_fig(fig, os.path.join(output_path, 'ratio_scatter.png'))


# ---------------------------------------------------------------------------
# Mineral SI plot
# ---------------------------------------------------------------------------

# All minerals that can precipitate, in display order. These mirror the sets used by
# Planet: clays in the pore space; carbonates, clays, silica, evaporites (+ reverse
# weathering clays when enabled) in the ocean.
_PORE_MINERALS  = list(clay_minerals)
_OCEAN_MINERALS = (carbonate_minerals + clay_minerals + silica_minerals +
                   evaporite_minerals + list(reverse_weathering_minerals))
_ALL_MINERALS   = list(dict.fromkeys(_PORE_MINERALS + _OCEAN_MINERALS))  # ordered, unique

_MINERAL_LABELS = {
    'Calcite':      'Calcite',
    'Siderite':     'Siderite (FeCO₃)',
    'Nahcolite':    'Nahcolite (NaHCO₃)',
    'Kaolinite':    'Kaolinite',
    'Goethite':     'Goethite',
    'SiO2(am)':     'Amorphous SiO₂',
    'Halite':       'Halite (NaCl)',
    'Sepiolite(d)': 'Sepiolite (Mg)',
    'Saponite-Na':  'Saponite-Na',
    'Greenalite':   'Greenalite (Fe)',
}


def plot_damkohler_contour(df, output_path, out_targets=(0.1, 1.0, 10.0)):
    """2D contourf of Damköhler number in (instellation × crust production) space.

    One panel per outgassing rate (vertically stacked), coloured by log10(Da).
    The Da = 1 boundary is drawn as a black contour.
    Uses the reference-crust, rw=True, depth=3000, f_HT=0 baseline.
    """
    subset = df[
        _ref_crust(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0)
    ]
    if subset.empty:
        print("No data for Da contour plot — skipping.")
        return

    all_out = sorted(subset['outgassing'].unique())
    out_values = list(dict.fromkeys(
        min(all_out, key=lambda x: abs(x - t)) for t in out_targets
    ))

    s_vals    = sorted(subset['instellation'].unique())
    crust_vals = sorted(subset['crust_production'].unique())

    sel = subset[subset['outgassing'].isin(out_values)]
    sel = _add_diag_columns(sel, output_path)

    s_arr     = np.array(s_vals)
    log_crust = np.log10(np.array(crust_vals))

    vmin, vmax = -3.0, 3.0
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cmap = cmr.prinsenvlag
    levels = np.linspace(vmin, vmax, 31)

    nrows = len(out_values)

    with plt.rc_context({'figure.constrained_layout.use': False}):
        fig, axes = plt.subplots(nrows, 1, figsize=figure_size('single', n_rows=nrows, row_height=2.0),
                                 sharex=True, sharey=True, squeeze=False)

        for idx, out in enumerate(out_values):
            ax = axes[idx, 0]
            sub = sel[np.isclose(sel['outgassing'], out)]

            pivot = sub.pivot_table(
                index='crust_production', columns='instellation',
                values='da', aggfunc='first',
            ).reindex(index=crust_vals, columns=s_vals)

            Z = np.log10(np.maximum(pivot.values.astype(float), 1e-10))
            Z = np.ma.masked_invalid(Z)

            ax.contourf(s_arr, log_crust, Z, levels=levels, cmap=cmap, norm=norm,
                        extend='both')
            if not np.all(Z.mask if np.ma.is_masked(Z) else False):
                ax.contour(s_arr, log_crust, Z, levels=[0.0],
                           colors='k', linewidths=1.5)

            ax.set_title(f'Outgassing = {out:g}×', fontsize=9)
            ax.set_yticks(log_crust)
            ax.set_yticklabels([f'{v:g}' for v in crust_vals], fontsize=7)
            ax.set_ylabel('Crust prod. (×Earth)')
            ax.grid(True, linestyle='--', alpha=0.3, color='k')

        axes[-1, 0].set_xlabel('Instellation (S/S₀)')

        fig.subplots_adjust(left=0.16, right=0.82, top=0.91, bottom=0.12, hspace=0.25)

        pos_top = axes[0, 0].get_position()
        pos_bot = axes[-1, 0].get_position()
        cbar_ax = fig.add_axes([0.85, pos_bot.y0, 0.03, pos_top.y1 - pos_bot.y0]) # type: ignore
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax, label=r'$\log_{10}(\mathrm{Damkohler Coefficient})$')
        cbar.ax.axhline(0, color='k', linewidth=1.5)

        # fig.suptitle('Damköhler number — basalt_49, rw=True, depth = 3000 m', fontsize=10)
        _save_fig(fig, os.path.join(output_path, 'da_contour.png'))


def plot_continental_baseline(df, output_path):
    """T, P_CO2, pH, salinity, and individual ion concentrations vs instellation
    for the Earth-like continental baseline.

    Runs with land_fraction=0.3, reference crust, rw=True, out=1×, crust=1×, depth=3000 m.
    """
    subset = df[
        (df['land_fraction'] == 0.3) &
        _ref_crust(df) &
        df['reverse_weathering'] &
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0) &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0)
    ]
    if subset.empty:
        print("No continental baseline data found — skipping.")
        return

    cols = ['T', 'P_CO2', 'pH', 'salinity']
    group_all = subset.sort_values('instellation')
    group_hab = group_all[
        (group_all['T'] > T_SNOWBALL) &
        (group_all['T'] < T_RUNAWAY)
    ]

    EARTH_S    = 1.0
    EARTH_T    = 288.0
    EARTH_PCO2 = 280e-6
    EARTH_PH   = 8.1
    EARTH_SAL  = (2.0e-3 * 61.0 + 0.1e-3 * 60.1 +
                  10.3e-3 * 40.1 + 52.8e-3 * 24.3 +
                  480e-3 * 23.0 + 550e-3 * 35.45)

    earth_vals = {'T': EARTH_T, 'P_CO2': EARTH_PCO2, 'pH': EARTH_PH, 'salinity': EARTH_SAL}

    # --- ion panel setup ---
    # (index, label, Earth mmol/kg or None) — Al(3), Fe(4), SO₄(9) excluded
    ION_SPEC = [
        (0, 'Alk',  2.3),
        (1, 'C',  2.0),
        (2, 'Si',   0.1),
        (5, 'Ca',  10.3),
        (6, 'Mg',  52.8),
        (7, 'Na', 480.0),
        (8, 'Cl', 550.0),
    ]
    ION_COLORS = [plt.cm.tab10(k / 10) for k in range(len(ION_SPEC))] # type: ignore

    # Load final b_ocean for each run from the JSON files
    ion_rows = []
    for _, row in group_hab.iterrows():
        fpath = os.path.join(output_path, f"{row['name']}.json")
        try:
            with open(fpath) as fh:
                d_json = json.load(fh)
            y = d_json['data']['y']
            n_el = len(y) - 3
            b = [max(float(y[2 + i][-1]), 1e-15) for i in range(n_el)]
            # Pad to 10 elements if an older file lacks some ions
            while len(b) < 10:
                b.append(1e-15)
            ion_rows.append((row['instellation'], b[:10]))
        except Exception:
            pass

    # --- figure 1: 4 summary panels ---
    if presentation:
        fig, _ax2d = plt.subplots(2, 2, figsize=figure_size('double', height=4.5),
                                   sharex=True)
        axes_summary = _ax2d.ravel()
    else:
        fig, axes_summary = plt.subplots(
            len(cols), 1, sharex=True,
            figsize=figure_size('single', n_rows=len(cols), row_height=2.0))

    # T panel shows all runs (including failed); other panels show only runs in habitable T range
    _plot_group_on_axes(axes_summary[:1], group_all, color='k', show_markers=False, cols=['T'])
    _plot_group_on_axes(axes_summary[1:], group_hab, color='k', show_markers=False, cols=['P_CO2', 'pH', 'salinity'])
    _style_axes(axes_summary, cols)
    if presentation:
        axes_summary[2].set_xlabel('Instellation (S/S₀)')  # bottom-left panel also needs x label

    for ax, col in zip(axes_summary, cols):
        ax.scatter(EARTH_S, earth_vals[col], marker='*', s=220, color='blue',
                   edgecolors='k', linewidths=0.7, zorder=6)
    axes_summary[0].annotate(
        'Earth', xy=(EARTH_S, EARTH_T), xytext=(EARTH_S + 0.06, EARTH_T - 6),
        fontsize=8, arrowprops=dict(arrowstyle='-', color='k', lw=0.8),
    )
    _save_fig(fig, os.path.join(output_path, 'continental_baseline.png'))

    # --- figure 2: ion ratio bar chart (model vs Earth seawater at S ≈ 1) ---
    figsize2 = figure_size('single', height=2.25)
    fig2, ax_ions = plt.subplots(1, 1, figsize=figsize2)

    if ion_rows:
        s_arr  = np.array([r[0] for r in ion_rows])
        b_mmol = np.array([r[1] for r in ion_rows]) * 1e3  # mol/kg → mmol/kg
        closest = int(np.argmin(np.abs(s_arr - EARTH_S)))
        b_model = b_mmol[closest]

        labels    = [spec[1] for spec in ION_SPEC]
        ratios    = [100 * (b_model[spec[0]] - spec[2]) / spec[2] for spec in ION_SPEC]
        x         = np.arange(len(labels))

        ax_ions.scatter(x, ratios, color=ION_COLORS, edgecolors='k', linewidths=0.6, s=80, zorder=3, alpha=0.6)
        ax_ions.set_xticks(x)
        ax_ions.set_xticklabels(labels)
        ax_ions.axvline(3.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1)
        trans = ax_ions.get_xaxis_transform()
        ax_ions.text(1.5, 1.02, 'Biotically controlled', transform=trans,
                     ha='center', va='bottom', fontsize=8)
        ax_ions.text(5.0, 1.02, 'Abiotically controlled', transform=trans,
                     ha='center', va='bottom', fontsize=8)

    ax_ions.spines['bottom'].set_position(('data', 0))
    ax_ions.spines['top'].set_visible(False)
    ax_ions.spines['right'].set_visible(False)
    ax_ions.set_ylabel('Model Difference (%)')
    # ax_ions.set_yscale('symlog', linthresh=0.1)
    ax_ions.grid(True, linestyle='--', alpha=0.4, axis='y')
    _save_fig(fig2, os.path.join(output_path, 'continental_baseline_ions.png'))


def _get_mineral_si(d):
    """Return {'pore': {mineral: SI}, 'ocean': {mineral: SI}, 'da': float, 'T': float}
    for all precipitating minerals in pore space and ocean."""
    from kamino.precipitation import get_precipitation
    from kamino.chemistry import ChemistryError
    nan_result = {'pore': {}, 'ocean': {}, 'da': np.nan, 'T': np.nan}
    try:
        y_list = d.get('data', {}).get('y', [])
        if not y_list or len(y_list) - 3 < 7:
            return nan_result

        b_ocean, P_pore, T_pore, T_seafloor, P_CO2, crust_rate, J_total = _pore_conditions(d)
        T_surface = float(d['T'])

        rw        = bool(d.get('reverse_weathering', False))
        pore_min  = list(clay_minerals)  # planet.pore_precipitating_minerals
        ocean_min = (carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals +
                     (list(reverse_weathering_minerals) if rw else []))

        _, diag = get_weathering_flux(
            P_pore, T_pore, P_CO2, b_ocean,
            alpha=float(d.get('alpha', 1.43)),
            rate=crust_rate, J=J_total,
            crust_composition=_crust_composition_of(d),
            sedimentation_rate=_sedimentation_rate(d, b_ocean, P_pore, T_seafloor),
            precipitating_minerals=pore_min
        )
        pore_si = diag.get('secondary_SI', {})

        try:
            _, _, ocean_si = get_precipitation(P_pore, T_seafloor, b_ocean, ocean_min,
                                               precipitation_timescale=1e6 * YR)
        except (ChemistryError, Exception):
            ocean_si = {}

        return {'pore': pore_si, 'ocean': ocean_si, 'da': float(diag['Da']), 'T': T_surface}
    except Exception:
        return nan_result


def plot_mineral_si(df, output_path):
    """SI of all precipitating minerals vs instellation, split into pore-space and ocean columns.

    One figure per crust production rate, lines coloured by outgassing rate.
    """
    base = _base(df)
    if base.empty:
        print("No base data for mineral SI plot — skipping.")
        return

    crust_rates     = sorted(base['crust_production'].unique())
    outgassing_vals = sorted(base['outgassing'].unique())
    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical

    pore_mineral_set  = set(_PORE_MINERALS)
    ocean_mineral_set = set(_OCEAN_MINERALS)

    n_min  = len(_ALL_MINERALS)
    panels = [('Pore space', 'pore', pore_mineral_set),
              ('Ocean',      'ocean', ocean_mineral_set)]

    for c in crust_rates:
        subset_c = base[base['crust_production'] == c]

        # Load mineral SI data for every run in this crust-rate slice
        si_by_name = {}
        for _, row in subset_c.iterrows():
            fpath = os.path.join(output_path, f'{row["name"]}.json')
            try:
                with open(fpath) as fh:
                    d = json.load(fh)
                rec = _get_mineral_si(d)
                si_by_name[row['name']] = {
                    'pore':         rec['pore'],
                    'ocean':        rec['ocean'],
                    'da':           rec['da'],
                    'T':            rec['T'],
                    'instellation': row['instellation'],
                    'outgassing':   row['outgassing'],
                }
            except Exception:
                pass

        # Diagnostic sheet (one per crust rate), read on screen -- not page-sized.
        fig, axes = plt.subplots(n_min, 2, figsize=diagnostic_size(n_min, 1, col_width=3.5,
                                                                   row_height=1.0, pad=0.0),
                                  sharex=True, squeeze=False)

        for o in outgassing_vals:
            color = cmap(norm(o))
            runs  = sorted(
                [v for v in si_by_name.values() if v['outgassing'] == o],
                key=lambda v: v['instellation'],
            )
            if len(runs) < 2:
                continue
            s_arr  = np.array([r['instellation'] for r in runs])
            da_arr = np.array([r['da']           for r in runs])
            T_arr  = np.array([r['T']            for r in runs])
            at_floor = 1.02 * T_arr - 16.7 <= 274.001

            for mi, mineral in enumerate(_ALL_MINERALS):
                for ci, (_, loc, min_set) in enumerate(panels):
                    ax = axes[mi, ci]
                    if mineral not in min_set:
                        continue
                    si_arr = np.array([r[loc].get(mineral, np.nan) for r in runs])
                    if np.all(np.isnan(si_arr)):
                        continue
                    _plot_line_da_style(ax, s_arr, si_arr, da_arr, color, at_floor=at_floor)

        # Style axes
        for mi, mineral in enumerate(_ALL_MINERALS):
            label = _MINERAL_LABELS.get(mineral, mineral)
            for ci, (col_title, loc, min_set) in enumerate(panels):
                ax = axes[mi, ci]
                ax.axhline(0, color='k', linewidth=0.6, linestyle='--', alpha=0.5)
                ax.set_xlim(0.25, 1.45)
                ax.grid(True, linestyle='--', alpha=0.4)
                if mineral in min_set:
                    ax.set_ylabel('SI')
                    ax.set_title(f'{label} — {col_title}', fontsize=7, pad=2)
                else:
                    ax.set_visible(False)

        axes[-1, 0].set_xlabel('Instellation (S/S₀)')
        axes[-1, 1].set_xlabel('Instellation (S/S₀)')

        _add_colorbar(fig, list(axes.ravel()), cmap, norm, 'Earth Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_min * 8)
        fig.legend(handles=DA_LEGEND, loc='outside lower center', ncol=_legend_ncol(DA_LEGEND, 3))
        fig.suptitle(f'Mineral saturation indices — Crust = {c}× Earth')
        _save_fig(fig, os.path.join(output_path, f'mineral_si_crust_{c}.png'))

def plot_habitability_phase_space(df, output_path):
    """
    Phase diagram mapping Instellation vs. Outgassing/Crust production ratio.
    Aggregates overlapping data points to determine the dominant macro-state
    at each coordinate, completely eliminating visual clutter.
    """
    base = _base(df).copy()
    if base.empty:
        print("No runs found — skipping phase space plot.")
        return

    base['ratio'] = base['outgassing'] / base['crust_production']

    # 1. Define Boolean masks to strictly enforce state conditions.
    #
    # Classification is driven by the FINAL STATE, not by which event stopped the run:
    # terminations no longer encode outcomes (see Planet.time_evolve). Final T decides
    # snowball vs hothouse vs habitable; only a run stopped at a CO2 wall has a genuinely
    # unknown fate, because there the model quit while the climate was still evolving.
    wall = base.get('domain_wall')
    if wall is None:
        wall = pd.Series([None] * len(base), index=base.index)

    # A run stopped at a CO2 wall quit while the climate was still evolving, so its fate is
    # genuinely unknown. One stopped at a temperature wall has already reached its outcome.
    cond_unknown = base['termination'].isin(OUT_OF_DOMAIN) & wall.isin(['co2_high', 'co2_low'])

    cond_snow = ((base['T'] <= T_SNOWBALL) | (wall == 'cold')) & ~cond_unknown
    cond_hot  = ((base['T'] >= T_RUNAWAY) | (wall == 'hot')) & ~cond_unknown & ~cond_snow
    cond_hab  = ~(cond_unknown | cond_snow | cond_hot)
    cond_acid = cond_unknown  # name kept: downstream code below still refers to it

    # 2. Assign a string label to each row based on its state
    base['macro_state'] = None
    base.loc[cond_snow, 'macro_state'] = 'Snowball'
    base.loc[cond_hab,  'macro_state'] = 'Habitable'
    base.loc[cond_hot,  'macro_state'] = 'Hothouse'
    base.loc[cond_acid, 'macro_state'] = 'Unknown'

    # Drop any edge cases that missed categorization
    base = base.dropna(subset=['macro_state'])

    # 3. AGGREGATE: Group by grid coordinate and find the most frequent state
    # This guarantees exactly ONE data point per (X, Y) coordinate
    def get_dominant_state(s):
        return s.value_counts().idxmax()
    
    agg_df = base.groupby(['ratio', 'instellation'])['macro_state'].agg(get_dominant_state).reset_index()

    # 4. Map these conditions to our visual styling
    macro_states = {
        'Snowball': {
            'name': 'Snowball',
            'color': '#4ea8ff', # Icy Blue
            'marker': 'v'
        },
        'Habitable': {
            'name': 'Habitable',
            'color': '#51c46f', # Lush Green
            'marker': 'o'
        },
        'Hothouse': {
            'name': 'Hothouse',
            'color': '#ff5e5e', # Hot Red
            'marker': '^'
        },
        'Unknown (stopped at CO₂ wall)': {
            'name': 'Unknown',
            'color': '#f39c12', # Acidic Yellow-Orange
            'marker': 's'
        }
    }

    # Use existing sizing conventions
    figsize = figure_size('single', height=3.0)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    legend_handles = []

    for label, config in macro_states.items():
        # Filter the AGGREGATED dataframe, not the base dataframe
        grp = agg_df[agg_df['macro_state'] == config['name']]
        if grp.empty:
            continue

        # Scatter the points for this state
        # Markers are slightly larger now that they aren't overlapping
        marker_size = 80 if presentation else 50
        
        ax.scatter(grp['ratio'], grp['instellation'],
                   c=config['color'], marker=config['marker'],
                   s=marker_size, edgecolors='k', linewidths=0.6, alpha=0.9, zorder=4)

        # Build custom legend handle
        legend_handles.append(Line2D([0], [0], marker=config['marker'], color='w',
                                     markerfacecolor=config['color'], markeredgecolor='k',
                                     markersize=9, label=label))

    # Format the axes
    ax.set_xscale('log')
    ax.set_xlim([1e-3, 1e3]) # type: ignore
    
    # Add a small margin to the Y-axis based on the data
    y_min = agg_df['instellation'].min()
    y_max = agg_df['instellation'].max()
    margin = (y_max - y_min) * 0.05 if y_max != y_min else 0.1
    ax.set_ylim([y_min - margin, y_max + margin]) # type: ignore

    ax.set_xlabel('Outgassing / Crust production rate')
    ax.set_ylabel('Instellation (S/S₀)')
    ax.grid(True, linestyle='--', alpha=0.4, zorder=0)

    # Add the legend below the plot
    fig.legend(handles=legend_handles, loc='outside lower center', ncol=2)

    _save_fig(fig, os.path.join(output_path, 'ratio_phase_space.png'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot kamino parameter sweep results.')
    parser.add_argument('--path', default=DEFAULT_OUTPUT_PATH,
                        help='Directory containing planet_*.json files.')
    for _knob in CHEM_KNOBS:
        parser.add_argument(f'--{_knob.replace("_", "-")}', type=float, default=None,
                            help=f'Pin {_knob} to this value in the main plots '
                                 f'(default: the most-run value).')
    parser.add_argument('--legacy', action='store_true',
                        help='Input predates the current run schema: upgrade it first '
                             '(see plot_legacy.py).')
    parser.add_argument('--depth', type=float, default=None,
                        help='Ocean depth (m) for the composition figures '
                             '(default: every depth that has a composition sweep).')
    args = parser.parse_args()

    for _knob in CHEM_KNOBS:
        _v = getattr(args, _knob)
        if _v is not None:
            CHEM_OVERRIDE[_knob] = _v

    df = load_data(args.path)

    if args.legacy:
        import plot_legacy
        df = plot_legacy.upgrade(df)
        plot_legacy.plot_named_compositions(df, args.path, split_panels=presentation)

    if df.empty:
        print("No data found. Check --path.")
        raise SystemExit(1)

    # Depths that actually carry a composition sweep. The composition figures were previously
    # hardcoded to 3000 m, so a deep composition sweep produced no plot at all.
    comp_depths = sorted(
        d for d in df['ocean_depth'].unique()
        if df[(df['ocean_depth'] == d)][['mg_si', 'delta_iw']].nunique().max() > 1
    ) or [3000.0]
    if args.depth is not None:
        comp_depths = [args.depth]
    print(f"Composition figures for depth(s): {[f'{d:g}' for d in comp_depths]}")

    plot_faceted_lines(df, args.path)
    plot_faceted_lines(df, args.path, all_results=False, split_panels=True)
    plot_faceted_lines(df, args.path, split_panels=True, all_results=False, sequence=True)
    plot_ocean_depth_effect(df, args.path, split_panels=presentation)
    plot_chemistry_constants(df, args.path, split_panels=presentation)
    plot_ratio_scatter(df, args.path)
    for _d in comp_depths:
        plot_crust_composition(df, args.path, show_markers=False, split_panels=presentation,
                               ocean_depth=_d)
        plot_composition_grid(df, args.path, show_markers=False, split_panels=presentation,
                              ocean_depth=_d)
    for _q in ('T', 'P_CO2'):
        plot_composition_map(df, args.path, ocean_depth=comp_depths[0], quantity=_q)
    plot_damkohler_contour(df, args.path)
    plot_habitability_phase_space(df, args.path)
    plot_continental_baseline(df, args.path)
    print("Done.")
