"""Figures for kamino parameter sweeps: python experiments/plot_results.py --path <sweep dir>."""
import os
import re
import sys
import glob
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import cmasher as cmr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
from kamino.constants import SOLAR_CONSTANT, STEFAN_BOLTZMANN, EARTH_MANTLE_MG_SI, EARTH_DELTA_IW, SEAFLOOR_T_FLOOR
from kamino.chemistry import elements
from kamino.planet import KD_MG_HT, K_NA_CONT_REMOVAL, PE_DEFAULT
from kamino.weathering import ALPHA_REF
from kamino.planet import OCEAN_ALBEDO
from kamino.climate.analytic import get_T_surface_analytic

import continental_baseline as cb

# --- Style and figure geometry ---
# KAMINO_PRESENTATION=1 switches to the presentation style, scaling every figure 2x at the paper aspect.
presentation = os.environ.get('KAMINO_PRESENTATION', '0').lower() in ('1', 'true', 'yes')
STYLE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          f'planetary-chem-{"presentation" if presentation else "paper"}.mplstyle')
if not os.path.exists(STYLE_FILE):
    raise SystemExit(f"missing style file {STYLE_FILE}")
plt.style.use(STYLE_FILE)

PAGE_WIDTHS = {'single': 240 / 72.27, 'double': 504 / 72.27}   # MNRAS column and text width (in)
ROW_HEIGHT_IN = 1.5
_PRES_SCALE = 2.0 if presentation else 1.0


def figure_size(width='single', height=None, n_rows=1, row_height=ROW_HEIGHT_IN):
    """(width, height) in inches for a page figure; height defaults to n_rows * row_height."""
    if width not in PAGE_WIDTHS:
        raise ValueError(f"width must be one of {sorted(PAGE_WIDTHS)}, not {width!r}")
    h = height if height is not None else n_rows * row_height
    return (PAGE_WIDTHS[width] * _PRES_SCALE, h * _PRES_SCALE)


def diagnostic_size(n_rows, n_cols, col_width=6.0, row_height=3.0, pad=1.5):
    """Fixed per-panel size for on-screen diagnostic grids that are not page-constrained."""
    return (col_width * n_cols + pad, row_height * n_rows)


# --- Colour maps: one per swept parameter, so a parameter looks the same in every figure ---
INSTELLATION_CMAP     = cmr.cosmic
OUTGASSING_CMAP       = cmr.get_sub_cmap('cmr.ember', 0.15, 0.90)
CRUST_CMAP            = cmr.get_sub_cmap('cmr.amber', 0.15, 0.90)
DEPTH_CMAP            = cmr.bubblegum_r
MG_SI_CMAP            = cmr.gem_r
DIW_CMAP              = cmr.get_sub_cmap('cmr.emerald_r', 0.25, 0.95)
PE_CMAP               = cmr.lavender
LAND_FRACTION_CMAP    = cmr.get_sub_cmap('cmr.savanna_r', 0.15, 0.85)
CHEM_KNOB_CMAP        = cmr.amber

PARAM_CMAPS = {
    'instellation':     INSTELLATION_CMAP,
    'outgassing':       OUTGASSING_CMAP,
    'crust_production': CRUST_CMAP,
    'ocean_depth':      DEPTH_CMAP,
    'mg_si':            MG_SI_CMAP,
    'delta_iw':         DIW_CMAP,
    'pe':               PE_CMAP,
    'land_fraction':    LAND_FRACTION_CMAP,
    'alpha':            CHEM_KNOB_CMAP,
    'kd_mg':            CHEM_KNOB_CMAP,
    'k_na':             CHEM_KNOB_CMAP,
}

DIVERGING_CMAP        = cmr.prinsenvlag   # log10(Da), centred on Da = 1
RELATIVE_CMAP         = cmr.fusion_r      # differences against the reference crust
WEATHERING_RATIO_CMAP = cmr.fusion        # log seafloor/continental: blue where seafloor dominates
QUANTITY_CMAP         = cmr.ember         # absolute value of a log-scaled quantity
TEMPERATURE_CMAP      = cmr.get_sub_cmap('cmr.ember', 0.12, 1.0)   # surface temperature contours

DEFAULT_OUTPUT_PATH = os.environ.get('KAMINO_SWEEP_OUTPUT', '/home/pt426/Code/kamino/sweep_output')

# --- Terminations and thresholds ---
TERM_LABELS = {
    'converged':      'Converged',
    'timeout':        'Timeout (2 Gyr)',
    'wall_timeout':   'Wall-clock cap',
    'out_of_domain':  'Outside model domain',
    'fallback_limit': 'Chemistry fallback cap',
    'chemistry_void': 'Chemistry solver void',
    'solver_failure': 'ODE solver failure',
}

# Where an out_of_domain run left the validity box.
WALL_LABELS = {
    'cold':     'Frozen (T → 181 K)',
    'hot':      'Runaway (T → 389 K)',
    'co2_high': 'CO₂ ceiling (10 bar)',
    'co2_low':  'CO₂ depleted (0.1 Pa)',
}

# Failed or incomplete integrations (wall_timeout, fallback_limit, ...) are in neither set.
OUT_OF_DOMAIN = {'out_of_domain'}
HABITABLE = {'converged', 'timeout'}
# A wall_timeout state was genuinely evaluated, so its Da can anchor a transition marker.
DA_TRUSTWORTHY = HABITABLE | {'wall_timeout'}

HAB_MARKERS    = {'converged': 'o', 'timeout': 's'}
FAILED_MARKERS = {                      # hollow, one per unreliable termination
    'out_of_domain':  's',
    'wall_timeout':   '^',
    'fallback_limit': 'D',
    'chemistry_void': 'v',
    'solver_failure': 'x',
}

T_SNOWBALL = 260.0
T_RUNAWAY  = 360.0
# The climate solver's validity box; a T pinned to either value is a sentinel, not a solution.
T_COLD_WALL = 181.0
T_HOT_WALL  = 389.0

DA_LEGEND = [
    Line2D([0], [0], color='k', linestyle='-',  linewidth=1.4, label='Da < 1 (kinetic)'),
    Line2D([0], [0], color='k', linestyle='--', linewidth=1.4, label='Da ≥ 1 (thermodynamic)'),
    Line2D([0], [0], color='k', linestyle='none', marker='o', markerfacecolor='none',
           markersize=5.5, markeredgewidth=1.2, label='Da = 1 transition'),
    Line2D([0], [0], color='k', linestyle=':',  linewidth=1.4, label=f'$T_\\mathrm{{seafloor}}$ at floor ({SEAFLOOR_T_FLOOR:g} K)'),
]

PANEL_COLS = ['T', 'P_CO2', 'pH', 'salinity']

# --- Reference planet ---
REF_MG_SI = float(EARTH_MANTLE_MG_SI)
REF_DIW   = float(EARTH_DELTA_IW)
REF_PE    = float(PE_DEFAULT)   # reducing ocean; every figure except plot_pe pins pe here (--pe overrides)
PE_LEGACY_DEFAULT = 4.0         # pe of runs written before pe was a parameter (PHREEQC's oxidising default)

# Composition-sweep Mg/Si values left off the Mg/Si colour-bar figures; 1.75 crowds the 1.8 end-member.
MG_SI_HIDDEN = (1.75,)

# Modern Earth; salinity sums the tracked seawater ions (g/kg), comparable with the model's column.
EARTH = {'S': 1.0, 'T': 288.0, 'P_CO2': 280e-6, 'pH': 8.1,
         'salinity': (2.0e-3 * 61.0 + 0.1e-3 * 60.1 + 10.3e-3 * 40.1 +
                      52.8e-3 * 24.3 + 480e-3 * 23.0 + 550e-3 * 35.45 +
                      28.2e-3 * 96.06 + 10.2e-3 * 39.10)}

# --- Continental habitable zone ---
# Edges of the Earth-like continental baseline (land 0.3, all else Earth), from continental_baseline.py.
# Outer: T crosses T_SNOWBALL between S = 0.45 and 0.50. Inner: runaway bracketed between S = 1.10 and 1.15.
CONTINENTAL_HZ_OUTER = 0.480
CONTINENTAL_HZ_INNER = 1.125
SHOW_HZ_EDGES = True   # global switch; figures also take show_hz to override it


def _draw_hz_edges(ax, show_hz=None):
    """Draw and label the continental HZ edges; show_hz=None defers to SHOW_HZ_EDGES."""
    if not (SHOW_HZ_EDGES if show_hz is None else show_hz):
        return
    trans = ax.get_xaxis_transform()
    for s, label, offset, ha in ((CONTINENTAL_HZ_OUTER, 'CWHZ outer edge', -3, 'right'),
                                 (CONTINENTAL_HZ_INNER, 'CWHZ inner edge', 3, 'left')):
        ax.axvline(s, color='0.35', linestyle=(0, (6, 3)), linewidth=1.0, alpha=0.85, zorder=1)
        ax.annotate(label, xy=(s, 0.97), xycoords=trans, xytext=(offset, 0),
                    textcoords='offset points', rotation=90, ha=ha, va='top',
                    fontsize=6, color='0.35', zorder=5,
                    bbox=dict(boxstyle='square,pad=0.1', facecolor='white', edgecolor='none',
                              alpha=0.75))


def equilibrium_temperature(instellation, albedo=0.3, greenhouse=0.5):
    """Grey-atmosphere equilibrium temperature (K) for an instellation in units of S0."""
    return (((1 - albedo) * instellation * SOLAR_CONSTANT)
            / (4 * STEFAN_BOLTZMANN * greenhouse)) ** 0.25


# --- Loading and diagnostics ---
# Molar masses (g/mol) for salinity; C as HCO3-, S as SO4 2-. y = [P_CO2, P_H2O, *elements] (older outputs add r_avg).
_ELEMENT_MASSES = {'C': 61.0, 'Si': 60.1, 'Al': 27.0, 'Fe': 55.8, 'Ca': 40.1,
                   'Mg': 24.3, 'Na': 23.0, 'Cl': 35.45, 'S': 96.06, 'K': 39.10}
_SAL_INDICES = [2 + i for i, e in enumerate(elements) if e in _ELEMENT_MASSES]
_SAL_MASSES  = [_ELEMENT_MASSES[e] for e in elements if e in _ELEMENT_MASSES]


def _salinity_from_y(y_list):
    """Salinity (g/kg) of a run's final ocean state."""
    try:
        return sum(
            float(y_list[i][-1]) * mass
            for i, mass in zip(_SAL_INDICES, _SAL_MASSES)
            if len(y_list) > i and len(y_list[i]) > 0
        )
    except Exception:
        return np.nan


# How runs whose pCO2 is still cycling at their end enter the figures (--oscillating):
#   'final' the last state, as before; 'mean' the time average over the tail (T, pCO2, salinity; pH and
#   Da stay final); 'exclude' dropped from every figure.
OSC_MODES = ('final', 'mean', 'exclude')
OSCILLATION_MODE = 'final'
OSC_TAIL = 0.4           # fraction of the run examined
OSC_MIN_TURNS = 3        # turning points in log pCO2 over the tail
OSC_MIN_RANGE = 0.2      # and a ln-range above this (~20 %)


def _oscillation(d, y_list):
    """(oscillating, mean T, mean pCO2 in bar, mean salinity) over the last OSC_TAIL of a run's time."""
    t = np.asarray((d.get('data') or {}).get('time') or [], dtype=float)
    if t.size < 10 or not y_list:
        return False, np.nan, np.nan, np.nan
    tq = np.linspace(t[-1] * (1 - OSC_TAIL), t[-1], 60)
    p = np.maximum(np.interp(tq, t, np.asarray(y_list[0], dtype=float)), 1.0)   # 1 Pa climate floor
    lp = np.log(p)
    dl = np.diff(lp)
    sgn = np.sign(dl[np.abs(dl) > 0.02])
    turns = int(np.sum(sgn[1:] != sgn[:-1])) if sgn.size > 1 else 0
    osc = turns >= OSC_MIN_TURNS and (lp.max() - lp.min()) > OSC_MIN_RANGE
    if not osc:
        return False, np.nan, np.nan, np.nan
    S = float(d['instellation']) * SOLAR_CONSTANT
    T_mean = float(np.mean([get_T_surface_analytic(S, pi, OCEAN_ALBEDO) for pi in p]))
    y_mean = [[float(np.mean(np.interp(tq, t, np.asarray(row, dtype=float))))] for row in y_list]
    return True, T_mean, float(np.mean(p)) / 1e5, _salinity_from_y(y_mean)


_DIAG_NAN = {'da': np.nan, 'calcite_si': np.nan, 'ocean_si': np.nan, 'alk_flux': np.nan, 'pH': np.nan}
_DIAG_KEYS = (('da', 'da'), ('calcite_si', 'calcite_si'), ('ocean_si', 'ocean_si'),
              ('alk_flux', 'alk_flux'), ('pH', 'pH_seafloor'))   # (column, key in the run's block)


def _diag_from_run(d):
    """Diagnostics Planet.time_evolve recorded for the final state, or NaNs for runs without them."""
    block = d.get('diagnostics')
    if not (isinstance(block, dict) and 'da' in block):
        return _DIAG_NAN.copy()
    return {key: np.nan if block.get(src) is None else float(block[src]) for key, src in _DIAG_KEYS}


_DIAG_CACHE = {}
RUN_PATH = None   # directory load_data read, so diagnostics come from there even when figures go elsewhere
_DIAG_VERSION = 1
_DIAG_CACHE_FILE = '.plot_diag_cache.json'
_diag_cache_loaded = set()
_diag_cache_dirty = set()


def _diag_cache_key(fpath):
    """Size and mtime of a run file, so a rewritten run invalidates its cached record."""
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
            continue
    if n:
        print(f"  reusing cached diagnostics for {n} run(s) from {_DIAG_CACHE_FILE}")


def _save_diag_cache(output_path):
    """Write the sidecar cache, merged with what is on disk, if anything new was computed."""
    if output_path not in _diag_cache_dirty:
        return
    _diag_cache_dirty.discard(output_path)
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
    """Add da, calcite_si, ocean_si, alk_flux and seafloor pH columns, cached per run file.

    Rows without a recorded pH keep the pH stored in the run.
    """
    run_path = RUN_PATH or output_path
    _load_diag_cache(run_path)
    records = []
    todo = sum(1 for n in df['name'] if os.path.join(run_path, f'{n}.json') not in _DIAG_CACHE)
    if todo:
        print(f"  reading diagnostics for {todo} run(s)...", flush=True)
    for name in df['name']:
        fpath = os.path.join(run_path, f'{name}.json')
        if fpath not in _DIAG_CACHE:
            try:
                with open(fpath) as fh:
                    _DIAG_CACHE[fpath] = _diag_from_run(json.load(fh))
            except Exception:
                _DIAG_CACHE[fpath] = _DIAG_NAN.copy()
            _diag_cache_dirty.add(run_path)
        records.append(_DIAG_CACHE[fpath])
    _save_diag_cache(run_path)
    diag_df = pd.DataFrame(records, index=df.index)
    df = df.assign(**diag_df.drop(columns=['pH']))
    df['pH'] = diag_df['pH'].where(diag_df['pH'].notna(), df['pH'])
    return df


def _t_end_gyr(d):
    """Integration limit in Gyr: recorded by time_evolve, else the run-name `_tend` tag, else 2."""
    if d.get('t_end_yr') is not None:
        return round(float(d['t_end_yr']) / 1e9, 6)
    m = re.search(r'_tend([0-9.]+)', d.get('name', ''))
    return float(m.group(1)) if m else 2.0


def load_data(output_path):
    """One row per finished run in `output_path`, with its sweep axes and final state."""
    global RUN_PATH
    RUN_PATH = output_path
    rows = []
    for f in sorted(glob.glob(os.path.join(output_path, 'planet_*.json'))):
        try:
            with open(f) as fh:
                d = json.load(fh)
        except (json.JSONDecodeError, OSError):
            print(f"  Skipping (unreadable, probably mid-write): {os.path.basename(f)}")
            continue
        if 'termination' not in d:
            print(f"  Skipping (no termination): {os.path.basename(f)}")
            continue
        y_list = d.get('data', {}).get('y', [])
        rows.append({
            'name':               d.get('name', ''),
            'instellation':       float(d['instellation']),
            'outgassing':         float(d['outgassing']),
            'crust_production':   float(d['crust_production_rate']),
            'reverse_weathering': bool(d.get('reverse_weathering', False)),
            'ocean_depth':        float(d['ocean_depth']),
            'mg_si':              float(d.get('mantle_mg_si', REF_MG_SI)),
            'delta_iw':           float(d.get('delta_iw', REF_DIW)),
            'pe':                 float(d['pe']) if d.get('pe') is not None else PE_LEGACY_DEFAULT,
            'f_HT':               float(d.get('f_HT', 0.0)),
            'alpha':              float(d.get('alpha', ALPHA_REF)),
            'kd_mg':              float(d.get('kd_mg_ht', KD_MG_HT)),
            'k_na':               float(d.get('k_na_cont_removal', K_NA_CONT_REMOVAL)),
            'land_fraction':      float(d.get('land_fraction', 0.0)),
            'termination':        d['termination'],
            'domain_wall':        d.get('domain_wall'),
            'end_time_yr':        d.get('end_time_yr', np.nan),
            'T':                  float(d.get('T', np.nan)),
            'P_CO2':              d.get('P_CO2', np.nan),
            'pH':                 d.get('pH', np.nan),   # replaced by the recorded seafloor pH in _add_diag_columns
            'salinity':           _salinity_from_y(y_list) if y_list else np.nan,
            'cl_ratio':           float(d.get('cl_outgassing_ratio') or 0.0),
            'oscillating':        False,
            't_end_gyr':          _t_end_gyr(d),
        })
        if d['termination'] == 'timeout':   # a cycling run never converges, so it ends at t_end
            osc, T_m, p_m, sal_m = _oscillation(d, y_list)
            if osc:
                if OSCILLATION_MODE == 'exclude':
                    rows.pop()
                    continue
                rows[-1]['oscillating'] = True
                if OSCILLATION_MODE == 'mean':
                    rows[-1].update(T=T_m, P_CO2=p_m, salinity=sal_m)
    df = pd.DataFrame(rows)
    # The alpha-outgassing plane once wrote out_1.0 beside the basic sweep's out_1: the same config twice.
    n0 = len(df)
    config = [c for c in df.columns if c not in ('name', 'termination', 'domain_wall', 'end_time_yr',
                                                   'T', 'P_CO2', 'pH', 'salinity', 'oscillating')]
    df = df.drop_duplicates(config).reset_index(drop=True)
    print(f"Loaded {len(df)} simulations" + (f" ({n0 - len(df)} duplicate configs dropped)."
                                             if len(df) < n0 else "."))
    n_osc = int(df['oscillating'].sum()) if 'oscillating' in df else 0
    if n_osc or OSCILLATION_MODE == 'exclude':
        print(f"  Oscillating runs: {n_osc} shown as their {'tail mean' if OSCILLATION_MODE == 'mean' else 'final state'}"
              if OSCILLATION_MODE != 'exclude' else "  Oscillating runs excluded (--oscillating exclude).")
    return df


# --- Run selection ---
CHEM_KNOBS = {
    'alpha': r'Reactive area scaling $\alpha$',
    'kd_mg': r'Mg$\rightarrow$Ca exchange $k_{Mg}$',
    'k_na':  r'Na sink $k_{Na}$',
}
CHEM_SHIPPED = {'alpha': ALPHA_REF, 'kd_mg': KD_MG_HT, 'k_na': K_NA_CONT_REMOVAL}
CHEM_OVERRIDE = {}     # set from the CLI to choose which chemistry the main plots show
_chem_pinned = set()   # (column, value) already reported


def _ref_crust(df):
    """Mask for Earth's mantle Mg/Si and core-formation dIW."""
    return np.isclose(df['mg_si'], REF_MG_SI) & np.isclose(df['delta_iw'], REF_DIW)


def _ref_redox(df):
    """Mask pinning pe to REF_PE; sweeps run both redox arms, so line figures must pick one."""
    return np.isclose(df['pe'], REF_PE)


def _chem_reference(df, col):
    """Most-run value of a chemistry constant, avoiding disabled (0) terms, then nearest the shipped one."""
    if col in CHEM_OVERRIDE:
        return CHEM_OVERRIDE[col]
    counts = df[col].value_counts()
    tied = [v for v, n in counts.items() if n == counts.max()]
    return min(tied, key=lambda v: (v == 0, abs(v - CHEM_SHIPPED[col])))


def _ref_chem(df):
    """Mask pinning each swept chemistry constant to its reference value."""
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


_setup_pinned = set()   # (column, value) already reported


def _ref_setup(df):
    """Mask pinning the Cl outgassing ratio and the integration limit to their most-run values.

    The paired Cl sweep (basic_cl: Cl on, 4.5 Gyr) shares every other axis with basic, so without
    this its runs land on the same lines and every crust = 1 figure zigzags between the two.
    """
    mask = pd.Series(True, index=df.index)
    for col in ('cl_ratio', 't_end_gyr'):
        if col not in df.columns or df[col].nunique() <= 1:
            continue
        ref = df[col].value_counts().idxmax()
        mask &= np.isclose(df[col], ref)
        if (col, ref) not in _setup_pinned:
            _setup_pinned.add((col, ref))
            others = sorted(v for v in df[col].unique() if not np.isclose(v, ref))
            print(f"  Pinning {col} = {ref:g} for the main plots "
                  f"(also present: {', '.join(f'{v:g}' for v in others)}).")
    return mask


def _sweep_mask(df, depth=3000, land=0.0, crust=True, chem=True, redox=True, f_ht=False, rw=True,
                setup=True, **pins):
    """Standard selection: reverse weathering on, one depth and land fraction, reference crust/chem/redox.

    Pass None or False to skip a filter (rw=False selects the no-reverse-weathering control, rw=None
    both arms, setup=False keeps every Cl ratio and run length); `pins` are further column values
    matched with np.isclose.
    """
    mask = (pd.Series(True, index=df.index) if rw is None
            else df['reverse_weathering'].astype(bool) == rw)
    if depth is not None:
        mask &= df['ocean_depth'] == depth
    if land is not None:
        mask &= np.isclose(df['land_fraction'], land)
    if crust:
        mask &= _ref_crust(df)
    if chem:
        mask &= _ref_chem(df)
    if redox:
        mask &= _ref_redox(df)
    if f_ht:
        mask &= df['f_HT'] == 0.0
    if setup:
        mask &= _ref_setup(df)
    for col, value in pins.items():
        mask &= np.isclose(df[col], value)
    return mask


def _base(df, mg_si=None):
    """Land-free basic-sweep runs at 3 km with reference chemistry and redox, at Earth's or a given Mg/Si."""
    pins = {} if mg_si is None else {'mg_si': mg_si, 'delta_iw': REF_DIW}
    return df[_sweep_mask(df, crust=mg_si is None, **pins) & (df['outgassing'] > 0)]


def basic_plane_mg_si(df, min_combos=4):
    """Mantle Mg/Si values with a full basic sweep, i.e. several (outgassing, crust) combinations."""
    pool = df[_sweep_mask(df, crust=False, delta_iw=REF_DIW) & (df['outgassing'] > 0)]
    if pool.empty:
        return []
    n_combos = pool.groupby('mg_si').apply(
        lambda g: g.groupby(['outgassing', 'crust_production']).ngroups, include_groups=False)
    return sorted(float(v) for v, n in n_combos.items() if n >= min_combos)


# --- Shared plot helpers ---
FIGURE_SUBDIR = 'figures'


def figure_path(output_path, name):
    """Path for a figure in the `figures/` subdirectory of the sweep directory."""
    directory = os.path.join(output_path, FIGURE_SUBDIR)
    os.makedirs(directory, exist_ok=True)
    return os.path.join(directory, name)


def _save_fig(fig, path, tight=False):
    """Save PNG and PDF at the created size; tight=True crops, for diagnostics only."""
    kw = {'bbox_inches': 'tight'} if tight else {}
    stem = os.path.splitext(path)[0]
    for ext in ('png', 'pdf'):
        fig.savefig(f'{stem}.{ext}', **kw)
    plt.close(fig)
    print(f"Saved {stem}.png / .pdf")


def _panel_groups(split):
    """[(columns, filename suffix), ...]: all four panels, or T/pCO2 and pH/salinity separately."""
    if split:
        return [(['T', 'P_CO2'], '_tp'), (['pH', 'salinity'], '_chem')]
    return [(PANEL_COLS, '')]


def _temperature_bands(ax, **kw):
    """Shade the snowball and runaway bands on a temperature axis."""
    ax.axhspan(T_SNOWBALL - 25, T_SNOWBALL, color='blue', alpha=0.12, **kw)
    ax.axhspan(T_RUNAWAY - 20, T_RUNAWAY, color='red', alpha=0.12, **kw)


def _label_along(ax, x, y, text, frac=0.88, offset=5):
    """Label a curve in place, rotated to its on-screen slope; needs the axis limits already set."""
    i = int(frac * (len(x) - 1))
    (x0, y0), (x1, y1) = ax.transData.transform([(x[max(i - 1, 0)], y[max(i - 1, 0)]), (x[i], y[i])])
    ax.annotate(text, xy=(x[i], y[i]), xytext=(0, offset), textcoords='offset points',
                rotation=np.degrees(np.arctan2(y1 - y0, x1 - x0)), rotation_mode='anchor',
                ha='center', va='bottom', fontsize=6, color='0.3', zorder=5,
                bbox=dict(boxstyle='square,pad=0.1', facecolor='white', edgecolor='none',
                          alpha=0.75))


def _style_axes(axes, cols, x_lims=(0.25, 1.45), show_hz=None, show_eq_temp=False):
    """Label and scale panels against instellation; HZ edges and the equilibrium curve go on T only."""
    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(*x_lims)
    for ax, col in zip(axes, cols):
        if col == 'T':
            _draw_hz_edges(ax, show_hz)
            ax.set_ylabel('Temperature (K)')
            _temperature_bands(ax)
            ax.set_ylim(235, 360)
            if show_eq_temp:
                s_eq = np.linspace(x_lims[0], x_lims[1], 300)
                t_eq = equilibrium_temperature(s_eq)
                ax.plot(s_eq, t_eq, color='k', linestyle='-.', linewidth=0.8, zorder=1, alpha=0.7)
                # Mid-curve: the planet's own line runs below it there, and the right end would
                # clash with the HZ inner-edge label.
                _label_along(ax, s_eq, t_eq, 'Equilibrium temperature', frac=0.45)
        elif col == 'P_CO2':
            ax.set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 20)
        elif col == 'pH':
            ax.set_ylabel('Ocean pH')
            ax.set_ylim(3, 12)   # blank oceans span pH 3.5-11.8
        elif col == 'salinity':
            ax.set_ylabel('Salinity (g/kg)')
            ax.set_yscale('log')
            ax.set_ylim(5e-2, 1e2)   # blank oceans reach 0.07 g/kg
        elif col == 'calcite_si':
            ax.set_ylabel('Calcite SI')
            ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
        elif col == 'alk_flux':
            ax.set_ylabel('Alk. flux (Tmol eq/yr)')
            ax.set_yscale('symlog', linthresh=0.001)
            ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
    axes[-1].set_xlabel('Instellation (S/S₀)')


def _style_combined_col(axes_c, ci, n_cols, title='', cols=None, show_hz=None):
    """Style one column of a multi-column grid; only the first column keeps y labels."""
    cols = PANEL_COLS if cols is None else cols
    col_axes = axes_c[:, ci]
    _style_axes(col_axes, cols, show_hz=show_hz)
    if title:
        axes_c[0, ci].set_title(title)
    if ci > 0:
        for ax in col_axes:
            ax.set_ylabel('')
            plt.setp(ax.get_yticklabels(), visible=False)
    if ci != n_cols // 2:
        col_axes[-1].set_xlabel('')


def _add_colorbar(fig, ax, cmap, norm, label, ticks=None, ticklabels=None, aspect=30):
    """Colour bar for a cmap/norm pair, with optional explicit ticks."""
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, location='right', pad=0.02, aspect=aspect)
    cbar.set_label(label)
    if ticks is not None:
        cbar.set_ticks(ticks)
        if ticklabels is not None:
            cbar.set_ticklabels(ticklabels)
    return cbar


def _outgassing_norm(values):
    return mcolors.LogNorm(vmin=min(values), vmax=max(values))


def _outgassing_colorbar(fig, ax, values, norm, aspect):
    """The 'Earth Outgassing' colour bar with a tick per value."""
    _add_colorbar(fig, ax, OUTGASSING_CMAP, norm, 'Earth Outgassing', ticks=values,
                  ticklabels=[f'{v}×' for v in values], aspect=aspect)


def _make_legend_handles(show_markers=True, prefix_handles=None):
    """Legend handles: DA_LEGEND (or `prefix_handles`) plus the termination markers."""
    handles = list(prefix_handles if prefix_handles is not None else DA_LEGEND)
    if show_markers:
        handles += [plt.scatter([], [], marker=m, s=28, color='k', label=TERM_LABELS[t])
                    for t, m in HAB_MARKERS.items()]
        handles += [plt.scatter([], [], marker=m, s=50, label=TERM_LABELS[t],
                                facecolors='none', edgecolors='k', linewidths=1.4)
                    for t, m in FAILED_MARKERS.items()]
    return handles


def _add_figure_legend(fig, axes, handles, loc='outside lower center', **kw):
    """Figure legend with as many columns as fit within the width of the panels."""
    axs = [a for a in np.ravel(np.asarray(axes, dtype=object)) if a is not None]
    fig.canvas.draw()
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


def _at_seafloor_floor(T):
    """True where the seafloor temperature implied by surface T sits at its floor (SEAFLOOR_T_FLOOR)."""
    return 1.02 * np.asarray(T, dtype=float) - 16.7 <= SEAFLOOR_T_FLOOR + 0.001


def _da_trustworthy(group):
    """Rows whose Da may set the line style: a real evaluated state, not a temperature-wall sentinel.

    A run stopped at the CO2 ceiling counts, since its temperature is a solution rather than a clamp.
    """
    wall = group['domain_wall'] if 'domain_wall' in group else pd.Series(None, index=group.index)
    return (group['termination'].isin(DA_TRUSTWORTHY)
            | (group['termination'].isin(OUT_OF_DOMAIN) & (wall == 'co2_high'))).to_numpy()


def _bracketed_crossing(group, cutoff, da_trust):
    """Last drawn row when Da crosses 1 in the step onto a temperature-wall run, else None.

    The wall run's Da comes from a clamped state, so the crossing can only be bracketed, not placed.
    """
    if not 0 < cutoff < len(group):
        return None
    nxt = group.iloc[cutoff]
    if nxt['termination'] not in OUT_OF_DOMAIN or nxt.get('domain_wall') not in ('hot', 'cold'):
        return None
    da = group['da'].to_numpy(dtype=float)[:cutoff]
    ok = np.flatnonzero(np.isfinite(da) & da_trust[:cutoff])
    if ok.size == 0 or da[ok[-1]] >= 1:
        return None
    return group.iloc[ok[-1]] if np.isfinite(nxt['da']) and nxt['da'] >= 1 else None


def _plot_line_da_style(ax, x, y, da, color, at_floor=None, trustworthy=None,
                        linewidth=1.4, alpha=0.8, zorder=3):
    """Line dotted where seafloor T is floored, dashed where Da >= 1, solid where Da < 1.

    Points whose Da is untrustworthy copy the nearest trustworthy style; kinetic/thermodynamic
    transitions get an open circle at the last trustworthy point before the change.
    """
    if len(x) < 2:
        return
    at_floor = np.zeros(len(x), dtype=bool) if at_floor is None else at_floor
    trustworthy = np.ones(len(x), dtype=bool) if trustworthy is None else trustworthy

    def _raw_ls(i):
        if at_floor[i]:
            return ':'
        if not (trustworthy[i] and np.isfinite(da[i])):
            return None
        return '-' if da[i] < 1 else '--'

    raw = [_raw_ls(i) for i in range(len(x))]
    resolved = list(raw)
    for i, s in enumerate(raw):
        if s is None:
            prv = next((resolved[j] for j in range(i - 1, -1, -1) if raw[j] is not None), None)
            nxt = next((raw[j] for j in range(i + 1, len(raw)) if raw[j] is not None), None)
            resolved[i] = prv or nxt or '-'

    seg_x, seg_y, current_ls = [x[0]], [y[0]], resolved[0]
    trans_x, trans_y = [], []
    for i in range(1, len(x)):
        ls = resolved[i]
        if ls == current_ls:
            seg_x.append(x[i])
            seg_y.append(y[i])
            continue
        ax.plot(seg_x, seg_y, color=color, linewidth=linewidth, alpha=alpha,
                linestyle=current_ls, zorder=zorder)
        if {current_ls, ls} == {'-', '--'}:
            j = i - 1
            while j > 0 and raw[j] is None:
                j -= 1
            trans_x.append(x[j])
            trans_y.append(y[j])
        seg_x, seg_y, current_ls = [seg_x[-1], x[i]], [seg_y[-1], y[i]], ls
    if len(seg_x) >= 2:
        ax.plot(seg_x, seg_y, color=color, linewidth=linewidth, alpha=alpha,
                linestyle=current_ls, zorder=zorder)
    if trans_x:
        ax.scatter(trans_x, trans_y, facecolors='none', edgecolors=color,
                   s=30, linewidths=1.2, zorder=5)


def _plot_group_on_axes(axes, group, color, linestyle='-', show_markers=True, cols=None):
    """Draw one instellation-sorted run group on each panel, Da-styled, with optional termination markers."""
    cols = PANEL_COLS if cols is None else cols
    hab = group[group['termination'].isin(HABITABLE)]
    if hab.empty:
        return
    non_hab = group[~group['termination'].isin(HABITABLE)]
    # End the line at the last converged/timeout/wall-clock run: a trailing tail is often a clamp sentinel.
    line_end = group['termination'].isin(DA_TRUSTWORTHY).values
    cutoff = (np.nonzero(line_end)[0].max() + 1) if line_end.any() else len(group)
    line_group = group.iloc[:cutoff]
    da_styled = 'da' in group.columns and linestyle == '-'
    if da_styled:
        da_trust = _da_trustworthy(group)
        bracket = _bracketed_crossing(group, cutoff, da_trust)

    for ax, col in zip(axes, cols):
        if len(line_group) > 1:
            if da_styled:
                _plot_line_da_style(ax, line_group['instellation'].values, line_group[col].values,
                                    line_group['da'].values, color,
                                    at_floor=_at_seafloor_floor(line_group['T'].values),
                                    trustworthy=da_trust[:cutoff])
                if bracket is not None and np.isfinite(bracket[col]):
                    # Same open circle as a measured crossing: the distinction does not need marking.
                    ax.scatter(bracket['instellation'], bracket[col], s=30, facecolors='none',
                               edgecolors=color, linewidths=1.2, zorder=5)
            else:
                ax.plot(line_group['instellation'], line_group[col], color=color, linewidth=1.4,
                        alpha=0.8, linestyle=linestyle, zorder=3)
        if not show_markers:
            continue
        for term, hab_marker in HAB_MARKERS.items():
            sub = hab[hab['termination'] == term]
            if not sub.empty:
                ax.scatter(sub['instellation'], sub[col], color=color,
                           marker=hab_marker, s=28, zorder=4)
        for _, row in non_hab.iterrows():
            if np.isfinite(row[col]):
                ax.scatter(row['instellation'], row[col],
                           marker=FAILED_MARKERS.get(row['termination'], 'x'), s=55,
                           facecolors='none', color=color, zorder=4, linewidths=1.4)


def _best_operating_point(pool, col, what):
    """(subset, outgassing, crust) for the pair spanning the most values of `col`, or None."""
    counts = pool.groupby(['outgassing', 'crust_production'])[col].nunique()
    if counts.empty or counts.max() < 2:
        print(f"No {what} variation at a fixed (outgassing, crust) -- skipping.")
        return None
    # Ties (e.g. the alpha x outgassing plane) go to the pair nearest Earth's (1x, 1x), not the first one.
    tied = [k for k, n in counts.items() if n == counts.max()]
    best_o, best_c = min(tied, key=lambda k: abs(np.log10(k[0])) + abs(np.log10(k[1])))
    subset = pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)]
    print(f"{what} plot: using outgassing={best_o:g}, crust={best_c:g} "
          f"({subset[col].nunique()} values).")
    return subset, best_o, best_c


def _x_limits(subset, default=(0.25, 1.45), frac=0.05):
    """Instellation limits with a small margin."""
    lo, hi = subset['instellation'].min(), subset['instellation'].max()
    if pd.isna(lo):
        return default
    margin = (hi - lo) * frac if hi != lo else 0.1
    return (lo - margin, hi + margin)


def _value_norm(values, pad=0.05):
    """Log norm when positive values span a decade or more, otherwise linear widened by `pad`."""
    lo, hi = min(values), max(values)
    if lo > 0 and hi / lo >= 10:
        return mcolors.LogNorm(vmin=lo, vmax=hi)
    span = hi - lo
    return mcolors.Normalize(vmin=lo - pad * span if span else lo - 1,
                             vmax=hi + pad * span if span else hi + 1)


def _colorbar_ticks(values, max_ticks=10):
    """Explicit ticks and labels only when there are few enough to label individually."""
    if len(values) > max_ticks:
        return None, None
    return list(values), [f'{v:g}' for v in values]


def _faceted_lines(subset, col, values, colours, cmap, norm, cbar_label, stem, output_path,
                   split_panels=False, show_markers=False, x_lims=None,
                   ticks=None, ticklabels=None, aspect_per_row=None,
                   width='single', height=None, show_hz=None, groups=None):
    """Stacked panels against instellation, one line per value of `col`, one file per panel group.

    `groups` overrides the panel groups, e.g. [(['pH'], '_ph')] for a single-panel version.
    """
    for cols, sfx in (groups if groups is not None else _panel_groups(split_panels)):
        n_rows = len(cols)
        fig, axes = plt.subplots(n_rows, 1, sharex=True, squeeze=False,
                                 figsize=figure_size(width, height, n_rows))
        axes = axes[:, 0]
        for value, colour in zip(values, colours):
            group = subset[subset[col] == value].sort_values('instellation')
            if not group.empty:
                _plot_group_on_axes(axes, group, colour, show_markers=show_markers, cols=cols)
        _style_axes(axes, cols, show_hz=show_hz, **({} if x_lims is None else {'x_lims': x_lims}))
        _add_colorbar(fig, list(axes), cmap, norm, cbar_label, ticks=ticks, ticklabels=ticklabels,
                      **({} if aspect_per_row is None else {'aspect': n_rows * aspect_per_row}))
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        _save_fig(fig, figure_path(output_path, f'{stem}{sfx}.png'))


# --- Basic and one-axis sweeps ---
def plot_basic(df, output_path, all_results=True, multiple_plots=False, split_panels=True,
               sequence=False, width='double', height=None, mg_si=None, show_hz=None):
    """Basic sweep against instellation: one column per crust rate, one line per outgassing rate.

    all_results draws every value on an on-screen grid; mg_si tags the filenames of end-member planes.
    """
    base = _base(df, mg_si=mg_si)
    if base.empty:
        print(f"No basic sweep data at Mg/Si = {mg_si:g} -- skipping." if mg_si is not None
              else "No basic sweep data -- skipping.")
        return
    base = _add_diag_columns(base, output_path)
    mg_tag = '' if mg_si is None else f'_mgsi{mg_si:g}'
    mg_title = '' if mg_si is None else f'   (mantle Mg/Si = {mg_si:g})'
    if all_results:
        crust_rates = sorted(base['crust_production'].unique())
        outgassing_vals = sorted(base['outgassing'].unique())
    else:
        crust_rates = [0.1, 1, 10]
        outgassing_vals = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    norm = _outgassing_norm(outgassing_vals)

    def draw(axes, rates, cols, show_markers, show=lambda c, o: True):
        for ci, c in enumerate(rates):
            sub = base[base['crust_production'] == c]
            for o in outgassing_vals:
                group = sub[sub['outgassing'] == o].sort_values('instellation')
                if show(c, o) and not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, OUTGASSING_CMAP(norm(o)),
                                        show_markers=show_markers, cols=cols)

    def finish(fig, axes, title, name, show_markers, aspect, tight=False):
        _outgassing_colorbar(fig, list(axes.ravel()), outgassing_vals, norm, aspect)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        fig.suptitle(title)
        _save_fig(fig, figure_path(output_path, name), tight=tight)

    # (filename, which lines to draw); the seq figures build the full figure up line by line.
    scenarios = [(f'sweep_basic{"_full" if all_results else ""}{mg_tag}', lambda c, o: True)]
    if sequence:
        scenarios += [(f'sweep_basic_seq1{mg_tag}', lambda c, o: np.isclose(c, 1) and np.isclose(o, 1)),
                      (f'sweep_basic_seq2{mg_tag}', lambda c, o: np.isclose(c, 1)),
                      (f'sweep_basic_seq3{mg_tag}', lambda c, o: np.isclose(o, 1))]

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        if multiple_plots:
            for c in crust_rates:
                fig, axes = plt.subplots(n_rows, 1, sharex=True, squeeze=False,
                                         figsize=diagnostic_size(n_rows, 1, col_width=7.0, pad=0.0))
                draw(axes, [c], cols, True)
                _style_axes(axes[:, 0], cols, show_hz=show_hz)
                finish(fig, axes, f'Crust production = {c}× Earth{mg_title}',
                       f'sweep_basic_crust{c}{mg_tag}{sfx}.png', True, n_rows * 7.5)

        n_cols = len(crust_rates)
        figsize = (diagnostic_size(n_rows, n_cols) if all_results
                   else figure_size(width, height, n_rows))
        for i, (stem, show) in enumerate(scenarios):
            fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True, sharey='row',
                                     squeeze=False)
            draw(axes, crust_rates, cols, all_results, show)
            for ci, c in enumerate(crust_rates):
                _style_combined_col(axes, ci, n_cols, title=f'{c}×', cols=cols, show_hz=show_hz)
            finish(fig, axes, f'Earth crust production rate{mg_title}', f'{stem}{sfx}.png',
                   all_results, n_rows * 10, tight=all_results and i == 0)


def plot_basic_ph(df, output_path, crust_production=1.0, width='single', height=2.2, show_hz=None):
    """Abridged basic sweep: the pH panel alone at one crust production rate."""
    base = _base(df)
    base = base[np.isclose(base['crust_production'], crust_production)]
    if base.empty:
        print(f"No basic sweep data at crust production {crust_production:g} -- skipping pH panel.")
        return
    base = _add_diag_columns(base, output_path)
    outgassing_vals = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    norm = _outgassing_norm(outgassing_vals)

    fig, ax = plt.subplots(1, 1, figsize=figure_size(width, height))
    for o in outgassing_vals:
        group = base[base['outgassing'] == o].sort_values('instellation')
        if not group.empty:
            _plot_group_on_axes([ax], group, OUTGASSING_CMAP(norm(o)), show_markers=False,
                                cols=['pH'])
    _style_axes([ax], ['pH'], show_hz=show_hz)
    ax.set_title(f'Crust production = {crust_production:g}× Earth')
    _outgassing_colorbar(fig, ax, outgassing_vals, norm, 15)
    _add_figure_legend(fig, [ax], _make_legend_handles(show_markers=False))
    _save_fig(fig, figure_path(output_path, f'sweep_basic_crust{crust_production:g}_ph.png'))


def plot_basic_mgsi_grid(df, output_path, split_panels=True, show_markers=False,
                         crust_production=None, mg_si_values=None, width='double', height=None,
                         show_hz=None):
    """Basic sweep at one crust rate with mantle Mg/Si as columns and outgassing as colour."""
    mg_vals = basic_plane_mg_si(df) if mg_si_values is None else sorted(mg_si_values)
    if len(mg_vals) < 2:
        print("Fewer than 2 Mg/Si values with a basic sweep -- skipping the Mg/Si grid.")
        return
    if crust_production is None:   # Earth's rate if every column has it, else the best-covered rate
        pool = pd.concat([_base(df, mg_si=m) for m in mg_vals])
        by_rate = pool.groupby('crust_production')['mg_si'].nunique()
        crust_production = (1.0 if by_rate.get(1.0, 0) == len(mg_vals)
                            else float(by_rate.idxmax()))

    subsets = {}
    for m in mg_vals:
        sub = _base(df, mg_si=m)
        sub = sub[sub['crust_production'] == crust_production]
        if not sub.empty:
            subsets[m] = _add_diag_columns(sub, output_path)
    mg_vals = [m for m in mg_vals if m in subsets]
    if len(mg_vals) < 2:
        print(f"Fewer than 2 Mg/Si values at crust production {crust_production:g}x -- skipping.")
        return
    outgassing_vals = sorted(set().union(*(set(v['outgassing']) for v in subsets.values())))
    print(f"Mg/Si basic grid: {[f'{m:g}' for m in mg_vals]} at crust production "
          f"{crust_production:g}x, {len(outgassing_vals)} outgassing values")
    norm = _outgassing_norm(outgassing_vals)

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        fig, axes = plt.subplots(n_rows, len(mg_vals), figsize=figure_size(width, height, n_rows),
                                 sharex=True, sharey='row', squeeze=False)
        for ci, m in enumerate(mg_vals):
            for o in outgassing_vals:
                group = subsets[m][subsets[m]['outgassing'] == o].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, OUTGASSING_CMAP(norm(o)),
                                        show_markers=show_markers, cols=cols)
            title = f'{m:g}' + (' (Earth)' if np.isclose(m, REF_MG_SI) else '')
            _style_combined_col(axes, ci, len(mg_vals), title=title, cols=cols, show_hz=show_hz)
        _outgassing_colorbar(fig, list(axes.ravel()), outgassing_vals, norm, n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        fig.suptitle(f'Mantle Mg/Si   (crust production = {crust_production:g}× Earth)')
        _save_fig(fig, figure_path(output_path, f'sweep_basic_mgsi_grid{sfx}.png'))


def _one_axis_sweep(pool, col, what, cmap, label, stem, output_path, select=None, norm_pad=0.05,
                    tick_fmt=None, **kw):
    """One line per value of `col` at the best (outgassing, crust) operating point, via _faceted_lines."""
    picked = _best_operating_point(pool, col, what)
    if picked is None:
        return
    subset = _add_diag_columns(picked[0] if select is None else select(picked[0]), output_path)
    values = sorted(subset[col].unique())
    norm = _value_norm(values, pad=norm_pad)
    ticks, ticklabels = _colorbar_ticks(values)
    if tick_fmt is not None and ticks is not None:
        ticklabels = [tick_fmt(v) for v in ticks]
    _faceted_lines(subset, col, values, [cmap(norm(v)) for v in values], cmap, norm, label, stem,
                   output_path, x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                   aspect_per_row=7.5, **kw)


DEPTHS_SHOWN = (300, 1000, 3000, 10000, 30000)   # m, the depth figure's lines when available


def plot_depth(df, output_path, show_markers=False, split_panels=True, width='single',
               height=None, show_hz=None):
    """Depth sweep: one line per ocean depth."""
    pool = df[_sweep_mask(df, depth=None)]
    if pool.empty:
        print("No data for ocean depth sweep — skipping.")
        return

    def thin(sub):
        # About five depths: the shown set if present (so a stray depth such as the continental
        # baseline's 3.7 km arm cannot displace the 3 km reference), else evenly spaced.
        depths = sorted(d for d in sub['ocean_depth'].unique() if d < 100000)
        shown = [d for d in DEPTHS_SHOWN if d in depths]
        if len(shown) >= 3:
            depths = shown
        elif len(depths) > 5:
            depths = sorted({depths[i] for i in np.linspace(0, len(depths) - 1, 5).round().astype(int)})
        print(f"  Ocean-depth plot: showing {[f'{d:g}' for d in depths]} m.")
        return sub[sub['ocean_depth'].isin(depths)]

    _one_axis_sweep(pool, 'ocean_depth', 'Ocean-depth', DEPTH_CMAP, 'Ocean Depth (km)', 'sweep_depth',
                    output_path, select=thin, norm_pad=0.0, tick_fmt=lambda v: f'{v / 1000:g}',
                    split_panels=split_panels, show_markers=show_markers, width=width,
                    height=height, show_hz=show_hz)


def plot_chemistry(df, output_path, show_markers=False, split_panels=True, width='single',
                   height=None, show_hz=None):
    """One figure per varying chemistry constant (alpha, kd_mg, k_na), the others held at reference."""
    pool_all = df[_sweep_mask(df, chem=False)]
    varying = [c for c in CHEM_KNOBS if not pool_all.empty and pool_all[c].nunique() > 1]
    if not varying:
        print("No chemistry-constant variation in the data — skipping.")
        return
    for col in varying:
        held = pd.Series(True, index=pool_all.index)
        for other in varying:
            if other != col:
                held &= (pool_all[other] == _chem_reference(pool_all, other))
        stem = 'sweep_alpha' if col == 'alpha' else f'sweep_chemistry_{col}'
        _one_axis_sweep(pool_all[held], col, col, CHEM_KNOB_CMAP, CHEM_KNOBS[col], stem, output_path,
                        split_panels=split_panels, show_markers=show_markers, width=width,
                        height=height, show_hz=show_hz)


def plot_pe(df, output_path, show_markers=False, split_panels=True, ocean_depth=3000,
            width='single', height=None, show_hz=None):
    """Redox sweep: one line per ocean pe."""
    pool = df[_sweep_mask(df, depth=ocean_depth, redox=False)]
    if pool.empty:
        print(f"No data for redox sweep at depth {ocean_depth:g} m — skipping.")
        return
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'
    _one_axis_sweep(pool, 'pe', 'pe', PE_CMAP, r'Ocean redox $p_e$', f'sweep_pe{tag}', output_path,
                    split_panels=split_panels, show_markers=show_markers, width=width,
                    height=height, show_hz=show_hz)


# --- Crust composition ---
def _composition_pool(df, ocean_depth=3000):
    """Runs usable for a composition figure: reference chemistry and redox, land-free, one depth."""
    return df[_sweep_mask(df, depth=ocean_depth, crust=False, f_ht=True)].copy()


def _composition_slice(pool):
    """(subset, outgassing, crust) for the pair spanning the most crust compositions, or None."""
    counts = pool.groupby(['outgassing', 'crust_production']).apply(
        lambda g: max(g['mg_si'].nunique(), g['delta_iw'].nunique()), include_groups=False)
    if counts.empty or counts.max() <= 1:
        return None
    best_o, best_c = counts.idxmax()
    return pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)], best_o, best_c


def _drop_hidden_mg_si(sub):
    return sub[~np.isclose(sub['mg_si'].to_numpy()[:, None], MG_SI_HIDDEN).any(axis=1)]


def plot_cross(df, output_path, split_panels=True, show_markers=False, ocean_depth=3000,
               width='single', height=None, show_hz=None):
    """Cross sweep: one figure per varying composition axis (Mg/Si, dIW), the other held at Earth's."""
    pool = _composition_pool(df, ocean_depth)
    sliced = None if pool.empty else _composition_slice(pool)
    if sliced is None:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping.")
        return
    subset, best_o, best_c = sliced
    subset = _drop_hidden_mg_si(subset)

    axes_spec = []
    if subset['mg_si'].nunique() > 1:
        axes_spec.append(('mg_si', 'Mantle Mg/Si', lambda v: f'{v:g}'))
    if subset['delta_iw'].nunique() > 1:
        axes_spec.append(('delta_iw', r'Core-formation $\Delta$IW', lambda v: f'{v:+g}'))
    if not axes_spec:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping.")
        return
    subset = _add_diag_columns(subset, output_path)
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'

    for key, label, fmt in axes_spec:
        cut, held = subset, []
        for other, ref in (('mg_si', REF_MG_SI), ('delta_iw', REF_DIW)):
            if other != key and subset[other].nunique() > 1:
                cut = cut[np.isclose(cut[other], ref)]
                held.append(f'{other}={ref:g}')
        if cut.empty or cut[key].nunique() < 2:
            print(f"  {key}: fewer than 2 values once {', '.join(held)} held -- skipping.")
            continue
        values = sorted(cut[key].unique())
        print(f"Crust composition [{key}] depth={ocean_depth:g}: {len(values)} values, "
              f"outgassing={best_o:g}, crust={best_c:g}"
              + (f", holding {', '.join(held)}" if held else ""))
        cmap = PARAM_CMAPS[key]
        numeric = [float(v) for v in values]
        norm = _value_norm(numeric)
        common = dict(show_markers=show_markers, ticks=numeric, ticklabels=[fmt(v) for v in values],
                      width=width, show_hz=show_hz)
        lines = (cut, key, values, [cmap(norm(n)) for n in numeric], cmap, norm, label,
                 f'sweep_cross_{key}{tag}', output_path)
        _faceted_lines(*lines, split_panels=split_panels, height=height, **common)
        if key == 'mg_si':   # abridged pH-only version
            _faceted_lines(*lines, height=2.2 if height is None else height,
                           groups=[(['pH'], '_ph')], aspect_per_row=15, **common)


def _diw_title(dw, shown):
    """Column title naming the redox end-members and the Earth reference."""
    if np.isclose(dw, REF_DIW):
        return f'{dw:+g} (Earth-like)'
    if len(shown) >= 2 and dw == min(shown):
        return f'{dw:+g} (reduced)'
    if len(shown) >= 2 and dw == max(shown):
        return f'{dw:+g} (oxidised)'
    return f'{dw:+g}'


def _three_columns(values, ref, n=3):
    """`n` representative values: the one nearest `ref`, then the extremes either side."""
    vals = sorted(values)
    if n >= len(vals):
        return vals
    centre = min(vals, key=lambda v: abs(v - ref))
    lo = [v for v in vals if v < centre]
    hi = [v for v in vals if v > centre]
    picked = [centre] + lo[:1] + hi[-1:]
    pool = [v for v in vals if v not in picked]
    while len(picked) < n and pool:
        picked.append(pool.pop(len(pool) // 2))
    return sorted(picked)[:n]


def _composition_grid(subset, col_key, col_vals, col_title, line_key, line_vals, cmap, cbar_label,
                      tick_fmt, suptitle, stem, output_path, split_panels, show_markers, width,
                      height, show_hz):
    """Grid with one column per `col_key` value and one line per `line_key` value."""
    norm = mcolors.Normalize(vmin=min(line_vals), vmax=max(line_vals))
    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        fig, axes = plt.subplots(n_rows, len(col_vals), figsize=figure_size(width, height, n_rows),
                                 sharex=True, sharey='row', squeeze=False)
        for ci, cv in enumerate(col_vals):
            col_df = subset[subset[col_key] == cv]
            for lv in line_vals:
                group = col_df[col_df[line_key] == lv].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, cmap(norm(lv)),
                                        show_markers=show_markers, cols=cols)
            _style_combined_col(axes, ci, len(col_vals), title=col_title(cv), cols=cols,
                                show_hz=show_hz)
        _add_colorbar(fig, list(axes.ravel()), cmap, norm, cbar_label, ticks=line_vals,
                      ticklabels=[tick_fmt(v) for v in line_vals], aspect=n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers))
        fig.suptitle(suptitle)
        _save_fig(fig, figure_path(output_path, f'{stem}{sfx}.png'))


def plot_composition(df, output_path, split_panels=True, show_markers=False, ocean_depth=3000,
                     min_lines=2, n_cols=3, width='double', height=None, show_hz=None):
    """Composition factorial: dIW columns coloured by Mg/Si, and the transpose (Mg/Si columns by dIW)."""
    pool = _composition_pool(df, ocean_depth)
    sliced = None if pool.empty else _composition_slice(pool)
    if sliced is None:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping grid.")
        return
    subset, best_o, best_c = sliced
    subset = _drop_hidden_mg_si(subset)
    if subset['mg_si'].nunique() < 2 or subset['delta_iw'].nunique() < 2:
        print("Composition grid needs both Mg/Si and dIW to vary -- skipping.")
        return
    subset = _add_diag_columns(subset, output_path)
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'
    mg_vals = sorted(subset['mg_si'].unique())
    diw_all = sorted(subset['delta_iw'].unique())
    common = dict(output_path=output_path, split_panels=split_panels, show_markers=show_markers,
                  width=width, height=height, show_hz=show_hz)

    usable = [d for d in diw_all if subset[subset['delta_iw'] == d]['mg_si'].nunique() >= min_lines]
    if not usable:
        print("No dIW column has enough Mg/Si values -- skipping grid.")
        return
    diw_vals = _three_columns(usable, REF_DIW, n_cols)
    dropped = sorted(set(diw_all) - set(diw_vals))
    print(f"Composition grid depth={ocean_depth:g}: dIW columns {diw_vals} x "
          f"{len(mg_vals)} Mg/Si lines, outgassing={best_o:g}, crust={best_c:g}"
          + (f" (not shown: dIW {dropped})" if dropped else ""))
    _composition_grid(subset, 'delta_iw', diw_vals, lambda dw: _diw_title(dw, diw_vals),
                      'mg_si', mg_vals, MG_SI_CMAP, 'Mantle Mg/Si', lambda v: f'{v:g}',
                      r'Core-formation $\Delta$IW', f'sweep_composition{tag}', **common)

    mg_usable = [m for m in mg_vals if subset[subset['mg_si'] == m]['delta_iw'].nunique() >= min_lines]
    if len(mg_usable) >= 2:
        _composition_grid(subset, 'mg_si', _three_columns(mg_usable, REF_MG_SI, n_cols),
                          lambda mg: f'{mg:g}', 'delta_iw', diw_all, DIW_CMAP,
                          r'Core-formation $\Delta$IW', lambda v: f'{v:+g}', 'Mantle Mg/Si',
                          f'sweep_composition_T{tag}', **common)


def plot_composition_map(df, output_path, s_vals=(0.7, 0.9, 1.0, 1.1), ocean_depth=3000,
                         quantity='T', relative=True, min_cells=3, width='double', height=2.6):
    """Composition sweep as a Mg/Si x dIW map, one panel per instellation.

    relative plots the difference from the Earth crust at the same instellation; rows or columns with
    fewer than min_cells in-domain runs are dropped, and out-of-domain cells are marked, not coloured.
    """
    if quantity not in ('T', 'P_CO2', 'pH'):
        raise ValueError(f'unsupported quantity {quantity!r}')
    pool = _composition_pool(df, ocean_depth)
    sliced = None if pool.empty else _composition_slice(pool)
    if sliced is None:
        print(f"No crust composition sweep data at depth {ocean_depth:g} m -- skipping map.")
        return
    subset, best_o, best_c = sliced
    if subset['mg_si'].nunique() < 2 or subset['delta_iw'].nunique() < 2:
        print("Composition map needs both Mg/Si and dIW to vary -- skipping.")
        return
    if quantity == 'pH':
        subset = _add_diag_columns(subset, output_path)

    live = subset[~subset['termination'].isin(OUT_OF_DOMAIN)]
    # Count distinct cells on the other axis, so a one-off cross-design value does not draw a stripe.
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
    log_q = quantity == 'P_CO2'

    def _cell(S, mg, dw):
        r = live[(live['instellation'] == S) & (live['mg_si'] == mg) & (live['delta_iw'] == dw)]
        if r.empty:
            return np.nan
        v = float(r[quantity].iloc[0])
        return np.log10(v) if log_q and v > 0 else (np.nan if log_q else v)

    grids = {}
    for S in s_present:
        g = np.array([[_cell(S, mg, dw) for mg in mg_vals] for dw in diw_vals])
        if relative:
            ref = _cell(S, REF_MG_SI, REF_DIW)
            g = g - ref if np.isfinite(ref) else g * np.nan
        grids[S] = g
    allv = np.concatenate([g[np.isfinite(g)].ravel() for g in grids.values()])
    if allv.size == 0:
        print("  no in-domain runs (or no reference cell) -- skipping map.")
        return

    unit = {'T': 'K', 'P_CO2': 'dex', 'pH': ''}[quantity]
    if relative:
        short = {'T': r'$\Delta T$', 'P_CO2': r'$\Delta\log P_{\mathrm{CO_2}}$',
                 'pH': r'$\Delta$pH'}[quantity]
        lim = float(np.nanpercentile(np.abs(allv), 98)) or 1.0
        norm = mcolors.TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
        cmap = RELATIVE_CMAP
        cbar_label = f'{short} vs Earth crust' + (f' ({unit})' if unit else '')
    else:
        qname = {'T': 'Temperature', 'P_CO2': '$P_{\\mathrm{CO_2}}$', 'pH': 'Ocean pH'}[quantity]
        norm = mcolors.Normalize(vmin=np.nanpercentile(allv, 2), vmax=np.nanpercentile(allv, 98))
        cmap = QUANTITY_CMAP if log_q else RELATIVE_CMAP
        cbar_label = qname + (f' ({unit})' if unit else '')
    print(f"Composition map [{quantity}, {'relative' if relative else 'absolute'}] "
          f"depth={ocean_depth:g}: {len(mg_vals)}x{len(diw_vals)} grid at S={s_present}"
          + (f"; dropped sparse Mg/Si {dropped[0]} dIW {dropped[1]}" if any(dropped) else ""))

    fig, axs = plt.subplots(1, len(s_present), figsize=figure_size(width, height), sharey=True,
                            squeeze=False)
    axs = axs[0]
    for ax, S in zip(axs, s_present):
        ax.pcolormesh(np.arange(len(mg_vals) + 1), np.arange(len(diw_vals) + 1), grids[S],
                      cmap=cmap, norm=norm, edgecolors='w', linewidth=0.4)
        for _, r in subset[subset['instellation'] == S].iterrows():
            if r['termination'] in OUT_OF_DOMAIN and r['mg_si'] in mg_vals and r['delta_iw'] in diw_vals:
                ax.plot(mg_vals.index(r['mg_si']) + 0.5, diw_vals.index(r['delta_iw']) + 0.5,
                        marker='x', color='0.4', markersize=5, mew=1.2)
        if relative and REF_MG_SI in mg_vals and REF_DIW in diw_vals:
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
    _save_fig(fig, figure_path(output_path,
                               f'sweep_composition_map_{quantity}{"" if relative else "_abs"}.png'))


def plot_damkohler_contour(df, output_path, out_targets=(0.1, 1.0, 10.0)):
    """log10(Da) over instellation x crust production, one panel per outgassing rate, with Da = 1 drawn."""
    subset = df[_sweep_mask(df, land=None, chem=False, f_ht=True)]
    if subset.empty:
        print("No data for Da contour plot — skipping.")
        return
    all_out = sorted(subset['outgassing'].unique())
    out_values = list(dict.fromkeys(min(all_out, key=lambda x: abs(x - t)) for t in out_targets))
    s_vals = sorted(subset['instellation'].unique())
    crust_vals = sorted(subset['crust_production'].unique())
    sel = _add_diag_columns(subset[subset['outgassing'].isin(out_values)], output_path)
    log_crust = np.log10(np.array(crust_vals))
    norm = mcolors.TwoSlopeNorm(vmin=-3.0, vcenter=0.0, vmax=3.0)
    levels = np.linspace(-3.0, 3.0, 31)

    with plt.rc_context({'figure.constrained_layout.use': False}):
        fig, axes = plt.subplots(len(out_values), 1, sharex=True, sharey=True, squeeze=False,
                                 figsize=figure_size('single', n_rows=len(out_values), row_height=2.0))
        for ax, out in zip(axes[:, 0], out_values):
            pivot = (sel[np.isclose(sel['outgassing'], out)]
                     .pivot_table(index='crust_production', columns='instellation', values='da',
                                  aggfunc='first')
                     .reindex(index=crust_vals, columns=s_vals))
            Z = np.ma.masked_invalid(np.log10(np.maximum(pivot.values.astype(float), 1e-10)))
            ax.contourf(s_vals, log_crust, Z, levels=levels, cmap=DIVERGING_CMAP, norm=norm,
                        extend='both')
            if not np.all(Z.mask if np.ma.is_masked(Z) else False):
                ax.contour(s_vals, log_crust, Z, levels=[0.0], colors='k', linewidths=1.5)
            ax.set_title(f'Outgassing = {out:g}×', fontsize=9)
            ax.set_yticks(log_crust)
            ax.set_yticklabels([f'{v:g}' for v in crust_vals], fontsize=7)
            ax.set_ylabel('Crust prod. (×Earth)')
            ax.grid(True, linestyle='--', alpha=0.3, color='k')
        axes[-1, 0].set_xlabel('Instellation (S/S₀)')

        fig.subplots_adjust(left=0.16, right=0.82, top=0.91, bottom=0.12, hspace=0.25)
        top, bottom = axes[0, 0].get_position(), axes[-1, 0].get_position()
        cbar_ax = fig.add_axes([0.85, bottom.y0, 0.03, top.y1 - bottom.y0])
        sm = plt.cm.ScalarMappable(cmap=DIVERGING_CMAP, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax, label=r'$\log_{10}(\mathrm{Damkohler Coefficient})$')
        cbar.ax.axhline(0, color='k', linewidth=1.5)
        _save_fig(fig, figure_path(output_path, 'da_contour.png'))


# --- Continental baseline, land-fraction series and weathering crossover ---
# The sweep design (land fractions, grid axes, Earth reference) is read from continental_baseline.py as cb.
ARM_COLOURS = {0.3: '#a4632a', 0.0: '#2a6fa4'}
ARM_LABELS = {0.3: 'Continental (land fraction 0.3)', 0.0: 'Ocean world (land free)'}
OLR_FIT_T_MAX = 350.0   # upper T limit of the Haqq-Misra et al. (2016) OLR fit the climate model uses

# (index into b_ocean, label, Earth seawater mmol/kg, calibrate_earth.py's targets). Al, Fe, SO4 omitted.
ION_SPEC = [(0, 'Alk', 2.3), (1, 'C', 2.1), (2, 'Si', 0.1), (5, 'Ca', 10.3),
            (6, 'Mg', 52.8), (7, 'Na', 469.0), (8, 'Cl', 546.0)]


def _arm(df, land, setup=None):
    """Runs at one land fraction with every other axis at the Earth reference.

    `setup` pins the Cl ratio and run length (e.g. the 'earth' sweep's); None takes the most-run values.
    """
    pins = setup or {}
    return df[_sweep_mask(df, depth=cb.OCEAN_DEPTH, land=land, f_ht=True, outgassing=cb.OUTGASSING,
                          crust_production=cb.CRUST_PRODUCTION, setup=not pins,
                          **pins)].sort_values('instellation')


def plot_continental_baseline(df, output_path, show_hz=None):
    """Earth-like continental baseline: T/pCO2 and pH/salinity against instellation, and ions against seawater.

    Drawn from continental_baseline's 'earth' sweep (seawater seed, Earth Cl ratio, 4 Gyr), the setup the
    constants are calibrated in, so the S = 1 run is the calibrated Earth.
    """
    group_all = _arm(df, cb.LAND_FRACTION, {'cl_ratio': cb.EARTH_CL_RATIO, 't_end_gyr': cb.EARTH_T_END_GYR})
    if group_all.empty:
        print("No 'earth' sweep runs (continental_baseline.py, SWEEP = 'earth') -- "
              "skipping the continental baseline figures.")
        return
    group_hab = group_all[(group_all['T'] > T_SNOWBALL) & (group_all['T'] < T_RUNAWAY)]

    for grp_cols, sfx in _panel_groups(True):
        fig, axes = plt.subplots(len(grp_cols), 1, sharex=True,
                                 figsize=figure_size('single', n_rows=len(grp_cols)))
        # The T panel shows every run; the others only runs inside the temperate band.
        if 'T' in grp_cols:
            _plot_group_on_axes([axes[grp_cols.index('T')]], group_all, color='k',
                                show_markers=False, cols=['T'])
        other = [c for c in grp_cols if c != 'T']
        if other:
            _plot_group_on_axes([axes[grp_cols.index(c)] for c in other], group_hab, color='k',
                                show_markers=False, cols=other)
        _style_axes(axes, grp_cols, show_hz=show_hz, show_eq_temp=True)
        for ax, col in zip(axes, grp_cols):
            ax.scatter(EARTH['S'], EARTH[col], marker='*', s=220, color='blue',
                       edgecolors='k', linewidths=0.7, zorder=6)
        if 'T' in grp_cols:
            axes[grp_cols.index('T')].annotate(
                'Earth', xy=(EARTH['S'], EARTH['T']), xytext=(EARTH['S'] + 0.06, EARTH['T'] - 6),
                fontsize=8, arrowprops=dict(arrowstyle='-', color='k', lw=0.8))
        _save_fig(fig, figure_path(output_path, f'continental_baseline{sfx}.png'))

    # Ions at the run nearest S = 1, model against Earth seawater, on a log axis.
    ion_rows = []
    for _, row in group_hab.iterrows():
        try:
            with open(os.path.join(RUN_PATH or output_path, f"{row['name']}.json")) as fh:
                y = json.load(fh)['data']['y']
            b = [max(float(y[2 + i][-1]), 1e-15) for i in range(len(elements))]
            ion_rows.append((row['instellation'], (b + [1e-15] * 10)[:10]))
        except Exception:
            pass

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.6))
    if ion_rows:
        s_arr = np.array([r[0] for r in ion_rows])
        b_model = np.array([r[1] for r in ion_rows])[int(np.argmin(np.abs(s_arr - EARTH['S'])))] * 1e3
        earth = np.array([spec[2] for spec in ION_SPEC])
        model = np.array([max(b_model[spec[0]], 1e-4) for spec in ION_SPEC])
        x = np.arange(len(ION_SPEC))
        ax.vlines(x, np.minimum(earth, model), np.maximum(earth, model), color='0.6',
                  linewidth=1.0, zorder=2)
        ax.scatter(x, earth, marker='o', s=46, facecolors='none', edgecolors='#4C72B0',
                   linewidths=1.4, zorder=3, label='Earth seawater')
        ax.scatter(x, model, marker='D', s=34, color='#C44E52', edgecolors='w',
                   linewidths=0.5, zorder=4, label='Model (calibrated)')
        ax.set_xticks(x)
        ax.set_xticklabels([spec[1] for spec in ION_SPEC])
        ax.axvline(2.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1)
        trans = ax.get_xaxis_transform()
        ax.text(1.0, 1.02, 'Biotically controlled', transform=trans, ha='center', va='bottom',
                fontsize=8)
        ax.text(4.5, 1.02, 'Abiotically controlled', transform=trans, ha='center', va='bottom',
                fontsize=8)
        ax.set_xlim(-0.6, len(ION_SPEC) - 0.4)
    ax.set_yscale('log')
    ax.set_ylabel('Concentration (mmol kg$^{-1}$)')
    ax.set_ylim(0.05, 1.0e3)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(True, linestyle='--', alpha=0.35, axis='y', which='major')
    ax.legend(frameon=False, fontsize=7, loc='lower right', handletextpad=0.4, borderaxespad=0.6)
    _save_fig(fig, figure_path(output_path, 'continental_baseline_ions.png'))


def _olr_limit(pco2_bar):
    """First local maximum of OLR(T): the Simpson-Nakajima radiation limit for this atmosphere."""
    from kamino.climate.analytic import OLR
    peak = OLR(180.0, pco2_bar)
    for T in np.arange(181.0, 391.0, 1.0):
        v = OLR(float(T), pco2_bar)
        if v < peak:
            break
        peak = v
    return peak


def _past_runaway(S, pco2_bar, albedo=0.3):
    """True when absorbed instellation exceeds the OLR limit, so the stored T is on the hot branch."""
    from kamino.climate.analytic import albedo_funtion
    pco2_bar = max(float(pco2_bar), 1e-5)      # the model's 1 Pa CO2 floor
    A = albedo_funtion(pco2_bar, albedo)
    return S * SOLAR_CONSTANT * (1 - A) * 0.25 > _olr_limit(pco2_bar)


def _draw_arm(axes, group, colour, cols):
    """Draw one arm and return whether any run is habitable; a never-habitable arm is drawn faint with hollow markers."""
    if group['termination'].isin(HABITABLE).any():
        _plot_group_on_axes(axes, group, colour, show_markers=False, cols=cols)
        return True
    clamped = (group['termination'].isin(OUT_OF_DOMAIN) &
               ((group['T'] >= T_HOT_WALL) | (group['T'] <= T_COLD_WALL)))
    shown = group[~clamped]
    for ax, col in zip(axes, cols):
        ax.plot(shown['instellation'], shown[col], color=colour, linewidth=1.2, alpha=0.5, zorder=2)
        for _, row in shown.iterrows():
            if np.isfinite(row[col]):
                ax.scatter(row['instellation'], row[col],
                           marker=FAILED_MARKERS.get(row['termination'], 'x'), s=22,
                           facecolors='none', edgecolors=colour, linewidths=1.0, zorder=4)
    return False


def _arm_handles(arms, habitable):
    """Legend handles for the land arms, faint with a marker when an arm is never habitable."""
    return [Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6, alpha=1.0 if habitable[l] else 0.5,
                   marker='' if habitable[l] else 's', markerfacecolor='none',
                   label=ARM_LABELS[l] + ('' if habitable[l] else ' — never habitable'))
            for l in arms]


def _alk_fluxes(group):
    """(continental, seafloor) alkalinity flux in Tmol eq/yr, both over the area the model applies them to.

    The stored seafloor flux is normalised on Earth's fixed seafloor area, so it is rescaled by
    (1 - land_fraction) / EARTH_OCEAN_FRACTION.
    """
    from kamino.weathering import get_continental_weathering_flux
    from kamino.chemistry import alk_idx
    from kamino.constants import YR, R_EARTH, EARTH_OCEAN_FRACTION
    surface = 4 * np.pi * R_EARTH ** 2
    cont = np.full(len(group), np.nan)
    for i, (_, r) in enumerate(group.iterrows()):
        T, p, land = r['T'], r['P_CO2'], r['land_fraction']
        if land <= 0:
            cont[i] = 0.0
        elif np.isfinite(T) and np.isfinite(p):
            f = get_continental_weathering_flux(float(T), float(p) * 1e5)   # pCO2 stored in bar
            cont[i] = float(f[alk_idx]) * land * surface * YR / 1e12
    sea = (group['alk_flux'].to_numpy(dtype=float)
           * (1.0 - group['land_fraction'].to_numpy(dtype=float)) / EARTH_OCEAN_FRACTION)
    return cont, sea


def _crossover_land_fraction(lands, ratios):
    """Land fraction where the flux ratio (either way up) crosses 1, log-interpolated; None if not bracketed."""
    pairs = sorted((float(l), float(r)) for l, r in zip(lands, ratios)
                   if l > 0 and np.isfinite(r) and r > 0)
    for (l0, r0), (l1, r1) in zip(pairs, pairs[1:]):
        if (r0 - 1.0) * (r1 - 1.0) <= 0 and r0 != r1:
            w = -np.log10(r0) / (np.log10(r1) - np.log10(r0))
            return float(10 ** (np.log10(l0) + w * (np.log10(l1) - np.log10(l0))))
    return None


def _crossings(ratio):
    """{S: crossover land fraction} from {(S, land): flux ratio}, for instellations that bracket one."""
    out = {}
    for s in sorted({s for s, _ in ratio}):
        lands = sorted(l for (s2, l) in ratio if s2 == s)
        got = _crossover_land_fraction(lands, [ratio[(s, l)] for l in lands])
        if got is not None:
            out[s] = got
    return out


def hz_edges(group):
    """(S_outer, S_inner, outer_kind, inner_kind) of the habitable band on one line, or None.

    kind is 'crossing' (interpolated between trustworthy runs), 'bracketed' (midpoint of the step to
    a run at the matching wall) or 'open' (last habitable point; the true edge lies further out).
    """
    g = group.sort_values('instellation')
    S = g['instellation'].to_numpy(dtype=float)
    T = g['T'].to_numpy(dtype=float)
    wall = (g['domain_wall'].to_numpy(dtype=object) if 'domain_wall' in g
            else np.full(len(g), None, dtype=object))
    trusted = g['termination'].isin(HABITABLE).to_numpy() & np.isfinite(T)
    # Past the runaway or above the OLR fit's range, a T inside the window is still not habitable.
    S_ok = np.array([not _past_runaway(s, p) for s, p in zip(S, g['P_CO2'].to_numpy(dtype=float))])
    hab = trusted & (T > T_SNOWBALL) & (T < T_RUNAWAY) & (T <= OLR_FIT_T_MAX) & S_ok
    if not hab.any():
        return None
    idx = np.flatnonzero(hab)
    i0, i1 = int(idx[0]), int(idx[-1])

    outer, outer_kind = S[i0], 'open'
    if i0 > 0:
        if trusted[i0 - 1] and T[i0 - 1] <= T_SNOWBALL:
            outer = float(np.interp(T_SNOWBALL, [T[i0 - 1], T[i0]], [S[i0 - 1], S[i0]]))
            outer_kind = 'crossing'
        elif wall[i0 - 1] == 'cold':
            outer, outer_kind = 0.5 * (S[i0 - 1] + S[i0]), 'bracketed'

    inner, inner_kind = S[i1], 'open'
    if i1 + 1 < len(S):
        if trusted[i1 + 1] and T[i1 + 1] >= T_RUNAWAY:
            inner = float(np.interp(T_RUNAWAY, [T[i1], T[i1 + 1]], [S[i1], S[i1 + 1]]))
            inner_kind = 'crossing'
        elif wall[i1 + 1] == 'hot' or not S_ok[i1 + 1]:
            inner, inner_kind = 0.5 * (S[i1] + S[i1 + 1]), 'bracketed'
    return float(outer), float(inner), outer_kind, inner_kind


def plot_baseline_vs_ocean(arms, output_path):
    """Continental arm against the ocean-world arm: T, pCO2, pH and salinity against instellation."""
    habitable = {land: group['termination'].isin(HABITABLE).any() for land, group in arms.items()}
    handles = _arm_handles(arms, habitable) + list(DA_LEGEND)
    for cols, sfx in _panel_groups(True):
        fig, axes = plt.subplots(len(cols), 1, sharex=True,
                                 figsize=figure_size('single', n_rows=len(cols), row_height=2.0))
        for land, group in arms.items():
            _draw_arm(axes, group, ARM_COLOURS[land], cols)
        _style_axes(axes, cols)
        for ax, col in zip(axes, cols):
            ax.scatter(EARTH['S'], EARTH[col], marker='*', s=180, color='gold',
                       edgecolors='k', linewidths=0.7, zorder=6)
        _add_figure_legend(fig, axes, handles)
        _save_fig(fig, figure_path(output_path, f'continental_vs_ocean{sfx}.png'))


def plot_habitable_zone(arms, output_path):
    """Temperature curves with the habitable zone of each arm as a bar; carets mark edges that are only bounds."""
    edges = {land: e for land, group in arms.items() if (e := hz_edges(group)) is not None}
    if cb.LAND_FRACTION not in edges:
        print("No habitable band on the continental arm -- skipping the habitable-zone figure.")
        return edges

    fig, (ax, ax_z) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 1],
                                   figsize=figure_size('single', height=4.0))
    habitable = {land: _draw_arm([ax], group, ARM_COLOURS[land], ['T'])
                 for land, group in arms.items()}
    _style_axes([ax], ['T'])
    ax.set_xlabel('')
    ax.scatter(EARTH['S'], EARTH['T'], marker='*', s=180, color='gold', edgecolors='k',
               linewidths=0.7, zorder=6)

    for row, land in enumerate(arms):   # every arm gets a row, even one with no habitable band
        colour, y = ARM_COLOURS[land], len(arms) - 1 - row
        if land not in edges:
            ax_z.text(0.5 * sum(ax.get_xlim()), y, 'no habitable zone', ha='center',
                      va='center', fontsize=7, color=colour, style='italic', zorder=6)
            continue
        lo, hi, lo_kind, hi_kind = edges[land]
        ax_z.barh(y, hi - lo, left=lo, height=0.5, color=colour, alpha=0.35,
                  edgecolor=colour, linewidth=1.4, zorder=3)
        for x, kind, marker in ((lo, lo_kind, '<'), (hi, hi_kind, '>')):
            if kind == 'open':
                ax_z.scatter(x, y, marker=marker, s=30, color=colour, zorder=5)
        ax_z.text(0.5 * (lo + hi), y, f'{lo:.2f}–{hi:.2f}', ha='center', va='center',
                  fontsize=7, zorder=6)
        ax.axvspan(lo, hi, color=colour, alpha=0.07, zorder=0)

    ax_z.set_yticks(range(len(arms)))
    ax_z.set_yticklabels([])
    ax_z.set_ylim(-0.6, len(arms) - 0.4)
    ax_z.set_ylabel('Habitable\nzone')
    ax_z.set_xlabel('Instellation (S/S₀)')
    ax_z.grid(True, axis='x', linestyle='--', alpha=0.4, zorder=0)
    ax_z.set_xlim(*ax.get_xlim())

    handles = _arm_handles(arms, habitable)
    if any(k == 'open' for e in edges.values() for k in e[2:]):
        handles.append(Line2D([0], [0], color='k', linestyle='none', marker='>', markersize=5,
                              label='Edge is a bound (sweep limit)'))
    _add_figure_legend(fig, [ax, ax_z], handles)
    _save_fig(fig, figure_path(output_path, 'continental_habitable_zone.png'))
    return edges


def _report(arms, edges):
    """Print the habitable-zone edges and warn if CONTINENTAL_HZ_* no longer match them."""
    print(f"\nHabitable zone (T between {T_SNOWBALL:.0f} K and {T_RUNAWAY:.0f} K), "
          f"{cb._pe_label(REF_PE)} ocean, {cb.OCEAN_DEPTH/1000:g} km, outgassing "
          f"{cb.OUTGASSING:g}x, crust {cb.CRUST_PRODUCTION:g}x Earth:")
    print(f"  {'arm':>32s} {'outer S':>8s} {'inner S':>8s} {'width':>7s}   how located")
    for land in arms:
        if land not in edges:
            walls = arms[land]['domain_wall'].dropna().value_counts()
            why = ', '.join(f'{n} {WALL_LABELS.get(w, w)}' for w, n in walls.items())
            print(f"  {ARM_LABELS[land]:>32s} {'--':>8} {'--':>8} {'none':>7}"
                  f"   no habitable run ({why})")
            continue
        lo, hi, lo_kind, hi_kind = edges[land]
        print(f"  {ARM_LABELS[land]:>32s} {lo:8.3f} {hi:8.3f} {hi - lo:7.3f}"
              f"   outer {lo_kind}, inner {hi_kind}")
    if len(edges) == 2:
        (a_lo, a_hi, _, _), (b_lo, b_hi, _, _) = edges[cb.LAND_FRACTION], edges[0.0]
        print(f"  continental zone is {(a_hi - a_lo) - (b_hi - b_lo):+.3f} S wide relative to "
              f"the ocean world ({(a_hi - a_lo) / (b_hi - b_lo):.2f}x)")
    if cb.LAND_FRACTION in edges:
        lo, hi = edges[cb.LAND_FRACTION][:2]
        for label, here, there in (('CONTINENTAL_HZ_OUTER', lo, CONTINENTAL_HZ_OUTER),
                                   ('CONTINENTAL_HZ_INNER', hi, CONTINENTAL_HZ_INNER)):
            if abs(here - there) > 5e-4:
                print(f"  NOTE plot_results.{label} = {there:g}, but this sweep measures "
                      f"{here:.3f}. Update it, or the HZ lines on every other figure are stale.")


def _land_series(df, output_path):
    """{land_fraction: runs with diagnostics} for the land-fraction series on disk."""
    series = {}
    for land in cb.LAND_FRACTIONS:
        sub = _arm(df, land)
        if not sub.empty:
            series[land] = _add_diag_columns(sub, output_path).sort_values('instellation')
    return series


def plot_land_fraction_series(series, output_path):
    """T, pCO2, pH and salinity against instellation, one line per land fraction (0 in the ocean-world blue)."""
    if len(series) < 2:
        print("Fewer than two land fractions on disk -- skipping the land-fraction series.")
        return
    positive = sorted(l for l in series if l > 0)
    norm = mcolors.LogNorm(vmin=min(positive), vmax=max(positive)) if len(positive) > 1 else None
    colours = {l: LAND_FRACTION_CMAP(norm(l)) if norm else ARM_COLOURS[cb.LAND_FRACTION]
               for l in positive}
    colours[0.0] = ARM_COLOURS[0.0]

    for cols, sfx in _panel_groups(True):
        fig, axes = plt.subplots(len(cols), 1, sharex=True,
                                 figsize=figure_size('single', n_rows=len(cols), row_height=2.0))
        for land, group in sorted(series.items(), reverse=True):
            _draw_arm(axes, group, colours[land], cols)
        _style_axes(axes, cols)
        if norm is not None:
            _add_colorbar(fig, list(axes), LAND_FRACTION_CMAP, norm, 'Land fraction', ticks=positive,
                          ticklabels=[f'{v:g}' for v in positive], aspect=len(cols) * 7.5)
        handles = ([Line2D([0], [0], color=ARM_COLOURS[0.0], linewidth=1.6, label='Land free (0)')]
                   if 0.0 in series else []) + list(DA_LEGEND)
        _add_figure_legend(fig, axes, handles)
        _save_fig(fig, figure_path(output_path, f'land_fraction_series{sfx}.png'))


def plot_weathering_crossover(series, output_path):
    """Continental and seafloor alkalinity flux against land fraction near S = 1, and the crossover against S.

    Only steady states (converged or 2 Gyr) are used; a run stopped at a wall is still evolving.
    """
    if len(series) < 2:
        print("Fewer than two land fractions on disk -- skipping the crossover figure.")
        return None
    rows, dropped = {}, 0   # (land, S) -> (continental, seafloor)
    for land, group in series.items():
        steady = group[group['termination'].isin(HABITABLE)]
        dropped += len(group) - len(steady)
        if steady.empty:
            continue
        for s, c, f in zip(steady['instellation'], *_alk_fluxes(steady)):
            if np.isfinite(c) and np.isfinite(f) and f > 0:
                rows[(float(land), float(s))] = (float(c), float(f))
    if not rows:
        print("No steady-state runs to compare fluxes on -- skipping the crossover figure.")
        return None
    if dropped:
        print(f"  crossover: ignoring {dropped} run(s) that never reached a steady state.")

    s_vals = sorted({s for _, s in rows})
    crossings = _crossings({(s, l): c / f for (l, s), (c, f) in rows.items()})
    s_ref = min(s_vals, key=lambda s: abs(s - EARTH['S']))
    fig, (ax, ax_c) = plt.subplots(1, 2, figsize=figure_size('double', height=2.8))

    lands_ref = sorted(l for (l, s) in rows if s == s_ref)
    if lands_ref:
        cont = [rows[(l, s_ref)][0] for l in lands_ref]
        sea = [rows[(l, s_ref)][1] for l in lands_ref]
        x = [max(l, 1e-4) for l in lands_ref]      # 0 has no place on a log axis
        ax.plot(x, cont, color=ARM_COLOURS[0.3], marker='o', markersize=3, linewidth=1.6,
                label='Continental')
        ax.plot(x, sea, color=ARM_COLOURS[0.0], marker='s', markersize=3, linewidth=1.6,
                label='Seafloor (LT)')
        if s_ref in crossings:
            ax.axvline(crossings[s_ref], color='0.35', linestyle=(0, (6, 3)), linewidth=1.0)
            ax.annotate(f'{crossings[s_ref]:.3g}', xy=(crossings[s_ref], max(cont)),
                        xytext=(3, -2), textcoords='offset points', fontsize=7)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Land fraction')
    ax.set_ylabel('Alkalinity flux (Tmol eq/yr)')
    ax.set_title(f'S = {s_ref:g}', fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=7, frameon=False)
    if crossings:
        ax_c.plot(list(crossings), list(crossings.values()), color='k', marker='o', markersize=3,
                  linewidth=1.4)
    ax_c.set_yscale('log')
    ax_c.set_xlabel('Instellation (S/S₀)')
    ax_c.set_ylabel('Crossover land fraction')
    ax_c.grid(True, linestyle='--', alpha=0.4)
    _save_fig(fig, figure_path(output_path, 'weathering_crossover.png'))

    print("\nContinental vs seafloor alkalinity flux (Tmol eq/yr), steady states only:")
    print(f"  {'land':>8} " + ' '.join(f'{s:>9.2f}' for s in s_vals))
    for land in sorted(series, reverse=True):
        cells = [f"{rows[(land, s)][0] / rows[(land, s)][1]:9.3g}" if (land, s) in rows
                 else f"{'--':>9}" for s in s_vals]
        print(f"  {land:8g} " + ' '.join(cells))
    print("  (continental / seafloor; < 1 means seafloor weathering dominates)")
    if crossings:
        print(f"  crossover land fraction: {min(crossings.values()):.3g} to "
              f"{max(crossings.values()):.3g} over S = {min(crossings):g}-{max(crossings):g}")
    else:
        ratios_all = [c / f for c, f in rows.values()]
        print(f"  no crossing inside the sampled land fractions -- ratio spans "
              f"{min(ratios_all):.3g} to {max(ratios_all):.3g}; "
              f"extend continental_baseline.LAND_FRACTIONS to bracket it.")
    return crossings


def _ratio_cells(series, output_path):
    """(ratio, not_steady, net_sink) for one panel: seafloor/continental flux keyed by (S, land).

    not_steady runs never reached a steady state; net_sink runs are steady with a negative seafloor flux.
    """
    ratio, not_steady, net_sink = {}, [], []
    for land, group in series.items():
        if land <= 0 or group.empty:
            continue
        group = _add_diag_columns(group, output_path)
        steady = group[group['termination'].isin(HABITABLE)]
        names = set(steady['name'])
        not_steady += [(float(r['instellation']), land) for _, r in group.iterrows()
                       if r['name'] not in names]
        if steady.empty:
            continue
        for s, c, f in zip(steady['instellation'], *_alk_fluxes(steady)):
            if np.isfinite(c) and np.isfinite(f) and c > 0 and f > 0:
                ratio[(float(s), land)] = f / c
            elif np.isfinite(f) and f <= 0:
                net_sink.append((float(s), land))
    return ratio, not_steady, net_sink


def _ratio_scale(log_values, step=0.5):
    """(bands, norm, ticks): half-decade bands and a norm diverging about ratio 1, not forced symmetric."""
    lo = np.floor(min(log_values) / step) * step
    hi = np.ceil(max(log_values) / step) * step
    norm = mcolors.TwoSlopeNorm(vmin=min(lo, -step), vcenter=0.0, vmax=max(hi, step))
    return np.arange(lo, hi + 0.5 * step, step), norm, [t for t in range(-9, 10) if lo <= t <= hi]


def _ratio_panel(ax, cells, bands, norm, lands=None, marker_size=3, line_width=1.2):
    """Contour one (ratio, not_steady, net_sink) panel; returns (contourf, ratio = 1 contour) or Nones."""
    ratio, not_steady, net_sink = cells
    s_vals = sorted({s for s, _ in ratio})
    lands = sorted({l for _, l in ratio}) if lands is None else lands
    cf = cs = None
    if len(s_vals) > 1 and len(lands) > 1:
        Z = np.full((len(lands), len(s_vals)), np.nan)
        for (s, l), v in ratio.items():
            Z[lands.index(l), s_vals.index(s)] = np.log10(v)
        Zm = np.ma.masked_invalid(Z)
        cf = ax.contourf(s_vals, lands, Zm, levels=bands, cmap=WEATHERING_RATIO_CMAP, norm=norm,
                         extend='both')
        if np.nanmin(Z) < 0 < np.nanmax(Z):
            cs = ax.contour(s_vals, lands, Zm, levels=[0.0], colors='k', linewidths=line_width)
    else:
        ax.text(0.5, 0.5, 'no steady state', transform=ax.transAxes, ha='center', va='center',
                fontsize=7, color='0.5', style='italic')
    for s, l in not_steady:
        ax.plot(s, l, marker='x', color='0.45', markersize=marker_size, mew=0.8, zorder=4)
    for s, l in net_sink:
        ax.plot(s, l, marker='o', markerfacecolor='none', markeredgecolor='0.25',
                markersize=marker_size + 1, mew=0.9, zorder=4)
    ax.set_yscale('log')
    return cf, cs


def _ratio_colorbar(fig, ax, cf, ticks, aspect):
    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=aspect, ticks=ticks)
    cbar.set_label('Seafloor / continental alkalinity flux')
    cbar.set_ticklabels(['1' if t == 0 else f'$10^{{{t}}}$' for t in ticks])
    if 0 in ticks:
        cbar.ax.axhline(0, color='k', linewidth=1.2)


def _ratio_facets(cells, rows, cols, row_label, col_label, scale, path, title=None, height=5.0):
    """Grid of ratio panels, largest row value on top, all on one shared colour scale."""
    bands, norm, ticks = scale
    fig, axes = plt.subplots(len(rows), len(cols), sharex=True, sharey=True, squeeze=False,
                             figsize=figure_size('double', height=height))
    cf = None
    for i, rv in enumerate(reversed(rows)):
        for j, cv in enumerate(cols):
            ax = axes[i, j]
            got, _ = _ratio_panel(ax, cells.get((rv, cv), ({}, [], [])), bands, norm)
            cf = got if got is not None else cf
            ax.grid(True, linestyle='--', alpha=0.3)
            if i == 0:
                ax.set_title(col_label(cv), fontsize=8)
            if j == 0:
                lbl = row_label(rv)
                ax.set_ylabel(f'{lbl}\nLand fraction' if lbl else 'Land fraction', fontsize=7)
            if i == len(rows) - 1 and j == len(cols) // 2:
                ax.set_xlabel('Instellation (S/S₀)')
    if cf is not None:
        _ratio_colorbar(fig, list(axes.ravel()), cf, ticks, aspect=30)
    if title:
        fig.suptitle(title, fontsize=9)
    _save_fig(fig, path)


def plot_weathering_ratio_map(series, output_path):
    """Seafloor / continental alkalinity flux over instellation x land fraction (steady states only)."""
    lands = sorted(l for l in series if l > 0)
    if len(lands) < 2:
        print("Fewer than two positive land fractions -- skipping the ratio map.")
        return None
    cells = _ratio_cells(series, output_path)
    ratio, not_steady = cells[0], cells[1]
    if not ratio:
        print("No steady-state runs with both fluxes positive -- skipping the ratio map.")
        return None
    s_vals = sorted({s for s, _ in ratio})
    bands, norm, ticks = _ratio_scale([np.log10(v) for v in ratio.values()])

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.9))
    cf, cs = _ratio_panel(ax, cells, bands, norm, lands=lands, marker_size=3.5, line_width=1.4)
    if cs is not None:
        ax.clabel(cs, fmt={0.0: 'equal'}, fontsize=7, inline=True)
    ax.set_xlabel('Instellation (S/S₀)')
    ax.set_ylabel('Land fraction')
    dx = 0.02 * (max(s_vals) - min(s_vals))   # margins so edge markers are not cut by the frame
    dy = 0.04 * (np.log10(max(lands)) - np.log10(min(lands)))
    ax.set_xlim(min(s_vals) - dx, max(s_vals) + dx)
    ax.set_ylim(10 ** (np.log10(min(lands)) - dy), 10 ** (np.log10(max(lands)) + dy))
    _ratio_colorbar(fig, ax, cf, ticks, aspect=22)
    _save_fig(fig, figure_path(output_path, 'weathering_ratio_map.png'))

    print("\nSeafloor / continental alkalinity flux (steady states only):")
    print(f"  {'land':>8} " + ' '.join(f'{s:>8.2f}' for s in s_vals))
    for land in reversed(lands):
        print(f"  {land:8g} " + ' '.join(f"{ratio[(s, land)]:8.2g}" if (s, land) in ratio
                                          else f"{'--':>8}" for s in s_vals))
    if not_steady:
        print(f"  ({len(not_steady)} cell(s) blank: no steady state)")
    return ratio


def _coarse_mask(df):
    """Runs on the coarse grid's instellation and land-fraction values."""
    return (df['instellation'].isin(cb.COARSE_INSTELLATION) &
            df['land_fraction'].apply(lambda v: any(np.isclose(v, l) for l in cb.COARSE_LAND_FRACTIONS)))


def _by_land(sub):
    return {float(l): sub[np.isclose(sub['land_fraction'], l)].sort_values('instellation')
            for l in sorted(sub['land_fraction'].unique())}


def _grid_slice(df, outgassing, crust, mg_si):
    """{land_fraction: runs} for one (outgassing, crust, Mg/Si) cell of the coarse grid."""
    return _by_land(df[_coarse_mask(df) & _sweep_mask(
        df, depth=cb.OCEAN_DEPTH, land=None, crust=False, f_ht=True, outgassing=outgassing,
        crust_production=crust, mg_si=mg_si, delta_iw=cb.DELTA_IW)])


def _alpha_slice(df, outgassing, alpha):
    """{land_fraction: runs} for one (outgassing, alpha) cell; alpha is pinned here, not by _ref_chem."""
    return _by_land(df[_coarse_mask(df) & _sweep_mask(
        df, depth=cb.OCEAN_DEPTH, land=None, crust=False, chem=False, f_ht=True,
        outgassing=outgassing, crust_production=cb.CRUST_PRODUCTION, mg_si=cb.MG_SI_EARTH,
        delta_iw=cb.DELTA_IW, alpha=alpha, kd_mg=cb.KD_MG_CALIB, k_na=cb.K_NA_CALIB)])


def plot_weathering_ratio_grid(df, output_path):
    """Ratio map faceted over crust production x outgassing, one figure per Mg/Si, on one shared scale.

    Mg/Si and crust production reach only the seafloor sink; continental weathering sees them via climate alone.
    """
    mg_vals = [m for m in cb.GRID_MG_SI if np.isclose(df['mg_si'], m).any()]
    if not mg_vals:
        print("No grid runs on disk -- skipping the faceted ratio map.")
        return None
    cells = {}
    for mg in mg_vals:
        for c in cb.GRID_CRUST:
            for o in cb.GRID_OUTGASSING:
                series = _grid_slice(df, o, c, mg)
                if series:
                    cells[(mg, c, o)] = _ratio_cells(series, output_path)
    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state grid runs with both fluxes positive -- skipping.")
        return None
    scale = _ratio_scale(allv)
    for mg in mg_vals:
        _ratio_facets({(c, o): v for (m, c, o), v in cells.items() if m == mg},
                      cb.GRID_CRUST, cb.GRID_OUTGASSING, lambda c: f'crust {c:g}x',
                      lambda o: f'outgassing {o:g}x', scale,
                      figure_path(output_path, f'weathering_ratio_grid_mgsi{mg:g}.png'),
                      f'Mantle Mg/Si = {mg:g}')

    print("\nCrossover land fraction across the grid (steady states only):")
    print(f"  {'Mg/Si':>6} {'crust':>7} {'out':>6}   crossover (by instellation)")
    for (mg, c, o), (ratio, *_) in sorted(cells.items()):
        pts = list(_crossings(ratio).values())
        span = (f"{min(pts):.2g}-{max(pts):.2g}" if pts else
                ("none in range" if ratio else "no steady state"))
        print(f"  {mg:6g} {c:7g} {o:6g}   {span}")
    return cells


def plot_alpha_scaling(df, output_path):
    """Crossover land fraction against alpha per outgassing rate, with the alpha^1 kinetic-limit reference."""
    rows = []
    for o in cb.GRID_OUTGASSING:
        for a in cb.GRID_ALPHA:
            series = _alpha_slice(df, o, a)
            if series:
                rows += [{'outgassing': o, 'alpha': a, 'instellation': s, 'f_star': f}
                         for s, f in _crossings(_ratio_cells(series, output_path)[0]).items()]
    if not rows:
        print("No crossovers found across the alpha grid -- skipping.")
        return None
    tab = pd.DataFrame(rows)

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=3.0))
    norm = mcolors.LogNorm(vmin=min(cb.GRID_OUTGASSING), vmax=max(cb.GRID_OUTGASSING))
    print("\nCrossover land fraction f* against alpha (geometric mean over instellation):")
    print(f"  {'outgassing':>10} " + ' '.join(f'{a:>10.4g}' for a in cb.GRID_ALPHA)
          + f" {'exponent':>9}")
    exponents, anchor = {}, None
    for o in cb.GRID_OUTGASSING:
        g = tab[tab.outgassing == o]
        if g.empty:
            continue
        # Geometric mean: f* spans decades.
        means = {a: float(np.exp(np.log(g[g.alpha == a].f_star).mean()))
                 for a in cb.GRID_ALPHA if (g.alpha == a).any()}
        slope = float('nan')
        if len(means) >= 2:
            slope = exponents[o] = float(np.polyfit(np.log10(list(means)),
                                                    np.log10(list(means.values())), 1)[0])
        print(f"  {o:10g} " + ' '.join(f'{means[a]:10.4g}' if a in means else f'{"--":>10}'
                                       for a in cb.GRID_ALPHA) + f" {slope:9.2f}")
        ax.plot(list(means), list(means.values()), marker='o', markersize=4,
                color=OUTGASSING_CMAP(norm(o)), linewidth=1.6, label=f'{o:g}x')
        if anchor is None and means:
            anchor = (min(means), means[min(means)])
    if anchor is not None:
        xs = np.array(cb.GRID_ALPHA, dtype=float)
        ax.plot(xs, anchor[1] * xs / anchor[0], color='0.4', linestyle=(0, (6, 3)), linewidth=1.2,
                label=r'$\propto \alpha$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Reactive area scaling $\alpha$')
    ax.set_ylabel(r'Crossover land fraction $f^*$')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=7, frameon=False, title='Outgassing', title_fontsize=7)
    _save_fig(fig, figure_path(output_path, 'alpha_scaling.png'))
    if exponents:
        print(f"  exponent d log f* / d log alpha: {min(exponents.values()):.2f} to "
              f"{max(exponents.values()):.2f} (alpha^1 would be 1.00)")
    tab.to_csv(os.path.join(output_path, 'alpha_crossover.csv'), index=False)
    return tab


def plot_alpha_ratio_grid(df, output_path, outgassing=1.0):
    """Ratio map at each reactive-area scaling alpha, one panel per value, at fixed outgassing.

    Outgassing is held fixed because it moves the crossover the same way alpha does, so sweeping
    both adds no information; the panels share axes and one colour scale.
    """
    cells = {}
    for a in cb.GRID_ALPHA:
        series = _alpha_slice(df, outgassing, a)
        if series:
            cells[(outgassing, a)] = _ratio_cells(series, output_path)
    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state alpha runs with both fluxes positive -- skipping.")
        return None
    _ratio_facets(cells, [outgassing], cb.GRID_ALPHA, lambda o: '',
                  lambda a: rf'$\alpha$ = {a:g}', _ratio_scale(allv),
                  figure_path(output_path, 'weathering_ratio_alpha_grid.png'), height=2.9)
    return cells


def plot_outgassing_ratio_grid(df, output_path, alpha=None):
    """The companion of plot_alpha_ratio_grid: one panel per outgassing rate at fixed alpha."""
    alpha = cb.GRID_ALPHA[0] if alpha is None else alpha
    cells = {}
    for o in cb.GRID_OUTGASSING:
        series = _alpha_slice(df, o, alpha)
        if series:
            cells[(alpha, o)] = _ratio_cells(series, output_path)
    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state outgassing runs with both fluxes positive -- skipping.")
        return None
    _ratio_facets(cells, [alpha], cb.GRID_OUTGASSING, lambda a: '',
                  lambda o: f'Outgassing = {o:g}×', _ratio_scale(allv),
                  figure_path(output_path, 'weathering_ratio_outgassing_grid.png'), height=2.9)
    return cells


def plot_continental(df, output_path):
    """Every continental figure: baseline vs ocean, habitable zone, land series and ratio maps."""
    arms = {}
    for land in cb.LAND_ARMS:
        sub = _arm(df, land)
        if not sub.empty:
            arms[land] = _add_diag_columns(sub, output_path).sort_values('instellation')
    if cb.LAND_FRACTION not in arms:
        print(f"No runs at land_fraction = {cb.LAND_FRACTION:g} -- skipping the continental figures.")
        return
    for land, group in arms.items():
        print(f"  {len(group)} run(s) on the {ARM_LABELS[land].lower()} arm.")

    plot_baseline_vs_ocean(arms, output_path)
    _report(arms, plot_habitable_zone(arms, output_path) or {})
    series = _land_series(df, output_path)
    if len(series) > len(cb.LAND_ARMS):
        print(f"\n  land fractions on disk: {sorted(series, reverse=True)}")
    plot_land_fraction_series(series, output_path)
    plot_weathering_crossover(series, output_path)
    plot_weathering_ratio_map(series, output_path)
    plot_weathering_ratio_grid(df, output_path)
    plot_alpha_scaling(df, output_path)
    plot_alpha_ratio_grid(df, output_path)
    plot_outgassing_ratio_grid(df, output_path)
    plot_continental_baseline(df, output_path)


# --- Reactive area (alpha), temperature planes, Da transition and climate-state map ---
ALPHA_MARKERS = ('o', 's', 'D', '^', 'v')   # alphas share the instellation colours, so shape tells them apart


def _alpha_collapse_exponent(sub, exponents=np.arange(-1.0, 2.01, 0.05)):
    """(best p, rms, rms by p) for collapsing T onto one quadratic in log10(outgassing / alpha**p)."""
    groups = [g for _, g in sub.groupby('instellation') if g['alpha'].nunique() > 1 and len(g) >= 4]
    if not groups:
        return np.nan, np.nan, {}
    rms = {}
    for p in exponents:
        total, n = 0.0, 0
        for g in groups:
            x = np.log10(g['outgassing'].to_numpy() / g['alpha'].to_numpy() ** p)
            y = g['T'].to_numpy()
            total += float(np.sum((y - np.polyval(np.polyfit(x, y, 2), x)) ** 2))
            n += len(y)
        rms[round(float(p), 3)] = np.sqrt(total / n)
    best = min(rms, key=rms.get)
    return best, rms[best], rms


def _alpha_combination_label(p):
    """Axis label for outgassing / alpha**p."""
    if np.isclose(p, 1.0):
        return r'Outgassing / $\alpha$ (×Earth)'
    if np.isclose(p, -1.0):
        return r'Outgassing $\times\ \alpha$ (×Earth)'
    if np.isclose(p, 0.0):
        return r'Outgassing (×Earth)'
    return rf'Outgassing / $\alpha^{{{p:g}}}$ (×Earth)'


def plot_alpha_outgassing(df, output_path, mg_si=REF_MG_SI, crust_production=1.0,
                          alpha_exponent=1.0, ocean_depth=3000, s_vals=(0.4, 0.6, 0.8, 1.0, 1.2)):
    """T against outgassing / alpha**p for kinetic (Da < 1) runs, coloured by instellation, marker per alpha.

    If alpha and outgassing trade off exactly, every alpha falls on one curve; the best-collapsing p is printed.
    """
    sub = df[_sweep_mask(df, depth=ocean_depth, crust=False, chem=False, mg_si=mg_si,
                         delta_iw=REF_DIW, crust_production=crust_production)
             & (df['outgassing'] > 0)].copy()
    sub = sub[np.isfinite(sub['T']) & (sub['T'] < T_HOT_WALL) & (sub['T'] > T_COLD_WALL)]
    if sub.empty:
        print("No in-domain runs for the alpha figure — skipping.")
        return
    # Past Da = 1 the sink no longer responds to alpha, so those runs say nothing about the trade-off.
    sub = _add_diag_columns(sub, output_path)
    sub = sub[np.isfinite(sub['da']) & (sub['da'] < 1.0)]
    if sub['alpha'].nunique() < 2:
        print("Fewer than two alpha values with Da < 1 — skipping alpha figure.")
        return

    best_p, best_rms, rms = _alpha_collapse_exponent(sub)
    print(f"Alpha vs outgassing [Mg/Si={mg_si:g}, crust={crust_production:g}]: {len(sub)} runs, "
          f"alpha = {', '.join(f'{a:g}' for a in sorted(sub['alpha'].unique()))}")
    print(f"  best collapse at p = {best_p:+.2f} (rms {best_rms:.1f} K); plotted p = "
          f"{alpha_exponent:+.2f} (rms {rms.get(round(float(alpha_exponent), 3), np.nan):.1f} K)")

    if s_vals is not None:   # a few well-separated instellation bands stay legible
        avail = sorted(sub['instellation'].unique())
        sub = sub[sub['instellation'].isin([min(avail, key=lambda v: abs(v - t)) for t in s_vals])]
    sub['combined'] = sub['outgassing'] / sub['alpha'] ** alpha_exponent
    s_vals = sorted(sub['instellation'].unique())
    norm = mcolors.Normalize(vmin=min(s_vals), vmax=max(s_vals))

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.6))
    alphas = sorted(sub['alpha'].unique())
    for a, marker in zip(alphas, ALPHA_MARKERS):
        grp = sub[np.isclose(sub['alpha'], a)]
        ax.scatter(grp['combined'], grp['T'], marker=marker, c=grp['instellation'],
                   cmap=INSTELLATION_CMAP, norm=norm, s=20, alpha=0.9, zorder=4,
                   linewidths=0.4, edgecolors='0.25')
    _temperature_bands(ax, zorder=0)
    ax.set_xscale('log')
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel(_alpha_combination_label(alpha_exponent))
    ax.grid(True, linestyle='--', alpha=0.4)
    _add_colorbar(fig, ax, INSTELLATION_CMAP, norm, 'Instellation (S/S₀)',
                  ticks=_colorbar_ticks(s_vals)[0])
    _add_figure_legend(fig, [ax], [Line2D([0], [0], marker=m, color='0.25', linestyle='none',
                                          markersize=4, markerfacecolor='0.6',
                                          label=rf'$\alpha$ = {a:g}')
                                   for a, m in zip(alphas, ALPHA_MARKERS)])
    tag = '' if np.isclose(alpha_exponent, 1.0) else f'_p{alpha_exponent:g}'
    _save_fig(fig, figure_path(output_path, f'alpha_outgassing{tag}.png'))


def _temperature_plane(sub, xcol, ycol, xlabel, ylabel, title, path, what):
    """Filled T contours over (log10 xcol, log10 ycol) from steady states; other runs are marked."""
    # Keep values shared across the plane, dropping neighbours from other sweeps' grids.
    x_cover = sub.groupby(xcol)[ycol].nunique()
    xs = sorted(x_cover[x_cover >= 0.5 * sub[ycol].nunique()].index)
    sub = sub[sub[xcol].isin(xs)]
    y_cover = sub.groupby(ycol)[xcol].nunique()
    ys = sorted(y_cover[y_cover >= 0.5 * len(xs)].index)
    sub = sub[sub[ycol].isin(ys)]
    if len(xs) < 2 or len(ys) < 2:
        print(f"Fewer than two shared {xcol} or {ycol} values — skipping {what}.")
        return
    n_dup = int(sub.duplicated([xcol, ycol]).sum())
    if n_dup:
        print(f"  {what}: {n_dup} duplicate ({xcol}, {ycol}) run(s); using the first.")

    steady = sub[sub['termination'].isin(HABITABLE) & np.isfinite(sub['T'])]
    Z = (steady.pivot_table(index=ycol, columns=xcol, values='T', aggfunc='first')
               .reindex(index=ys, columns=xs))
    if Z.notna().sum().sum() < 4:
        print(f"Too few steady states — skipping {what}.")
        return
    x, y = np.log10(xs), np.log10(ys)
    Zm = np.ma.masked_invalid(Z.to_numpy(dtype=float))
    t_lo, t_hi = float(np.nanmin(Z.to_numpy())), float(np.nanmax(Z.to_numpy()))
    levels = np.arange(np.floor(t_lo / 10) * 10, np.ceil(t_hi / 10) * 10 + 10, 10)
    if len(levels) < 3:
        levels = np.linspace(t_lo - 1, t_hi + 1, 5)

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.8))
    cf = ax.contourf(x, y, Zm, levels=levels, cmap=TEMPERATURE_CMAP)
    cs = ax.contour(x, y, Zm, levels=levels, colors='k', linewidths=0.4, alpha=0.5)
    ax.clabel(cs, levels=levels[::2], fmt='%d K', fontsize=6, inline=True)
    unsteady = sub.drop(steady.index)
    if not unsteady.empty:
        ax.plot(np.log10(unsteady[xcol]), np.log10(unsteady[ycol]), linestyle='none',
                marker='x', color='0.45', markersize=3.5, mew=0.9, zorder=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=8)
    dx, dy = 0.03 * (x[-1] - x[0]), 0.03 * (y[-1] - y[0])   # keep edge markers inside the frame
    ax.set_xlim(x[0] - dx, x[-1] + dx)
    ax.set_ylim(y[0] - dy, y[-1] + dy)
    fig.colorbar(cf, ax=ax, pad=0.02, aspect=22).set_label('Surface temperature (K)')
    if not unsteady.empty:
        _add_figure_legend(fig, [ax], [Line2D([0], [0], linestyle='none', marker='x', color='0.45',
                                              markersize=3.5, mew=0.9, label='No steady state')])
    print(f"{what}: {len(xs)} {xcol} x {len(ys)} {ycol}, {len(steady)} steady, "
          f"{len(unsteady)} not, T {t_lo:.1f}-{t_hi:.1f} K")
    _save_fig(fig, path)


def plot_alpha_outgassing_plane(df, output_path, instellation=0.8, mg_si=REF_MG_SI,
                                crust_production=1.0, ocean_depth=3000):
    """Surface temperature over the (log alpha, log outgassing) plane at one instellation."""
    sub = df[_sweep_mask(df, depth=ocean_depth, crust=False, chem=False, instellation=instellation,
                         mg_si=mg_si, delta_iw=REF_DIW, crust_production=crust_production,
                         kd_mg=_chem_reference(df, 'kd_mg'), k_na=_chem_reference(df, 'k_na'))
             & (df['outgassing'] > 0)]
    _temperature_plane(sub, 'alpha', 'outgassing', r'$\log_{10}\,\alpha$ (reactive area scaling)',
                       r'$\log_{10}$ outgassing (×Earth)', f'S = {instellation:g}',
                       figure_path(output_path, 'alpha_outgassing_plane.png'),
                       f'Alpha-outgassing plane at S = {instellation:g}')


def plot_outgassing_crust_plane(df, output_path, instellation=0.8):
    """Surface temperature over the (log outgassing, log crust production) plane; constant ratio is a diagonal."""
    sub = _base(df)
    _temperature_plane(sub[np.isclose(sub['instellation'], instellation)], 'outgassing',
                       'crust_production', r'$\log_{10}$ outgassing (×Earth)',
                       r'$\log_{10}$ crust production (×Earth)', f'S = {instellation:g}',
                       figure_path(output_path, 'outgassing_crust_plane.png'),
                       f'Outgassing-crust plane at S = {instellation:g}')


RW_RATIO_NEUTRAL = 5e-4   # |T_RW / T_noRW - 1| below this (~0.15 K) is drawn as no change


def plot_reverse_weathering_ratio(df, output_path, crust_values=(0.01, 0.1, 1.0, 10.0)):
    """T(reverse weathering) / T(none) over instellation x outgassing, one panel per crust production.

    Only cells steady in both arms are contoured; runs on a domain wall would give a ratio of exactly 1.
    """
    arms = {rw: df[_sweep_mask(df, rw=rw) & (df['outgassing'] > 0)] for rw in (True, False)}
    if arms[False].empty:
        print("No runs without reverse weathering — skipping reverse weathering ratio plot.")
        return
    keys = ['instellation', 'outgassing', 'crust_production']
    for rw, arm in arms.items():
        n_dup = int(arm.duplicated(keys).sum())
        if n_dup:   # e.g. the alpha-outgassing plane writes out_1.0 beside the basic sweep's out_1
            print(f"  RW ratio: {n_dup} duplicate {'RW' if rw else 'no-RW'} run(s); using the first.")
            arms[rw] = arm.drop_duplicates(keys)
    pairs = arms[True].merge(arms[False], on=keys, suffixes=('_rw', '_norw'))
    pairs = pairs[pairs['crust_production'].apply(lambda c: any(np.isclose(c, v) for v in crust_values))]
    if pairs.empty:
        print("No paired runs at the requested crust production rates — skipping.")
        return
    steady = {arm: pairs[f'termination_{arm}'].isin(HABITABLE) & np.isfinite(pairs[f'T_{arm}'])
              for arm in ('rw', 'norw')}
    pairs['ratio'] = np.where(steady['rw'] & steady['norw'], pairs['T_rw'] / pairs['T_norw'], np.nan)
    one_arm = steady['rw'] ^ steady['norw']   # reverse weathering changes whether a steady state exists

    s_vals = sorted(pairs['instellation'].unique())
    out_vals = sorted(pairs['outgassing'].unique())
    dev = np.nanmax(np.abs(pairs['ratio'] - 1.0))
    if not np.isfinite(dev):
        print("No cell is steady in both arms — skipping reverse weathering ratio plot.")
        return
    # Bands log-spaced in |ratio - 1| (1-2-5), symmetric about 1, so 0.1% and 10% effects share one
    # scale; the middle band spans 1 so no-change cells (ratio 1 +- float noise) are neutral.
    edges = [m * 10.0 ** e for e in range(-4, 1) for m in (1, 2, 5)]
    edges = [x for x in edges if x >= RW_RATIO_NEUTRAL]
    edges = edges[:next((i for i, x in enumerate(edges) if x >= dev), len(edges) - 1) + 1]
    levels = 1.0 + np.array([-x for x in reversed(edges)] + edges)
    norm = mcolors.BoundaryNorm(levels, ncolors=RELATIVE_CMAP.N, extend='both')
    cbar_ticks = levels

    fig, axes = plt.subplots(1, len(crust_values), sharex=True, sharey=True, squeeze=False,
                             figsize=figure_size('double', height=2.4))
    cf = None
    for ax, crust in zip(axes[0], crust_values):
        sel = pairs[np.isclose(pairs['crust_production'], crust)]
        Z = (sel.pivot_table(index='outgassing', columns='instellation', values='ratio', aggfunc='first')
                .reindex(index=out_vals, columns=s_vals))
        Zm = np.ma.masked_invalid(Z.to_numpy(dtype=float))
        ax.set_facecolor('0.88')   # masked cells grey, so they differ from the white ratio-1 band
        if Zm.count() >= 4:
            cf = ax.contourf(s_vals, out_vals, Zm, levels=levels, cmap=RELATIVE_CMAP, norm=norm,
                             extend='both')
        for mask, marker in ((~(steady['rw'] | steady['norw']), 'x'), (one_arm, 'o')):
            m = sel[mask.loc[sel.index]]
            if marker == 'x':
                ax.plot(m['instellation'], m['outgassing'], linestyle='none', marker='x', color='0.45',
                        markersize=3, mew=0.8, zorder=4)
            else:
                ax.plot(m['instellation'], m['outgassing'], linestyle='none', marker='o',
                        markerfacecolor='none', markeredgecolor='0.25', markersize=4, mew=0.9, zorder=4)
        ax.set_yscale('log')
        ax.set_title(f'Crust prod. = {crust:g}×', fontsize=8)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_xlabel('Instellation (S/S₀)')
        n_ratio = int(np.isfinite(sel['ratio']).sum())
        print(f"  RW ratio, crust {crust:g}: {n_ratio} steady pairs, "
              f"{int(one_arm.loc[sel.index].sum())} steady in one arm only, "
              f"ratio {np.nanmin(sel['ratio']):.4f}-{np.nanmax(sel['ratio']):.4f}" if n_ratio else
              f"  RW ratio, crust {crust:g}: no steady pairs")
    axes[0, 0].set_ylabel('Outgassing (×Earth)')
    if cf is not None:
        cbar = fig.colorbar(cf, ax=list(axes[0]), pad=0.02, aspect=25, ticks=cbar_ticks,
                            spacing='uniform')
        cbar.set_label(r'$T_\mathrm{RW}\,/\,T_\mathrm{no\;RW}$')
        cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{round(v, 6):g}'))
        cbar.ax.tick_params(labelsize=6)
    _add_figure_legend(fig, list(axes[0]), [
        Line2D([0], [0], linestyle='none', marker='x', color='0.45', markersize=3, mew=0.8,
               label='No steady state in either'),
        Line2D([0], [0], linestyle='none', marker='o', markerfacecolor='none', markeredgecolor='0.25',
               markersize=4, mew=0.9, label='Steady in one only'),
    ])
    _save_fig(fig, figure_path(output_path, 'reverse_weathering_ratio.png'))


DA_TRANSITION_LABELS = {
    'crust_production': 'Crust production (×Earth)',
    'outgassing':       'Outgassing (×Earth)',
}
DA_TRANSITION_MG_SI = (0.8, 1.25, 1.8)   # panels, intersected with the basic sweeps on disk


def _da_transition_instellation(group):
    """(S where Da first rises through 1, interpolated in log Da, whether a bracketing run is wall-clamped).

    Returns (nan, False) when the line never goes from kinetic to thermodynamic.
    """
    g = group.sort_values('instellation')
    da = g['da'].to_numpy(dtype=float)
    ok = np.isfinite(da) & (da > 0)
    if ok.sum() < 2:
        return np.nan, False
    s = g['instellation'].to_numpy(dtype=float)[ok]
    T = g['T'].to_numpy(dtype=float)[ok]
    ld = np.log10(da[ok])
    # Upward crossings only: a frozen, CO2-starved state can sit at Da > 1 at low S.
    idx = np.flatnonzero((ld[:-1] < 0) & (ld[1:] >= 0))
    if idx.size == 0:
        return np.nan, False
    i = idx[0]
    s_crit = s[i] - ld[i] * (s[i + 1] - s[i]) / (ld[i + 1] - ld[i])
    clamped = bool(np.nanmax(T[i:i + 2]) >= T_HOT_WALL - 0.1 or
                   np.nanmin(T[i:i + 2]) <= T_COLD_WALL + 0.1)
    return float(s_crit), clamped


def plot_da_transition(df, output_path, mg_si_values=DA_TRANSITION_MG_SI,
                       line_by='crust_production', cmap=None, show_hz=None, max_outgassing=3.0):
    """Instellation of the kinetic-to-thermodynamic transition, one panel per Mg/Si.

    line_by picks which tectonic axis carries the lines and colour bar; the other is the y-axis.
    """
    if line_by not in DA_TRANSITION_LABELS:
        raise ValueError(f"line_by must be one of {sorted(DA_TRANSITION_LABELS)}, not {line_by!r}")
    y_by = 'outgassing' if line_by == 'crust_production' else 'crust_production'
    cmap = PARAM_CMAPS[line_by] if cmap is None else cmap
    available = basic_plane_mg_si(df)
    mg_vals = [m for m in mg_si_values if any(np.isclose(m, a) for a in available)]
    if not mg_vals:
        print("No basic Mg/Si plane with a full sweep — skipping Da transition figure.")
        return

    subset = pd.concat([_add_diag_columns(sel, output_path) for mg in mg_vals
                        if not (sel := _base(df, mg_si=None if np.isclose(mg, REF_MG_SI) else mg)).empty])
    # Main-sweep outgassing values only; other sweeps' grids exist at crust 1x alone.
    main_out = [o for o in cb.ps.outgassing if o <= max_outgassing]
    subset = subset[subset['outgassing'].apply(lambda v: any(np.isclose(v, o) for o in main_out))]

    records = []
    for (mg, out, crust), group in subset.groupby(['mg_si', 'outgassing', 'crust_production']):
        s_crit, clamped = _da_transition_instellation(group)
        records.append({'mg_si': mg, 'outgassing': out, 'crust_production': crust,
                        's_crit': s_crit, 'clamped': clamped})
    trans = pd.DataFrame(records).dropna(subset=['s_crit'])
    if trans.empty:
        print("No Da = 1 crossings found — skipping Da transition figure.")
        return
    n_clamped = int(trans['clamped'].sum())
    if n_clamped:
        print(f"  Da transition: {n_clamped} of {len(trans)} crossings bracketed by a T-wall state.")

    line_vals = sorted(trans[line_by].unique())
    norm = _value_norm(line_vals, pad=0.0)
    ticks, ticklabels = _colorbar_ticks(line_vals)
    fig, axes = plt.subplots(1, len(mg_vals), sharex=True, sharey=True, squeeze=False,
                             figsize=figure_size('double', height=2.6))
    axes = axes[0]
    for ax, mg in zip(axes, mg_vals):
        panel = trans[np.isclose(trans['mg_si'], mg)]
        for c in line_vals:
            line = panel[np.isclose(panel[line_by], c)].sort_values(y_by)
            if not line.empty:
                ax.plot(line['s_crit'], line[y_by], color=cmap(norm(c)), linewidth=1.4, zorder=3)
        ax.set_yscale('log')
        ax.set_title(f'Mg/Si = {mg:g}')
        ax.grid(True, linestyle='--', alpha=0.3)
        _draw_hz_edges(ax, show_hz)
    axes[0].set_ylabel(DA_TRANSITION_LABELS[y_by])
    axes[len(axes) // 2].set_xlabel('Instellation (S/S$_0$)')

    marks = list(trans['s_crit'])
    if SHOW_HZ_EDGES if show_hz is None else show_hz:
        marks += [CONTINENTAL_HZ_OUTER, CONTINENTAL_HZ_INNER]
    lo, hi = min(marks), max(marks)
    pad = 0.08 * (hi - lo) if hi > lo else 0.1
    axes[0].set_xlim(lo - 4 * pad, hi + 4 * pad)   # each margin holds a regime label
    for ax in axes:
        trans_ax = ax.get_xaxis_transform()
        # Below mid-height so the boxes clear the HZ edge labels.
        for x, text, face, edge in ((lo - 2 * pad, 'Kinetic\n(stable) regime', '#d4eddb', '#3f8f55'),
                                    (hi + 2 * pad, 'Thermodynamic\n(unstable) regime', '#f8d6d1', '#b8483a')):
            ax.text(x, 0.35, text, transform=trans_ax, rotation=90, ha='center', va='center',
                    fontsize=7, color='0.15', zorder=5,
                    bbox=dict(boxstyle='round,pad=0.35', facecolor=face, edgecolor=edge, linewidth=0.6))
    _add_colorbar(fig, list(axes), cmap, norm, DA_TRANSITION_LABELS[line_by], ticks=ticks,
                  ticklabels=ticklabels)
    stem = 'da_transition' if line_by == 'crust_production' else 'da_transition_by_outgassing'
    _save_fig(fig, figure_path(output_path, f'{stem}.png'))


def _largest_rectangle(mask):
    """Inclusive (row0, row1, col0, col1) of the largest all-True rectangle in a 2-D mask, or None."""
    best, best_area = None, 0
    heights = np.zeros(mask.shape[1], dtype=int)
    for j in range(mask.shape[0]):
        heights = np.where(mask[j], heights + 1, 0)
        for i0 in range(mask.shape[1]):
            h = heights[i0]
            for i1 in range(i0, mask.shape[1]):
                h = min(h, heights[i1])
                if h == 0:
                    break
                if h * (i1 - i0 + 1) > best_area:
                    best_area, best = h * (i1 - i0 + 1), (j - h + 1, j, i0, i1)
    return best


def plot_habitability_phase_space(df, output_path, s_max=1.2):
    """Map of climate state over outgassing / crust production and instellation, labelled by region.

    Each cell takes the majority state of its runs, drawn as a hatched block outlined where the state changes.
    """
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Rectangle

    base = _base(df).copy()
    if base.empty:
        print("No runs found — skipping phase space plot.")
        return
    base['ratio'] = base['outgassing'] / base['crust_production']
    # Only a run stopped at a CO2 wall has an unknown fate; the final T decides the rest.
    wall = base['domain_wall']
    unknown = base['termination'].isin(OUT_OF_DOMAIN) & wall.isin(['co2_high', 'co2_low'])
    snow = ((base['T'] <= T_SNOWBALL) | (wall == 'cold')) & ~unknown
    hot = ((base['T'] >= T_RUNAWAY) | (wall == 'hot')) & ~unknown & ~snow
    base['state'] = np.select([unknown, snow, hot], [3, 0, 2], 1)   # 0 snowball, 1 habitable, 2 hothouse, 3 unknown

    # Merge near-identical ratios from different sweep grids (e.g. 0.03, 0.032, 0.0333) into one column.
    logr = np.log10(np.sort(base['ratio'].unique()))
    groups = np.concatenate([[0], np.cumsum(np.diff(logr) > 0.05)])
    centres = np.array([logr[groups == g].mean() for g in range(groups[-1] + 1)])
    base['col'] = groups[np.searchsorted(logr, np.log10(base['ratio']).clip(logr[0], logr[-1]))]
    s_vals = np.sort(base['instellation'].unique())
    cover = base.groupby('col')['instellation'].nunique()   # thin columns would leave holes
    cols = [c for c in range(len(centres)) if cover.get(c, 0) >= 0.5 * len(s_vals)]
    base = base[base['col'].isin(cols)]

    grid = np.full((len(s_vals), len(cols)), -1)
    n_mixed = 0
    for (c, s), g in base.groupby(['col', 'instellation']):
        counts = g['state'].value_counts()
        grid[np.searchsorted(s_vals, s), cols.index(c)] = counts.idxmax()
        n_mixed += len(counts) > 1

    def _edges(c):
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])
    xe, ye = 10 ** _edges(centres[cols]), _edges(s_vals)

    style = {0: ('#9cc5ea', '#3a78b5', '\\\\\\\\\\\\'),
             1: ('#a8dbb4', '#3f8f55', None),
             2: ('#f2a99f', '#b8483a', '//////'),
             3: ('#d9d9d9', '#8a8a8a', '...')}
    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=3.0))
    for k, (face, edge, hatch) in style.items():
        cells = [Rectangle((xe[i], ye[j]), xe[i + 1] - xe[i], ye[j + 1] - ye[j])
                 for j, i in zip(*np.nonzero(grid == k))]
        if cells:
            ax.add_collection(PatchCollection(cells, facecolor=face, edgecolor=edge, hatch=hatch,
                                              linewidth=0, zorder=2))

    padded = np.pad(grid, 1, constant_values=-1)   # outline state changes, including against empty cells
    segs = []
    for j in range(1, padded.shape[0] - 1):
        for i in range(padded.shape[1] - 1):
            if padded[j, i] != padded[j, i + 1]:
                segs.append([(xe[i], ye[j - 1]), (xe[i], ye[j])])
    for j in range(padded.shape[0] - 1):
        for i in range(1, padded.shape[1] - 1):
            if padded[j, i] != padded[j + 1, i]:
                segs.append([(xe[i - 1], ye[j]), (xe[i], ye[j])])
    ax.add_collection(LineCollection(segs, colors='0.15', linewidths=0.9, zorder=3))

    ax.set_xscale('log')
    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[0], min(ye[-1], s_max))
    ax.set_xlabel('Outgassing / Crust production rate')
    ax.set_ylabel('Instellation (S/S₀)')
    # Label each state at the centre of its largest visible block of cells.
    visible = ye[:-1] < s_max
    for k, name in enumerate(['Snowball', 'Habitable', 'Hothouse', 'Unknown\n(CO₂ wall)']):
        rect = _largest_rectangle((grid == k) & visible[:, None])
        if rect is None:
            continue
        j0, j1, i0, i1 = rect
        ax.text(10 ** (0.5 * (np.log10(xe[i0]) + np.log10(xe[i1 + 1]))),
                0.5 * (ye[j0] + min(ye[j1 + 1], s_max)), name, ha='center', va='center',
                fontsize=7, zorder=5,
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', edgecolor='none', alpha=0.85))
    print(f"Ratio phase space: {len(cols)} ratio columns x {len(s_vals)} instellations, "
          f"{n_mixed} cell(s) where runs disagree (majority shown).")
    _save_fig(fig, figure_path(output_path, 'ratio_phase_space.png'))


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
    parser.add_argument('--oscillating', choices=OSC_MODES, default=OSCILLATION_MODE,
                        help="How runs still cycling at their end are shown: their final state, the mean over "
                             "the last 40%% of the run, or not at all (default: %(default)s).")
    parser.add_argument('--pe', type=float, default=None,
                        help='Pin ocean pe to this value in the main plots '
                             '(default: the model reference, planet.PE_DEFAULT).')
    args = parser.parse_args()
    for _knob in CHEM_KNOBS:
        if getattr(args, _knob) is not None:
            CHEM_OVERRIDE[_knob] = getattr(args, _knob)
    if args.pe is not None:
        REF_PE = args.pe

    OSCILLATION_MODE = args.oscillating
    df = load_data(args.path)
    if args.legacy:
        import plot_legacy
        df = plot_legacy.upgrade(df)
        plot_legacy.plot_named_compositions(df, args.path, split_panels=True)
    if df.empty:
        print("No data found. Check --path.")
        raise SystemExit(1)

    # Depths carrying a composition sweep, and depths with a resolved redox sweep (more than two pe values).
    comp_depths = [args.depth] if args.depth is not None else (sorted(
        d for d in df['ocean_depth'].unique()
        if df[df['ocean_depth'] == d][['mg_si', 'delta_iw']].nunique().max() > 1) or [3000.0])
    redox_depths = sorted(d for d in df['ocean_depth'].unique()
                          if df[df['ocean_depth'] == d]['pe'].nunique() > 2) or [3000.0]
    print(f"Composition figures for depth(s): {[f'{d:g}' for d in comp_depths]}")
    print(f"Redox figures for depth(s): {[f'{d:g}' for d in redox_depths]}")
    basic_mg = basic_plane_mg_si(df)
    end_members = [m for m in basic_mg if not np.isclose(m, REF_MG_SI)]
    if end_members:
        print(f"Basic sweep Mg/Si planes: {[f'{m:g}' for m in basic_mg]}")

    path = args.path
    # Basic sweep, at Earth's Mg/Si and each end-member.
    for _mg in [None] + end_members:
        plot_basic(df, path, split_panels=True, mg_si=_mg)
        plot_basic(df, path, all_results=False, split_panels=True, mg_si=_mg)
    plot_basic_ph(df, path)
    plot_basic_mgsi_grid(df, path, split_panels=True)
    # One-axis sweeps.
    plot_depth(df, path, split_panels=False)
    plot_chemistry(df, path, split_panels=True)
    for _d in redox_depths:
        plot_pe(df, path, split_panels=True, ocean_depth=_d)
    # Crust composition.
    for _d in comp_depths:
        plot_cross(df, path, split_panels=True, ocean_depth=_d)
        plot_composition(df, path, split_panels=True, ocean_depth=_d)
    for _q in ('T', 'P_CO2'):
        plot_composition_map(df, path, ocean_depth=comp_depths[0], quantity=_q)
    # Tectonics, alpha and the Da transition.
    plot_damkohler_contour(df, path)
    plot_habitability_phase_space(df, path)
    plot_outgassing_crust_plane(df, path)
    plot_reverse_weathering_ratio(df, path)
    plot_da_transition(df, path)
    plot_da_transition(df, path, line_by='outgassing')
    plot_alpha_outgassing(df, path)
    plot_alpha_outgassing(df, path, alpha_exponent=-1.0)
    plot_alpha_outgassing_plane(df, path)
    # Continental figures.
    plot_continental(df, path)
    print("Done.")
