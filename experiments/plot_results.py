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
from kamino.planet import KD_MG_HT, K_NA_CONT_REMOVAL, PE_DEFAULT
from kamino.weathering import ALPHA_REF
from kamino.planet import _S_TERR_EARTH

import continental_baseline as cb

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


# ---------------------------------------------------------------------------
# Colour maps
# ---------------------------------------------------------------------------
# ONE colour map per swept input parameter, defined here and nowhere else, so that a parameter
# is recognisable by colour across every figure that facets on it -- outgassing is always
# tropical, ocean depth always bubblegum, and so on. Change a name here to recolour every figure
# that uses that parameter in one go. Anything cmasher or matplotlib exposes works.
#
# Several are sub-mapped: the full range of most sequential maps starts at black and ends at
# near-white, and with only 5-7 grid values per axis the end members were landing in both, where
# they are indistinguishable from the axis text and from the panel background respectively.

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

# Maps for OUTPUT quantities, which are not parameters and so are deliberately kept out of
# PARAM_CMAPS: a diverging map for signed quantities about a meaningful centre (log Da about
# 1, a difference against the Earth crust) and a sequential pair for absolute values.
DIVERGING_CMAP        = cmr.prinsenvlag   # log10(Da), centred on Da = 1
RELATIVE_CMAP         = cmr.fusion_r      # differences against the reference crust
WEATHERING_RATIO_CMAP = cmr.fusion        # log seafloor/continental: blue where seafloor dominates
QUANTITY_CMAP         = cmr.ember         # absolute value of a log-scaled quantity
TEMPERATURE_CMAP      = cmr.get_sub_cmap('cmr.ember', 0.12, 1.0)   # surface temperature contours

DEFAULT_OUTPUT_PATH = os.environ.get('KAMINO_SWEEP_OUTPUT', '/home/pt426/Code/kamino/sweep_output')

TERM_LABELS = {
    'converged':      'Converged',
    'timeout':        'Timeout (2 Gyr)',
    'wall_timeout':   'Wall-clock cap',
    'out_of_domain':  'Outside model domain',
    'fallback_limit': 'Chemistry fallback cap',
    'chemistry_void': 'Chemistry solver void',
    'solver_failure': 'ODE solver failure',
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
# NOTE 'wall_timeout', 'fallback_limit', 'chemistry_void' and 'solver_failure' are deliberately
# in NEITHER set: each is an incomplete or failed integration, so it is neither known-habitable
# nor known-out-of-domain. Every one of them gets its own FAILED_MARKERS entry (see below) so a
# run's specific unreliability is visible rather than lumped into a single generic marker.
HABITABLE = {'converged', 'timeout'}

# Terminations whose Da is trustworthy enough to anchor a kinetic<->thermodynamic transition
# marker (see _plot_group_on_axes / _plot_line_da_style). Wider than HABITABLE: a wall_timeout
# run's ABORT state is a real state dY_dt already evaluated once, and Planet._final_diagnostics
# now re-evaluates it a second time (deadline suspended) to record real diagnostics there -- the
# run simply didn't reach 2 Gyr or an event, which says nothing about whether ITS Da is trustworthy.
# The others do not get this benefit of the doubt: an out_of_domain state can be a fabricated
# ceiling/floor sentinel (T pinned to an exact 389.00 K clamp is not a converged temperature),
# and fallback_limit/chemistry_void/solver_failure all mean the chemistry solver or integrator
# was actively struggling at that exact point, not merely stopped early on an otherwise-fine one.
#
# NOTE this only helps sweeps re-run after Planet._final_diagnostics started covering the abort
# state (see planet.py) -- older output still has da=NaN for wall_timeout points regardless of
# this set, since np.isfinite(da[i]) gates them out first either way.
DA_TRUSTWORTHY = HABITABLE | {'wall_timeout'}

T_SNOWBALL = 260.0
T_RUNAWAY  = 360.0

# The validity box the climate solver is scanned over (planet.py, T_LO/T_HI). A run pinned to
# either value did not converge onto that temperature -- it is a sentinel -- so any diagnostic
# evaluated there is a property of the clamp, not of the planet.
T_COLD_WALL = 181.0
T_HOT_WALL  = 389.0

# ---------------------------------------------------------------------------
# Habitable-zone edges of the Earth-like continental baseline
# ---------------------------------------------------------------------------
# Instellation limits of the temperate band on the reference Earth-like planet, measured by
# experiments/continental_baseline.py: land fraction 0.3 with every other axis at Earth --
# 1x outgassing, 1x crust production, 3 km ocean, Earth crust (Mg/Si 1.25, dIW -2), reverse
# weathering on, reducing ocean. Drawn as vertical lines on TEMPERATURE panels, so any figure of
# T against instellation can be read against the same reference planet. Temperature only: the
# edges say where the temperate band is entered and left, which is not a statement about pCO2,
# pH, salinity or mineral saturation, so marking those panels would imply a threshold in a
# quantity the edges say nothing about.
#
# OUTER is a crossing: T interpolated onto T_SNOWBALL between the S = 0.45 and S = 0.50 runs.
# It is CO2-SUPPLY limited, not radiation limited. The WHAK continental sink stays near 2x the
# modern Earth rate all the way out (beta = 0.3 makes the pCO2 term cancel the temperature
# term), so pCO2 never reaches the several bars a maximum greenhouse needs -- the analytic
# climate model's own maximum-greenhouse edge is nearer S = 0.41, and Kopparapu et al. (2013)
# put it at 0.35.
#
# INNER is bracketed between the S = 1.10 and S = 1.15 runs. At 1.15 absorbed instellation
# exceeds the OLR (Simpson-Nakajima) limit, so no cool-branch solution exists and the planet is
# in a runaway greenhouse; the climate model puts that threshold at S = 1.141 once CO2 is
# exhausted. Do not read the inner edge to better than the 0.05 grid spacing.
#
# These are properties of ONE reference planet. A figure that varies ocean depth, crust
# composition or outgassing is being compared against Earth-like continents, not against its
# own habitable zone.
CONTINENTAL_HZ_OUTER = 0.480
CONTINENTAL_HZ_INNER = 1.125

# Whether the HZ lines are drawn. THIS is the switch: set it to True to put the edges on every
# instellation figure in one go. Off by default so existing figures are unchanged. Every
# plotting function also takes show_hz=True/False to override it for one figure.
SHOW_HZ_EDGES = True


def _draw_hz_edges(ax, show_hz=None):
    """Mark the continental baseline's habitable-zone edges on a TEMPERATURE axis.

    `show_hz` of None defers to the module default, so a caller that does not care need not
    thread the flag; True or False decides for this axis alone.
    """
    if not (SHOW_HZ_EDGES if show_hz is None else show_hz):
        return
    trans = ax.get_xaxis_transform()
    # Each label sits outside the habitable band, beside its line.
    for s, label, offset, ha in ((CONTINENTAL_HZ_OUTER, 'HZ outer edge', -3, 'right'),
                                 (CONTINENTAL_HZ_INNER, 'HZ inner edge', 3, 'left')):
        ax.axvline(s, color='0.35', linestyle=(0, (6, 3)), linewidth=1.0, alpha=0.85, zorder=1)
        ax.annotate(label, xy=(s, 0.97), xycoords=trans, xytext=(offset, 0),
                    textcoords='offset points', rotation=90, ha=ha, va='top',
                    fontsize=6, color='0.35', zorder=5,
                    bbox=dict(boxstyle='square,pad=0.1', facecolor='white', edgecolor='none',
                              alpha=0.75))


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

# Composition-sweep Mg/Si values left off the Mg/Si colour-bar figures; 1.75 crowds the 1.8 end-member.
MG_SI_HIDDEN = (1.75,)

# Reference ocean redox: the model's own default (abiotic, reducing -- see planet.PE_DEFAULT and
# development_history.md section 28). Every sweep since 2026-08-27 runs both redox arms
# (parameter_sweep.PE_STATES), so this is what every figure EXCEPT plot_pe pins pe to,
# exactly as REF_MG_SI/REF_DIW pin the composition axes for every figure except the crust ones.
# Overridable from the CLI (--pe), the same way the CHEM_KNOBS are.
REF_PE = float(PE_DEFAULT)

# What a run predating the `pe` parameter (before 2026-08-27) actually ran at: PHREEQC's own
# implicit default, which section 28.3 found is firmly oxidising. NOT the same as REF_PE --
# defaulting missing pe to REF_PE would silently relabel every pre-pe run as reducing.
PE_LEGACY_DEFAULT = 4.0

COMP_COLORS = {
    'komatiite_42': '#7b2d8b',
    'komatiite_44': '#c44bc4',
    'basalt_47':    '#e08040',
    'basalt_49':    '#d4b000',
    'basalt_51':    '#6abf69',
}
HAB_MARKERS    = {'converged': 'o', 'timeout': 's'}
# One distinct (hollow) marker per unreliable termination, so a run's specific failure mode is
# visible on the plot instead of every non-habitable point collapsing onto a single 'x'.
FAILED_MARKERS = {
    'out_of_domain':  's',
    'wall_timeout':   '^',
    'fallback_limit': 'D',
    'chemistry_void': 'v',
    'solver_failure': 'x',
}

DA_LEGEND = [
    Line2D([0], [0], color='k', linestyle='-',  linewidth=1.4, label='Da < 1 (kinetic)'),
    Line2D([0], [0], color='k', linestyle='--', linewidth=1.4, label='Da ≥ 1 (thermodynamic)'),
    Line2D([0], [0], color='k', linestyle=':',  linewidth=1.4, label='$T_\\mathrm{sf}$ at floor (274 K)'),
]

PANEL_COLS = ['T', 'P_CO2', 'pH', 'salinity']



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
    pe = float(d['pe']) if d.get('pe') is not None else None  # None -> library default (oxidising, see PE_LEGACY_DEFAULT)
    fast_minerals = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals
    try:
        F_prec, _, _ = get_precipitation(P_pore, T_seafloor, b_ocean, fast_minerals,
                                         precipitation_timescale=float(d.get('tau_prec', 1e5 * YR)),
                                         pe=pe)
        if rw:
            F_rw, _, _ = get_precipitation(P_pore, T_seafloor, b_ocean, list(reverse_weathering_minerals),
                                           precipitation_timescale=float(d.get('tau_rw', 5e6 * YR)),
                                           pe=pe)
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
        pe = float(d['pe']) if d.get('pe') is not None else None  # None -> library default (oxidising, see PE_LEGACY_DEFAULT)

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
            pe=pe,
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
                precipitation_timescale=float(d.get('tau_prec', 1e5 * YR)), pe=pe)
            pH_recomputed = float(pH_recomputed)
        except (ChemistryError, Exception):
            pH_recomputed = np.nan

        # Pore SI: the pore space now precipitates clays only, so Calcite no longer
        # appears in secondary_SI. Evaluate it directly on the post-weathering pore
        # fluid (b_pore) — the driving force for pore-space calcite.
        try:
            b_pore = np.asarray(diag['b_pore'], dtype=float)
            _, _, si_p = get_precipitation(P_pore, T_pore, np.maximum(b_pore, 0.0), ['Calcite'],
                                           precipitation_timescale=1e6 * YR, pe=pe)
            calcite_si = float(si_p.get('Calcite', np.nan))
        except (ChemistryError, Exception):
            calcite_si = np.nan

        # Ocean SI: only reliable when Ca_ocean > 0; set to NaN otherwise to avoid
        # the spurious -∞ that results when ocean Ca has been depleted to the ODE floor.
        if b_ocean[_ca_idx] > 1e-6:
            try:
                _, _, si_o = get_precipitation(P_pore, T_seafloor, b_ocean, ['Calcite'],
                                               precipitation_timescale=1e6 * YR, pe=pe)
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
            'pe':                 float(d['pe']) if d.get('pe') is not None else PE_LEGACY_DEFAULT,
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


def _ref_redox(df):
    """Mask pinning ocean pe to the reference (reducing, abiotic) redox state.

    Every sweep since 2026-08-27 runs both redox arms (parameter_sweep.PE_STATES), doubling
    every combination of the other axes. Without this, every line-based figure would silently
    draw one point per instellation for the oxidising arm and one for the reducing arm and join
    them as if they were a single trajectory. pe gets its own figure (plot_pe), which
    is the one caller that must NOT apply this mask.
    """
    return np.isclose(df['pe'], REF_PE)


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

    Everything except plot_chemistry shows ONE chemistry. Without this, runs that
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


def _base(df, mg_si=None):
    """Sweep 1: reference crust and chemistry, rw=True, depth=3000, outgassing>0, ocean world.

    `mg_si` selects a non-Earth mantle Mg/Si instead of the reference one, for the basic sweep
    repeated at the composition end-members (parameter_sweep.sweep_basic_low_mgsi /
    sweep_basic_high_mgsi). dIW stays at its reference either way, so the selection is still a
    single crust rather than a mixture.
    """
    crust = (_ref_crust(df) if mg_si is None else
             np.isclose(df['mg_si'], mg_si) & np.isclose(df['delta_iw'], REF_DIW))
    return df[
        df['reverse_weathering'] &
        crust &
        _ref_chem(df) &
        _ref_redox(df) &
        (df['ocean_depth'] == 3000) &
        (df['land_fraction'] == 0.0) &
        (df['outgassing'] > 0)
    ]


def basic_plane_mg_si(df, min_combos=4):
    """Mantle Mg/Si values whose basic sweep -- the outgassing x crust-production plane -- was run.

    The composition and cross sweeps also write rows at non-Earth Mg/Si, but only at the ONE
    (outgassing, crust) pair those designs hold fixed; treating those as a basic sweep would draw
    a one-point "plane". Requiring several combinations restricts this to the Mg/Si values
    sweep_basic / sweep_basic_low_mgsi / sweep_basic_high_mgsi actually cover, so the figures
    below appear when those sweeps have been run and not otherwise.
    """
    pool = df[
        df['reverse_weathering'] &
        _ref_chem(df) &
        _ref_redox(df) &
        np.isclose(df['delta_iw'], REF_DIW) &
        (df['ocean_depth'] == 3000) &
        (df['land_fraction'] == 0.0) &
        (df['outgassing'] > 0)
    ]
    if pool.empty:
        return []
    n_combos = pool.groupby('mg_si').apply(
        lambda g: g.groupby(['outgassing', 'crust_production']).ngroups, include_groups=False)
    return sorted(float(v) for v, n in n_combos.items() if n >= min_combos)


# ---------------------------------------------------------------------------
# Shared plot helpers
# ---------------------------------------------------------------------------

def _panel_groups(split):
    """Return [(cols, filename_suffix), ...] — one entry normally, two when split."""
    if split:
        return [(['T', 'P_CO2'], '_tp'), (['pH', 'salinity'], '_chem')]
    return [(PANEL_COLS, '')]


def _style_axes(axes, cols, x_lims=(0.25, 1.45), show_hz=None, show_eq_temp=False):
    """Style a set of axes given the column names they represent.

    `show_hz` draws the continental baseline's habitable-zone edges on the TEMPERATURE panel
    only. The edges are a statement about temperature -- where the band is entered and left --
    so putting them on a pCO2 or salinity panel would assert a threshold in a quantity they say
    nothing about. `show_eq_temp` adds the equilibrium-temperature curve to that panel.
    """
    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(*x_lims)
    for ax, col in zip(axes, cols):
        if col == 'T':
            _draw_hz_edges(ax, show_hz)
            ax.set_ylabel('Temperature (K)')
            ax.axhspan(T_SNOWBALL - 25, T_SNOWBALL, color='blue', alpha=0.12)
            ax.axhspan(T_RUNAWAY - 20,  T_RUNAWAY,  color='red',  alpha=0.12)
            ax.set_ylim(235, 360)
            if show_eq_temp:
                s_eq = np.linspace(x_lims[0], x_lims[1], 300)
                ax.plot(s_eq, equilbrium_temperature(s_eq), color='k', linestyle='-.',
                        linewidth=0.8, zorder=1, alpha=0.7)
        elif col == 'P_CO2':
            ax.set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
            ax.set_yscale('log')
            ax.set_ylim(1e-5, 20)
        elif col == 'pH':
            ax.set_ylabel('Ocean pH')
            ax.set_ylim(5, 9)
        elif col == 'salinity':
            ax.set_ylabel('Salinity (g/kg)')
            ax.set_yscale('log')
            ax.set_ylim(1e-1, 1e2)
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


def _make_legend_handles(show_markers=True, prefix_handles=None, show_hz=None, cols=None):
    """Build legend handle list: prefix (default DA_LEGEND) + optional marker entries.

    The HZ edges are labelled on the axes themselves (`_draw_hz_edges`), so they get no entry.
    """
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


FIGURE_SUBDIR = 'figures'


def figure_path(output_path, name):
    """Path for a figure: the `figures/` subdirectory of the sweep directory, made on demand.

    Keeps the sweep directory to run JSON and the diagnostics cache -- a full plot run drops ~70
    files, which previously buried the runs they were made from.
    """
    directory = os.path.join(output_path, FIGURE_SUBDIR)
    os.makedirs(directory, exist_ok=True)
    return os.path.join(directory, name)


def _save_fig(fig, path, tight=False):
    r"""Write a figure at EXACTLY the size it was created with.

    `bbox_inches='tight'` is deliberately off. The style sets constrained_layout, which already
    fits every artist inside the canvas, and 'tight' then re-crops to the content box -- which
    for these figures came out LARGER than the requested width (3.43 in against the 3.32 in
    column). A figure wider than the column gets scaled down by \includegraphics, shrinking the
    type below the size the style chose. Pass tight=True for diagnostics, where exact width does
    not matter.
    """
    kw = {'bbox_inches': 'tight'} if tight else {}
    stem = os.path.splitext(path)[0]
    for ext in ('png', 'pdf'):
        fig.savefig(f'{stem}.{ext}', **kw)
    plt.close(fig)
    print(f"Saved {stem}.png / .pdf")


def _style_combined_col(axes_c, ci, n_cols, title='', cols=None, show_hz=None):
    """Style one column of a multi-column grid: axis labels, tick visibility, x-label."""
    if cols is None:
        cols = PANEL_COLS
    col_axes = axes_c[:, ci]
    _style_axes(col_axes, cols, show_hz=show_hz)
    if title:
        axes_c[0, ci].set_title(title)
    if ci > 0:
        for ax in col_axes:
            ax.set_ylabel('')
        for row in range(axes_c.shape[0]):
            plt.setp(axes_c[row, ci].get_yticklabels(), visible=False)
    if ci != n_cols // 2:
        col_axes[-1].set_xlabel('')


def _plot_line_da_style(ax, x, y, da, color, at_floor=None, trustworthy=None,
                        linewidth=1.4, alpha=0.8, zorder=3):
    """Draw a line: dotted where seafloor T is floored, dashed where Da≥1, solid where Da<1.

    Places an open circle at each kinetic↔thermodynamic transition (solid↔dashed only;
    transitions into/out of the dotted floor regime are not marked).

    A point's Da is only trusted (usable to CONFIRM a transition) when `trustworthy[i]` --
    normally the run's termination being genuinely HABITABLE (converged/timeout). Two other
    cases produce a Da that must not anchor a marker even when it is a finite number:

    - `NaN`: a run that hit the wall-clock cap (`wall_timeout`, sometimes `fallback_limit`) is
      cut off before `time_evolve`'s final-state re-evaluation runs, so no diagnostics were ever
      recorded.
    - A finite but FABRICATED Da: an `out_of_domain` run's stored state can be a sentinel the
      model clamped to when it gave up (T pinned to an exact ceiling -- e.g. 389.00 K, identical
      to 15 significant figures across otherwise-unrelated runs -- is not a converged physical
      temperature), and Da computed from that sentinel is enormous but meaningless. At crust
      production >= a few x Earth combined with even modest outgassing, many lines leave the
      valid domain (or stall in wall_timeout) before ever reaching a confirmed Da >= 1, so the
      only "dashed" point available is this sentinel -- using it fabricated a transition that was
      never actually observed within the model's valid domain.

    Untrustworthy points are still drawn (resolved to whichever TRUSTWORTHY neighbour is
    nearest, so the line stays continuous -- a run with no trustworthy Da anywhere stays solid),
    but several consecutive untrustworthy points (a real occurrence near the runaway transition,
    where the integration is stiffest and most likely to hit the wall or leave the domain) all
    copy the same nearest trustworthy style, which can carry a "kinetic" classification several
    grid points past where Da actually left that regime -- so a transition circle is placed at
    the LAST point on the old-regime side that is itself trustworthy, not at the raw style-flip
    index.
    """
    if len(x) < 2:
        return

    if at_floor is None:
        at_floor = np.zeros(len(x), dtype=bool)
    if trustworthy is None:
        trustworthy = np.ones(len(x), dtype=bool)

    def _raw_ls(i):
        """None means "Da not trustworthy at this point" -- resolved below, not drawn directly."""
        if at_floor[i]:
            return ':'
        if not (trustworthy[i] and np.isfinite(da[i])):
            return None
        return '-' if da[i] < 1 else '--'

    raw = [_raw_ls(i) for i in range(len(x))]
    resolved = list(raw)
    for i, s in enumerate(raw):
        if s is not None:
            continue
        prv = next((resolved[j] for j in range(i - 1, -1, -1) if raw[j] is not None), None)
        nxt = next((raw[j] for j in range(i + 1, len(raw)) if raw[j] is not None), None)
        resolved[i] = prv or nxt or '-'

    seg_x = [x[0]]
    seg_y = [y[0]]
    current_ls = resolved[0]
    trans_x, trans_y = [], []

    for i in range(1, len(x)):
        ls = resolved[i]
        if ls == current_ls:
            seg_x.append(x[i])
            seg_y.append(y[i])
        else:
            ax.plot(seg_x, seg_y, color=color, linewidth=linewidth,
                    alpha=alpha, linestyle=current_ls, zorder=zorder)
            # Mark only solid↔dashed transitions (skip dotted floor segments), at the last
            # point on the old-regime side with a directly measured Da -- walking back past
            # any unknown (wall_timeout) points rather than mis-marking one of them.
            if {current_ls, ls} == {'-', '--'}:
                j = i - 1
                while j > 0 and raw[j] is None:
                    j -= 1
                trans_x.append(x[j])
                trans_y.append(y[j])
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
                   width='single', height=None, show_hz=None, groups=None):
    """The standard one-variable faceted line figure, one file per panel group.

    `groups` overrides the panel groups, e.g. [(['pH'], '_ph')] for an abridged single panel.

    Every figure that plots T / P_CO2 / pH / salinity against instellation with one line per value
    of a single variable -- ocean depth, a chemistry constant, mantle Mg/Si, dIW -- shares this
    body. They differed only in which rows they select and how the lines are labelled, so the
    selection stays with the caller and the assembly lives here.

    `colours` is parallel to `values` rather than derived from `norm`, so a caller can colour by
    something other than the facet value itself (plot_legacy does, by melt SiO2).

    `aspect_per_row` is multiplied by the panel count rather than passed through, because the
    number of panels differs between the split and combined panel groups.
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
        _style_axes(axes, cols, show_hz=show_hz,
                    **({} if x_lims is None else {'x_lims': x_lims}))
        _add_colorbar(fig, list(axes), cmap, norm, cbar_label, ticks=ticks,
                      ticklabels=ticklabels,
                      **({} if aspect_per_row is None else {'aspect': n_rows * aspect_per_row}))
        _add_figure_legend(fig, axes,
                           _make_legend_handles(show_markers=show_markers, show_hz=show_hz,
                                                cols=cols))
        _save_fig(fig, figure_path(output_path, f'{stem}{sfx}.png'))


def _plot_group_on_axes(axes, group, color, linestyle='-', show_markers=True, cols=None):
    if cols is None:
        cols = PANEL_COLS
    hab = group[group['termination'].isin(HABITABLE)]
    if hab.empty:
        return
    non_hab = group[~group['termination'].isin(HABITABLE)]
    has_da = 'da' in group.columns

    # Truncate the drawn LINE (never the markers below, which still cover the full group) at
    # the last point trustworthy enough to draw a confident trend through. `group` is always
    # instellation-sorted by its caller, so a trailing run of untrustworthy points is almost
    # always the tail end of the sweep leaving the model domain -- an out_of_domain state there
    # can be a fabricated ceiling/floor sentinel (T pinned to an exact 389.00 K clamp), so
    # drawing the line through it manufactured a cliff that read as a real, sharp trend rather
    # than "the model gave up here." A trustworthy point in the interior (bridging a
    # wall_timeout gap) is untouched -- only the trailing tail is cut.
    trustworthy_all = group['termination'].isin(DA_TRUSTWORTHY).values
    cutoff = (np.nonzero(trustworthy_all)[0].max() + 1) if trustworthy_all.any() else len(group)
    line_group = group.iloc[:cutoff]
    trustworthy = trustworthy_all[:cutoff]

    for ax, col in zip(axes, cols):
        if len(line_group) > 1:
            if has_da and linestyle == '-':
                T_sf = 1.02 * line_group['T'].values - 16.7
                at_floor = T_sf <= 274.001
                _plot_line_da_style(ax, line_group['instellation'].values, line_group[col].values,
                                    line_group['da'].values, color, at_floor=at_floor,
                                    trustworthy=trustworthy)
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
            marker = FAILED_MARKERS.get(row['termination'], 'x')
            val = row[col]
            if np.isfinite(val):
                ax.scatter(row['instellation'], val, marker=marker, s=55,
                           facecolors='none', color=color, zorder=4, linewidths=1.4)


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

def plot_basic(df, output_path, all_results=True, multiple_plots=False,
                       split_panels=True, sequence=False, width='double', height=None,
                       mg_si=None, show_hz=None):
    """The basic sweep: T, P_CO2, pH, salinity vs instellation per crust rate, coloured by outgassing.

    `mg_si` draws the same plane at a non-reference mantle Mg/Si; the value is tagged into every
    filename so the end-member figures cannot overwrite the reference ones.
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
        crust_rates     = sorted(base['crust_production'].unique())
        outgassing_vals = sorted(base['outgassing'].unique())
    else:
        crust_rates = [0.1, 1, 10]
        outgassing_vals = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = OUTGASSING_CMAP

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
                _h = _make_legend_handles(show_hz=show_hz, cols=cols)
                fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
                fig.suptitle(f'Crust production = {c}× Earth{mg_title}')
                _save_fig(fig, figure_path(output_path,
                                            f'sweep_basic_crust{c}{mg_tag}{sfx}.png'))

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
            _style_combined_col(axes_c, ci, n_cols, title=f'{c}×', cols=cols,
                                show_hz=show_hz)

        _h = _make_legend_handles(show_markers=all_results, show_hz=show_hz, cols=cols)
        fig_c.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
        _add_colorbar(fig_c, list(axes_c.ravel()), cmap, norm, 'Earth Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_rows * 10)
        fig_c.suptitle(f'Earth crust production rate{mg_title}')
        fname = f'sweep_basic{"_full" if all_results else ""}{mg_tag}{sfx}.png'
        _save_fig(fig_c, figure_path(output_path, fname), tight=all_results)

        if sequence:
            ref_crust, ref_out = 1.0, 1.0
            seq_scenarios = [
                ('single',      f'sweep_basic_seq1{mg_tag}{sfx}'),
                ('out_sweep',   f'sweep_basic_seq2{mg_tag}{sfx}'),
                ('crust_sweep', f'sweep_basic_seq3{mg_tag}{sfx}'),
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
                    _style_combined_col(axes_s, ci, n_cols, title=f'{c}×', cols=cols,
                                        show_hz=show_hz)
                _add_colorbar(fig_s, list(axes_s.ravel()), cmap, norm, 'Earth Outgassing',
                              ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                              aspect=n_rows * 10)
                _h = _make_legend_handles(show_markers=all_results, show_hz=show_hz, cols=cols)
                fig_s.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 4))
                fig_s.suptitle(f'Earth crust production rate{mg_title}')
                _save_fig(fig_s, figure_path(output_path, seq_fname + '.png'))


def plot_basic_ph(df, output_path, crust_production=1.0, width='single', height=2.2,
                  show_hz=None):
    """Abridged basic sweep: the pH panel alone at one crust production rate, coloured by outgassing."""
    base = _base(df)
    base = base[np.isclose(base['crust_production'], crust_production)]
    if base.empty:
        print(f"No basic sweep data at crust production {crust_production:g} -- skipping pH panel.")
        return
    base = _add_diag_columns(base, output_path)
    outgassing_vals = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))

    fig, ax = plt.subplots(1, 1, figsize=figure_size(width, height))
    for o in outgassing_vals:
        group = base[base['outgassing'] == o].sort_values('instellation')
        if not group.empty:
            _plot_group_on_axes([ax], group, OUTGASSING_CMAP(norm(o)), show_markers=False,
                                cols=['pH'])
    _style_axes([ax], ['pH'], show_hz=show_hz)
    ax.set_title(f'Crust production = {crust_production:g}× Earth')
    _add_colorbar(fig, ax, OUTGASSING_CMAP, norm, 'Earth Outgassing', ticks=outgassing_vals,
                  ticklabels=[f'{v}×' for v in outgassing_vals], aspect=15)
    _add_figure_legend(fig, [ax], _make_legend_handles(show_markers=False, show_hz=show_hz,
                                                       cols=['pH']))
    _save_fig(fig, figure_path(output_path, f'sweep_basic_crust{crust_production:g}_ph.png'))


def plot_basic_mgsi_grid(df, output_path, split_panels=True, show_markers=False,
                        crust_production=None, mg_si_values=None, width='double', height=None,
                        show_hz=None):
    """The basic_low_mgsi / basic_high_mgsi sweeps: Mg/Si as columns, outgassing as colour.

    plot_basic draws the outgassing x crust-production plane at ONE crust composition, so
    running it per Mg/Si gives one figure each and the compositions can only be compared by
    flipping between files. This holds crust production fixed and puts Mg/Si along the columns
    instead, which is what the low/high Mg/Si basic sweeps were run to show: whether the
    outgassing family -- the strength of the weathering feedback -- depends on crust chemistry.

    Distinct from plot_cross, which colours BY Mg/Si at a single (outgassing, crust)
    pair: that shows the composition effect for one tectonic state, this shows whether the
    composition effect survives across the outgassing axis.
    """
    mg_vals = basic_plane_mg_si(df) if mg_si_values is None else sorted(mg_si_values)
    if len(mg_vals) < 2:
        print("Fewer than 2 Mg/Si values with a basic sweep -- skipping the Mg/Si grid.")
        return

    # Earth crust production if it was run, else whichever rate covers the most Mg/Si columns.
    if crust_production is None:
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

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = OUTGASSING_CMAP
    n_cols = len(mg_vals)

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figure_size(width, height, n_rows),
                                 sharex=True, sharey='row', squeeze=False)
        for ci, m in enumerate(mg_vals):
            sub = subsets[m]
            for o in outgassing_vals:
                group = sub[sub['outgassing'] == o].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes[:, ci], group, cmap(norm(o)),
                                        show_markers=show_markers, cols=cols)
            title = f'{m:g}' + (' (Earth)' if np.isclose(m, REF_MG_SI) else '')
            _style_combined_col(axes, ci, n_cols, title=title, cols=cols, show_hz=show_hz)

        _add_colorbar(fig, list(axes.ravel()), cmap, norm, 'Earth Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers,
                                                      show_hz=show_hz, cols=cols))
        fig.suptitle(f'Mantle Mg/Si   (crust production = {crust_production:g}× Earth)')
        _save_fig(fig, figure_path(output_path, f'sweep_basic_mgsi_grid{sfx}.png'))


def plot_depth(df, output_path, show_markers=False, split_panels=True,
                            width='single', height=None, show_hz=None):
    """The depth sweep: T, P_CO2, pH, salinity vs instellation, coloured by ocean depth."""
    # The depth sweep fixes (outgassing, crust) at their sweep defaults and varies ocean_depth.
    # That default changed over time (outgassing 1.0 -> 0.1), so instead of hardcoding a value
    # (which silently produced an empty/one-depth plot), pick whichever (outgassing, crust) pair
    # actually spans the most distinct depths.
    pool = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_chem(df) &
        _ref_redox(df) &
        (df['land_fraction'] == 0.0)
    ]
    if pool.empty:
        print("No data for ocean depth sweep — skipping.")
        return
    picked = _best_operating_point(pool, 'ocean_depth', 'Ocean-depth')
    if picked is None:
        return
    raw_subset = picked[0]

    # All 10 swept depths makes this figure busy without adding much. Drop 50 km outright --
    # it's the least trustworthy point in the sweep (development_history.md 29.3: 50 km is
    # non-monotonic against 30 km at S=1.0, and its runs concentrate the wall_timeout tail) --
    # then thin the rest to ~5 evenly log-spaced depths so the remaining lines stay legible.
    all_depths = sorted(d for d in raw_subset['ocean_depth'].unique() if d < 50000)
    if len(all_depths) > 5:
        idx = np.linspace(0, len(all_depths) - 1, 5).round().astype(int)
        show_depths = sorted({all_depths[i] for i in idx})
    else:
        show_depths = all_depths
    print(f"  Ocean-depth plot: showing {[f'{d:g}' for d in show_depths]} m "
          f"(dropped {[f'{d:g}' for d in all_depths if d not in show_depths]} m, "
          f"and 50000 m outright).")
    subset = _add_diag_columns(raw_subset[raw_subset['ocean_depth'].isin(show_depths)],
                               output_path)

    depths = sorted(subset['ocean_depth'].unique())
    # pad=0: depth is a physical range that should map to the colour map exactly.
    norm = _value_norm(depths, pad=0.0)
    cmap = DEPTH_CMAP
    # Ticks stay in metres (the norm's native units); only the printed labels convert to km --
    # 300-50000 m as ten tick labels was the widest, most cluttered colorbar of any figure here.
    ticks, _ = _colorbar_ticks(depths)
    ticklabels = [f'{v / 1000:g}' for v in ticks] if ticks is not None else None
    _faceted_lines(subset, 'ocean_depth', depths, [cmap(norm(d)) for d in depths],
                   cmap, norm, 'Ocean Depth (km)', 'sweep_depth', output_path,
                   split_panels=split_panels, show_markers=show_markers,
                   x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                   aspect_per_row=7.5, width=width, height=height, show_hz=show_hz)


def plot_chemistry(df, output_path, show_markers=False, split_panels=True,
                             width='single', height=None, show_hz=None):
    """The alpha and chemistry sweeps: T, P_CO2, pH, salinity vs instellation per constant value.

    One figure per constant that actually varies (alpha, kd_mg, k_na); the other two are held
    at their most common value so each figure isolates a single knob.
    """
    pool_all = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_redox(df) &
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
        cmap = CHEM_KNOB_CMAP
        ticks, ticklabels = _colorbar_ticks(values)
        # alpha is its own named sweep (parameter_sweep.sweep_alpha); kd_mg and k_na are the
        # two arms of sweep_chemistry, so they are filed under that name.
        stem = 'sweep_alpha' if col == 'alpha' else f'sweep_chemistry_{col}'
        _faceted_lines(subset, col, values, [cmap(norm(v)) for v in values],
                       cmap, norm, CHEM_KNOBS[col], stem, output_path,
                       split_panels=split_panels, show_markers=show_markers,
                       x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                       aspect_per_row=7.5, width=width, height=height, show_hz=show_hz)


def plot_pe(df, output_path, show_markers=False, split_panels=True,
                      ocean_depth=3000, width='single', height=None, show_hz=None):
    """The pe sweep: T, P_CO2, pH, salinity vs instellation for each ocean redox state (pe).

    `pe` (an abiotic, reducing ocean vs an oxidised one) is not a continuous knob like the
    CHEM_KNOBS -- it switches the iron sink between Siderite and Goethite -- but it is swept the
    same way (parameter_sweep.sweep_pe): every other axis held at the Earth reference while pe
    runs from oxic seawater to below the Goethite saturation boundary. See
    development_history.md section 28/31 for why the model treats ocean redox as a free
    parameter at all.
    """
    pool = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_chem(df) &
        (df['ocean_depth'] == ocean_depth) &
        (df['land_fraction'] == 0.0)
    ]
    if pool.empty:
        print(f"No data for redox sweep at depth {ocean_depth:g} m — skipping.")
        return

    picked = _best_operating_point(pool, 'pe', 'pe')
    if picked is None:
        return
    subset = _add_diag_columns(picked[0], output_path)

    values = sorted(subset['pe'].unique())
    norm = _value_norm(values)
    cmap = PE_CMAP
    ticks, ticklabels = _colorbar_ticks(values)
    tag = '' if ocean_depth == 3000 else f'_d{ocean_depth:g}'
    _faceted_lines(subset, 'pe', values, [cmap(norm(v)) for v in values],
                   cmap, norm, r'Ocean redox $p_e$', f'sweep_pe{tag}', output_path,
                   split_panels=split_panels, show_markers=show_markers,
                   x_lims=_x_limits(subset), ticks=ticks, ticklabels=ticklabels,
                   aspect_per_row=7.5, width=width, height=height, show_hz=show_hz)


def _composition_pool(df, ocean_depth=3000):
    """Runs usable for a composition figure: reference chemistry, land-free, one depth."""
    return df[
        df['reverse_weathering'] &
        _ref_chem(df) &
        _ref_redox(df) &
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


def plot_cross(df, output_path, split_panels=True, show_markers=False,
                           ocean_depth=3000, width='single', height=None, show_hz=None):
    """The cross sweep: T, P_CO2, pH, salinity vs instellation, one figure per composition axis.

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
    subset = subset[~np.isclose(subset['mg_si'].to_numpy()[:, None], MG_SI_HIDDEN).any(axis=1)]

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

        cmap = PARAM_CMAPS[key]
        numeric = [float(v) for v in values]
        norm = _value_norm(numeric)
        cbar_label, ticklabels = label, [fmt(v) for v in values]

        _faceted_lines(cut, key, values, [cmap(norm(n)) for n in numeric],
                       cmap, norm, cbar_label, f'sweep_cross_{key}{tag}', output_path,
                       split_panels=split_panels, show_markers=show_markers,
                       ticks=numeric, ticklabels=ticklabels, width=width, height=height,
                       show_hz=show_hz)
        if key == 'mg_si':
            # Abridged pH-only version.
            _faceted_lines(cut, key, values, [cmap(norm(n)) for n in numeric],
                           cmap, norm, cbar_label, f'sweep_cross_{key}{tag}', output_path,
                           show_markers=show_markers, ticks=numeric, ticklabels=ticklabels,
                           width=width, height=(height if height is not None else 2.2),
                           show_hz=show_hz, groups=[(['pH'], '_ph')], aspect_per_row=15)


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


def plot_composition(df, output_path, split_panels=True, show_markers=False,
                          ocean_depth=3000, min_lines=2, n_cols=3,
                          width='double', height=None, show_hz=None):
    """The composition sweep, in the style of the basic one: dIW as columns, Mg/Si as colour.

    Same layout as `sweep_basic_full` (crust rate as columns, outgassing as colour), with the
    two crust axes substituted. This uses the WHOLE factorial rather than the one-axis cuts
    `plot_cross` draws: every panel is a fixed dIW, and every line within it a fixed
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
    subset = subset[~np.isclose(subset['mg_si'].to_numpy()[:, None], MG_SI_HIDDEN).any(axis=1)]
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
    cmap = MG_SI_CMAP
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
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers,
                                                      show_hz=show_hz, cols=cols))
        fig.suptitle(r'Core-formation $\Delta$IW')
        _save_fig(fig, figure_path(output_path, f'sweep_composition{tag}{sfx}.png'))

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
            _style_combined_col(axes, ci, len(mg_cols), title=f'{mg:g}', cols=cols,
                               show_hz=show_hz)
        _add_colorbar(fig, list(axes.ravel()), cmap_d, norm_d, r'Core-formation $\Delta$IW',
                      ticks=diw_all, ticklabels=[f'{v:+g}' for v in diw_all],
                      aspect=n_rows * 10)
        _add_figure_legend(fig, axes, _make_legend_handles(show_markers=show_markers,
                                                      show_hz=show_hz, cols=cols))
        fig.suptitle('Mantle Mg/Si')
        _save_fig(fig, figure_path(output_path, f'sweep_composition_T{tag}{sfx}.png'))


def plot_composition_map(df, output_path, s_vals=(0.7, 0.9, 1.0, 1.1), ocean_depth=3000,
                         quantity='T', relative=True, min_cells=3,
                         width='double', height=2.6):
    """The composition sweep as a map over the Mg/Si x dIW plane, one panel per instellation.

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
        cmap = RELATIVE_CMAP
        cbar_label = f'{short} vs Earth crust' + (f' ({unit})' if unit else '')
    else:
        norm = mcolors.Normalize(vmin=np.nanpercentile(allv, 2), vmax=np.nanpercentile(allv, 98))
        cmap = QUANTITY_CMAP if log_q else RELATIVE_CMAP
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
    _save_fig(fig, figure_path(output_path,
                                f'sweep_composition_map_{quantity}{"" if relative else "_abs"}.png'))


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
    cmap = INSTELLATION_CMAP

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
    _save_fig(fig, figure_path(output_path, 'ratio_scatter.png'))


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
        _ref_redox(df) &
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
    cmap = DIVERGING_CMAP
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
        _save_fig(fig, figure_path(output_path, 'da_contour.png'))


def plot_continental_baseline(df, output_path, show_hz=None):
    """T, P_CO2, pH, salinity, and individual ion concentrations vs instellation
    for the Earth-like continental baseline.

    Runs with land_fraction=0.3, reference crust, rw=True, out=1×, crust=1×, depth=3000 m.
    """
    subset = df[
        (df['land_fraction'] == 0.3) &
        _ref_crust(df) &
        _ref_chem(df) &
        df['reverse_weathering'] &
        _ref_redox(df) &
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0) &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0)
    ]
    if subset.empty:
        print("No continental baseline data found — skipping.")
        return

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
    # (index, label, Earth mmol/kg). Al(3) and Fe(4) are excluded: both are trace and the model
    # holds them at ~0, so a log axis has nothing to show. SO4(9) is excluded because it is PINNED
    # (planet.py F_net[so4_idx] = 0) -- whatever the seed carries is what it reports for ever, so
    # plotting it would claim credit for an input. The `pinned` flag and its italic-label handling
    # are kept wired so it can be put back in one line if the sulfur cycle is ever closed.
    #
    # Earth values are calibrate_earth.py's own fit targets, so the figure is scored against the
    # same numbers the calibration was. They are the standard S = 35 seawater values per kg of
    # solution; Na 469 / Cl 546 rather than the 480 / 550 previously hardcoded here.
    ION_SPEC = [
        (0, 'Alk',   2.3,   False),
        (1, 'C',     2.1,   False),
        (2, 'Si',    0.1,   False),
        (5, 'Ca',   10.3,   False),
        (6, 'Mg',   52.8,   False),
        (7, 'Na',  469.0,   False),
        (8, 'Cl',  546.0,   False),
    ]
    # Colour separates the two SERIES (Earth vs model), not the ions -- the ion is already given
    # by x position, so per-ion colour was encoding nothing. Marker shape repeats the distinction
    # so identity never rests on colour alone.
    C_EARTH, C_MODEL = '#4C72B0', '#C44E52'

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

    # --- figure 1: summary panels, split T/P_CO2 from pH/salinity like every other figure ---
    for grp_cols, sfx in _panel_groups(True):
        n_rows = len(grp_cols)
        fig, axes = plt.subplots(n_rows, 1, sharex=True,
                                 figsize=figure_size('single', n_rows=n_rows, row_height=1.5))
        if n_rows == 1:
            axes = [axes]

        # T panel shows all runs (including failed); the rest show only runs in habitable T range
        if 'T' in grp_cols:
            _plot_group_on_axes([axes[grp_cols.index('T')]], group_all, color='k',
                                show_markers=False, cols=['T'])
        other = [c for c in grp_cols if c != 'T']
        if other:
            _plot_group_on_axes([axes[grp_cols.index(c)] for c in other], group_hab, color='k',
                                show_markers=False, cols=other)
        _style_axes(axes, grp_cols, show_hz=show_hz, show_eq_temp=True)

        for ax, col in zip(axes, grp_cols):
            ax.scatter(EARTH_S, earth_vals[col], marker='*', s=220, color='blue',
                       edgecolors='k', linewidths=0.7, zorder=6)
        if 'T' in grp_cols:
            axes[grp_cols.index('T')].annotate(
                'Earth', xy=(EARTH_S, EARTH_T), xytext=(EARTH_S + 0.06, EARTH_T - 6),
                fontsize=8, arrowprops=dict(arrowstyle='-', color='k', lw=0.8),
            )
        _save_fig(fig, figure_path(output_path, f'continental_baseline{sfx}.png'))

    # --- figure 2: absolute ion concentrations, model vs Earth seawater at S ~ 1 ---
    # A dumbbell rather than the previous percent-difference scatter: the quantity of interest is
    # the concentration itself, and the ions span five decades (Si 0.1 to Cl 546 mM), so a log
    # axis showing both values directly is more honest than a ratio -- a "+40%" and a "+4000%"
    # tell you nothing about whether either number is physically reasonable.
    figsize2 = figure_size('single', height=2.6)
    fig2, ax_ions = plt.subplots(1, 1, figsize=figsize2)

    if ion_rows:
        s_arr  = np.array([r[0] for r in ion_rows])
        b_mmol = np.array([r[1] for r in ion_rows]) * 1e3  # mol/kg -> mmol/kg
        closest = int(np.argmin(np.abs(s_arr - EARTH_S)))
        b_model = b_mmol[closest]

        labels = [spec[1] for spec in ION_SPEC]
        earth  = np.array([spec[2] for spec in ION_SPEC])
        model  = np.array([max(b_model[spec[0]], 1e-4) for spec in ION_SPEC])
        pinned = np.array([spec[3] for spec in ION_SPEC])
        x      = np.arange(len(labels))

        # Connector first, so the markers sit on top of it.
        ax_ions.vlines(x, np.minimum(earth, model), np.maximum(earth, model),
                       color='0.6', linewidth=1.0, zorder=2)
        ax_ions.scatter(x, earth, marker='o', s=46, facecolors='none', edgecolors=C_EARTH,
                        linewidths=1.4, zorder=3, label='Earth seawater')
        ax_ions.scatter(x, model, marker='D', s=34, color=C_MODEL, edgecolors='w',
                        linewidths=0.5, zorder=4, label='Model (calibrated)')

        ax_ions.set_xticks(x)
        ax_ions.set_xticklabels(labels)
        # Mark the pinned species: SO4 has no source term, so its agreement is an input.
        for xi, is_pin in zip(x, pinned):
            if is_pin:
                ax_ions.get_xticklabels()[xi].set_style('italic')
                ax_ions.get_xticklabels()[xi].set_color('0.45')

        ax_ions.axvline(2.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6, zorder=1)
        trans = ax_ions.get_xaxis_transform()
        ax_ions.text(1.0, 1.02, 'Biotically controlled', transform=trans,
                     ha='center', va='bottom', fontsize=8)
        ax_ions.text(4.5, 1.02, 'Abiotically controlled', transform=trans,
                     ha='center', va='bottom', fontsize=8)
        ax_ions.set_xlim(-0.6, len(labels) - 0.4)

    ax_ions.set_yscale('log')
    ax_ions.set_ylabel('Concentration (mmol kg$^{-1}$)')
    ax_ions.set_ylim(0.05, 1.0e3)   # data spans Si 0.1 to Cl 546
    ax_ions.spines['top'].set_visible(False)
    ax_ions.spines['right'].set_visible(False)
    ax_ions.grid(True, linestyle='--', alpha=0.35, axis='y', which='major')
    ax_ions.legend(frameon=False, fontsize=7, loc='lower right', handletextpad=0.4,
                   borderaxespad=0.6)
    _save_fig(fig2, figure_path(output_path, 'continental_baseline_ions.png'))


# ---------------------------------------------------------------------------
# Continental baseline, land-fraction series and crossover figures
# ---------------------------------------------------------------------------
# Figures for the runs experiments/continental_baseline.py produces. The sweep design (land
# fractions, grid axes, Earth reference values) is read from that module as `cb`.

# Land-free blue against continental brown, the one colour decision these figures make.
ARM_COLOURS = {0.3: '#a4632a', 0.0: '#2a6fa4'}
ARM_LABELS = {0.3: 'Continental (land fraction 0.3)', 0.0: 'Ocean world (land free)'}

# Modern Earth, for the reference marker. Salinity is the sum of the model's tracked ions at
# their seawater concentrations, so it is comparable with the model's own salinity column.
EARTH = {'S': 1.0, 'T': 288.0, 'P_CO2': 280e-6, 'pH': 8.1,
         'salinity': (2.0e-3 * 61.0 + 0.1e-3 * 60.1 + 10.3e-3 * 40.1 +
                      52.8e-3 * 24.3 + 480e-3 * 23.0 + 550e-3 * 35.45)}


# When a run leaves the validity box the model records the box's own limit rather than a computed
# temperature -- an exact 389 K or 400 K at the hot end, 181 K at the cold end. plot_results says
# the same of the Da those states carry. Drawing a curve through them manufactures a plateau that
# reads as physics, so they are dropped from any curve drawn here.
T_CLAMP_HOT = 389.0
T_CLAMP_COLD = 181.0

# Validity of the climate model's OLR parameterisation (kamino.climate.analytic), which is the
# Haqq-Misra et al. (2016) polynomial fit to the Kopparapu et al. (2013, 2014) 1-D
# radiative-convective columns: 1e-5 bar < pCO2 < 10 bar and 150 K < T < 350 K, error <= 3.3 W/m2.
# A state outside that box is an extrapolation of the fit, not a prediction of the model.
OLR_FIT_T_MAX = 350.0


def _olr_limit(pco2_bar):
    """First local maximum of OLR(T) -- the Simpson-Nakajima radiation limit for this atmosphere.

    OLR is NOT monotonic in T: water vapour makes it plateau near 271 W/m2 (at low CO2) and then
    fall before the hot branch climbs again. Instellation above that plateau admits no cool-branch
    solution, which is the runaway greenhouse.
    """
    from kamino.climate.analytic import OLR
    peak = OLR(180.0, pco2_bar)
    for T in np.arange(181.0, 391.0, 1.0):
        v = OLR(float(T), pco2_bar)
        if v < peak:
            break
        peak = v
    return peak


def _past_runaway(S, pco2_bar, albedo=0.3):
    """True when absorbed instellation exceeds the OLR limit, i.e. the planet is in runaway.

    This is the check `get_T_surface_analytic` does NOT make. When no cool-branch root exists it
    returns the first sign change it finds, which lies on the HOT branch beyond the runaway --
    a number near 357 K that a plain `T < 360` habitability test happily accepts. Measured on
    this grid that put the continental inner edge at S = 1.15, one grid point too far.
    """
    from kamino.climate.analytic import albedo_funtion
    from kamino.constants import SOLAR_CONSTANT
    pco2_bar = max(float(pco2_bar), 1e-5)      # the model's own 1 Pa CO2 floor
    A = albedo_funtion(pco2_bar, albedo)
    return S * SOLAR_CONSTANT * (1 - A) * 0.25 > _olr_limit(pco2_bar)


def _drop_clamped(group):
    """Drop out-of-domain rows whose stored T is a box limit rather than a computed value."""
    clamped = (group['termination'].isin(OUT_OF_DOMAIN) &
               ((group['T'] >= T_CLAMP_HOT) | (group['T'] <= T_CLAMP_COLD)))
    return group[~clamped]


def _draw_arm(axes, group, colour, cols):
    """Draw one arm, and report whether any of it is habitable.

    `plot_results._plot_group_on_axes` returns without drawing when a group contains no habitable
    run at all, which is the right call for a facet of a larger figure but wrong here: at Earth's
    outgassing rate EVERY land-free run leaves the domain, and silently omitting the line would
    hide the very comparison this figure exists to make.

    Such an arm is drawn faint, with its clamp sentinels removed and plot_results' hollow
    per-termination markers on every point, so it reads as "measured, but not habitable". The
    line stays SOLID deliberately: the line styles are already spoken for by DA_LEGEND, so a
    dashed fallback would read as "Da >= 1" rather than "not habitable". Colour separates the
    arms; the hollow markers say these are not habitable states.
    """
    if group['termination'].isin(HABITABLE).any():
        _plot_group_on_axes(axes, group, colour, show_markers=False, cols=cols)
        return True
    shown = _drop_clamped(group)
    for ax, col in zip(axes, cols):
        ax.plot(shown['instellation'], shown[col], color=colour, linewidth=1.2, alpha=0.5,
                zorder=2)
        for _, row in shown.iterrows():
            if np.isfinite(row[col]):
                ax.scatter(row['instellation'], row[col],
                           marker=FAILED_MARKERS.get(row['termination'], 'x'), s=22,
                           facecolors='none', edgecolors=colour, linewidths=1.0, zorder=4)
    return False


def _arm(df, land):
    """Rows for one land-fraction arm of the baseline, at the Earth reference on every other axis."""
    return df[
        _ref_crust(df) &
        _ref_redox(df) &
        _ref_chem(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == cb.OCEAN_DEPTH) &
        (df['outgassing'] == cb.OUTGASSING) &
        (df['crust_production'] == cb.CRUST_PRODUCTION) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['land_fraction'], land)
    ].sort_values('instellation')


def _alk_fluxes(group):
    """Continental and seafloor alkalinity flux at each run's final state, Tmol eq/yr.

    Both are put on the SAME basis -- the flux the ODE actually applies -- so the ratio means
    what it looks like:

    * continental is `get_continental_weathering_flux(T, pCO2)` over `land_fraction * surface`,
      which is how planet.py applies it. On modern Earth this is 8 Tmol eq/yr by calibration
      (constants.EARTH_CONTINENTAL_WEATHERING_REF), so the number is readable on sight.
    * seafloor is the recorded `alk_flux` diagnostic rescaled from the FIXED reference area it is
      stored on to the area the model actually integrates over. planet.py normalises the
      diagnostic on A_SEAFLOOR_EARTH (0.7 of the surface, a constant) so that it always agrees
      with plot_results, but `F_diss` is applied over `seafloor_area = (1 - land_fraction) * A`.
      The conversion is therefore x (1 - land_fraction) / EARTH_OCEAN_FRACTION -- exactly 1 at
      Earth's land fraction, and 1.43x on a land-free world.

    Returns (continental, seafloor) arrays aligned with `group`.
    """
    from kamino.weathering import get_continental_weathering_flux
    from kamino.chemistry import alk_idx
    from kamino.constants import YR, R_EARTH, EARTH_OCEAN_FRACTION

    surface = 4 * np.pi * R_EARTH ** 2
    cont = np.full(len(group), np.nan)
    for i, (_, r) in enumerate(group.iterrows()):
        T, p, land = r['T'], r['P_CO2'], r['land_fraction']
        if not (np.isfinite(T) and np.isfinite(p)) or land <= 0:
            cont[i] = 0.0 if land <= 0 else np.nan
            continue
        f = get_continental_weathering_flux(float(T), float(p) * 1e5)   # pCO2 stored in bar
        cont[i] = float(f[alk_idx]) * land * surface * YR / 1e12
    sea = (group['alk_flux'].to_numpy(dtype=float)
           * (1.0 - group['land_fraction'].to_numpy(dtype=float)) / EARTH_OCEAN_FRACTION)
    return cont, sea


def _crossover_land_fraction(lands, ratios):
    """Land fraction where the two alkalinity fluxes are equal, log-interpolated.

    `ratios` may be given either way up. The crossing sits where log10(ratio) = 0, and inverting
    every ratio flips the sign of both the numerator and the denominator of the interpolation
    weight, so the land fraction it returns is identical. Callers here pass continental/seafloor
    in one place and seafloor/continental in the other; both are correct.

    Returns None when the sampled land fractions do not bracket a crossing -- the sweep then
    bounds the crossover rather than locating it, which the caller must say rather than
    extrapolate off the end of the grid.
    """
    pairs = sorted((float(l), float(r)) for l, r in zip(lands, ratios)
                   if l > 0 and np.isfinite(r) and r > 0)
    for (l0, r0), (l1, r1) in zip(pairs, pairs[1:]):
        if (r0 - 1.0) * (r1 - 1.0) <= 0 and r0 != r1:
            w = (0.0 - np.log10(r0)) / (np.log10(r1) - np.log10(r0))
            return float(10 ** (np.log10(l0) + w * (np.log10(l1) - np.log10(l0))))
    return None


def hz_edges(group):
    """Instellation limits of the habitable band along one instellation line.

    Returns ``(S_outer, S_inner, outer_kind, inner_kind)``, or None where no run on the line is
    habitable. A run counts as habitable when its integration is trustworthy (converged, or ran
    to 2 Gyr) AND its final surface temperature lies between the snowball and runaway
    thresholds -- the same two numbers `plot_results._style_axes` draws as walls.

    How each edge was located is reported rather than assumed, because the three cases are not
    equally good and the difference matters when the number is quoted:

    ``crossing``   both bracketing runs are trustworthy, so T(S) is interpolated onto the
                   threshold. This is a measurement.
    ``bracketed``  the neighbour left the model domain at the matching wall -- frozen below the
                   outer edge, runaway above the inner one. That is a real outcome, but its
                   stored T is a clamp sentinel (an exact 181 K or 389 K), so interpolating
                   through it would invent a slope. The edge is placed at the midpoint of the
                   grid interval and is uncertain by half a grid step.
    ``open``       there is no neighbour (the sweep ran out of range) or the neighbour's
                   integration is not trustworthy. The edge is the last habitable grid point and
                   is a BOUND -- the true edge is at least this far out.
    """
    g = group.sort_values('instellation')
    S = g['instellation'].to_numpy(dtype=float)
    T = g['T'].to_numpy(dtype=float)
    wall = (g['domain_wall'].to_numpy(dtype=object) if 'domain_wall' in g
            else np.full(len(g), None, dtype=object))
    trusted = g['termination'].isin(HABITABLE).to_numpy() & np.isfinite(T)

    # A temperature window alone is not a habitability test in this model. Two states pass
    # `T < T_RUNAWAY` without being habitable at all: one past the runaway greenhouse, whose T is
    # read off the hot branch (see `_past_runaway`), and one above the OLR fit's 350 K ceiling,
    # where the climate model is extrapolating. Both are excluded here rather than by moving
    # T_RUNAWAY, which is a plot_results convention shared with every other figure.
    S_ok = np.array([not _past_runaway(s, p) for s, p in
                     zip(S, g['P_CO2'].to_numpy(dtype=float))])
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
            # `not S_ok` is the runaway greenhouse: the neighbour has no cool-branch solution, so
            # the edge lies in this interval. Not interpolated -- the neighbour's T is on the hot
            # branch, so a line drawn through it has no meaning.
            inner, inner_kind = 0.5 * (S[i1] + S[i1 + 1]), 'bracketed'

    return float(outer), float(inner), outer_kind, inner_kind


def plot_baseline_vs_ocean(arms, output_path):
    """T, pCO2, pH and salinity against instellation, continental arm against ocean arm.

    Both lines carry plot_results' Damkohler styling, so the comparison also shows whether the
    two arms sit in the same weathering regime -- which they do not: continental weathering is
    transport-limited on Earth (Da >> 1) while the land-free worlds are kinetically limited.
    """
    habitable = {land: group['termination'].isin(HABITABLE).any()
                 for land, group in arms.items()}
    handles = [Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6,
                         alpha=1.0 if habitable[l] else 0.5,
                         marker='' if habitable[l] else 's', markerfacecolor='none',
                         label=ARM_LABELS[l] + ('' if habitable[l]
                                                else ' — never habitable'))
               for l in arms] + list(DA_LEGEND)

    for cols, sfx in _panel_groups(True):
        fig, axes = plt.subplots(len(cols), 1, sharex=True,
                                    figsize=figure_size('single', n_rows=len(cols),
                                                           row_height=2.0))
        for land, group in arms.items():
            _draw_arm(axes, group, ARM_COLOURS[land], cols)
        _style_axes(axes, cols)
        for ax, col in zip(axes, cols):
            ax.scatter(EARTH['S'], EARTH[col], marker='*', s=180, color='gold',
                       edgecolors='k', linewidths=0.7, zorder=6)
        _add_figure_legend(fig, axes, handles)
        _save_fig(fig, figure_path(output_path, f'continental_vs_ocean{sfx}.png'))


def plot_habitable_zone(arms, output_path):
    """The headline figure: where the model keeps a planet temperate, with land and without.

    Upper panel is the temperature curve each zone is read off; lower panel is the zone itself,
    one bar per arm on the same instellation axis. Edges located only as bounds (see `hz_edges`)
    carry a caret pointing the way the true edge lies, so a bar that is merely wider than the
    sweep could resolve cannot be read as a measured one.
    """
    edges = {}
    for land, group in arms.items():
        got = hz_edges(group)
        if got is not None:
            edges[land] = got
    if cb.LAND_FRACTION not in edges:
        print("No habitable band on the continental arm -- skipping the habitable-zone figure.")
        return edges

    fig, (ax, ax_z) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 1],
                                      figsize=figure_size('single', height=4.0))

    habitable = {}
    for land, group in arms.items():
        habitable[land] = _draw_arm([ax], group, ARM_COLOURS[land], ['T'])
    _style_axes([ax], ['T'])
    ax.set_xlabel('')
    ax.scatter(EARTH['S'], EARTH['T'], marker='*', s=180, color='gold', edgecolors='k',
               linewidths=0.7, zorder=6)

    # Every arm gets a row, including one with no habitable band at all -- that is the result at
    # Earth outgassing, and a missing row would read as a missing run rather than an empty zone.
    for row, land in enumerate(arms):
        colour = ARM_COLOURS[land]
        y = len(arms) - 1 - row
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

    handles = [Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6,
                         alpha=1.0 if habitable[l] else 0.5,
                         marker='' if habitable[l] else 's', markerfacecolor='none',
                         label=ARM_LABELS[l] + ('' if habitable[l]
                                                else ' — never habitable'))
               for l in arms]
    if any(k == 'open' for e in edges.values() for k in e[2:]):
        handles.append(Line2D([0], [0], color='k', linestyle='none', marker='>', markersize=5,
                                 label='Edge is a bound (sweep limit)'))
    _add_figure_legend(fig, [ax, ax_z], handles)
    _save_fig(fig, figure_path(output_path, 'continental_habitable_zone.png'))
    return edges


def _report(arms, edges):
    """Print the habitable-zone edges, with how each was located. See `hz_edges`."""
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

    # plot_results draws these edges as vertical lines on every other instellation figure, from
    # its own hardcoded copy. Say so loudly when the two disagree rather than let every figure
    # in the paper quote a stale zone -- the same reason parameter_sweep._warn_constant_drift
    # exists for the chemistry constants.
    if cb.LAND_FRACTION in edges:
        lo, hi = edges[cb.LAND_FRACTION][:2]
        for label, here, there in (('CONTINENTAL_HZ_OUTER', lo, CONTINENTAL_HZ_OUTER),
                                   ('CONTINENTAL_HZ_INNER', hi, CONTINENTAL_HZ_INNER)):
            if abs(here - there) > 5e-4:
                print(f"  NOTE plot_results.{label} = {there:g}, but this sweep measures "
                      f"{here:.3f}. Update it, or the HZ lines on every other figure are stale.")


def _land_series(df, output_path):
    """{land_fraction: instellation-sorted rows with diagnostics}, for the land fractions present."""
    series = {}
    for land in cb.LAND_FRACTIONS:
        sub = _arm(df, land)
        if not sub.empty:
            series[land] = _add_diag_columns(sub, output_path).sort_values('instellation')
    return series


def _land_colours(lands):
    """Colour per land fraction, log-scaled over the positive ones.

    0 cannot sit on a log scale and is not just 'a bit less land' -- it is the land-free ocean
    world every other sweep runs. It keeps the baseline figures' blue and its own legend entry.
    """
    positive = sorted(l for l in lands if l > 0)
    cmap = LAND_FRACTION_CMAP
    # A single positive value gives LogNorm(vmin == vmax), which cannot be normalised. That
    # happens whenever only the baseline's own land fraction is on disk -- i.e. before the
    # land-fraction sweep has been run -- so it is the ordinary case, not an error.
    norm = (mcolors.LogNorm(vmin=min(positive), vmax=max(positive))
            if len(positive) > 1 else None)
    if norm is not None:
        colours = {l: cmap(norm(l)) for l in positive}
    else:
        colours = {l: ARM_COLOURS[cb.LAND_FRACTION] for l in positive}
    if 0.0 in lands:
        colours[0.0] = ARM_COLOURS[0.0]
    return colours, cmap, norm


def plot_land_fraction_series(series, output_path):
    """T, pCO2, pH and salinity against instellation, one line per land fraction."""
    if len(series) < 2:
        print("Fewer than two land fractions on disk -- skipping the land-fraction series.")
        return
    colours, cmap, norm = _land_colours(series)

    for cols, sfx in _panel_groups(True):
        fig, axes = plt.subplots(len(cols), 1, sharex=True,
                                    figsize=figure_size('single', n_rows=len(cols),
                                                           row_height=2.0))
        for land, group in sorted(series.items(), reverse=True):
            _draw_arm(axes, group, colours[land], cols)
        _style_axes(axes, cols)
        if norm is not None:
            positive = sorted(l for l in series if l > 0)
            _add_colorbar(fig, list(axes), cmap, norm, 'Land fraction',
                             ticks=positive, ticklabels=[f'{v:g}' for v in positive],
                             aspect=len(cols) * 7.5)
        handles = []
        if 0.0 in series:
            handles.append(Line2D([0], [0], color=ARM_COLOURS[0.0], linewidth=1.6,
                                     label='Land free (0)'))
        handles += list(DA_LEGEND)
        _add_figure_legend(fig, axes, handles)
        _save_fig(fig, figure_path(output_path, f'land_fraction_series{sfx}.png'))


def plot_weathering_crossover(series, output_path):
    """Where seafloor weathering overtakes continental weathering as land fraction falls.

    Left: both alkalinity fluxes against land fraction at the instellation nearest Earth's.
    Right: the crossover land fraction across instellation.

    Only runs that reached a steady state (converged, or integrated to 2 Gyr) are used. A run
    stopped at a domain wall is a planet still evolving when the model gave up, and its fluxes
    are not a balance of anything -- reading a crossover off one would be reading it off a
    transient.
    """
    if len(series) < 2:
        print("Fewer than two land fractions on disk -- skipping the crossover figure.")
        return None

    # (land, S) -> (continental, seafloor), steady states only.
    rows, dropped = {}, 0
    for land, group in series.items():
        steady = group[group['termination'].isin(HABITABLE)]
        dropped += len(group) - len(steady)
        if steady.empty:
            continue
        cont, sea = _alk_fluxes(steady)
        for s, c, f in zip(steady['instellation'], cont, sea):
            if np.isfinite(c) and np.isfinite(f) and f > 0:
                rows[(float(land), float(s))] = (float(c), float(f))
    if not rows:
        print("No steady-state runs to compare fluxes on -- skipping the crossover figure.")
        return None
    if dropped:
        print(f"  crossover: ignoring {dropped} run(s) that never reached a steady state.")

    s_vals = sorted({s for _, s in rows})
    crossings = {}
    for s in s_vals:
        lands = [l for (l, ss) in rows if ss == s]
        ratios = [rows[(l, s)][0] / rows[(l, s)][1] for l in lands]
        got = _crossover_land_fraction(lands, ratios)
        if got is not None:
            crossings[s] = got

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
        ax_c.plot(list(crossings), [crossings[s] for s in crossings], color='k',
                  marker='o', markersize=3, linewidth=1.4)
    ax_c.set_yscale('log')
    ax_c.set_xlabel('Instellation (S/S₀)')
    ax_c.set_ylabel('Crossover land fraction')
    ax_c.grid(True, linestyle='--', alpha=0.4)
    _save_fig(fig, figure_path(output_path, 'weathering_crossover.png'))

    print("\nContinental vs seafloor alkalinity flux (Tmol eq/yr), steady states only:")
    print(f"  {'land':>8} " + ' '.join(f'{s:>9.2f}' for s in s_vals))
    for land in sorted(series, reverse=True):
        cells = []
        for s in s_vals:
            v = rows.get((land, s))
            cells.append(f"{v[0] / v[1]:9.3g}" if v else f"{'--':>9}")
        print(f"  {land:8g} " + ' '.join(cells))
    print("  (continental / seafloor; < 1 means seafloor weathering dominates)")
    if crossings:
        lo, hi = min(crossings.values()), max(crossings.values())
        print(f"  crossover land fraction: {lo:.3g} to {hi:.3g} over S = "
              f"{min(crossings):g}-{max(crossings):g}")
    else:
        ratios_all = [c / f for c, f in rows.values()]
        print(f"  no crossing inside the sampled land fractions -- ratio spans "
              f"{min(ratios_all):.3g} to {max(ratios_all):.3g}; extend continental_baseline.LAND_FRACTIONS to bracket it.")
    return crossings


def plot_weathering_ratio_map(series, output_path, levels=13):
    """Contour map of seafloor / continental alkalinity flux over instellation x land fraction.

    The quantity is a POLARITY -- which of the two sinks is winning -- so it is contoured as
    log10(seafloor / continental) on a diverging scale with a neutral midpoint pinned to a ratio
    of 1, and the ratio = 1 contour is drawn as a solid line. That line is the answer to "where
    does seafloor weathering take over": everything above it (toward less land) is
    seafloor-dominated, everything below is continental-dominated.

    Land fraction 0 is NOT on the map. Continental weathering there is exactly zero, so the ratio
    is infinite rather than large -- it is the limit the map runs toward, not a row in it.

    Only steady states (converged, or integrated to 2 Gyr) are contoured. A run stopped at a
    domain wall was still evolving when the model gave up, so its two fluxes are not a balance of
    anything; those cells are left blank and marked, rather than interpolated through silently.
    """
    lands = sorted(l for l in series if l > 0)
    if len(lands) < 2:
        print("Fewer than two positive land fractions -- skipping the ratio map.")
        return None

    ratio, dropped_pts = {}, []
    for land in lands:
        group = series[land]
        steady = group[group['termination'].isin(HABITABLE)]
        for _, r in group.iterrows():
            if r['name'] not in set(steady['name']):
                dropped_pts.append((float(r['instellation']), land))
        if steady.empty:
            continue
        cont, sea = _alk_fluxes(steady)
        for s, c, f in zip(steady['instellation'], cont, sea):
            if np.isfinite(c) and np.isfinite(f) and c > 0 and f > 0:
                ratio[(float(s), land)] = f / c

    if not ratio:
        print("No steady-state runs with both fluxes positive -- skipping the ratio map.")
        return None

    s_vals = sorted({s for s, _ in ratio})
    Z = np.full((len(lands), len(s_vals)), np.nan)
    for i, land in enumerate(lands):
        for j, s in enumerate(s_vals):
            v = ratio.get((s, land))
            if v is not None:
                Z[i, j] = np.log10(v)
    Zm = np.ma.masked_invalid(Z)

    # Diverging about ratio = 1, but NOT forced symmetric. The ratio runs from ~1e-3.5 to only
    # ~1e0.5, so a symmetric range would reserve half the ramp for values that do not occur and
    # squeeze every real contrast into one end. TwoSlopeNorm scales the two sides independently,
    # which keeps the neutral midpoint pinned to 1 -- the only value that means anything here --
    # while both halves still use their full colour range.
    zmin, zmax = float(np.nanmin(Z)), float(np.nanmax(Z))
    step = 0.5                                   # half-decade bands, so 0 is always a boundary
    lo = np.floor(zmin / step) * step
    hi = np.ceil(zmax / step) * step
    bands = np.arange(lo, hi + 0.5 * step, step)
    norm = mcolors.TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=max(hi, step))
    cmap = WEATHERING_RATIO_CMAP   # diverging, neutral (white) midpoint at ratio = 1

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.9))
    cf = ax.contourf(s_vals, lands, Zm, levels=bands, cmap=cmap, norm=norm, extend='both')
    # The crossover itself, drawn on top of the fill.
    if np.nanmin(Z) < 0 < np.nanmax(Z):
        cs = ax.contour(s_vals, lands, Zm, levels=[0.0], colors='k', linewidths=1.4)
        ax.clabel(cs, fmt={0.0: 'equal'}, fontsize=7, inline=True)

    for s, land in dropped_pts:
        ax.plot(s, land, marker='x', color='0.45', markersize=3.5, mew=0.9, zorder=4)

    ax.set_yscale('log')
    ax.set_xlabel('Instellation (S/S₀)')
    ax.set_ylabel('Land fraction')
    # A small margin on both axes so the markers on the edge rows and columns are not sliced in
    # half by the frame; the y margin is taken in log space, where that axis lives.
    dx = 0.02 * (max(s_vals) - min(s_vals))
    dy = 0.04 * (np.log10(max(lands)) - np.log10(min(lands)))
    ax.set_xlim(min(s_vals) - dx, max(s_vals) + dx)
    ax.set_ylim(10 ** (np.log10(min(lands)) - dy), 10 ** (np.log10(max(lands)) + dy))

    # Ticks as ratios, not decades of a log ratio -- the reader wants "10x", not "1 dex".
    ticks = [t for t in range(-9, 10) if lo <= t <= hi]
    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=22, ticks=ticks)
    cbar.set_label('Seafloor / continental alkalinity flux')
    cbar.set_ticklabels([('1' if t == 0 else f'$10^{{{t}}}$') for t in ticks])
    if np.nanmin(Z) < 0 < np.nanmax(Z):
        cbar.ax.axhline(0, color='k', linewidth=1.2)

    _save_fig(fig, figure_path(output_path, 'weathering_ratio_map.png'))

    print("\nSeafloor / continental alkalinity flux (steady states only):")
    print(f"  {'land':>8} " + ' '.join(f'{s:>8.2f}' for s in s_vals))
    for i, land in enumerate(reversed(lands)):
        row = Z[len(lands) - 1 - i]
        cells = [(f"{10 ** v:8.2g}" if np.isfinite(v) else f"{'--':>8}") for v in row]
        print(f"  {land:8g} " + ' '.join(cells))
    if dropped_pts:
        print(f"  ({len(dropped_pts)} cell(s) blank: no steady state)")
    return ratio


def _grid_slice(df, outgassing, crust, mg_si):
    """{land_fraction: rows} for one (outgassing, crust production, Mg/Si) cell of the grid.

    Restricted to the COARSE axes even where finer runs exist. The Earth-reference cell
    (out 1x, crust 1x, Mg/Si 1.25) is also where the land-fraction series ran, so without this it
    would carry ~3x the instellation samples and two extra land fractions. Contour interpolation
    depends on sampling density, so that one panel would be smoother and reach further in
    instellation than its neighbours -- a difference in the sampling, read as a difference in the
    physics. The fine runs keep their own figure (plot_weathering_ratio_map).
    """
    sub = df[
        df['instellation'].isin(cb.COARSE_INSTELLATION) &
        df['land_fraction'].apply(
            lambda v: any(np.isclose(v, l) for l in cb.COARSE_LAND_FRACTIONS)) &
        _ref_redox(df) &
        _ref_chem(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == cb.OCEAN_DEPTH) &
        (df['outgassing'] == outgassing) &
        (df['crust_production'] == crust) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], cb.DELTA_IW)
    ]
    return {float(l): sub[np.isclose(sub['land_fraction'], l)].sort_values('instellation')
            for l in sorted(sub['land_fraction'].unique())}


def _ratio_cells(series, output_path):
    """Ratio cells for one panel: ``(ratio, not_steady, net_sink)``.

    A cell can be missing for two quite different reasons, and they are returned separately so a
    figure can mark them differently rather than leaving identical blanks:

    ``not_steady``  the run never reached a steady state (left the model domain, or hit the
                    wall-clock cap), so its fluxes are a transient, not a balance.
    ``net_sink``    the run IS a steady state but its seafloor alkalinity flux is NEGATIVE -- the
                    pore space precipitates more than the basalt dissolves, so the seafloor is a
                    net alkalinity sink. That is a real outcome, not a failure; it just has no
                    place on a log ratio. It shows up where continental weathering is enormous
                    (~150 Tmol/yr at 10x outgassing with land), which floods the ocean with
                    cations until pore precipitation overwhelms dissolution.
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
        cont, sea = _alk_fluxes(steady)
        for s, c, f in zip(steady['instellation'], cont, sea):
            if np.isfinite(c) and np.isfinite(f) and c > 0 and f > 0:
                ratio[(float(s), land)] = f / c
            elif np.isfinite(f) and f <= 0:
                net_sink.append((float(s), land))
    return ratio, not_steady, net_sink


def plot_weathering_ratio_grid(df, output_path, step=0.5):
    """The ratio map faceted over the coarse grid: crust production x outgassing, per Mg/Si.

    One figure per mantle Mg/Si, so the two compositions are compared panel-for-panel rather than
    by colour. Every panel shares ONE colour scale, computed across BOTH figures -- otherwise each
    panel would renormalise to its own range and the question the grid exists to answer (does the
    crossover move with tectonics or crust chemistry?) would be invisible, because every panel
    would look alike whatever its numbers were.

    Same conventions as the single map: diverging about a ratio of 1, land fraction 0 excluded
    (continental weathering is exactly zero there, so the ratio is infinite), steady states only,
    and cells without one left blank and marked.

    READ THE Mg/Si AND CRUST-PRODUCTION AXES WITH CARE. Both feed the SEAFLOOR side only:
    `mantle_mg_si` reaches the model through `Planet.crust_composition`, and
    `crust_production_rate` through `J_total`, and neither is an argument to
    `get_continental_weathering_flux`, which sees only T and pCO2 against a `F_alk_ref` pinned to
    modern Earth and a cation split fixed to modern river chemistry. So continental weathering
    cannot respond to crust chemistry or tectonic rate except through the shared climate.

    That is not a small correction. Measured at land 0.003, S = 0.8, going Mg/Si 1.25 -> 1.8:
    the seafloor flux rises only 1.3-1.8x while the continental flux FALLS to 0.48-0.68x, because
    the stronger seafloor sink draws pCO2 down (1.45 -> 0.76 bar) and cools the planet
    (330.5 -> 321.6 K), weakening WHAK. In 4 of 5 cells the climate-mediated continental change
    is the larger of the two. The Mg/Si signal here is therefore mostly an indirect response of a
    composition-blind continental law; a continental crust that tracked mantle Mg/Si would weather
    faster too and cancel part of it, so treat the shift as an UPPER BOUND.

    The crust-production axis is likewise one-sided. On a real planet tectonic vigour also drives
    orogeny and uplift, hence physical erosion and the supply of fresh silicate to continental
    weathering -- the supply-limited regime of West et al. (2005) and Maher & Chamberlain (2014).
    The seafloor law here carries transport/supply limitation (sedimentation, Damkohler number)
    but the continental law is pure kinetic WHAK with no runoff and no supply term, so that
    coupling has no route into the model at all.
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
    lo = np.floor(min(allv) / step) * step
    hi = np.ceil(max(allv) / step) * step
    bands = np.arange(lo, hi + 0.5 * step, step)
    norm = mcolors.TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=max(hi, step))
    cmap = WEATHERING_RATIO_CMAP
    ticks = [t for t in range(-9, 10) if lo <= t <= hi]

    for mg in mg_vals:
        fig, axes = plt.subplots(len(cb.GRID_CRUST), len(cb.GRID_OUTGASSING), sharex=True,
                                    sharey=True, squeeze=False,
                                    figsize=figure_size('double', height=5.0))
        cf = None
        for i, c in enumerate(reversed(cb.GRID_CRUST)):
            for j, o in enumerate(cb.GRID_OUTGASSING):
                ax = axes[i, j]
                got = cells.get((mg, c, o))
                ratio, skipped, sinks = got if got else ({}, [], [])
                s_vals = sorted({s for s, _ in ratio})
                lands = sorted({l for _, l in ratio})
                if len(s_vals) > 1 and len(lands) > 1:
                    Z = np.full((len(lands), len(s_vals)), np.nan)
                    for a, land in enumerate(lands):
                        for b, sv in enumerate(s_vals):
                            v = ratio.get((sv, land))
                            if v is not None:
                                Z[a, b] = np.log10(v)
                    Zm = np.ma.masked_invalid(Z)
                    cf = ax.contourf(s_vals, lands, Zm, levels=bands, cmap=cmap, norm=norm,
                                     extend='both')
                    if np.nanmin(Z) < 0 < np.nanmax(Z):
                        ax.contour(s_vals, lands, Zm, levels=[0.0], colors='k', linewidths=1.2)
                else:
                    ax.text(0.5, 0.5, 'no steady state', transform=ax.transAxes, ha='center',
                            va='center', fontsize=7, color='0.5', style='italic')
                for sv, land in skipped:
                    ax.plot(sv, land, marker='x', color='0.45', markersize=3, mew=0.8)
                for sv, land in sinks:
                    ax.plot(sv, land, marker='o', markerfacecolor='none', markeredgecolor='0.25',
                            markersize=4, mew=0.9)
                ax.set_yscale('log')
                ax.grid(True, linestyle='--', alpha=0.3)
                if i == 0:
                    ax.set_title(f'outgassing {o:g}x', fontsize=8)
                if j == 0:
                    ax.set_ylabel(f'crust {c:g}x\nLand fraction', fontsize=7)
                if i == len(cb.GRID_CRUST) - 1:
                    ax.set_xlabel('Instellation (S/S0)')

        if cf is not None:
            cbar = fig.colorbar(cf, ax=list(axes.ravel()), pad=0.02, aspect=30, ticks=ticks)
            cbar.set_label('Seafloor / continental alkalinity flux')
            cbar.set_ticklabels([('1' if t == 0 else f'$10^{{{t}}}$') for t in ticks])
            cbar.ax.axhline(0, color='k', linewidth=1.2)
        fig.suptitle(f'Mantle Mg/Si = {mg:g}', fontsize=9)
        _save_fig(fig, figure_path(output_path, f'weathering_ratio_grid_mgsi{mg:g}.png'))

    print("\nCrossover land fraction across the grid (steady states only):")
    print(f"  {'Mg/Si':>6} {'crust':>7} {'out':>6}   crossover (by instellation)")
    for (mg, c, o), (ratio, *_) in sorted(cells.items()):
        s_vals = sorted({s for s, _ in ratio})
        pts = []
        for sv in s_vals:
            lands = sorted({l for (s2, l) in ratio if s2 == sv})
            got = _crossover_land_fraction(lands, [ratio[(sv, l)] for l in lands])
            if got is not None:
                pts.append(got)
        span = (f"{min(pts):.2g}-{max(pts):.2g}" if pts else
                ("none in range" if ratio else "no steady state"))
        print(f"  {mg:6g} {c:7g} {o:6g}   {span}")
    return cells


def _alpha_slice(df, outgassing, alpha, crust=None, mg_si=None):
    """{land_fraction: rows} for one (outgassing, alpha) cell.

    Deliberately does NOT use _ref_chem: that helper pins alpha to the most-run value, which
    would silently discard the alpha = 10 and 50 arms and leave a "sweep" of one point. kd_mg and
    k_na are still pinned, explicitly, to the calibrated values.
    """
    crust = cb.CRUST_PRODUCTION if crust is None else crust
    mg_si = cb.MG_SI_EARTH if mg_si is None else mg_si
    sub = df[
        df['instellation'].isin(cb.COARSE_INSTELLATION) &
        df['land_fraction'].apply(
            lambda v: any(np.isclose(v, l) for l in cb.COARSE_LAND_FRACTIONS)) &
        _ref_redox(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == cb.OCEAN_DEPTH) &
        (df['outgassing'] == outgassing) &
        (df['crust_production'] == crust) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], cb.DELTA_IW) &
        np.isclose(df['alpha'], alpha) &
        np.isclose(df['kd_mg'], cb.KD_MG_CALIB) &
        np.isclose(df['k_na'], cb.K_NA_CALIB)
    ]
    return {float(l): sub[np.isclose(sub['land_fraction'], l)].sort_values('instellation')
            for l in sorted(sub['land_fraction'].unique())}


def plot_alpha_scaling(df, output_path):
    """Does the crossover land fraction really go as alpha^1?

    In the kinetic limit the seafloor flux is linear in alpha while continental weathering does
    not see alpha at all, so f* should scale as alpha^1. The climate feedback should DAMP that:
    raising alpha strengthens the sink, which cools the planet and draws CO2 down, weakening both
    fluxes again. An exponent below 1 is therefore the expected outcome, and its size is the
    result -- it says how much of alpha's nominal leverage survives the feedback.
    """
    rows = []
    for o in cb.GRID_OUTGASSING:
        for a in cb.GRID_ALPHA:
            series = _alpha_slice(df, o, a)
            if not series:
                continue
            ratio, _, _ = _ratio_cells(series, output_path)
            if not ratio:
                continue
            for sv in sorted({s for s, _ in ratio}):
                lands = sorted({l for (s2, l) in ratio if s2 == sv})
                got = _crossover_land_fraction(lands, [ratio[(sv, l)] for l in lands])
                if got is not None:
                    rows.append({'outgassing': o, 'alpha': a, 'instellation': sv, 'f_star': got})
    if not rows:
        print("No crossovers found across the alpha grid -- skipping.")
        return None
    tab = pd.DataFrame(rows)

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=3.0))
    cmap = OUTGASSING_CMAP
    norm = mcolors.LogNorm(vmin=min(cb.GRID_OUTGASSING), vmax=max(cb.GRID_OUTGASSING))

    print("\nCrossover land fraction f* against alpha (geometric mean over instellation):")
    print(f"  {'outgassing':>10} " + ' '.join(f'{a:>10.4g}' for a in cb.GRID_ALPHA)
          + f" {'exponent':>9}")
    exponents, anchor = {}, None
    for o in cb.GRID_OUTGASSING:
        g = tab[tab.outgassing == o]
        if g.empty:
            continue
        # Geometric mean over instellation: f* spans decades, so an arithmetic mean would be
        # dominated by whichever instellation sits nearest the runaway.
        means = {a: float(np.exp(np.log(g[g.alpha == a].f_star).mean()))
                 for a in cb.GRID_ALPHA if (g.alpha == a).any()}
        cells = [(f'{means[a]:10.4g}' if a in means else f'{"--":>10}') for a in cb.GRID_ALPHA]
        slope = float('nan')
        if len(means) >= 2:
            slope = float(np.polyfit(np.log10(list(means)),
                                     np.log10(list(means.values())), 1)[0])
            exponents[o] = slope
        print(f"  {o:10g} " + ' '.join(cells) + f" {slope:9.2f}")
        ax.plot(list(means), list(means.values()), marker='o', markersize=4,
                color=cmap(norm(o)), linewidth=1.6, label=f'{o:g}x')
        if anchor is None and means:
            anchor = (min(means), means[min(means)])

    if anchor is not None:
        a0, f0 = anchor
        xs = np.array(cb.GRID_ALPHA, dtype=float)
        ax.plot(xs, f0 * xs / a0, color='0.4', linestyle=(0, (6, 3)), linewidth=1.2,
                label=r'$\propto \alpha$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Reactive area scaling $\alpha$')
    ax.set_ylabel(r'Crossover land fraction $f^*$')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=7, frameon=False, title='Outgassing', title_fontsize=7)
    _save_fig(fig, figure_path(output_path, 'alpha_scaling.png'))

    if exponents:
        v = list(exponents.values())
        print(f"  exponent d log f* / d log alpha: {min(v):.2f} to {max(v):.2f} "
              f"(alpha^1 would be 1.00)")
    tab.to_csv(os.path.join(output_path, 'alpha_crossover.csv'), index=False)
    return tab


def plot_alpha_ratio_grid(df, output_path, step=0.5):
    """Ratio map faceted over alpha x outgassing: columns are outgassing, rows are alpha.

    Same conventions as the other ratio maps -- diverging about a ratio of 1 with the neutral
    midpoint pinned there, land fraction 0 excluded (continental weathering is exactly zero, so
    the ratio is infinite), steady states only, cells without one left blank and marked.

    ONE colour scale across all nine panels. Per-panel normalisation would make every panel look
    alike whatever its numbers were, which would hide the whole point: alpha slides the ratio
    bodily up the land-fraction axis (measured d log f*/d log alpha = 0.78-1.16) while outgassing
    decides whether a crossover exists at all.

    To transpose the layout, swap ROWS and COLS below.
    """
    ROWS, COLS = cb.GRID_ALPHA, cb.GRID_OUTGASSING          # rows: alpha, columns: outgassing
    row_label, col_label = r'$\alpha$', 'outgassing'

    cells = {}
    for a in ROWS:
        for o in COLS:
            series = _alpha_slice(df, o, a)
            if series:
                cells[(a, o)] = _ratio_cells(series, output_path)

    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state alpha runs with both fluxes positive -- skipping.")
        return None
    lo = np.floor(min(allv) / step) * step
    hi = np.ceil(max(allv) / step) * step
    bands = np.arange(lo, hi + 0.5 * step, step)
    norm = mcolors.TwoSlopeNorm(vmin=min(lo, -step), vcenter=0.0, vmax=max(hi, step))
    cmap = WEATHERING_RATIO_CMAP

    fig, axes = plt.subplots(len(ROWS), len(COLS), sharex=True, sharey=True, squeeze=False,
                                figsize=figure_size('double', height=5.0))
    cf = None
    for i, rv in enumerate(reversed(ROWS)):          # largest alpha at the top
        for j, cv in enumerate(COLS):
            ax = axes[i, j]
            ratio, skipped, sinks = cells.get((rv, cv), ({}, [], []))
            s_vals = sorted({s for s, _ in ratio})
            lands = sorted({l for _, l in ratio})
            if len(s_vals) > 1 and len(lands) > 1:
                Z = np.full((len(lands), len(s_vals)), np.nan)
                for a_, land in enumerate(lands):
                    for b_, sv in enumerate(s_vals):
                        v = ratio.get((sv, land))
                        if v is not None:
                            Z[a_, b_] = np.log10(v)
                Zm = np.ma.masked_invalid(Z)
                cf = ax.contourf(s_vals, lands, Zm, levels=bands, cmap=cmap, norm=norm,
                                 extend='both')
                if np.nanmin(Z) < 0 < np.nanmax(Z):
                    ax.contour(s_vals, lands, Zm, levels=[0.0], colors='k', linewidths=1.2)
            else:
                ax.text(0.5, 0.5, 'no steady state', transform=ax.transAxes, ha='center',
                        va='center', fontsize=7, color='0.5', style='italic')
            for sv, land in skipped:
                ax.plot(sv, land, marker='x', color='0.45', markersize=3, mew=0.8)
            for sv, land in sinks:
                ax.plot(sv, land, marker='o', markerfacecolor='none', markeredgecolor='0.25',
                        markersize=4, mew=0.9)
            ax.set_yscale('log')
            ax.grid(True, linestyle='--', alpha=0.3)
            if i == 0:
                ax.set_title(f'{col_label} {cv:g}x', fontsize=8)
            if j == 0:
                ax.set_ylabel(f'{row_label} = {rv:g}\nLand fraction', fontsize=7)
            if i == len(ROWS) - 1:
                ax.set_xlabel('Instellation (S/S0)')

    if cf is not None:
        ticks = [t for t in range(-9, 10) if lo <= t <= hi]
        cbar = fig.colorbar(cf, ax=list(axes.ravel()), pad=0.02, aspect=30, ticks=ticks)
        cbar.set_label('Seafloor / continental alkalinity flux')
        cbar.set_ticklabels([('1' if t == 0 else f'$10^{{{t}}}$') for t in ticks])
        cbar.ax.axhline(0, color='k', linewidth=1.2)
    _save_fig(fig, figure_path(output_path, 'weathering_ratio_alpha_grid.png'))
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
    edges = plot_habitable_zone(arms, output_path)
    _report(arms, edges or {})

    # The land-fraction series only draws once intermediate land fractions are on disk.
    series = _land_series(df, output_path)
    if len(series) > len(cb.LAND_ARMS):
        print(f"\n  land fractions on disk: {sorted(series, reverse=True)}")
    plot_land_fraction_series(series, output_path)
    plot_weathering_crossover(series, output_path)
    plot_weathering_ratio_map(series, output_path)
    plot_weathering_ratio_grid(df, output_path)
    plot_alpha_scaling(df, output_path)
    plot_alpha_ratio_grid(df, output_path)
    plot_continental_baseline(df, output_path)


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
    cmap = OUTGASSING_CMAP

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
        _save_fig(fig, figure_path(output_path, f'mineral_si_crust_{c}.png'))

# ---------------------------------------------------------------------------
# Reactive area (alpha) against outgassing
# ---------------------------------------------------------------------------
# Marker per alpha value for plot_alpha_outgassing. The point of the figure is that the three
# alphas land on ONE curve, so they must be distinguishable by shape while sharing the
# instellation colour scale -- colouring by alpha instead would make the collapse unreadable.
ALPHA_MARKERS = ('o', 's', 'D', '^', 'v')


def _alpha_collapse_exponent(sub, exponents=np.arange(-1.0, 2.01, 0.05)):
    """The p that best collapses T onto a single curve in outgassing / alpha**p.

    Measured, not assumed: at each instellation the runs are fitted with a quadratic in
    log10(outgassing / alpha**p) and the residuals pooled across instellations. p = 1 is exact
    trade-off (one decade of alpha cancels one decade of outgassing); p = -1 is the product.
    Returns (best_p, rms_at_best, rms_by_exponent) or (nan, nan, {}) if nothing is comparable.
    """
    groups = [g for _, g in sub.groupby('instellation')
              if g['alpha'].nunique() > 1 and len(g) >= 4]
    if not groups:
        return np.nan, np.nan, {}
    rms = {}
    for p in exponents:
        total, n = 0.0, 0
        for g in groups:
            x = np.log10(g['outgassing'].to_numpy() / g['alpha'].to_numpy() ** p)
            y = g['T'].to_numpy()
            resid = y - np.polyval(np.polyfit(x, y, 2), x)
            total += float(np.sum(resid ** 2))
            n += len(resid)
        rms[round(float(p), 3)] = np.sqrt(total / n)
    best = min(rms, key=lambda k: rms[k])
    return best, rms[best], rms


def plot_alpha_outgassing(df, output_path, mg_si=REF_MG_SI, crust_production=1.0,
                          alpha_exponent=1.0, ocean_depth=3000,
                          s_vals=(0.4, 0.6, 0.8, 1.0, 1.2)):
    """Temperature against outgassing and reactive area combined into one variable.

    Raising alpha raises the seafloor weathering sink; raising outgassing raises the CO2 source.
    If the two are exactly interchangeable, T depends only on their RATIO, and the runs at every
    alpha fall on one curve at fixed instellation. `alpha_exponent` p sets the combination
    plotted, outgassing / alpha**p: p = 1 is exact trade-off, p = -1 is the product, p = 0 drops
    alpha entirely. The measured best-collapsing p is printed for comparison.

    Only Da < 1 runs are drawn: past Da = 1 the sink is thermodynamically limited and stops
    responding to alpha, so those runs carry no information about the trade-off and only add
    scatter. Held fixed at one mantle Mg/Si and one crust production rate, since both change the
    sink for reasons that have nothing to do with alpha. Coloured by instellation, with a marker per
    alpha: a collapse is the three marker shapes interleaving along a single colour band.
    `s_vals` thins the instellation grid to a few well-separated bands (as plot_ratio_scatter
    does) -- the full 19-value grid renders as a colour continuum in which no band is legible.
    The exponent is measured on every run, not only the ones drawn.
    """
    sub = df[
        df['reverse_weathering'] &
        _ref_redox(df) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], REF_DIW) &
        (df['ocean_depth'] == ocean_depth) &
        (df['land_fraction'] == 0.0) &
        (df['outgassing'] > 0) &
        np.isclose(df['crust_production'], crust_production)
    ].copy()
    # T clamped to a validity wall is a sentinel, not a temperature the climate solver found,
    # and a band of runs pinned to exactly 389 K would read as a real plateau in the collapse.
    sub = sub[np.isfinite(sub['T']) & (sub['T'] < T_HOT_WALL) & (sub['T'] > T_COLD_WALL)]
    if sub.empty:
        print("No in-domain runs for the alpha figure — skipping.")
        return
    # Kinetic runs only. Past Da = 1 the sink is thermodynamically limited, so it no longer
    # responds to alpha at all -- those runs pile up on whatever temperature the saturated sink
    # allows and scatter across the whole x range, which is what made the figure unreadable.
    # The trade-off alpha vs outgassing is a statement about the KINETIC regime.
    sub = _add_diag_columns(sub, output_path)
    sub = sub[np.isfinite(sub['da']) & (sub['da'] < 1.0)]
    if sub['alpha'].nunique() < 2:
        print("Fewer than two alpha values with Da < 1 — skipping alpha figure.")
        return

    best_p, best_rms, rms = _alpha_collapse_exponent(sub)
    shown_rms = rms.get(round(float(alpha_exponent), 3), np.nan)
    print(f"Alpha vs outgassing [Mg/Si={mg_si:g}, crust={crust_production:g}]: {len(sub)} runs, "
          f"alpha = {', '.join(f'{a:g}' for a in sorted(sub['alpha'].unique()))}")
    print(f"  best collapse at p = {best_p:+.2f} (rms {best_rms:.1f} K); "
          f"plotted p = {alpha_exponent:+.2f} (rms {shown_rms:.1f} K)")

    if s_vals is not None:
        keep = [min(sorted(sub['instellation'].unique()), key=lambda v: abs(v - t))
                for t in s_vals]
        sub = sub[sub['instellation'].isin(keep)]
    sub['combined'] = sub['outgassing'] / sub['alpha'] ** alpha_exponent

    s_vals = sorted(sub['instellation'].unique())
    norm = mcolors.Normalize(vmin=min(s_vals), vmax=max(s_vals))
    cmap = INSTELLATION_CMAP

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=2.6))

    alphas = sorted(sub['alpha'].unique())
    for a, marker in zip(alphas, ALPHA_MARKERS):
        grp = sub[np.isclose(sub['alpha'], a)]
        if grp.empty:
            continue
        ax.scatter(grp['combined'], grp['T'], marker=marker, c=grp['instellation'],
                   cmap=cmap, norm=norm, s=20, alpha=0.9, zorder=4,
                   linewidths=0.4, edgecolors='0.25')

    ax.axhspan(T_SNOWBALL - 25, T_SNOWBALL, color='blue', alpha=0.12, zorder=0)
    ax.axhspan(T_RUNAWAY - 20,  T_RUNAWAY,  color='red',  alpha=0.12, zorder=0)
    ax.set_xscale('log')
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel(_alpha_combination_label(alpha_exponent))
    ax.grid(True, linestyle='--', alpha=0.4)

    _add_colorbar(fig, ax, cmap, norm, 'Instellation (S/S₀)',
                  ticks=_colorbar_ticks(s_vals)[0])

    handles = [Line2D([0], [0], marker=m, color='0.25', linestyle='none', markersize=4,
                      markerfacecolor='0.6', label=rf'$\alpha$ = {a:g}')
               for a, m in zip(alphas, ALPHA_MARKERS)]
    _add_figure_legend(fig, [ax], handles)

    tag = '' if np.isclose(alpha_exponent, 1.0) else f'_p{alpha_exponent:g}'
    _save_fig(fig, figure_path(output_path, f'alpha_outgassing{tag}.png'))


def _alpha_combination_label(p):
    """Axis label for outgassing / alpha**p, written the way that p is normally read."""
    if np.isclose(p, 1.0):
        return r'Outgassing / $\alpha$ (×Earth)'
    if np.isclose(p, -1.0):
        return r'Outgassing $\times\ \alpha$ (×Earth)'
    if np.isclose(p, 0.0):
        return r'Outgassing (×Earth)'
    return rf'Outgassing / $\alpha^{{{p:g}}}$ (×Earth)'


def _temperature_plane(sub, xcol, ycol, xlabel, ylabel, title, path, what):
    """Filled surface-temperature contours over (log10 xcol, log10 ycol), saved to `path`.

    Only steady states are contoured; runs that never reached one are marked, since their stored
    T is often a domain-wall sentinel rather than a computed temperature.
    """
    # Values from other sweeps' grids sit beside this plane's own; keep the ones shared across it.
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
    # Small margin so markers on the edge rows and columns are not cut by the frame.
    dx, dy = 0.03 * (x[-1] - x[0]), 0.03 * (y[-1] - y[0])
    ax.set_xlim(x[0] - dx, x[-1] + dx)
    ax.set_ylim(y[0] - dy, y[-1] + dy)
    cbar = fig.colorbar(cf, ax=ax, pad=0.02, aspect=22)
    cbar.set_label('Surface temperature (K)')
    if not unsteady.empty:
        _add_figure_legend(fig, [ax], [Line2D([0], [0], linestyle='none', marker='x', color='0.45',
                                              markersize=3.5, mew=0.9, label='No steady state')])
    print(f"{what}: {len(xs)} {xcol} x {len(ys)} {ycol}, {len(steady)} steady, "
          f"{len(unsteady)} not, T {t_lo:.1f}-{t_hi:.1f} K")
    _save_fig(fig, path)


def plot_alpha_outgassing_plane(df, output_path, instellation=0.8, mg_si=REF_MG_SI,
                                crust_production=1.0, ocean_depth=3000):
    """Surface temperature over the (log alpha, log outgassing) plane at one instellation."""
    sub = df[
        df['reverse_weathering'] &
        _ref_redox(df) &
        np.isclose(df['instellation'], instellation) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], REF_DIW) &
        (df['ocean_depth'] == ocean_depth) &
        (df['land_fraction'] == 0.0) &
        (df['outgassing'] > 0) &
        np.isclose(df['crust_production'], crust_production) &
        np.isclose(df['kd_mg'], _chem_reference(df, 'kd_mg')) &
        np.isclose(df['k_na'], _chem_reference(df, 'k_na'))
    ]
    _temperature_plane(sub, 'alpha', 'outgassing',
                       r'$\log_{10}\,\alpha$ (reactive area scaling)',
                       r'$\log_{10}$ outgassing (×Earth)', f'S = {instellation:g}',
                       figure_path(output_path, 'alpha_outgassing_plane.png'),
                       f'Alpha-outgassing plane at S = {instellation:g}')


def plot_outgassing_crust_plane(df, output_path, instellation=0.8):
    """Surface temperature over the (log outgassing, log crust production) plane at one instellation.

    Lines of constant outgassing / crust production ratio are diagonals of slope 1 here, so this
    shows directly whether temperature follows the ratio or the two rates separately.
    """
    sub = _base(df)
    sub = sub[np.isclose(sub['instellation'], instellation)]
    _temperature_plane(sub, 'outgassing', 'crust_production',
                       r'$\log_{10}$ outgassing (×Earth)',
                       r'$\log_{10}$ crust production (×Earth)', f'S = {instellation:g}',
                       figure_path(output_path, 'outgassing_crust_plane.png'),
                       f'Outgassing-crust plane at S = {instellation:g}')


# ---------------------------------------------------------------------------
# The kinetic -> thermodynamic transition in phase space
# ---------------------------------------------------------------------------
# Axis labels for the two tectonic axes, so neither arrangement of plot_da_transition
# hardcodes them. The colours come from PARAM_CMAPS, like every other faceted figure.
DA_TRANSITION_LABELS = {
    'crust_production': 'Crust production (×Earth)',
    'outgassing':       'Outgassing (×Earth)',
}

# Mg/Si planes shown as the three panels, left to right. These are the values
# parameter_sweep.sweep_basic_low_mgsi / sweep_basic / sweep_basic_high_mgsi cover; a value with
# no basic sweep behind it would draw an empty panel, so the caller's list is intersected with
# basic_plane_mg_si(df).
DA_TRANSITION_MG_SI = (0.8, 1.25, 1.8)


def _da_transition_instellation(group):
    """Instellation at which one (Mg/Si, outgassing, crust) line crosses Da = 1.

    The crossing is interpolated linearly in log10(Da), which is the quantity that is actually
    smooth in instellation -- Da spans ~6 decades along a single line, so interpolating Da
    itself would put the root almost on top of the first supra-unity point.

    Returns (S_crit, clamped): `clamped` flags a crossing whose two bracketing runs include a
    state pinned to a temperature wall. Da there was evaluated on a sentinel state rather than
    on a temperature the climate solver actually found, so the crossing is drawn hollow.
    Returns (nan, False) when the line never crosses -- entirely kinetic or entirely
    thermodynamic over the swept instellation range.
    """
    g = group.sort_values('instellation')
    da = g['da'].to_numpy(dtype=float)
    ok = np.isfinite(da) & (da > 0)
    if ok.sum() < 2:
        return np.nan, False
    s   = g['instellation'].to_numpy(dtype=float)[ok]
    T   = g['T'].to_numpy(dtype=float)[ok]
    ld  = np.log10(da[ok])

    # Kinetic -> thermodynamic only; a frozen, CO2-starved state can sit at Da > 1 at low S.
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
    """Where the ocean world stops being kinetically limited, in (instellation x outgassing).

    THE summary figure for the transition. Each instellation sweep -- one fixed (Mg/Si,
    outgassing, crust production) -- runs from a kinetically limited seafloor at low
    instellation to a thermodynamically limited one at high instellation; the Da = 1 crossing is
    reduced to a single instellation by _da_transition_instellation. Joining those crossings
    across outgassing gives one line per crust-production rate: the boundary of the kinetic
    regime. Everything LEFT of a line is kinetic (the weathering sink can still respond to
    warming, so the carbon cycle stabilises); everything RIGHT is thermodynamic (the sink is
    saturated, cannot respond, and the planet runs away).

    One panel per mantle Mg/Si. `line_by` picks which of the two tectonic axes carries the
    lines and the colourbar; the other one becomes the y-axis. Both arrangements show the same
    147 crossings -- which axis reads as the stronger control is exactly what they compare.
    """
    if line_by not in DA_TRANSITION_LABELS:
        raise ValueError(f"line_by must be one of {sorted(DA_TRANSITION_LABELS)}, not {line_by!r}")
    y_by = 'outgassing' if line_by == 'crust_production' else 'crust_production'
    if cmap is None:
        cmap = PARAM_CMAPS[line_by]
    available = basic_plane_mg_si(df)
    mg_vals = [m for m in mg_si_values
               if any(np.isclose(m, a) for a in available)]
    if not mg_vals:
        print("No basic Mg/Si plane with a full sweep — skipping Da transition figure.")
        return

    frames = []
    for mg in mg_vals:
        sel = _base(df, mg_si=(None if np.isclose(mg, REF_MG_SI) else mg))
        if not sel.empty:
            frames.append(_add_diag_columns(sel, output_path))
    subset = pd.concat(frames)
    # Main-sweep outgassing only: other sweeps' grids (3.2, 5.6, ...) exist at crust 1x alone.
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

    # Not drawn any more, but a crossing bracketed by a temperature-wall sentinel is a weaker
    # measurement than the rest -- report the count rather than losing it silently.
    n_clamped = int(trans['clamped'].sum())
    if n_clamped:
        print(f"  Da transition: {n_clamped} of {len(trans)} crossings bracketed by a T-wall state.")

    hz_on = SHOW_HZ_EDGES if show_hz is None else show_hz
    line_vals = sorted(trans[line_by].unique())
    norm   = _value_norm(line_vals, pad=0.0)
    colours = {c: cmap(norm(c)) for c in line_vals}
    ticks, ticklabels = _colorbar_ticks(line_vals)

    n_cols = len(mg_vals)
    fig, axes = plt.subplots(1, n_cols, sharex=True, sharey=True, squeeze=False,
                             figsize=figure_size('double', height=2.6))
    axes = axes[0]

    for ax, mg in zip(axes, mg_vals):
        panel = trans[np.isclose(trans['mg_si'], mg)]
        for c in line_vals:
            line = panel[np.isclose(panel[line_by], c)].sort_values(y_by)
            if line.empty:
                continue
            ax.plot(line['s_crit'], line[y_by], color=colours[c],
                    linestyle='-', linewidth=1.4, zorder=3)
        ax.set_yscale('log')
        ax.set_title(f'Mg/Si = {mg:g}')
        ax.grid(True, linestyle='--', alpha=0.3)
        _draw_hz_edges(ax, show_hz)

    axes[0].set_ylabel(DA_TRANSITION_LABELS[y_by])
    axes[len(axes) // 2].set_xlabel('Instellation (S/S$_0$)')
    # Limits from the crossings themselves, not from the swept instellation range: the sweep
    # runs out to 1.45 but no line transitions past ~1.15, which left a third of the axis empty.
    # The HZ edges are included when drawn, since the outer one sits well left of every crossing
    # and would otherwise fall outside the axis.
    marks = list(trans['s_crit'])
    if hz_on:
        marks += [CONTINENTAL_HZ_OUTER, CONTINENTAL_HZ_INNER]
    lo, hi = min(marks), max(marks)
    pad = 0.08 * (hi - lo) if hi > lo else 0.1
    # Wide margins so the lines are not crowded; each side's margin holds its regime label.
    x_kinetic, x_thermo = lo - 2 * pad, hi + 2 * pad
    axes[0].set_xlim(lo - 4 * pad, hi + 4 * pad)

    label_kw = dict(rotation=90, ha='center', va='center', fontsize=7, color='0.15', zorder=5)
    kinetic_box = dict(boxstyle='round,pad=0.35', facecolor='#d4eddb', edgecolor='#3f8f55',
                       linewidth=0.6)
    thermo_box = dict(boxstyle='round,pad=0.35', facecolor='#f8d6d1', edgecolor='#b8483a',
                      linewidth=0.6)
    for ax in axes:
        trans_ax = ax.get_xaxis_transform()
        # Below mid-height so the boxes clear the HZ edge labels hanging from the top.
        ax.text(x_kinetic, 0.35, 'Kinetic\n(stable) regime', transform=trans_ax,
                bbox=kinetic_box, **label_kw)
        ax.text(x_thermo, 0.35, 'Thermodynamic\n(unstable) regime', transform=trans_ax,
                bbox=thermo_box, **label_kw)

    _add_colorbar(fig, list(axes), cmap, norm, DA_TRANSITION_LABELS[line_by],
                  ticks=ticks, ticklabels=ticklabels)

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
    """Region map of climate state over outgassing / crust production ratio and instellation.

    Each (ratio, instellation) cell takes the majority state of its runs, drawn as a filled block
    with a per-state hatch, outlined where the state changes and labelled inside each region.
    """
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Rectangle

    base = _base(df).copy()
    if base.empty:
        print("No runs found — skipping phase space plot.")
        return
    base['ratio'] = base['outgassing'] / base['crust_production']

    # Final state decides snowball/hothouse; only a run stopped at a CO2 wall has an unknown fate.
    wall = base.get('domain_wall')
    if wall is None:
        wall = pd.Series([None] * len(base), index=base.index)
    cond_unknown = base['termination'].isin(OUT_OF_DOMAIN) & wall.isin(['co2_high', 'co2_low'])
    cond_snow = ((base['T'] <= T_SNOWBALL) | (wall == 'cold')) & ~cond_unknown
    cond_hot = ((base['T'] >= T_RUNAWAY) | (wall == 'hot')) & ~cond_unknown & ~cond_snow
    states = ['Snowball', 'Habitable', 'Hothouse', 'Unknown']
    base['state'] = np.select([cond_unknown, cond_snow, cond_hot], [3, 0, 2], 1)

    # Merge near-duplicate ratios from different sweep grids (e.g. 0.03, 0.032, 0.0333) into one column.
    logr = np.log10(np.sort(base['ratio'].unique()))
    groups = np.concatenate([[0], np.cumsum(np.diff(logr) > 0.05)])
    centres = np.array([logr[groups == g].mean() for g in range(groups[-1] + 1)])
    base['col'] = groups[np.searchsorted(logr, np.log10(base['ratio']).clip(logr[0], logr[-1]))]

    s_vals = np.sort(base['instellation'].unique())
    # Keep columns covering at least half the instellation range, or thin columns leave holes.
    cover = base.groupby('col')['instellation'].nunique()
    cols = [c for c in range(len(centres)) if cover.get(c, 0) >= 0.5 * len(s_vals)]
    base = base[base['col'].isin(cols)]
    x_c = centres[cols]

    grid = np.full((len(s_vals), len(cols)), -1)
    n_mixed = 0
    for (c, s), g in base.groupby(['col', 'instellation']):
        counts = g['state'].value_counts()
        grid[np.searchsorted(s_vals, s), cols.index(c)] = counts.idxmax()
        n_mixed += len(counts) > 1

    def _edges(c):
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])
    xe, ye = 10 ** _edges(x_c), _edges(s_vals)

    style = {0: ('#9cc5ea', '#3a78b5', '\\\\\\\\\\\\'),
             1: ('#a8dbb4', '#3f8f55', None),
             2: ('#f2a99f', '#b8483a', '//////'),
             3: ('#d9d9d9', '#8a8a8a', '...')}

    fig, ax = plt.subplots(1, 1, figsize=figure_size('single', height=3.0))
    for k, (face, edge, hatch) in style.items():
        cells = [Rectangle((xe[i], ye[j]), xe[i + 1] - xe[i], ye[j + 1] - ye[j])
                 for j, i in zip(*np.nonzero(grid == k))]
        if not cells:
            continue
        ax.add_collection(PatchCollection(cells, facecolor=face, edgecolor=edge, hatch=hatch,
                                          linewidth=0, zorder=2))

    # Outline every cell edge where the state changes, including against empty cells.
    padded = np.pad(grid, 1, constant_values=-1)
    segs = []
    for j in range(padded.shape[0]):
        for i in range(padded.shape[1] - 1):
            if padded[j, i] != padded[j, i + 1] and 0 < j < padded.shape[0] - 1:
                segs.append([(xe[i], ye[j - 1]), (xe[i], ye[j])])
    for j in range(padded.shape[0] - 1):
        for i in range(padded.shape[1]):
            if padded[j, i] != padded[j + 1, i] and 0 < i < padded.shape[1] - 1:
                segs.append([(xe[i - 1], ye[j]), (xe[i], ye[j])])
    ax.add_collection(LineCollection(segs, colors='0.15', linewidths=0.9, zorder=3))

    ax.set_xscale('log')
    ax.set_xlim(xe[0], xe[-1])
    ax.set_ylim(ye[0], min(ye[-1], s_max))
    ax.set_xlabel('Outgassing / Crust production rate')
    ax.set_ylabel('Instellation (S/S₀)')
    # Label each state at the centre of its largest visible block of cells, so the text sits inside it.
    names = ['Snowball', 'Habitable', 'Hothouse', 'Unknown\n(CO₂ wall)']
    visible = ye[:-1] < s_max
    for k, name in enumerate(names):
        rect = _largest_rectangle((grid == k) & visible[:, None])
        if rect is None:
            continue
        j0, j1, i0, i1 = rect
        x_mid = 10 ** (0.5 * (np.log10(xe[i0]) + np.log10(xe[i1 + 1])))
        y_mid = 0.5 * (ye[j0] + min(ye[j1 + 1], s_max))
        ax.text(x_mid, y_mid, name, ha='center', va='center', fontsize=7, zorder=5,
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
    parser.add_argument('--pe', type=float, default=None,
                        help='Pin ocean pe to this value in the main plots '
                             '(default: the model reference, planet.PE_DEFAULT).')
    args = parser.parse_args()

    for _knob in CHEM_KNOBS:
        _v = getattr(args, _knob)
        if _v is not None:
            CHEM_OVERRIDE[_knob] = _v
    if args.pe is not None:
        REF_PE = args.pe

    df = load_data(args.path)

    if args.legacy:
        import plot_legacy
        df = plot_legacy.upgrade(df)
        plot_legacy.plot_named_compositions(df, args.path, split_panels=True)

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

    # Depths that carry a resolved redox (pe) sweep -- more than the bracketing reducing/
    # oxidising pair every other sweep runs (parameter_sweep.sweep_pe / sweep_pe_deep).
    redox_depths = sorted(
        d for d in df['ocean_depth'].unique()
        if df[df['ocean_depth'] == d]['pe'].nunique() > 2
    ) or [3000.0]
    print(f"Redox figures for depth(s): {[f'{d:g}' for d in redox_depths]}")

    # Mg/Si values with a full basic sweep. The reference (Earth) plane is drawn untagged, as
    # before; each end-member gets the same pair of figures with the value in the filename.
    basic_mg = basic_plane_mg_si(df)
    end_members = [m for m in basic_mg if not np.isclose(m, REF_MG_SI)]
    if end_members:
        print(f"Basic sweep Mg/Si planes: {[f'{m:g}' for m in basic_mg]}")

    plot_basic(df, args.path, split_panels=True)
    plot_basic(df, args.path, all_results=False, split_panels=True)
    plot_basic_ph(df, args.path)
    for _mg in end_members:
        plot_basic(df, args.path, split_panels=True, mg_si=_mg)
        plot_basic(df, args.path, all_results=False, split_panels=True, mg_si=_mg)
    plot_basic_mgsi_grid(df, args.path, split_panels=True)
    # plot_basic(df, args.path, split_panels=True, all_results=False, sequence=True)
    plot_depth(df, args.path, split_panels=False)
    plot_chemistry(df, args.path, split_panels=True)
    for _d in redox_depths:
        plot_pe(df, args.path, split_panels=True, ocean_depth=_d)
    plot_ratio_scatter(df, args.path)
    for _d in comp_depths:
        plot_cross(df, args.path, show_markers=False, split_panels=True,
                               ocean_depth=_d)
        plot_composition(df, args.path, show_markers=False, split_panels=True,
                              ocean_depth=_d)
    for _q in ('T', 'P_CO2'):
        plot_composition_map(df, args.path, ocean_depth=comp_depths[0], quantity=_q)
    plot_damkohler_contour(df, args.path)
    plot_habitability_phase_space(df, args.path)
    plot_outgassing_crust_plane(df, args.path)
    plot_da_transition(df, args.path)
    plot_da_transition(df, args.path, line_by='outgassing')
    plot_alpha_outgassing(df, args.path)
    plot_alpha_outgassing(df, args.path, alpha_exponent=-1.0)
    plot_alpha_outgassing_plane(df, args.path)
    plot_continental(df, args.path)
    print("Done.")
