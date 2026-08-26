import os
import sys
import glob
import json
import re
import argparse
import functools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
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
from kamino.climate.analytic import get_T_surface_analytic
from kamino.planet import OCEAN_ALBEDO, LAND_ALBEDO, KD_MG_HT, K_NA_CONT_REMOVAL
from kamino.weathering import ALPHA_REF
from kamino.planet import _S_TERR_EARTH

fig_width_half = 3.5
fig_subplot_height = 1.5

presentation = True

if presentation:
    plt.style.use('experiments/planetary-chem-presentation.mplstyle')
else:
    plt.style.use('experiments/planetary-chem-paper.mplstyle')

DEFAULT_OUTPUT_PATH = '/data/pt426/kamino_experiments_fast_18/'

TERM_LABELS = {
    'converged':     'Converged',
    'timeout':       'Timeout (2 Gyr)',
    'out_of_domain': 'Outside model domain',
    # Legacy terminations, kept so pre-domain-event runs in output/ still plot.
    'snowball':      'Snowball (legacy)',
    'hothouse':      'Hothouse (legacy)',
    'co2_ceiling':   'Outside model domain (legacy)',
    'acid_ocean':    'Outside model domain (legacy)',
    'co2_floor':     'CO₂ floor (legacy)',
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
OUT_OF_DOMAIN = {'out_of_domain', 'co2_ceiling', 'acid_ocean'}
HABITABLE = {'converged', 'timeout', 'co2_floor'}
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
COMP_LABELS = {
    'komatiite_42': 'Komatiite (42% SiO₂)',
    'komatiite_44': 'Komatiite (44% SiO₂)',
    'basalt_47':    'Basalt (47% SiO₂)',
    'basalt_49':    'Basalt (49% SiO₂)',
    'basalt_51':    'Basalt (51% SiO₂)',
}

HAB_MARKERS    = {'converged': 'o', 'timeout': 's', 'co2_floor': 'P'}
FAILED_MARKERS = {'out_of_domain': 's', 'snowball': 'v', 'hothouse': '^',
                  'co2_ceiling': 's', 'acid_ocean': 's'}

DA_LEGEND = [
    Line2D([0], [0], color='k', linestyle='-',  linewidth=1.8, label='Da < 1 (kinetic)'),
    Line2D([0], [0], color='k', linestyle='--', linewidth=1.8, label='Da ≥ 1 (thermodynamic)'),
    Line2D([0], [0], color='k', linestyle=':',  linewidth=1.8, label='$T_\\mathrm{sf}$ at floor (274 K)'),
    Line2D([0], [0], color='k', linestyle='-.', linewidth=0.8, label='Equilibrium temperature'),
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


def _recompute_T(instellation, P_CO2_bar, land_fraction=0.0):
    """Surface temperature from (instellation, final P_CO2), NOT from the JSON's stored 'T'.

    The stored 'T' is planet.py's self._T, a side effect set on EVERY call to dY_dt --
    including Jacobian finite-difference probes and solve_ivp's internal event root-finding
    trials -- so whichever call happened to run last is not guaranteed to correspond to the
    accepted final state that P_CO2 (also in the JSON) is taken from. Near a domain wall this
    can matter a lot: the climate response can be a genuine cliff, so a routine 1% Jacobian
    perturbation is enough to flip between a real state and the analytic model's literal 400.0
    "no equilibrium found" sentinel, producing a T that has nothing to do with the reported
    P_CO2 -- visible as a spurious temperature drop right at high-instellation domain-wall
    terminations. planet.py now recomputes T from the true final state before writing new
    output (see Planet.time_evolve), but this fixes it for JSONs already on disk, and is cheap
    (pure function, no PHREEQC) so it is applied to every row unconditionally.

    P_CO2 itself does NOT need this treatment: it is read directly from sol.y[0, -1], not from
    a side-effect channel, so it was never corrupted by the same mechanism.
    """
    albedo = LAND_ALBEDO * land_fraction + OCEAN_ALBEDO * (1 - land_fraction)
    P_CO2_Pa = float(np.clip(P_CO2_bar * 1e5, 0.0, 1e6))
    return float(get_T_surface_analytic(instellation * SOLAR_CONSTANT, P_CO2_Pa, albedo, False))


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
    T_surface = _recompute_T(float(d['instellation']), float(d['P_CO2']),
                             float(d.get('land_fraction', 0.0)))
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

        # Recompute pH the same way T is recomputed (see _recompute_T): self._pH in the JSON
        # is planet.py's side-effect value, subject to the identical last-call-wins corruption.
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


_DIAG_CACHE = {}


def _add_diag_columns(df, output_path):
    """Add da, calcite_si, alk_flux and (corrected) pH columns by re-reading each JSON file.

    Results are cached per run because several plots request diagnostics for
    overlapping subsets, and each one costs a couple of PHREEQC solves.

    Overwrites the 'pH' column (loaded from the JSON's possibly-corrupted stored value, see
    _recompute_T) with the recomputed one wherever _diag_from_json succeeded -- callers of this
    function already pay the PHREEQC cost the recompute needs, so the correction is free here.
    Rows where the recompute itself failed keep the original stored 'pH' rather than becoming
    NaN, since a stale-but-present value is more useful than none for a plot.
    """
    records = []
    for name in df['name']:
        fpath = os.path.join(output_path, f'{name}.json')
        if fpath in _DIAG_CACHE:
            records.append(_DIAG_CACHE[fpath])
            continue
        try:
            with open(fpath) as fh:
                d = json.load(fh)
            rec = _diag_from_json(d)
        except Exception:
            rec = _DIAG_NAN.copy()
        _DIAG_CACHE[fpath] = rec
        records.append(rec)
    diag_df = pd.DataFrame(records, index=df.index)
    df = df.assign(**diag_df.drop(columns=['pH']))
    df['pH'] = diag_df['pH'].where(diag_df['pH'].notna(), df['pH'])
    return df


def load_data(output_path):
    files = sorted(glob.glob(os.path.join(output_path, 'planet_*.json')))
    rows = []
    for f in files:
        with open(f) as fh:
            d = json.load(fh)
        if 'termination' not in d:
            print(f"  Skipping (no termination): {os.path.basename(f)}")
            continue

        name = d.get('name', '')
        comp_match = re.search(r'_comp_(.+?)(?:_rw|_fht_|$)', name)
        comp_name = comp_match.group(1) if comp_match else ''

        y_list = d.get('data', {}).get('y', [])
        salinity = _salinity_from_y(y_list) if y_list else np.nan

        rows.append({
            'name':               name,
            'instellation':       float(d['instellation']),
            'outgassing':         float(d['outgassing']),
            'crust_production':   float(d['crust_production_rate']),
            'reverse_weathering': bool(d.get('reverse_weathering', False)),
            'crust_carbonate':    float(d.get('crust_carbonate_content', 0.0)),
            'ocean_depth':        float(d['ocean_depth']),
            'comp_name':          comp_name,
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
            # 'T' recomputed from (instellation, P_CO2) rather than trusting the JSON's
            # stored 'T' -- see _recompute_T's docstring. Falls back to the stored value
            # when P_CO2 is missing (e.g. a run that errored before any state existed).
            'T':                  (_recompute_T(float(d['instellation']), float(d['P_CO2']),
                                                float(d.get('land_fraction', 0.0)))
                                    if d.get('P_CO2') is not None else d.get('T', np.nan)),
            'P_CO2':              d.get('P_CO2', np.nan),
            'pH':                 d.get('pH', np.nan),  # corrected in _diag_from_json when available
            'salinity':           salinity,
        })
    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} simulations.")
    return df


def _ref_crust(df):
    """Mask for the reference crust.

    Old sweeps tagged a named composition in the run name (`_comp_basalt_49`); the
    current model derives the mineralogy from mantle Mg/Si and dIW, so select Earth values.
    """
    named = df['comp_name'].astype(bool)
    legacy_ref  = named & (df['comp_name'] == 'basalt_49')
    derived_ref = (~named) & np.isclose(df['mg_si'], REF_MG_SI) & np.isclose(df['delta_iw'], REF_DIW)
    return legacy_ref | derived_ref


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
        (df['crust_carbonate'] == 0.0) &
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
    return len(handles) if presentation else fallback


def _save_fig(fig, path):
    fig.savefig(path, bbox_inches='tight')
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

def plot_faceted_lines(df, output_path, all_results=True, multiple_plots=False, split_panels=False, sequence=False):
    """T, P_CO2, pH, salinity vs instellation per crust rate, coloured by outgassing."""
    base = _base(df)
    base = _add_diag_columns(base, output_path)

    if all_results:
        crust_rates     = sorted(base['crust_production'].unique())
        outgassing_vals = sorted(base['outgassing'].unique())
    else:
        crust_rates = [0.1, 1, 10]
        outgassing_vals = [0.1, 0.3, 1, 3, 10]

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)

        if multiple_plots:
            for c in crust_rates:
                subset_c = base[base['crust_production'] == c]
                fig, axes = plt.subplots(n_rows, 1, figsize=(7, n_rows * 3), sharex=True)
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
        figsize = (fig_width_half * 2 * 2, n_rows * fig_subplot_height * 2) if presentation else (fig_width_half * 2, n_rows * fig_subplot_height)
        full_figsize = (6 * n_cols + 1.5, n_rows * 3) if all_results else figsize
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
        _save_fig(fig_c, os.path.join(output_path, fname))

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


def plot_ocean_depth_effect(df, output_path, show_markers=False, split_panels=False):
    """T, P_CO2, pH, salinity vs instellation for Earth-like tectonics, coloured by ocean depth."""
    # The depth sweep fixes (outgassing, crust) at their sweep defaults and varies ocean_depth.
    # That default changed over time (outgassing 1.0 -> 0.1), so instead of hardcoding a value
    # (which silently produced an empty/one-depth plot), pick whichever (outgassing, crust) pair
    # actually spans the most distinct depths.
    pool = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        _ref_chem(df) &
        (df['crust_carbonate'] == 0.0) &
        (df['land_fraction'] == 0.0)
    ]
    if pool.empty:
        print("No data for ocean depth sweep — skipping.")
        return
    ndepth = pool.groupby(['outgassing', 'crust_production'])['ocean_depth'].nunique()
    if ndepth.max() < 2:
        print("No ocean-depth variation in the data — skipping depth plot.")
        return
    best_o, best_c = ndepth.idxmax()
    subset = pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)]
    print(f"Ocean-depth plot: using outgassing={best_o:g}, crust={best_c:g} "
          f"({subset['ocean_depth'].nunique()} depths).")
    subset = _add_diag_columns(subset, output_path)

    depths = sorted(subset['ocean_depth'].unique())
    if len(depths) > 1 and (max(depths) / max(1, min(depths))) >= 10:
        norm = mcolors.LogNorm(vmin=min(depths), vmax=max(depths))
    else:
        norm = mcolors.Normalize(vmin=min(depths), vmax=max(depths))
    cmap = cmr.bubblegum_r

    min_x = subset['instellation'].min()
    max_x = subset['instellation'].max()
    margin = (max_x - min_x) * 0.05 if not pd.isna(min_x) and max_x != min_x else 0.1
    x_lims = (min_x - margin, max_x + margin) if not pd.isna(min_x) else (0.25, 1.45)

    ticks = depths if len(depths) <= 10 else None
    ticklabels = [f'{v:g}' for v in depths] if len(depths) <= 10 else None

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        figsize = (fig_width_half * 2, n_rows * fig_subplot_height * 2) if presentation else (fig_width_half, n_rows * fig_subplot_height)
        fig, axes = plt.subplots(n_rows, 1, figsize=figsize, sharex=True)
        for d in depths:
            group = subset[subset['ocean_depth'] == d].sort_values('instellation')
            if not group.empty:
                _plot_group_on_axes(axes, group, cmap(norm(d)), show_markers=show_markers, cols=cols)
        _style_axes(axes, cols, x_lims=x_lims)
        _add_colorbar(fig, list(axes), cmap, norm, 'Ocean Depth (m)',
                      ticks=ticks, ticklabels=ticklabels, aspect=n_rows * 7.5) # type: ignore
        _h = _make_legend_handles(show_markers=show_markers)
        fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 2 if show_markers else 1))
        _save_fig(fig, os.path.join(output_path, f'lines_ocean_depth{sfx}.png'))


def plot_chemistry_constants(df, output_path, show_markers=False, split_panels=False):
    """Sweeps 4/5: T, P_CO2, pH, salinity vs instellation for each chemistry-constant value.

    One figure per constant that actually varies (alpha, kd_mg, k_na); the other two are held
    at their most common value so each figure isolates a single knob.
    """
    pool_all = df[
        df['reverse_weathering'] &
        _ref_crust(df) &
        (df['crust_carbonate'] == 0.0) &
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

        # Don't hardcode the (outgassing, crust) the constants were swept at -- pick whichever
        # pair actually spans the most values, as the depth and crust plots do.
        nvals = pool.groupby(['outgassing', 'crust_production'])[col].nunique()
        if nvals.empty or nvals.max() < 2:
            print(f"No {col} variation at a fixed (outgassing, crust) — skipping {col}.")
            continue
        best_o, best_c = nvals.idxmax()
        subset = pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)]
        print(f"{col} plot: using outgassing={best_o:g}, crust={best_c:g} "
              f"({subset[col].nunique()} values).")
        subset = _add_diag_columns(subset, output_path)

        values = sorted(subset[col].unique())
        # Log scale where the values span decades, but only if none is zero (k_mg/k_na can be 0,
        # which is the 'term switched off' ablation).
        if min(values) > 0 and max(values) / min(values) >= 10:
            norm = mcolors.LogNorm(vmin=min(values), vmax=max(values))
        else:
            span = max(values) - min(values)
            norm = mcolors.Normalize(vmin=min(values) - 0.05 * span if span else min(values) - 1,
                                     vmax=max(values) + 0.05 * span if span else max(values) + 1)
        cmap = cmr.ember

        min_x, max_x = subset['instellation'].min(), subset['instellation'].max()
        margin = (max_x - min_x) * 0.05 if not pd.isna(min_x) and max_x != min_x else 0.1
        x_lims = (min_x - margin, max_x + margin) if not pd.isna(min_x) else (0.25, 1.45)

        ticks = values if len(values) <= 10 else None
        ticklabels = [f'{v:g}' for v in values] if len(values) <= 10 else None

        for cols, sfx in _panel_groups(split_panels):
            n_rows = len(cols)
            figsize = (fig_width_half * 2, n_rows * fig_subplot_height * 2) if presentation else (fig_width_half, n_rows * fig_subplot_height)
            fig, axes = plt.subplots(n_rows, 1, figsize=figsize, sharex=True)
            for v in values:
                group = subset[subset[col] == v].sort_values('instellation')
                if not group.empty:
                    _plot_group_on_axes(axes, group, cmap(norm(v)), show_markers=show_markers, cols=cols)
            _style_axes(axes, cols, x_lims=x_lims)
            _add_colorbar(fig, list(axes), cmap, norm, CHEM_KNOBS[col],
                          ticks=ticks, ticklabels=ticklabels, aspect=n_rows * 7.5) # type: ignore
            _h = _make_legend_handles(show_markers=show_markers)
            fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 2 if show_markers else 1))
            _save_fig(fig, os.path.join(output_path, f'lines_{col}{sfx}.png'))


def plot_crust_composition(df, output_path, split_panels=False, show_markers=False):
    """Sweep 3: T, P_CO2, pH, salinity vs instellation for each crust composition.

    Crust mineralogy is set by the mantle potential temperature (hotter mantle → more
    olivine-rich, lower-SiO₂ melt), so runs are grouped and coloured by T_p; at fixed
    T_p the mantle Mg/Si is used instead. Legacy sweeps that varied a named composition
    are grouped by that name.
    """
    pool = df[
        df['reverse_weathering'] &
        _ref_chem(df) &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0) &
        (df['land_fraction'] == 0.0)
    ].copy()

    # The crust sweep fixes (outgassing, crust) at their sweep defaults and varies the crust
    # knob (Mg-Si / dIW / named composition). That default drifted over time (outgassing 1.0 ->
    # 0.1), so don't hardcode it: pick the (outgassing, crust) pair that spans the most distinct
    # crust values. crust_knob is whichever of comp_name/mg_si/delta_iw actually varies.
    def _knobcol(g):
        for k in ['comp_name', 'mg_si', 'delta_iw']:
            col = g[k].astype(str) if k == 'comp_name' else g[k]
            if col.nunique() > 1:
                return k
        return None
    if _knobcol(pool) is None:
        print("No crust composition sweep data found — skipping.")
        return
    kc = _knobcol(pool)
    nvals = pool.groupby(['outgassing', 'crust_production'])[kc].nunique()
    best_o, best_c = nvals.idxmax()
    subset = pool[(pool['outgassing'] == best_o) & (pool['crust_production'] == best_c)]
    print(f"Crust-composition plot: knob={kc}, using outgassing={best_o:g}, crust={best_c:g} "
          f"({subset[kc].nunique()} values).")

    # Pick whichever crust knob actually varies in this sweep.
    if subset['comp_name'].astype(bool).any() and subset['comp_name'].nunique() > 1:
        key, label, fmt = 'comp_name', 'Crust composition', lambda v: COMP_LABELS.get(v, v)
        values = [c for c in ['komatiite_42', 'komatiite_44', 'basalt_47', 'basalt_49', 'basalt_51']
                  if c in set(subset['comp_name'])]
    elif subset['mg_si'].nunique() > 1:
        key, label, fmt = 'mg_si', 'Mantle Mg/Si', lambda v: f'{v:g}'
        values = sorted(subset['mg_si'].unique())
    elif subset['delta_iw'].nunique() > 1:
        key, label, fmt = 'delta_iw', r'Core-formation $\Delta$IW', lambda v: f'{v:+g}'
        values = sorted(subset['delta_iw'].unique())
    else:
        print("No crust composition sweep data found — skipping.")
        return

    subset = _add_diag_columns(subset, output_path)

    cmap = cmr.gem
    if key == 'comp_name':
        # Colour by the SiO₂ content encoded in the name
        numeric = [int(c.split('_')[-1]) for c in values]
        norm = mcolors.Normalize(vmin=42, vmax=53)
        cbar_label = 'SiO₂ content (%)'
        ticklabels = [f'{v}%' for v in numeric]
    else:
        numeric = [float(v) for v in values]
        span = max(numeric) - min(numeric)
        norm = mcolors.Normalize(vmin=min(numeric) - 0.05 * span if span else min(numeric) - 1,
                                 vmax=max(numeric) + 0.05 * span if span else max(numeric) + 1)
        cbar_label = label
        ticklabels = [fmt(v) for v in values]

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        figsize = (fig_width_half * 2, n_rows * fig_subplot_height * 2) if presentation else (fig_width_half, n_rows * fig_subplot_height)
        fig, axes = plt.subplots(n_rows, 1, figsize=figsize, sharex=True)
        for value, num in zip(values, numeric):
            group = subset[subset[key] == value].sort_values('instellation')
            if not group.empty:
                _plot_group_on_axes(axes, group, cmap(norm(num)), cols=cols, show_markers=show_markers)
        _style_axes(axes, cols)
        _add_colorbar(fig, list(axes), cmap, norm, cbar_label,
                      ticks=numeric, ticklabels=ticklabels)
        _h = _make_legend_handles(show_markers=show_markers)
        fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 2 if show_markers else 1))
        _save_fig(fig, os.path.join(output_path, f'lines_crust_composition{sfx}.png'))


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

    figsize = (fig_width_half * 2, fig_subplot_height * 2 * 2) if presentation else (fig_width_half, fig_subplot_height * 2)
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
        fig, axes = plt.subplots(nrows, 1, figsize=(3.5, nrows * 2),
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
        fig, _ax2d = plt.subplots(2, 2, figsize=(fig_width_half * 2, fig_subplot_height * 2 * 2),
                                   sharex=True)
        axes_summary = _ax2d.ravel()
    else:
        fig, axes_summary = plt.subplots(len(cols), 1, figsize=(3.5, len(cols) * 2), sharex=True)

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
    figsize2 = (fig_width_half * 2, fig_subplot_height * 2) if presentation else (fig_width_half, fig_subplot_height * 1.5)
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
        # T_seafloor/T_pore above are already derived from the recomputed T_surface (see
        # _pore_conditions); recover it rather than re-reading the JSON's raw (possibly
        # corrupted) 'T', which _pore_conditions no longer trusts. T_seafloor = max(1.02*T-16.7, 274).
        T_surface = _recompute_T(float(d['instellation']), float(d['P_CO2']),
                                 float(d.get('land_fraction', 0.0)))

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

        fig, axes = plt.subplots(n_min, 2, figsize=(3.5, n_min * 1),
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

    cond_unknown = (
        (base['termination'].isin(OUT_OF_DOMAIN) & wall.isin(['co2_high', 'co2_low']))
        # Legacy runs: pre-domain-event output has no domain_wall, and the old ceiling /
        # acid_ocean cutoffs stopped mid-evolution, so their fate is likewise unknown.
        | base['termination'].isin(['co2_ceiling', 'acid_ocean'])
    )

    cond_snow = ((base['T'] <= T_SNOWBALL) | base['termination'].isin(['snowball', 'co2_floor'])) & ~cond_unknown
    cond_hot  = ((base['T'] >= T_RUNAWAY)  | (base['termination'] == 'hothouse')) & ~cond_unknown & ~cond_snow
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
    figsize = (fig_width_half * 2, fig_subplot_height * 3) if presentation else (fig_width_half, fig_subplot_height * 2)
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
    args = parser.parse_args()

    for _knob in CHEM_KNOBS:
        _v = getattr(args, _knob)
        if _v is not None:
            CHEM_OVERRIDE[_knob] = _v

    df = load_data(args.path)

    if df.empty:
        print("No data found. Check --path.")
    else:
        plot_faceted_lines(df, args.path)
        plot_faceted_lines(df, args.path, all_results=False, split_panels=True)
        plot_faceted_lines(df, args.path, split_panels=True, all_results=False, sequence=True)
        plot_ocean_depth_effect(df, args.path, split_panels=presentation)
        plot_chemistry_constants(df, args.path, split_panels=presentation)
        plot_ratio_scatter(df, args.path)
        plot_crust_composition(df, args.path, show_markers=False, split_panels=presentation)
        plot_damkohler_contour(df, args.path)
        plot_habitability_phase_space(df, args.path)
        plot_continental_baseline(df, args.path)
        print("Done.")
