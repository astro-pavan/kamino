import os
import sys
import glob
import json
import re
import argparse
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
from kamino.weathering import get_weathering_flux, J_ref_normalised, rate_ref, A_seafloor
from kamino.constants import G, EARTH_CRUST_PRODUCTION_RATE_PER_AREA, YR
from kamino.chemistry import alk_idx
from kamino.mineral_info import carbonate_minerals, clay_minerals, silica_minerals, reverse_weathering_minerals

fig_width_half = 3.5
fig_subplot_height = 1.5

presentation = True

if presentation:
    plt.style.use('experiments/planetary-chem-presentation.mplstyle')
else:
    plt.style.use('experiments/planetary-chem-paper.mplstyle')

DEFAULT_OUTPUT_PATH = '/data/pt426/kamino_experiments_fast_2/'

TERM_LABELS = {
    'snowball':   'Snowball',
    'hothouse':   'Hothouse',
    'acid_ocean': 'Acid Ocean (>5 bar CO₂)',
    'converged':  'Converged',
    'co2_floor':  'CO₂ floor (depleted)',
    'timeout':    'Timeout (2 Gyr)',
}
HABITABLE = {'converged', 'timeout', 'co2_floor'}
T_SNOWBALL = 260.0
T_RUNAWAY  = 360.0

# Molar masses for b_ocean elements (g/mol), skipping Alkalinity (index 2 of y):
# y[3]=C×61, y[4]=Si×60.1, y[5]=Al×27, y[6]=Fe×55.8, y[7]=Ca×40.1, y[8]=Mg×24.3,
# y[9]=Na×23.0, y[10]=Cl×35.45  (added with NaCl chemistry)
_SAL_INDICES = [3, 4, 5, 6, 7, 8, 9, 10]
_SAL_MASSES  = [61.0, 60.1, 27.0, 55.8, 40.1, 24.3, 23.0, 35.45]

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
FAILED_MARKERS = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

DA_LEGEND = [
    Line2D([0], [0], color='k', linestyle='-',  linewidth=1.8, label='Da < 1 (kinetic)'),
    Line2D([0], [0], color='k', linestyle='--', linewidth=1.8, label='Da ≥ 1 (thermodynamic)'),
    Line2D([0], [0], color='k', linestyle=':',  linewidth=1.8, label='$T_\\mathrm{sf}$ at floor (274 K)'),
]

PANEL_COLS = ['T', 'P_CO2', 'pH', 'salinity']


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


_DIAG_NAN = {'da': np.nan, 'calcite_si': np.nan, 'ocean_si': np.nan, 'alk_flux': np.nan}


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
        b_ocean = np.maximum(np.array([float(y_list[i][-1]) for i in range(2, 2 + n_elements)]), 0.0)

        mass      = float(d.get('mass',   5.972e24))
        radius    = float(d.get('radius', 6.371e6))
        gravity   = G * mass / radius**2
        ocean_depth  = float(d['ocean_depth'])
        P_background = float(d['background_pressure'])
        T_surface    = float(d['T'])
        P_CO2        = float(d['P_CO2']) * 1e5          # bar → Pa
        P_H2O        = float(y_list[1][-1]) if y_list[1] else 0.0

        P_surface  = P_background + P_CO2 + P_H2O
        T_seafloor = max(1.02 * T_surface - 16.7, 274.0)
        T_pore     = T_seafloor + 9
        P_pore     = P_surface + 1000 * gravity * ocean_depth

        crust_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * float(d['crust_production_rate'])
        f_HT       = float(d.get('f_HT', 0.0))
        J_LT       = J_ref_normalised * (crust_rate / rate_ref) * (1 - f_HT)

        crust_composition = d.get('crust_composition', {})
        rw = bool(d.get('reverse_weathering', False))
        pore_minerals = carbonate_minerals + clay_minerals  # reverse weathering is in ocean sediments, not pore space

        flux, diag = get_weathering_flux(
            P_pore, T_pore, P_CO2, b_ocean,
            rate=crust_rate, J=J_LT,
            crust_composition=crust_composition,
            precipitating_minerals=pore_minerals,
        )

        # Pore SI: taken from secondary_SI which holds the pre-precipitation SI.
        # This is always the relevant metric for pore-space calcite precipitation.
        calcite_si = float(diag.get('secondary_SI', {}).get('Calcite', np.nan))

        # Ocean SI: only reliable when Ca_ocean > 0; set to NaN otherwise to avoid
        # the spurious -∞ that results when ocean Ca has been depleted to the ODE floor.
        from kamino.precipitation import get_precipitation
        from kamino.chemistry import ChemistryError, ca_idx as _ca_idx
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
                'ocean_si': ocean_si, 'alk_flux': alk_flux}
    except Exception:
        return _DIAG_NAN.copy()


def _add_diag_columns(df, output_path):
    """Add da, calcite_si, and alk_flux columns by re-reading each JSON file."""
    records = []
    for name in df['name']:
        fpath = os.path.join(output_path, f'{name}.json')
        try:
            with open(fpath) as fh:
                d = json.load(fh)
            records.append(_diag_from_json(d))
        except Exception:
            records.append(_DIAG_NAN.copy())
    diag_df = pd.DataFrame(records, index=df.index)
    return df.assign(**diag_df)


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
            'f_HT':               float(d.get('f_HT', 0.0)),
            'land_fraction':      float(d.get('land_fraction', 0.0)),
            'termination':        d['termination'],
            'end_time_yr':        d.get('end_time_yr', np.nan),
            'T':                  d.get('T', np.nan),
            'P_CO2':              d.get('P_CO2', np.nan),
            'pH':                 d.get('pH', np.nan),
            'salinity':           salinity,
        })
    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} simulations.")
    return df


def _base(df):
    """Sweep 1: basalt_49, rw=True, depth=3000, f_HT=0, outgassing>0, ocean world."""
    return df[
        df['reverse_weathering'] &
        (df['comp_name'] == 'basalt_49') &
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
        elif col == 'P_CO2':
            ax.set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
            ax.set_yscale('log')
            ax.set_ylim(1e-8, 20)
        elif col == 'pH':
            ax.set_ylabel('Ocean pH')
            ax.set_ylim(4.5, 12)
        elif col == 'salinity':
            ax.set_ylabel('Dissolved Ions (g/kg)')
            ax.set_yscale('log')
            ax.set_ylim(1e-2, 1e3)
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
                           facecolors='none', edgecolors=color, zorder=4, linewidths=1.4)


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

def plot_faceted_lines(df, output_path, all_results=True, split_panels=False):
    """T, P_CO2, pH, salinity vs instellation per crust rate, coloured by outgassing."""
    base = _base(df)
    base = _add_diag_columns(base, output_path)

    if all_results:
        crust_rates     = sorted(base['crust_production'].unique())
        outgassing_vals = sorted(base['outgassing'].unique())
    else:
        crust_rates = [0.1, 1, 10]
        outgassing_vals = sorted(base['outgassing'].unique())

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)

        if all_results:
            for c in crust_rates:
                subset_c = base[base['crust_production'] == c]
                fig, axes = plt.subplots(n_rows, 1, figsize=(7, n_rows * 3), sharex=True)
                for o in outgassing_vals:
                    group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
                    if not group.empty:
                        _plot_group_on_axes(axes, group, cmap(norm(o)), cols=cols)
                _style_axes(axes, cols)
                _add_colorbar(fig, list(axes), cmap, norm, 'Outgassing',
                              ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                              aspect=n_rows * 7.5)
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
        _add_colorbar(fig_c, list(axes_c.ravel()), cmap, norm, 'Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_rows * 10)
        fig_c.suptitle('Crust production rate')
        fname = f'lines_combined{"_full" if all_results else ""}{sfx}.png'
        _save_fig(fig_c, os.path.join(output_path, fname))


def plot_ocean_depth_effect(df, output_path, show_markers=False, split_panels=False):
    """T, P_CO2, pH, salinity vs instellation for Earth-like tectonics, coloured by ocean depth."""
    subset = df[
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0) &
        df['reverse_weathering'] &
        (df['comp_name'] == 'basalt_49') &
        (df['crust_carbonate'] == 0.0)
    ]
    if subset.empty:
        print("No data for ocean depth sweep — skipping.")
        return
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
                      ticks=ticks, ticklabels=ticklabels, aspect=n_rows * 7.5)
        _h = _make_legend_handles(show_markers=show_markers)
        fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 2 if show_markers else 1))
        _save_fig(fig, os.path.join(output_path, f'lines_ocean_depth{sfx}.png'))


def plot_crust_composition(df, output_path, split_panels=False, show_markers=False):
    """Sweep 3: T, P_CO2, pH, salinity vs instellation for each crust composition (Earth-like baseline)."""
    compositions = ['komatiite_44', 'basalt_47', 'basalt_49', 'basalt_51']

    subset = df[
        df['comp_name'].isin(compositions) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0) &
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0)
    ]
    if subset.empty:
        print("No crust composition sweep data found — skipping.")
        return
    subset = _add_diag_columns(subset, output_path)

    sio2_vals = [int(c.split('_')[-1]) for c in compositions]  # [44, 47, 49, 51]
    cmap = cmr.gem
    norm = mcolors.Normalize(vmin=42, vmax=53)

    for cols, sfx in _panel_groups(split_panels):
        n_rows = len(cols)
        figsize = (fig_width_half * 2, n_rows * fig_subplot_height * 2) if presentation else (fig_width_half, n_rows * fig_subplot_height)
        fig, axes = plt.subplots(n_rows, 1, figsize=figsize, sharex=True)
        for comp, sio2 in zip(compositions, sio2_vals):
            group = subset[subset['comp_name'] == comp].sort_values('instellation')
            if not group.empty:
                _plot_group_on_axes(axes, group, cmap(norm(sio2)), cols=cols, show_markers=show_markers)
        _style_axes(axes, cols)
        _add_colorbar(fig, list(axes), cmap, norm, 'SiO₂ content (%)',
                      ticks=sio2_vals, ticklabels=[f'{v}%' for v in sio2_vals])
        _h = _make_legend_handles(show_markers=show_markers)
        fig.legend(handles=_h, loc='outside lower center', ncol=_legend_ncol(_h, 2 if show_markers else 1))
        # fig.suptitle('Crust Composition Effect (out=1×, crust=1×, depth=3000 m, rw=True)')
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

# All minerals that can precipitate, in display order
_PORE_MINERALS  = carbonate_minerals + clay_minerals + reverse_weathering_minerals
_OCEAN_MINERALS = carbonate_minerals + clay_minerals + silica_minerals + reverse_weathering_minerals
_ALL_MINERALS   = list(dict.fromkeys(_PORE_MINERALS + _OCEAN_MINERALS))  # ordered, unique

_MINERAL_LABELS = {
    'Calcite':     'Calcite',
    'Siderite':    'Siderite (FeCO₃)',
    'Kaolinite':   'Kaolinite',
    'Goethite':    'Goethite',
    'SiO2(am)':   'Amorphous SiO₂',
    'Saponite-Mg': 'Saponite-Mg',
}


def plot_damkohler_contour(df, output_path, out_targets=(0.1, 1.0, 10.0)):
    """2D contourf of Damköhler number in (instellation × crust production) space.

    One panel per outgassing rate (vertically stacked), coloured by log10(Da).
    The Da = 1 boundary is drawn as a black contour.
    Uses the basalt_49, rw=True, depth=3000, f_HT=0 baseline.
    """
    subset = df[
        (df['comp_name'] == 'basalt_49') &
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
        cbar_ax = fig.add_axes([0.85, pos_bot.y0, 0.03, pos_top.y1 - pos_bot.y0])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax, label=r'$\log_{10}(\mathrm{Damkohler Coefficient})$')
        cbar.ax.axhline(0, color='k', linewidth=1.5)

        # fig.suptitle('Damköhler number — basalt_49, rw=True, depth = 3000 m', fontsize=10)
        _save_fig(fig, os.path.join(output_path, 'da_contour.png'))


def plot_continental_baseline(df, output_path):
    """T, P_CO2, pH, salinity, and individual ion concentrations vs instellation
    for the Earth-like continental baseline.

    Runs with land_fraction=0.3, basalt_49, rw=True, out=1×, crust=1×, depth=3000 m.
    """
    subset = df[
        (df['land_fraction'] == 0.3) &
        (df['comp_name'] == 'basalt_49') &
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
    group = subset.sort_values('instellation')

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
        (1, 'DIC',  2.0),
        (2, 'Si',   0.1),
        (5, 'Ca',  10.3),
        (6, 'Mg',  52.8),
        (7, 'Na', 480.0),
        (8, 'Cl', 550.0),
    ]
    ION_COLORS = [plt.cm.tab10(k / 10) for k in range(len(ION_SPEC))]

    # Load final b_ocean for each run from the JSON files
    ion_rows = []
    for _, row in group.iterrows():
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
    fig, axes_summary = plt.subplots(len(cols), 1, figsize=(3.5, len(cols) * 2), sharex=True)
    _plot_group_on_axes(axes_summary, group, color='k', show_markers=False, cols=cols)
    _style_axes(axes_summary, cols)

    for ax, col in zip(axes_summary, cols):
        ax.scatter(EARTH_S, earth_vals[col], marker='*', s=220, color='blue',
                   edgecolors='k', linewidths=0.7, zorder=6)
    axes_summary[0].annotate(
        'Earth', xy=(EARTH_S, EARTH_T), xytext=(EARTH_S + 0.06, EARTH_T - 6),
        fontsize=8, arrowprops=dict(arrowstyle='-', color='k', lw=0.8),
    )
    _save_fig(fig, os.path.join(output_path, 'continental_baseline.png'))

    # --- figure 2: ion concentrations ---
    fig2, ax_ions = plt.subplots(1, 1, figsize=(3.5, 3))

    if ion_rows:
        s_vals = np.array([r[0] for r in ion_rows])
        b_mmol = np.array([r[1] for r in ion_rows]) * 1e3  # mol/kg → mmol/kg

        for (idx, label, earth_val), color in zip(ION_SPEC, ION_COLORS):
            ax_ions.plot(s_vals, b_mmol[:, idx], color=color, linewidth=1.5, label=label)
            if earth_val is not None:
                ax_ions.scatter(EARTH_S, earth_val, marker='*', s=150,
                                color=color, edgecolors='k', linewidths=0.5, zorder=6)

    ax_ions.set_ylabel('Concentration (mmol/kg)')
    ax_ions.set_yscale('log')
    ax_ions.set_xlabel('Instellation ($S/S_0$)')
    ax_ions.grid(True, linestyle='--', alpha=0.4)
    ax_ions.set_xlim(0.25, 1.45)
    ax_ions.legend(ncols=2, fontsize=7, loc='lower left',
                   handlelength=1.2, columnspacing=0.8, labelspacing=0.3)
    _save_fig(fig2, os.path.join(output_path, 'continental_baseline_ions.png'))


def _get_mineral_si(d):
    """Return {'pore': {mineral: SI}, 'ocean': {mineral: SI}, 'da': float, 'T': float}
    for all precipitating minerals in pore space and ocean."""
    from kamino.precipitation import get_precipitation
    from kamino.chemistry import ChemistryError
    nan_result = {'pore': {}, 'ocean': {}, 'da': np.nan, 'T': np.nan}
    try:
        y_list = d.get('data', {}).get('y', [])
        if not y_list or len(y_list) < 9:
            return nan_result
        b_ocean = np.maximum(np.array([float(y_list[i][-1]) for i in range(2, 9)]), 0.0)

        mass      = float(d.get('mass',   5.972e24))
        radius    = float(d.get('radius', 6.371e6))
        gravity   = G * mass / radius**2
        P_bg      = float(d['background_pressure'])
        T_surface = float(d['T'])
        P_CO2     = float(d['P_CO2']) * 1e5
        P_H2O     = float(y_list[1][-1]) if y_list[1] else 0.0
        P_surface = P_bg + P_CO2 + P_H2O
        T_seafloor = max(1.02 * T_surface - 16.7, 274.0)
        T_pore     = T_seafloor + 9
        P_pore     = P_surface + 1000 * gravity * float(d['ocean_depth'])

        crust_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * float(d['crust_production_rate'])
        J_LT       = J_ref_normalised * (crust_rate / rate_ref) * (1 - float(d.get('f_HT', 0.0)))
        cc         = d.get('crust_composition', {})
        rw         = bool(d.get('reverse_weathering', False))
        pore_min   = carbonate_minerals + clay_minerals + (reverse_weathering_minerals if rw else [])
        ocean_min  = carbonate_minerals + clay_minerals + silica_minerals + (reverse_weathering_minerals if rw else [])

        _, diag = get_weathering_flux(
            P_pore, T_pore, P_CO2, b_ocean,
            rate=crust_rate, J=J_LT,
            crust_composition=cc,
            precipitating_minerals=pore_min,
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

        _add_colorbar(fig, list(axes.ravel()), cmap, norm, 'Outgassing',
                      ticks=outgassing_vals, ticklabels=[f'{v}×' for v in outgassing_vals],
                      aspect=n_min * 8)
        fig.legend(handles=DA_LEGEND, loc='outside lower center', ncol=_legend_ncol(DA_LEGEND, 3))
        fig.suptitle(f'Mineral saturation indices — Crust = {c}× Earth')
        _save_fig(fig, os.path.join(output_path, f'mineral_si_crust_{c}.png'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot kamino parameter sweep results.')
    parser.add_argument('--path', default=DEFAULT_OUTPUT_PATH,
                        help='Directory containing planet_*.json files.')
    args = parser.parse_args()

    df = load_data(args.path)

    if df.empty:
        print("No data found. Check --path.")
    else:
        plot_faceted_lines(df, args.path)
        plot_faceted_lines(df, args.path, all_results=False, split_panels=True)
        plot_ocean_depth_effect(df, args.path, split_panels=presentation)
        plot_ratio_scatter(df, args.path)
        plot_crust_composition(df, args.path, show_markers=False, split_panels=presentation)
        plot_damkohler_contour(df, args.path)
        # plot_continental_baseline(df, args.path)
        #plot_magnesium_chemistry(df, args.path)
        print("Done.")
