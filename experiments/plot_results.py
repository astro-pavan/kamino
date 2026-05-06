import os
import glob
import json
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import cmasher as cmr

DEFAULT_OUTPUT_PATH = '/data/pt426/kamino_experiments/'

TERM_COLORS = {
    'snowball':   '#777777',
    'hothouse':   '#777777',
    'acid_ocean': '#777777',
    'converged':  '#6abf69',
    'co2_floor':  '#a0c878',
    'timeout':    '#b8b8b8',
}
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
# y[3]=C×61, y[4]=Si×60.1, y[5]=Al×27, y[6]=Fe×55.8, y[7]=Ca×40.1, y[8]=Mg×24.3
_SAL_INDICES = [3, 4, 5, 6, 7, 8]
_SAL_MASSES  = [61.0, 60.1, 27.0, 55.8, 40.1, 24.3]

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


def _salinity_from_y(y_list):
    try:
        return sum(
            float(y_list[i][-1]) * mass
            for i, mass in zip(_SAL_INDICES, _SAL_MASSES)
            if len(y_list) > i and len(y_list[i]) > 0
        )
    except Exception:
        return np.nan


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
    """Sweep 1: basalt_49, rw=True, depth=3000, f_HT=0, outgassing>0."""
    return df[
        df['reverse_weathering'] &
        (df['comp_name'] == 'basalt_49') &
        (df['crust_carbonate'] == 0.0) &
        (df['ocean_depth'] == 3000) &
        (df['f_HT'] == 0.0) &
        (df['outgassing'] > 0)
    ]


def _style_4panel_axes(axes, x_lims=(0.25, 1.45)):
    axes[0].set_ylabel('Temperature (K)')
    axes[0].axhspan(T_SNOWBALL - 20, T_SNOWBALL, color='red', alpha=0.12)
    axes[0].axhspan(T_RUNAWAY - 10,  T_RUNAWAY,  color='blue',  alpha=0.12)
    axes[0].set_ylim(235, 360)

    axes[1].set_ylabel('P_CO₂ (bar)')
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-8, 20)

    axes[2].set_ylabel('Ocean pH')
    axes[2].set_ylim(4.5, 12)

    axes[3].set_ylabel('Dissolved Ions (g/kg)')
    axes[3].set_yscale('log')
    axes[3].set_ylim(1e-4, 1e2)
    axes[3].set_xlabel('Instellation (S/S₀)')

    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(*x_lims)


def _plot_group_on_axes(axes, group, color, marker_map, linestyle='-'):
    cols = ['T', 'P_CO2', 'pH', 'salinity']
    hab = group[group['termination'].isin(HABITABLE)]
    for ax, col in zip(axes, cols):
        if not hab.empty:
            ax.plot(hab['instellation'], hab[col], color=color, linewidth=1.8,
                    linestyle=linestyle, zorder=3)
            ax.scatter(hab['instellation'], hab[col], color=color, s=28, zorder=4)

    # Non-habitable markers only on the temperature panel
    for _, row in group[~group['termination'].isin(HABITABLE)].iterrows():
        marker = marker_map.get(row['termination'], 'x')
        ec = TERM_COLORS.get(row['termination'], 'k')
        val = row['T']
        if np.isfinite(val):
            axes[0].scatter(row['instellation'], val, marker=marker,
                            s=55, color=ec, zorder=4, linewidths=1.2)


def plot_termination_map(df, output_path):
    """Grid of termination types: rows=outgassing, cols=crust_production, panels=instellation."""
    base = _base(df)
    instellations = sorted(base['instellation'].unique())
    outgassings   = sorted(base['outgassing'].unique())
    crust_rates   = sorted(base['crust_production'].unique())

    ncols = 3
    nrows = -(-len(instellations) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, nrows * 3.6))
    axes = np.array(axes).flatten()

    for ax, s in zip(axes, instellations):
        subset = base[base['instellation'] == s]
        termZ = {(row['outgassing'], row['crust_production']): row['termination']
                 for _, row in subset.iterrows()}

        for yi, o in enumerate(outgassings):
            for xi, c in enumerate(crust_rates):
                term = termZ.get((o, c), None)
                color = TERM_COLORS.get(term, '#e8e8e8')
                ax.add_patch(plt.Rectangle((xi - 0.5, yi - 0.5), 1, 1, color=color, zorder=1))

        ax.set_xlim(-0.5, len(crust_rates) - 0.5)
        ax.set_ylim(-0.5, len(outgassings) - 0.5)
        ax.set_xticks(range(len(crust_rates)))
        ax.set_yticks(range(len(outgassings)))
        ax.set_xticklabels([f'{v}×' for v in crust_rates])
        ax.set_yticklabels([f'{v}×' for v in outgassings])
        ax.set_xlabel('Crust production rate')
        ax.set_ylabel('Outgassing')
        ax.set_title(f'S/S₀ = {s}', fontsize=11)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    for ax in axes[len(instellations):]:
        ax.set_visible(False)

    patches = [mpatches.Patch(color=TERM_COLORS[t], label=TERM_LABELS[t])
               for t in TERM_COLORS]
    fig.legend(handles=patches, title='Outcome', fontsize=9,
               loc='lower right', bbox_to_anchor=(0.98, 0.01))
    fig.suptitle('Long-run Outcome Phase Space', fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout()
    path = os.path.join(output_path, 'termination_map.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {path}")


def plot_faceted_lines(df, output_path):
    """T, P_CO2, pH, salinity vs instellation per crust rate, coloured by outgassing."""
    base = _base(df)
    crust_rates     = sorted(base['crust_production'].unique())
    outgassing_vals = sorted(base['outgassing'].unique())

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical
    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    for c in crust_rates:
        subset_c = base[base['crust_production'] == c]
        fig, axes = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
        fig.subplots_adjust(hspace=0.05, right=0.84)

        for o in outgassing_vals:
            group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
            if group.empty:
                continue
            _plot_group_on_axes(axes, group, cmap(norm(o)), marker_map)

        _style_4panel_axes(axes)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ax = fig.add_axes([0.86, 0.08, 0.03, 0.84])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Outgassing factor')
        cbar.set_ticks(outgassing_vals)
        cbar.set_ticklabels([f'{v}×' for v in outgassing_vals])

        for term, marker in marker_map.items():
            axes[0].scatter([], [], marker=marker, color=TERM_COLORS[term],
                            s=50, label=TERM_LABELS[term])
        axes[0].legend(fontsize=7, loc='upper left')

        fig.suptitle(f'Crust production = {c}× Earth', fontsize=13, fontweight='bold')
        path = os.path.join(output_path, f'lines_crust_{c}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {path}")

    # Combined plot: crust rates as columns
    n_cols = len(crust_rates)
    fig_c, axes_c = plt.subplots(4, n_cols, figsize=(6 * n_cols + 1.5, 12),
                                  sharex=True, sharey='row')
    fig_c.subplots_adjust(hspace=0.05, wspace=0.05, right=0.88)

    for ci, c in enumerate(crust_rates):
        subset_c = base[base['crust_production'] == c]
        col_axes = axes_c[:, ci]

        for o in outgassing_vals:
            group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
            if group.empty:
                continue
            _plot_group_on_axes(col_axes, group, cmap(norm(o)), marker_map)

        _style_4panel_axes(col_axes)
        axes_c[0, ci].set_title(f'{c}×', fontsize=10)

        if ci > 0:
            for ax in col_axes:
                ax.set_ylabel('')
            for row in range(4):
                plt.setp(axes_c[row, ci].get_yticklabels(), visible=False)

        if ci != n_cols // 2:
            col_axes[3].set_xlabel('')

    for term, marker in marker_map.items():
        axes_c[0, 0].scatter([], [], marker=marker, color=TERM_COLORS[term],
                             s=50, label=TERM_LABELS[term])
    axes_c[0, 0].legend(fontsize=6, loc='upper left')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig_c.add_axes([0.90, 0.08, 0.02, 0.84])
    cbar = fig_c.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Outgassing factor')
    cbar.set_ticks(outgassing_vals)
    cbar.set_ticklabels([f'{v}×' for v in outgassing_vals])

    fig_c.suptitle('Sweep 1: Crust production rate →', fontsize=13, fontweight='bold')
    path = os.path.join(output_path, 'lines_combined.png')
    fig_c.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig_c)
    print(f"Saved {path}")


def plot_heatmaps(df, output_path):
    """T, P_CO2, pH, salinity as heatmaps over (crust_production, outgassing) at each instellation."""
    base = _base(df)
    instellations = sorted(base['instellation'].unique())
    outgassings   = sorted(base['outgassing'].unique())
    crust_rates   = sorted(base['crust_production'].unique())

    nc = len(crust_rates)
    no = len(outgassings)
    x_edges = np.arange(-0.5, nc)
    y_edges = np.arange(-0.5, no)

    vars_cfg = [
        ('T',        'Temperature (K)',       'RdBu_r',  mcolors.Normalize(vmin=270, vmax=350)),
        ('P_CO2',    'P_CO₂ (bar)',           'inferno', mcolors.LogNorm(vmin=1e-8, vmax=10)),
        ('pH',       'Ocean pH',              'cividis', mcolors.Normalize(vmin=5, vmax=9)),
        ('salinity', 'Dissolved Ions (g/kg)', 'YlOrRd',  mcolors.LogNorm(vmin=1e-3, vmax=10)),
    ]

    for var, label, cmap_name, norm in vars_cfg:
        fig, axes = plt.subplots(1, len(instellations),
                                 figsize=(len(instellations) * 3.0, 3.6), sharey=True)

        mesh = None
        for ax, s in zip(axes, instellations):
            subset = base[base['instellation'] == s]
            Z     = np.full((no, nc), np.nan)
            termZ = np.full((no, nc), '', dtype=object)

            for _, row in subset.iterrows():
                xi = crust_rates.index(row['crust_production'])
                yi = outgassings.index(row['outgassing'])
                termZ[yi, xi] = row['termination']
                if row['termination'] in HABITABLE:
                    Z[yi, xi] = row[var]

            mesh = ax.pcolormesh(x_edges, y_edges, Z,
                                 cmap=cmap_name, norm=norm, zorder=1)

            for yi in range(no):
                for xi in range(nc):
                    term = termZ[yi, xi]
                    if term and term not in HABITABLE:
                        color = TERM_COLORS.get(term, '#aaaaaa')
                        ax.add_patch(plt.Rectangle(
                            (xi - 0.5, yi - 0.5), 1, 1,
                            color=color, alpha=0.85, zorder=2))

            ax.set_xticks(range(nc))
            ax.set_yticks(range(no))
            ax.set_xticklabels([f'{v}×' for v in crust_rates], fontsize=8)
            ax.set_yticklabels([f'{v}×' for v in outgassings], fontsize=8)
            ax.set_xlabel('Crust production', fontsize=9)
            ax.set_title(f'S/S₀ = {s}', fontsize=10)
            ax.tick_params(length=0)

        axes[0].set_ylabel('Outgassing', fontsize=9)

        if mesh is not None:
            cbar = fig.colorbar(mesh, ax=axes[-1], fraction=0.046, pad=0.04)
            cbar.set_label(label, fontsize=9)

        patches = [mpatches.Patch(color=TERM_COLORS[t], label=TERM_LABELS[t])
                   for t in ['snowball', 'hothouse', 'acid_ocean']]
        fig.legend(handles=patches, fontsize=7, title='Non-habitable',
                   loc='lower right', bbox_to_anchor=(1.22, 0.0))

        clean = var.replace('_', '').replace('/', '')
        fig.suptitle(f'Final {label}', fontsize=11, fontweight='bold')
        fig.tight_layout()
        path = os.path.join(output_path, f'heatmap_{clean}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {path}")


def plot_zero_outgassing_carb01(df, output_path):
    """T, P_CO2, pH, salinity vs instellation for outgassing=0, crust_carbonate=0.1, coloured by crust rate."""
    subset = df[
        (df['outgassing'] == 0.0) &
        (df['crust_carbonate'] == 0.1) &
        (df['ocean_depth'] == 3000)
    ]
    if subset.empty:
        print("No data for out=0, carb=0.1 — skipping.")
        return

    crust_rates = sorted(subset['crust_production'].unique())
    norm = mcolors.LogNorm(vmin=min(crust_rates), vmax=max(crust_rates))
    cmap = plt.get_cmap('viridis')
    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    fig, axes = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
    fig.subplots_adjust(hspace=0.05, right=0.84)

    for c in crust_rates:
        group = subset[subset['crust_production'] == c].sort_values('instellation')
        if group.empty:
            continue
        _plot_group_on_axes(axes, group, cmap(norm(c)), marker_map)

    _style_4panel_axes(axes, x_lims=(0.2, 1.5))

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.86, 0.08, 0.03, 0.84])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Crust production rate (×Earth)')
    cbar.set_ticks(crust_rates)
    cbar.set_ticklabels([f'{v}×' for v in crust_rates])

    for term, marker in marker_map.items():
        axes[0].scatter([], [], marker=marker, color=TERM_COLORS[term],
                        s=50, label=TERM_LABELS[term])
    axes[0].legend(fontsize=7, loc='upper left')

    fig.suptitle('Zero Outgassing, 10% Carbonate Crust (ocean depth 3000 m)',
                 fontsize=13, fontweight='bold')
    path = os.path.join(output_path, 'lines_out0_carb01.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {path}")


def plot_ocean_depth_effect(df, output_path):
    """T, P_CO2, pH, salinity vs instellation for Earth-like tectonics, coloured by ocean depth."""
    subset = df[
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0) &
        df['reverse_weathering'] &
        (df['comp_name'] == 'basalt_49') &
        (df['crust_carbonate'] == 0.0) &
        (df['f_HT'] == 0.0)
    ]

    if subset.empty:
        print("No data for ocean depth sweep — skipping.")
        return

    depths = sorted(subset['ocean_depth'].unique())

    if len(depths) > 1 and (max(depths) / max(1, min(depths))) >= 10:
        norm = mcolors.LogNorm(vmin=min(depths), vmax=max(depths))
    else:
        norm = mcolors.Normalize(vmin=min(depths), vmax=max(depths))

    cmap = plt.get_cmap('viridis')
    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    fig, axes = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
    fig.subplots_adjust(hspace=0.05, right=0.84)

    for d in depths:
        group = subset[subset['ocean_depth'] == d].sort_values('instellation')
        if group.empty:
            continue
        _plot_group_on_axes(axes, group, cmap(norm(d)), marker_map)

    min_x = subset['instellation'].min()
    max_x = subset['instellation'].max()
    margin = (max_x - min_x) * 0.05 if not pd.isna(min_x) and max_x != min_x else 0.1
    x_lims = (min_x - margin, max_x + margin) if not pd.isna(min_x) else (0.25, 1.45)

    _style_4panel_axes(axes, x_lims=x_lims)
    axes[2].set_ylim(4.5, 9.5)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.86, 0.08, 0.03, 0.84])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Ocean Depth (m)')

    if len(depths) <= 10:
        cbar.set_ticks(depths)
        cbar.set_ticklabels([f'{v:g}' for v in depths])

    for term, marker in marker_map.items():
        axes[0].scatter([], [], marker=marker, color=TERM_COLORS[term],
                        s=50, label=TERM_LABELS[term])
    axes[0].legend(fontsize=7, loc='upper left')

    fig.suptitle('Effect of Ocean Depth (Earth Tectonics Baseline)',
                 fontsize=13, fontweight='bold')
    path = os.path.join(output_path, 'lines_ocean_depth.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {path}")


def plot_crust_composition(df, output_path):
    """Sweep 3: T, P_CO2, pH, salinity vs instellation for each crust composition (Earth-like baseline)."""
    compositions = ['komatiite_42', 'komatiite_44', 'basalt_47', 'basalt_49', 'basalt_51']

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

    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}
    cols = ['T', 'P_CO2', 'pH', 'salinity']

    fig, axes = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
    fig.subplots_adjust(hspace=0.05, right=0.78)

    for comp in compositions:
        group = subset[subset['comp_name'] == comp].sort_values('instellation')
        if group.empty:
            continue
        color = COMP_COLORS.get(comp, 'k')
        label = COMP_LABELS.get(comp, comp)

        hab = group[group['termination'].isin(HABITABLE)]
        for ax, col in zip(axes, cols):
            if not hab.empty:
                line, = ax.plot(hab['instellation'], hab[col], color=color, linewidth=1.8, zorder=3)
                ax.scatter(hab['instellation'], hab[col], color=color, s=28, zorder=4)
                if ax is axes[0]:
                    line.set_label(label)

        for _, row in group[~group['termination'].isin(HABITABLE)].iterrows():
            marker = marker_map.get(row['termination'], 'x')
            val = row['T']
            if np.isfinite(val):
                axes[0].scatter(row['instellation'], val, marker=marker,
                                s=55, color=color, zorder=4, linewidths=1.2)

    _style_4panel_axes(axes)
    axes[0].legend(fontsize=7, loc='upper left', title='Composition')

    fig.suptitle('Crust Composition Effect\n(out=1×, crust=1×, depth=3000 m, rw=True)',
                 fontsize=12, fontweight='bold')
    path = os.path.join(output_path, 'lines_crust_composition.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {path}")


def plot_magnesium_chemistry(df, output_path):
    """Sweep 4: Compare mg chemistry (basalt_49_no_iron + rw=True) vs no_mg (Ca-only + rw=False)."""
    combined = df[df['comp_name'].isin(['mg', 'no_mg'])]

    if combined.empty:
        print("No magnesium chemistry sweep data found — skipping.")
        return

    crust_rates     = sorted(combined['crust_production'].unique())
    outgassing_vals = sorted(combined['outgassing'].unique())

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = cmr.tropical
    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    for c in crust_rates:
        fig, axes = plt.subplots(4, 1, figsize=(7, 12), sharex=True)
        fig.subplots_adjust(hspace=0.05, right=0.84)

        for o in outgassing_vals:
            color = cmap(norm(o))

            mg_group = combined[
                (combined['comp_name'] == 'mg') &
                (combined['crust_production'] == c) &
                (combined['outgassing'] == o)
            ].sort_values('instellation')

            nomg_group = combined[
                (combined['comp_name'] == 'no_mg') &
                (combined['crust_production'] == c) &
                (combined['outgassing'] == o)
            ].sort_values('instellation')

            for group, ls in [(mg_group, '-'), (nomg_group, '--')]:
                if not group.empty:
                    _plot_group_on_axes(axes, group, color, marker_map, linestyle=ls)

        _style_4panel_axes(axes)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ax = fig.add_axes([0.86, 0.08, 0.03, 0.84])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Outgassing factor')
        cbar.set_ticks(outgassing_vals)
        cbar.set_ticklabels([f'{v}×' for v in outgassing_vals])

        legend_handles = [
            Line2D([0], [0], color='k', linestyle='-',  linewidth=1.8, label='With Mg chemistry (rw=True)'),
            Line2D([0], [0], color='k', linestyle='--', linewidth=1.8, label='No Mg (Ca-only, rw=False)'),
        ]
        for term, marker in marker_map.items():
            legend_handles.append(
                plt.scatter([], [], marker=marker, color=TERM_COLORS[term], s=50, label=TERM_LABELS[term])
            )
        axes[0].legend(handles=legend_handles, fontsize=7, loc='upper left')

        fig.suptitle(f'Magnesium Chemistry Comparison — Crust = {c}× Earth',
                     fontsize=12, fontweight='bold')
        path = os.path.join(output_path, f'lines_mg_crust_{c}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {path}")

    # Combined plot: crust rates as columns
    n_cols = len(crust_rates)
    fig_c, axes_c = plt.subplots(4, n_cols, figsize=(6 * n_cols + 1.5, 12),
                                  sharex=True, sharey='row')
    fig_c.subplots_adjust(hspace=0.05, wspace=0.05, right=0.88)

    for ci, c in enumerate(crust_rates):
        col_axes = axes_c[:, ci]

        for o in outgassing_vals:
            color = cmap(norm(o))
            mg_group = combined[
                (combined['comp_name'] == 'mg') &
                (combined['crust_production'] == c) &
                (combined['outgassing'] == o)
            ].sort_values('instellation')
            nomg_group = combined[
                (combined['comp_name'] == 'no_mg') &
                (combined['crust_production'] == c) &
                (combined['outgassing'] == o)
            ].sort_values('instellation')
            for group, ls in [(mg_group, '-'), (nomg_group, ':')]:
                if not group.empty:
                    _plot_group_on_axes(col_axes, group, color, marker_map, linestyle=ls)

        _style_4panel_axes(col_axes)
        axes_c[0, ci].set_title(f'{c}×', fontsize=10)

        if ci > 0:
            for ax in col_axes:
                ax.set_ylabel('')
            for row in range(4):
                plt.setp(axes_c[row, ci].get_yticklabels(), visible=False)

        if ci != n_cols // 2:
            col_axes[3].set_xlabel('')

    legend_handles = [
        Line2D([0], [0], color='k', linestyle='-',  linewidth=1.8, label='With Mg chemistry (rw=True)'),
        Line2D([0], [0], color='k', linestyle=':', linewidth=1.8, label='No Mg (Ca-only, rw=False)'),
    ]
    for term, marker in marker_map.items():
        legend_handles.append(
            plt.scatter([], [], marker=marker, color=TERM_COLORS[term], s=50, label=TERM_LABELS[term])
        )
    axes_c[0, 0].legend(handles=legend_handles, fontsize=6, loc='upper left')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig_c.add_axes([0.90, 0.08, 0.02, 0.84])
    cbar = fig_c.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Outgassing factor')
    cbar.set_ticks(outgassing_vals)
    cbar.set_ticklabels([f'{v}×' for v in outgassing_vals])

    fig_c.suptitle('Magnesium Chemistry Comparison — Crust production rate →',
                   fontsize=12, fontweight='bold')
    path = os.path.join(output_path, 'lines_mg_combined.png')
    fig_c.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig_c)
    print(f"Saved {path}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot kamino parameter sweep results.')
    parser.add_argument('--path', default=DEFAULT_OUTPUT_PATH,
                        help='Directory containing planet_*.json files.')
    args = parser.parse_args()

    df = load_data(args.path)
    if df.empty:
        print("No data found. Check --path.")
    else:
        plot_termination_map(df, args.path)
        plot_faceted_lines(df, args.path)
        plot_heatmaps(df, args.path)
        plot_zero_outgassing_carb01(df, args.path)
        plot_ocean_depth_effect(df, args.path)
        plot_crust_composition(df, args.path)
        plot_magnesium_chemistry(df, args.path)
        print("Done.")
