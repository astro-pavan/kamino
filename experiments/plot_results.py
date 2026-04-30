import os
import glob
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

DEFAULT_OUTPUT_PATH = '/home/pavan/PhD/kamino_experiments/'

TERM_COLORS = {
    'snowball':   '#5b9bd5',
    'hothouse':   '#e05c5c',
    'acid_ocean': '#e8935e',
    'converged':  '#6abf69',
    'timeout':    '#b8b8b8',
}
TERM_LABELS = {
    'snowball':   'Snowball',
    'hothouse':   'Hothouse',
    'acid_ocean': 'Acid Ocean (>5 bar CO₂)',
    'converged':  'Converged',
    'timeout':    'Timeout (2 Gyr)',
}
HABITABLE = {'converged', 'timeout'}
T_SNOWBALL = 260.0
T_RUNAWAY  = 350.0


def load_data(output_path):
    files = sorted(glob.glob(os.path.join(output_path, 'planet_*.json')))
    rows = []
    for f in files:
        with open(f) as fh:
            d = json.load(fh)
        if 'termination' not in d:
            print(f"  Skipping (no termination): {os.path.basename(f)}")
            continue
        rows.append({
            'name':             d.get('name', ''),
            'instellation':     float(d['instellation']),
            'outgassing':       float(d['outgassing']),
            'crust_production': float(d['crust_production_rate']),
            'reverse_weathering': bool(d.get('reverse_weathering', False)),
            'crust_carbonate':  float(d.get('crust_carbonate_content', 0.0)),
            'ocean_depth':      float(d["ocean_depth"]),
            'termination':      d['termination'],
            'end_time_yr':      d.get('end_time_yr', np.nan),
            'T':                d.get('T', np.nan),
            'P_CO2':            d.get('P_CO2', np.nan),
            'pH':               d.get('pH', np.nan),
        })
    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} simulations.")
    return df


def _base(df):
    return df[~df['reverse_weathering'] & (df['crust_carbonate'] == 0.0) & (df['ocean_depth'] == 3000) & (df['outgassing'] > 0)]


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
    """T, P_CO2, pH vs instellation per crust rate, coloured by outgassing."""
    base = _base(df)
    crust_rates     = sorted(base['crust_production'].unique())
    outgassing_vals = sorted(base['outgassing'].unique())

    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = plt.get_cmap('plasma')

    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    for c in crust_rates:
        subset_c = base[base['crust_production'] == c]
        fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
        fig.subplots_adjust(hspace=0.05, right=0.84)

        for o in outgassing_vals:
            group = subset_c[subset_c['outgassing'] == o].sort_values('instellation')
            if group.empty:
                continue
            color = cmap(norm(o))

            hab = group[group['termination'].isin(HABITABLE)]
            for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
                if not hab.empty:
                    ax.plot(hab['instellation'], hab[col],
                            color=color, linewidth=1.8, zorder=3)
                    ax.scatter(hab['instellation'], hab[col],
                               color=color, s=28, zorder=4)

            for _, row in group[~group['termination'].isin(HABITABLE)].iterrows():
                marker = marker_map.get(row['termination'], 'x')
                ec = TERM_COLORS.get(row['termination'], 'k')
                for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
                    val = row[col]
                    if np.isfinite(val):
                        ax.scatter(row['instellation'], val, marker=marker,
                                   s=55, color=ec, zorder=4, linewidths=1.2)

        axes[0].set_ylabel('Temperature (K)')
        axes[0].axhspan(T_SNOWBALL - 20, T_SNOWBALL, color=TERM_COLORS['snowball'], alpha=0.12)
        axes[0].axhspan(T_RUNAWAY - 10,  T_RUNAWAY,  color=TERM_COLORS['hothouse'],  alpha=0.12)
        axes[0].set_ylim(235, 360)

        axes[1].set_ylabel('P_CO₂ (bar)')
        axes[1].set_yscale('log')
        axes[1].set_ylim(1e-8, 20)

        axes[2].set_ylabel('Ocean pH')
        axes[2].set_xlabel('Instellation (S/S₀)')
        axes[2].set_ylim(4.5, 9.5)

        for ax in axes:
            ax.grid(True, linestyle='--', alpha=0.4)
            ax.set_xlim(0.2, 1.4)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ax = fig.add_axes([0.86, 0.12, 0.03, 0.76])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Outgassing factor')
        cbar.set_ticks(outgassing_vals)
        cbar.set_ticklabels([f'{v}×' for v in outgassing_vals])

        # Marker legend for non-habitable outcomes
        for term, marker in marker_map.items():
            axes[0].scatter([], [], marker=marker, color=TERM_COLORS[term],
                            s=50, label=TERM_LABELS[term])
        axes[0].legend(fontsize=7, loc='upper left')

        fig.suptitle(f'Crust production = {c}× Earth', fontsize=13, fontweight='bold')
        path = os.path.join(output_path, f'lines_crust_{c}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {path}")


def plot_heatmaps(df, output_path):
    """T, P_CO2, pH as heatmaps over (crust_production, outgassing) at each instellation."""
    base = _base(df)
    instellations = sorted(base['instellation'].unique())
    outgassings   = sorted(base['outgassing'].unique())
    crust_rates   = sorted(base['crust_production'].unique())

    nc = len(crust_rates)
    no = len(outgassings)
    x_edges = np.arange(-0.5, nc)
    y_edges = np.arange(-0.5, no)

    vars_cfg = [
        ('T',     'Temperature (K)',   'RdBu_r',  mcolors.Normalize(vmin=270, vmax=350)),
        ('P_CO2', 'P_CO₂ (bar)',       'inferno', mcolors.LogNorm(vmin=1e-8, vmax=10)),
        ('pH',    'Ocean pH',           'cividis', mcolors.Normalize(vmin=5, vmax=9)),
    ]

    for var, label, cmap_name, norm in vars_cfg:
        fig, axes = plt.subplots(1, len(instellations),
                                 figsize=(len(instellations) * 3.0, 3.6), sharey=True)

        mesh = None
        for ax, s in zip(axes, instellations):
            subset = base[base['instellation'] == s]
            Z    = np.full((no, nc), np.nan)
            termZ = np.full((no, nc), '', dtype=object)

            for _, row in subset.iterrows():
                xi = crust_rates.index(row['crust_production'])
                yi = outgassings.index(row['outgassing'])
                termZ[yi, xi] = row['termination']
                if row['termination'] in HABITABLE:
                    Z[yi, xi] = row[var]

            mesh = ax.pcolormesh(x_edges, y_edges, Z,
                                 cmap=cmap_name, norm=norm, zorder=1)

            # Overlay non-habitable cells
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
    """T, P_CO2, pH vs instellation for outgassing=0, crust_carbonate=0.1, coloured by crust rate."""
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

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    fig.subplots_adjust(hspace=0.05, right=0.84)

    for c in crust_rates:
        group = subset[subset['crust_production'] == c].sort_values('instellation')
        if group.empty:
            continue
        color = cmap(norm(c))

        hab = group[group['termination'].isin(HABITABLE)]
        for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
            if not hab.empty:
                ax.plot(hab['instellation'], hab[col], color=color, linewidth=1.8, zorder=3)
                ax.scatter(hab['instellation'], hab[col], color=color, s=28, zorder=4)

        for _, row in group[~group['termination'].isin(HABITABLE)].iterrows():
            marker = marker_map.get(row['termination'], 'x')
            ec = TERM_COLORS.get(row['termination'], 'k')
            for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
                val = row[col]
                if np.isfinite(val):
                    ax.scatter(row['instellation'], val, marker=marker,
                               s=55, color=ec, zorder=4, linewidths=1.2)

    axes[0].set_ylabel('Temperature (K)')
    axes[0].axhspan(T_SNOWBALL - 20, T_SNOWBALL, color=TERM_COLORS['snowball'], alpha=0.12)
    axes[0].axhspan(T_RUNAWAY - 10,  T_RUNAWAY,  color=TERM_COLORS['hothouse'],  alpha=0.12)
    axes[0].set_ylim(235, 360)

    axes[1].set_ylabel('P_CO₂ (bar)')
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-8, 20)

    axes[2].set_ylabel('Ocean pH')
    axes[2].set_xlabel('Instellation (S/S₀)')
    axes[2].set_ylim(4.5, 9.5)

    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(0.2, 1.5)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.86, 0.12, 0.03, 0.76])
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
    """T, P_CO2, pH vs instellation for Earth-like tectonics, coloured by ocean depth."""
    # Filter for standard outgassing and crust production to isolate the ocean depth effect
    subset = df[
        (df['outgassing'] == 1.0) &
        (df['crust_production'] == 1.0) &
        (~df['reverse_weathering']) &
        (df['crust_carbonate'] == 0.0)
    ]
    
    if subset.empty:
        print("No standard Earth-like runs found to plot the ocean depth effect — skipping.")
        return

    depths = sorted(subset['ocean_depth'].unique())
    
    # Use LogNorm if there's a wide range of depths, otherwise use linear Normalize
    if len(depths) > 1 and (max(depths) / max(1, min(depths))) >= 10:
        norm = mcolors.LogNorm(vmin=min(depths), vmax=max(depths))
    else:
        norm = mcolors.Normalize(vmin=min(depths), vmax=max(depths))
        
    cmap = plt.get_cmap('viridis')
    marker_map = {'snowball': 'v', 'hothouse': '^', 'acid_ocean': 's'}

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    fig.subplots_adjust(hspace=0.05, right=0.84)

    for d in depths:
        group = subset[subset['ocean_depth'] == d].sort_values('instellation')
        if group.empty:
            continue
        color = cmap(norm(d))

        # Plot habitable continuous lines
        hab = group[group['termination'].isin(HABITABLE)]
        for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
            if not hab.empty:
                ax.plot(hab['instellation'], hab[col], color=color, linewidth=1.8, zorder=3)
                ax.scatter(hab['instellation'], hab[col], color=color, s=28, zorder=4)

        # Plot non-habitable specific markers
        for _, row in group[~group['termination'].isin(HABITABLE)].iterrows():
            marker = marker_map.get(row['termination'], 'x')
            ec = TERM_COLORS.get(row['termination'], 'k')
            for ax, col in zip(axes, ['T', 'P_CO2', 'pH']):
                val = row[col]
                if np.isfinite(val):
                    ax.scatter(row['instellation'], val, marker=marker,
                               s=55, color=ec, zorder=4, linewidths=1.2)

    # Styling Axes (Mirroring existing structure)
    axes[0].set_ylabel('Temperature (K)')
    axes[0].axhspan(T_SNOWBALL - 20, T_SNOWBALL, color=TERM_COLORS['snowball'], alpha=0.12)
    axes[0].axhspan(T_RUNAWAY - 10,  T_RUNAWAY,  color=TERM_COLORS['hothouse'],  alpha=0.12)
    axes[0].set_ylim(235, 360)

    axes[1].set_ylabel('P_CO₂ (bar)')
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-8, 20)

    axes[2].set_ylabel('Ocean pH')
    axes[2].set_xlabel('Instellation (S/S₀)')
    axes[2].set_ylim(4.5, 9.5)

    # Dynamic X-axis depending on the subset range
    min_x = subset['instellation'].min()
    max_x = subset['instellation'].max()
    margin = (max_x - min_x) * 0.05 if not pd.isna(min_x) and max_x != min_x else 0.1
    x_lims = (min_x - margin, max_x + margin) if not pd.isna(min_x) else (0.2, 1.5)

    for ax in axes:
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(*x_lims)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.86, 0.12, 0.03, 0.76])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Ocean Depth (m)')
    
    # If there are only a few depths being swept, show them explicitly on the colorbar
    if len(depths) <= 10:
        cbar.set_ticks(depths)
        cbar.set_ticklabels([f'{v:g}' for v in depths])

    # Legend for termination markers
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
        print("Done.")
