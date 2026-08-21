"""Ocean-only instellation sweep.

Sweeps instellation × outgassing with no continents (land_fraction=0).
Fixed: basalt_49, rw=True, depth=3000 m, crust=1×.
f_HT and tau_rw use Earth-calibrated values (iter-6 abiotic calibration).
"""

import itertools
import multiprocessing as mp
import os
import json
import glob
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cmasher as cmr

import kamino.planet as p_mod
from kamino.planet import Planet
from kamino.constants import M_EARTH, R_EARTH
from kamino.mineral_info import basalt_49

_STYLE = os.path.join(os.path.dirname(__file__), 'planetary-chem-paper.mplstyle')
if os.path.exists(_STYLE):
    plt.style.use(_STYLE)

OUTPUT_PATH = '/data/pt426/kamino_experiments/ocean_sweep/'

INSTELLATION = [
    0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
    0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4,
]
OUTGASSING = [0.3, 1.0, 3.0]
F_HT   = 0.039   # Earth-calibrated (iter-6 abiotic)
TAU_RW = 5e6     # Myr — Earth-calibrated reverse weathering timescale

# (b_ocean index, label, display unit, Earth reference value or None)
# b_ocean: 0=Alk, 1=DIC, 2=Si, 3=Al, 4=Fe, 5=Ca, 6=Mg, 7=Na, 8=Cl, 9=S
ION_SPEC = [
    (0, 'Alk',  'meq/kg', 2.3e-3),
    (1, 'DIC',  'mmol/kg', 2.0e-3),
    (2, 'Si',   'mmol/kg', 0.1e-3),
    (5, 'Ca',   'mmol/kg', 10.3e-3),
    (6, 'Mg',   'mmol/kg', 52.8e-3),
    (7, 'Na',   'mmol/kg', 480e-3),
    (8, 'Cl',   'mmol/kg', 550e-3),
]

# b_ocean index → molar mass (g/mol) for salinity proxy
_SAL_SPEC = [(1, 61.0), (2, 60.1), (3, 27.0), (4, 55.8),
             (5, 40.1), (6, 24.3), (7, 23.0), (8, 35.45)]

HABITABLE = {'converged', 'timeout', 'co2_floor'}
TERM_LABELS = {
    'converged':     'Converged',
    'timeout':       'Timeout (2 Gyr)',
    'out_of_domain': 'Outside model domain',
    # Legacy terminations, kept so pre-domain-event runs still plot.
    'co2_floor':     'CO₂ floor (legacy)',
    'snowball':      'Snowball (legacy)',
    'hothouse':      'Hothouse (legacy)',
    'co2_ceiling':   'Outside model domain (legacy)',
    'acid_ocean':    'Outside model domain (legacy)',
}
HAB_MARKERS    = {'converged': 'o', 'timeout': 's', 'co2_floor': 'P'}
FAILED_MARKERS = {'out_of_domain': 's', 'snowball': 'v', 'hothouse': '^',
                  'co2_ceiling': 's', 'acid_ocean': 's'}


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

def run_simulation(s, o, output_path):
    p_mod.output_path = output_path
    run_name = f'ocean_s_{s}_out_{o}'
    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=1.0,
            outgassing=o,
            ocean_depth=3000,
            land_fraction=0.0,
            crust_composition=basalt_49,
            reverse_weathering=True,
            f_HT=F_HT,
            tau_rw=TAU_RW * p_mod.YR,
            name=run_name,
        )
        p.time_evolve()
        return run_name, None
    except Exception as e:
        return run_name, str(e)


def run_sweep():
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    p_mod.output_path = OUTPUT_PATH

    combos = list(itertools.product(INSTELLATION, OUTGASSING))
    total  = len(combos)
    workers = min(total, 26)
    print(f"Running {total} ocean-only simulations ({workers} workers)...")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {executor.submit(run_simulation, s, o, OUTPUT_PATH): (s, o)
                   for s, o in combos}
        done = 0
        for future in as_completed(futures):
            done += 1
            run_name, error = future.result()
            if error:
                print(f"[{done}/{total}] FAILED {run_name}: {error}")
            else:
                print(f"[{done}/{total}] Done: {run_name}")

    print("All simulations complete.")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_results():
    files = sorted(glob.glob(os.path.join(OUTPUT_PATH, 'ocean_s_*.json')))
    rows = []
    for f in files:
        with open(f) as fh:
            d = json.load(fh)
        if 'termination' not in d:
            continue
        y = d.get('data', {}).get('y', [])
        if not y:
            continue
        n_el = len(y) - 3   # y = [P_CO2, P_H2O, *b_ocean, r_avg]
        b = [max(float(y[2 + i][-1]), 1e-15) for i in range(n_el)]
        while len(b) < 10:
            b.append(1e-15)

        salinity = sum(b[idx] * mass for idx, mass in _SAL_SPEC if idx < len(b))

        rows.append({
            'instellation': float(d['instellation']),
            'outgassing':   float(d['outgassing']),
            'termination':  d['termination'],
            'T':            float(d.get('T', np.nan)),
            'P_CO2':        float(d.get('P_CO2', np.nan)),
            'pH':           float(d.get('pH', np.nan)),
            'salinity':     salinity,
            'b':            b[:10],
        })

    rows.sort(key=lambda r: (r['outgassing'], r['instellation']))
    print(f"Loaded {len(rows)} results.")
    return rows


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _add_markers(ax, group, y_col, color):
    hab   = [r for r in group if r['termination'] in HABITABLE]
    nhab  = [r for r in group if r['termination'] not in HABITABLE]
    for r in hab:
        m = HAB_MARKERS.get(r['termination'], 'o')
        ax.scatter(r['instellation'], r[y_col], marker=m, s=28, color=color, zorder=4)
    for r in nhab:
        m = FAILED_MARKERS.get(r['termination'], 'x')
        if np.isfinite(r[y_col]):
            ax.scatter(r['instellation'], r[y_col], marker=m, s=50,
                       facecolors='none', edgecolors=color, linewidths=1.4, zorder=4)


def _add_colorbar(fig, axes, cmap, norm, out_vals):
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, location='right', pad=0.02, aspect=30)
    cbar.set_label('Outgassing (× Earth)')
    cbar.set_ticks(out_vals)
    cbar.set_ticklabels([f'{v:g}×' for v in out_vals])
    return cbar


def _save(fig, name):
    path = os.path.join(OUTPUT_PATH, name)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {path}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_climate_summary(rows):
    out_vals = sorted(set(r['outgassing'] for r in rows))
    norm  = mcolors.LogNorm(vmin=min(out_vals), vmax=max(out_vals))
    cmap  = cmr.tropical

    fig, axes = plt.subplots(4, 1, figsize=(7, 10), sharex=True)
    ax_t, ax_pco2, ax_ph, ax_sal = axes

    for o in out_vals:
        color = cmap(norm(o))
        group = [r for r in rows if r['outgassing'] == o]

        for ax, col in [(ax_t, 'T'), (ax_pco2, 'P_CO2'), (ax_ph, 'pH'), (ax_sal, 'salinity')]:
            xs = [r['instellation'] for r in group]
            ys = [r[col] for r in group]
            ax.plot(xs, ys, color=color, linewidth=1.8, alpha=0.85)
            _add_markers(ax, group, col, color)

    ax_t.set_ylabel('Temperature (K)')
    ax_t.axhspan(235, 260, color='blue', alpha=0.12, label='Snowball')
    ax_t.axhspan(340, 360, color='red',  alpha=0.12, label='Hothouse')
    ax_t.set_ylim(235, 370)

    ax_pco2.set_ylabel('$P_{\\mathrm{CO_2}}$ (bar)')
    ax_pco2.set_yscale('log')
    ax_pco2.set_ylim(1e-8, 20)

    ax_ph.set_ylabel('Ocean pH')
    ax_ph.set_ylim(4.5, 12)

    ax_sal.set_ylabel('Dissolved ions (g/kg)')
    ax_sal.set_yscale('log')

    for ax in axes:
        ax.set_xlim(0.25, 1.45)
        ax.grid(True, linestyle='--', alpha=0.4)
    axes[-1].set_xlabel('Instellation ($S/S_0$)')

    _add_colorbar(fig, list(axes), cmap, norm, out_vals)
    _save(fig, 'ocean_climate_summary.png')


def plot_ion_concentrations(rows):
    out_vals = sorted(set(r['outgassing'] for r in rows))
    norm  = mcolors.LogNorm(vmin=min(out_vals), vmax=max(out_vals))
    cmap  = cmr.tropical

    n = len(ION_SPEC)
    fig, axes = plt.subplots(n, 1, figsize=(7, n * 2.2), sharex=True)

    for o in out_vals:
        color = cmap(norm(o))
        group = [r for r in rows if r['outgassing'] == o]
        hab   = [r for r in group if r['termination'] in HABITABLE]
        if not hab:
            continue
        s_arr = np.array([r['instellation'] for r in hab])
        b_arr = np.array([r['b'] for r in hab])

        for ax, (idx, label, unit, _) in zip(axes, ION_SPEC):
            ax.plot(s_arr, b_arr[:, idx] * 1e3, color=color, linewidth=1.8, alpha=0.85)

    for ax, (idx, label, unit, earth_ref) in zip(axes, ION_SPEC):
        ax.set_ylabel(f'{label} ({unit})')
        ax.set_yscale('log')
        ax.set_xlim(0.25, 1.45)
        ax.grid(True, linestyle='--', alpha=0.4)
        if earth_ref is not None:
            ax.axhline(earth_ref * 1e3, color='k', linestyle=':', linewidth=0.9,
                       alpha=0.6, label='Earth')

    axes[-1].set_xlabel('Instellation ($S/S_0$)')
    axes[0].legend(fontsize=7, loc='upper right')

    _add_colorbar(fig, list(axes), cmap, norm, out_vals)
    _save(fig, 'ocean_ion_concentrations.png')


def plot_ions_combined(rows):
    """All ions on one panel, outgassing=1× only, with Earth reference stars."""
    o_ref = min(set(r['outgassing'] for r in rows), key=lambda x: abs(x - 1.0))
    hab   = [r for r in rows if r['outgassing'] == o_ref and r['termination'] in HABITABLE]
    if not hab:
        print("No habitable results for combined ion plot — skipping.")
        return

    hab.sort(key=lambda r: r['instellation'])
    s_arr = np.array([r['instellation'] for r in hab])
    b_arr = np.array([r['b'] for r in hab])

    ion_colors = [plt.cm.tab10(k / 10) for k in range(len(ION_SPEC))]

    fig, ax = plt.subplots(figsize=(7, 4))
    for (idx, label, unit, earth_ref), color in zip(ION_SPEC, ion_colors):
        ax.plot(s_arr, b_arr[:, idx] * 1e3, color=color, linewidth=1.8, label=label)
        if earth_ref is not None:
            ax.scatter(1.0, earth_ref * 1e3, marker='*', s=180, color=color,
                       edgecolors='k', linewidths=0.5, zorder=6)

    ax.set_ylabel('Concentration (mmol/kg or meq/kg)')
    ax.set_yscale('log')
    ax.set_xlabel('Instellation ($S/S_0$)')
    ax.set_xlim(0.25, 1.45)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(ncols=2, fontsize=8, loc='lower left', handlelength=1.2)
    ax.set_title(f'Ocean-only, outgassing={o_ref:g}× — ★ = Earth ocean')

    _save(fig, 'ocean_ions_combined.png')


def plot_salinity_breakdown(rows):
    """Stacked-area chart showing how each ion contributes to total salinity vs instellation."""
    o_ref = min(set(r['outgassing'] for r in rows), key=lambda x: abs(x - 1.0))
    hab   = sorted(
        [r for r in rows if r['outgassing'] == o_ref and r['termination'] in HABITABLE],
        key=lambda r: r['instellation'],
    )
    if not hab:
        return

    s_arr = np.array([r['instellation'] for r in hab])
    b_arr = np.array([r['b'] for r in hab])

    # Contributions to salinity (g/kg)
    labels_sal = ['DIC', 'Si', 'Al', 'Fe', 'Ca', 'Mg', 'Na', 'Cl']
    contribs = np.column_stack([b_arr[:, idx] * mass for idx, mass in _SAL_SPEC])
    colors_sal = [plt.cm.tab10(k / 10) for k in range(len(labels_sal))]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.stackplot(s_arr, contribs.T, labels=labels_sal, colors=colors_sal, alpha=0.85)
    ax.set_ylabel('Dissolved ions (g/kg)')
    ax.set_xlabel('Instellation ($S/S_0$)')
    ax.set_xlim(0.25, 1.45)
    ax.set_ylim(bottom=0)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(ncols=4, fontsize=7, loc='upper left', handlelength=1.0)
    ax.set_title(f'Salinity composition — ocean-only, outgassing={o_ref:g}×')

    _save(fig, 'ocean_salinity_breakdown.png')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Ocean-only instellation sweep.')
    parser.add_argument('--plot-only', action='store_true',
                        help='Skip simulations and replot from existing results.')
    args = parser.parse_args()

    if not args.plot_only:
        run_sweep()

    rows = load_results()
    if rows:
        plot_climate_summary(rows)
        plot_ion_concentrations(rows)
        plot_ions_combined(rows)
        plot_salinity_breakdown(rows)
    else:
        print("No results found. Run without --plot-only first.")
