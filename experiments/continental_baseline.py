"""The Earth-like continental baseline: an instellation sweep at land_fraction = 0.3.

Every sweep in `parameter_sweep.py` runs at `land_fraction = 0` -- the land-free ocean worlds the
paper is about. This script runs one instellation line at `land_fraction = 0.3`, Earth's, with
everything else held at Earth values, so the model's habitable zone can be quoted for an
Earth-like planet and compared like for like against the ocean worlds.

With land, `Planet.dY_dt` switches on two extra terms: the kinetic Walker-type continental
silicate flux `F_cont` and shelf carbonate burial `F_shelf_prec`. Nothing else about the planet
changes -- `LAND_ALBEDO == OCEAN_ALBEDO` in planet.py, so the two arms receive identical
instellation and differ only in their chemistry.

The land-free arm is part of the sweep rather than assumed present: its run names are IDENTICAL
to the ones `sweep_basic` already wrote (the `_land` tag is suppressed at 0, see
`parameter_sweep._run_name`), so those runs are reused for free and only the missing ones -- the
extension past S = 1.2 -- cost anything.

    # run the sweep, then draw the figures (resumes; runs already on disk are reused)
    /data/pt426/big-venv/bin/python experiments/continental_baseline.py

    # re-draw the figures from runs already on disk
    /data/pt426/big-venv/bin/python experiments/continental_baseline.py --plot-only
"""

import argparse
import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import numpy as np

import parameter_sweep as ps
from parameter_sweep import (ALPHA_CALIB, KD_MG_CALIB, K_NA_CALIB, OUTPUT_PATH, WORKERS,
                             PE_REDUCING, PE_OXIDISING, _pe_label, _run_name, run_simulation)
from kamino.constants import EARTH_MANTLE_MG_SI, EARTH_DELTA_IW

# ── The baseline planet: Earth, on every axis ─────────────────────────────────────────────────
# Earth's land fraction. `Planet` scales the continental flux linearly off this
# (`_s_terr = _S_TERR_EARTH * land_fraction / 0.3`), so 0.3 is the value the terrestrial
# denudation rate was calibrated at, not an arbitrary point on a continuum.
LAND_FRACTION = 0.3

# Both arms of the comparison. 0.0 is the ocean world every other sweep runs; it is listed second
# so the continental runs -- the ones that do not exist yet -- are submitted first.
LAND_ARMS = [LAND_FRACTION, 0.0]

# These are ints on purpose. `_run_name` interpolates them with plain str(), so 1 and 1.0 give
# 'crust_1' and 'crust_1.0' -- two names for one config, and the ocean arm would stop matching
# the runs `sweep_basic` already wrote. Match the types parameter_sweep uses.
OUTGASSING = 1                # x Earth
CRUST_PRODUCTION = 1          # x Earth
OCEAN_DEPTH = 3000            # m

REVERSE_WEATHERING = True
MG_SI = float(EARTH_MANTLE_MG_SI)      # mantle molar Mg/Si, Earth's 1.25
DELTA_IW = float(EARTH_DELTA_IW)       # core-formation oxygen fugacity, Earth's -2

# The reducing arm only. Every parameter_sweep sweep runs both redox states because the model has
# no basis for preferring one; here the abiotic (reducing) state is the model's own default and
# the one every figure is drawn at, so the oxidising arm would double the cost of a sweep nothing
# plots. Pass --both-redox to run it anyway.
PE_STATES = [PE_REDUCING]

# Matches parameter_sweep's grid out to 1.2, so the land-free arm is already on disk, then
# extends to 1.45. The extension is the point: continental weathering is a far stronger CO2 sink
# than seafloor weathering alone, so the Earth-like arm is expected to stay temperate past the
# instellation at which the ocean worlds run away -- and an inner edge outside the swept range is
# a bound, not a measurement.
INSTELLATION = list(ps.instellation) + [1.25, 1.3, 1.35, 1.4, 1.45]


def _combos(instellation=None, land_arms=None, pe_states=None):
    """Combos in `run_simulation` argument order, with land fraction last.

    Ordered through `parameter_sweep._cost_rank`, and within a cost tier continental runs before
    ocean ones, since the ocean arm is mostly already on disk.

    Treat that ordering as a rough guide rather than cheap-first. `_cost_rank` calls S > 1.12
    free because on a land-free world those runs leave the model domain within seconds -- which
    is exactly what continental weathering prevents. Measured on this grid, S = 1.15 stays
    temperate and integrates the full 2 Gyr, so the runs the heuristic puts first are among the
    EXPENSIVE ones on the land arm.
    """
    combos = [
        (s, OUTGASSING, CRUST_PRODUCTION, OCEAN_DEPTH, REVERSE_WEATHERING, MG_SI, DELTA_IW,
         ALPHA_CALIB, KD_MG_CALIB, K_NA_CALIB, pe, land)
        for s, land, pe in itertools.product(instellation or INSTELLATION,
                                             land_arms or LAND_ARMS,
                                             pe_states or PE_STATES)
    ]
    combos.sort(key=lambda c: (ps._cost_rank(c), c[11] == 0.0))
    return combos


def run(combos, output_path=OUTPUT_PATH):
    """Execute a combo list. Mirrors `parameter_sweep.run_combos`, carrying land fraction through."""
    if not output_path.endswith('/'):
        output_path += '/'
    ps.p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    # Distinct configs must map to distinct filenames or one silently overwrites the other, and
    # the resume path then hands back the survivor's result for both (the fast_13 resume trap).
    names = [_run_name(*combo) for combo in combos]
    if len(set(names)) != len(names):
        duplicated = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(f"{len(names) - len(set(names))} run name collision(s), e.g. "
                         f"{duplicated[:3]}. Two configs would share an output file.")

    on_disk = sum(1 for n in names if os.path.exists(os.path.join(output_path, f'{n}.json')))
    total = len(combos)
    print(f"Running {total} simulations with {WORKERS} worker processes "
          f"({on_disk} already on disk, reused unless parameter_sweep.RERUN)...")
    print(f"Output: {output_path}")
    print(f"  instellation: {min(INSTELLATION):g} to {max(INSTELLATION):g}")
    print(f"  land_fraction: {sorted({c[11] for c in combos})}")
    print(f"  pe: {[f'{v:g} ({_pe_label(v)})' for v in sorted({c[10] for c in combos})]}")
    print(f"  outgassing={OUTGASSING:g}x  crust={CRUST_PRODUCTION:g}x  "
          f"depth={OCEAN_DEPTH:g} m  Mg/Si={MG_SI:g}  dIW={DELTA_IW:g}")
    print(f"  alpha={ALPHA_CALIB:g}  kd_mg_ht={KD_MG_CALIB:g}  k_na={K_NA_CALIB:g}")
    ps._warn_constant_drift()

    completed = aborted = 0
    with ProcessPoolExecutor(max_workers=WORKERS, mp_context=mp.get_context('spawn')) as executor:
        futures = [executor.submit(run_simulation, *combo[:11], output_path, combo[11])
                   for combo in combos]
        for future in as_completed(futures):
            completed += 1
            run_name, error, T, termination = future.result()
            if error:
                print(f"[{completed}/{total}] FAILED {run_name}: {error}", flush=True)
                continue
            if termination == 'fallback_limit':
                aborted += 1
            T_str = f"{T:.1f} K" if T is not None else "T unknown"
            print(f"[{completed}/{total}] Done: {run_name} ({T_str}, "
                  f"{termination or 'unknown'})", flush=True)

    print("All simulations complete.")
    if aborted:
        print(f"{aborted}/{total} run(s) hit the fallback cap and were recorded "
              f"as 'fallback_limit'.")


# ══ Figures ═══════════════════════════════════════════════════════════════════════════════════
# Drawn through plot_results, so these share its style file, page geometry, panel groups,
# Damkohler line styling and termination markers with every other figure in the paper.
#
# plot_results is imported lazily: `spawn` re-imports this module in every worker process, and
# pulling matplotlib and pandas into 24 workers that will never draw anything is pure cost.

# Land-free blue against continental brown, the one colour decision these figures make.
ARM_COLOURS = {0.3: '#a4632a', 0.0: '#2a6fa4'}
ARM_LABELS = {0.3: 'Continental (land fraction 0.3)', 0.0: 'Ocean world (land free)'}

# Modern Earth, for the reference marker. Salinity is the sum of the model's tracked ions at
# their seawater concentrations, so it is comparable with the model's own salinity column.
EARTH = {'S': 1.0, 'T': 288.0, 'P_CO2': 280e-6, 'pH': 8.1,
         'salinity': (2.0e-3 * 61.0 + 0.1e-3 * 60.1 + 10.3e-3 * 40.1 +
                      52.8e-3 * 24.3 + 480e-3 * 23.0 + 550e-3 * 35.45)}


def _plot_results():
    import plot_results
    return plot_results


def _arm(df, land, pr):
    """Rows for one land-fraction arm of the baseline, at the Earth reference on every other axis."""
    return df[
        pr._ref_crust(df) &
        pr._ref_redox(df) &
        pr._ref_chem(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == OCEAN_DEPTH) &
        (df['outgassing'] == OUTGASSING) &
        (df['crust_production'] == CRUST_PRODUCTION) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['land_fraction'], land)
    ].sort_values('instellation')


def hz_edges(group, pr):
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
    trusted = g['termination'].isin(pr.HABITABLE).to_numpy() & np.isfinite(T)
    hab = trusted & (T > pr.T_SNOWBALL) & (T < pr.T_RUNAWAY)
    if not hab.any():
        return None

    idx = np.flatnonzero(hab)
    i0, i1 = int(idx[0]), int(idx[-1])

    outer, outer_kind = S[i0], 'open'
    if i0 > 0:
        if trusted[i0 - 1] and T[i0 - 1] <= pr.T_SNOWBALL:
            outer = float(np.interp(pr.T_SNOWBALL, [T[i0 - 1], T[i0]], [S[i0 - 1], S[i0]]))
            outer_kind = 'crossing'
        elif wall[i0 - 1] == 'cold':
            outer, outer_kind = 0.5 * (S[i0 - 1] + S[i0]), 'bracketed'

    inner, inner_kind = S[i1], 'open'
    if i1 + 1 < len(S):
        if trusted[i1 + 1] and T[i1 + 1] >= pr.T_RUNAWAY:
            inner = float(np.interp(pr.T_RUNAWAY, [T[i1], T[i1 + 1]], [S[i1], S[i1 + 1]]))
            inner_kind = 'crossing'
        elif wall[i1 + 1] == 'hot':
            inner, inner_kind = 0.5 * (S[i1] + S[i1 + 1]), 'bracketed'

    return float(outer), float(inner), outer_kind, inner_kind


def plot_baseline_vs_ocean(arms, output_path, pr):
    """T, pCO2, pH and salinity against instellation, continental arm against ocean arm.

    Both lines carry plot_results' Damkohler styling, so the comparison also shows whether the
    two arms sit in the same weathering regime -- which they do not: continental weathering is
    transport-limited on Earth (Da >> 1) while the land-free worlds are kinetically limited.
    """
    handles = [pr.Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6, label=ARM_LABELS[l])
               for l in arms] + list(pr.DA_LEGEND)

    for cols, sfx in pr._panel_groups(True):
        fig, axes = pr.plt.subplots(len(cols), 1, sharex=True,
                                    figsize=pr.figure_size('single', n_rows=len(cols),
                                                           row_height=2.0))
        for land, group in arms.items():
            pr._plot_group_on_axes(axes, group, ARM_COLOURS[land], show_markers=False, cols=cols)
        pr._style_axes(axes, cols)
        for ax, col in zip(axes, cols):
            ax.scatter(EARTH['S'], EARTH[col], marker='*', s=180, color='gold',
                       edgecolors='k', linewidths=0.7, zorder=6)
        pr._add_figure_legend(fig, axes, handles)
        pr._save_fig(fig, os.path.join(output_path, f'continental_vs_ocean{sfx}.png'))


def plot_habitable_zone(arms, output_path, pr):
    """The headline figure: where the model keeps a planet temperate, with land and without.

    Upper panel is the temperature curve each zone is read off; lower panel is the zone itself,
    one bar per arm on the same instellation axis. Edges located only as bounds (see `hz_edges`)
    carry a caret pointing the way the true edge lies, so a bar that is merely wider than the
    sweep could resolve cannot be read as a measured one.
    """
    edges = {}
    for land, group in arms.items():
        got = hz_edges(group, pr)
        if got is not None:
            edges[land] = got
    if LAND_FRACTION not in edges:
        print("No habitable band on the continental arm -- skipping the habitable-zone figure.")
        return edges

    fig, (ax, ax_z) = pr.plt.subplots(2, 1, sharex=True, height_ratios=[3, 1],
                                      figsize=pr.figure_size('single', height=4.0))

    for land, group in arms.items():
        pr._plot_group_on_axes([ax], group, ARM_COLOURS[land], show_markers=False, cols=['T'])
    pr._style_axes([ax], ['T'])
    ax.set_xlabel('')
    ax.scatter(EARTH['S'], EARTH['T'], marker='*', s=180, color='gold', edgecolors='k',
               linewidths=0.7, zorder=6)

    for row, (land, (lo, hi, lo_kind, hi_kind)) in enumerate(edges.items()):
        colour = ARM_COLOURS[land]
        y = len(edges) - 1 - row
        ax_z.barh(y, hi - lo, left=lo, height=0.5, color=colour, alpha=0.35,
                  edgecolor=colour, linewidth=1.4, zorder=3)
        for x, kind, marker in ((lo, lo_kind, '<'), (hi, hi_kind, '>')):
            if kind == 'open':
                ax_z.scatter(x, y, marker=marker, s=30, color=colour, zorder=5)
        ax_z.text(0.5 * (lo + hi), y, f'{lo:.2f}–{hi:.2f}', ha='center', va='center',
                  fontsize=7, zorder=6)
        ax.axvspan(lo, hi, color=colour, alpha=0.07, zorder=0)

    ax_z.set_yticks(range(len(edges)))
    ax_z.set_yticklabels([])
    ax_z.set_ylim(-0.6, len(edges) - 0.4)
    ax_z.set_ylabel('Habitable\nzone')
    ax_z.set_xlabel('Instellation (S/S₀)')
    ax_z.grid(True, axis='x', linestyle='--', alpha=0.4, zorder=0)
    ax_z.set_xlim(*ax.get_xlim())

    handles = [pr.Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6, label=ARM_LABELS[l])
               for l in arms]
    handles.append(pr.Line2D([0], [0], color='k', linestyle='none', marker='>', markersize=5,
                             label='Edge is a bound (sweep limit)'))
    pr._add_figure_legend(fig, [ax, ax_z], handles)
    pr._save_fig(fig, os.path.join(output_path, 'continental_habitable_zone.png'))
    return edges


def _report(edges, pr):
    """Print the habitable-zone edges, with how each was located. See `hz_edges`."""
    if not edges:
        return
    print(f"\nHabitable zone (T between {pr.T_SNOWBALL:.0f} K and {pr.T_RUNAWAY:.0f} K), "
          f"{_pe_label(pr.REF_PE)} ocean, {OCEAN_DEPTH/1000:g} km, outgassing "
          f"{OUTGASSING:g}x, crust {CRUST_PRODUCTION:g}x Earth:")
    print(f"  {'arm':>32s} {'outer S':>8s} {'inner S':>8s} {'width':>7s}   how located")
    for land, (lo, hi, lo_kind, hi_kind) in edges.items():
        print(f"  {ARM_LABELS[land]:>32s} {lo:8.3f} {hi:8.3f} {hi - lo:7.3f}"
              f"   outer {lo_kind}, inner {hi_kind}")
    if len(edges) == 2:
        (a_lo, a_hi, _, _), (b_lo, b_hi, _, _) = edges[LAND_FRACTION], edges[0.0]
        print(f"  continental zone is {(a_hi - a_lo) - (b_hi - b_lo):+.3f} S wide relative to "
              f"the ocean world ({(a_hi - a_lo) / (b_hi - b_lo):.2f}x)")


def make_plots(output_path=OUTPUT_PATH, pe=None):
    pr = _plot_results()
    if pe is not None:
        pr.REF_PE = pe
    df = pr.load_data(output_path)
    if df.empty:
        print(f"No runs found in {output_path}.")
        return

    arms = {}
    for land in LAND_ARMS:
        sub = _arm(df, land, pr)
        if not sub.empty:
            arms[land] = pr._add_diag_columns(sub, output_path).sort_values('instellation')
    if LAND_FRACTION not in arms:
        print(f"No runs at land_fraction = {LAND_FRACTION:g} -- run the sweep first.")
        return
    for land, group in arms.items():
        print(f"  {len(group)} run(s) on the {ARM_LABELS[land].lower()} arm.")

    plot_baseline_vs_ocean(arms, output_path, pr)
    edges = plot_habitable_zone(arms, output_path, pr)
    _report(edges, pr)
    # plot_results' own continental figures: the four-panel baseline and the ion-ratio chart
    # against modern seawater. They select on this same reference, so they read these runs.
    pr.plot_continental_baseline(df, output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--path', default=OUTPUT_PATH,
                        help='Directory for the run JSONs and figures.')
    parser.add_argument('--plot-only', action='store_true',
                        help='Draw the figures from runs already on disk; run nothing.')
    parser.add_argument('--no-plots', action='store_true',
                        help='Run the sweep and stop, without drawing anything.')
    parser.add_argument('--both-redox', action='store_true',
                        help=f'Also run the oxidising arm (pe = {PE_OXIDISING:g}). The figures '
                             f'are drawn at one pe either way -- see --pe.')
    parser.add_argument('--pe', type=float, default=None,
                        help='Ocean pe the figures are drawn at (default: the model reference, '
                             f'{PE_REDUCING:g}, reducing).')
    args = parser.parse_args()

    if not args.plot_only:
        pe_states = [PE_REDUCING, PE_OXIDISING] if args.both_redox else PE_STATES
        run(_combos(pe_states=pe_states), output_path=args.path)

    if not args.no_plots:
        make_plots(output_path=args.path, pe=args.pe)
    print("Done.")
