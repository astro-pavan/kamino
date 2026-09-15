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
import pandas as pd

import parameter_sweep as ps
from parameter_sweep import (ALPHA_CALIB, KD_MG_CALIB, K_NA_CALIB, OUTPUT_PATH, WORKERS,
                             PE_REDUCING, PE_OXIDISING, _pe_label, _run_name, run_simulation)
from kamino.constants import EARTH_MANTLE_MG_SI, EARTH_DELTA_IW

# ── The baseline planet: Earth, on every axis ─────────────────────────────────────────────────
# Earth's land fraction. `Planet` scales the continental flux linearly off this
# (`_s_terr = _S_TERR_EARTH * land_fraction / 0.3`), so 0.3 is the value the terrestrial
# denudation rate was calibrated at, not an arbitrary point on a continuum.
LAND_FRACTION = 0.3

# Earth's mantle molar Mg/Si (1.25). Named apart from the GRID_MG_SI axis below so the two
# cannot be confused: this one is the reference every baseline figure pins to.
MG_SI_EARTH = float(EARTH_MANTLE_MG_SI)

# Both arms of the baseline comparison. 0.0 is the ocean world every other sweep runs; it is
# listed second so the continental runs -- the ones that do not exist yet -- are submitted first.
LAND_ARMS = [LAND_FRACTION, 0.0]

# ── The land-fraction series ──────────────────────────────────────────────────────────────────
# Turns the continental sink down from Earth's to nothing, to find where seafloor weathering
# takes over as the dominant alkalinity source.
#
# This is a clean one-variable experiment because of how planet.py scales the two sinks against
# the ocean's mass, which is itself set by the seafloor area (2026-09-04 area fix):
#
#   continental  F_sil * (f * A)     / (d * (1-f) * A * 1000)   ->  scales as f / (1-f)
#   seafloor     flux_LT * ((1-f)*A) / (d * (1-f) * A * 1000)   ->  independent of f
#
# So land fraction turns the continental sink up and down and leaves the seafloor sink's
# CONCENTRATION rate geometrically untouched. The seafloor flux still responds, but only through
# the climate and ocean chemistry it shares with the continents, never through geometry.
#
# LOG spacing, because the crossover is nowhere near the middle of a linear range. Measured on
# the land 0.3 runs, continental alkalinity is ~21 Tmol eq/yr against a seafloor ~0.018 -- a
# ratio near 1200, and near 1700 once the area fix makes continental 1.43x stronger. Three and a
# half decades of land fraction are needed to close that, and half-decade steps locate the
# crossover to within a factor of ~3.
#
# 0.2 and 0.1 are kept at full resolution because that is the range real terrestrial planets
# plausibly occupy; below 0.03 the grid only has to bracket a crossing. Exactly 0.0 is included
# as the end member -- it is the ocean world, already on disk, and costs nothing.
LAND_FRACTIONS = [0.3, 0.2, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0]

# ── The coarse multi-axis grid ────────────────────────────────────────────────────────────────
# Trades resolution on the two axes above for the three the series holds fixed, to ask whether
# the continental/seafloor crossover MOVES with tectonics and crust chemistry or just sits where
# the Earth-reference series put it (~1e-3 land fraction).
#
# Instellation coarsens to 0.1 steps over 0.4-1.2: outside that every run left the domain in the
# fine series, so the trimmed range costs nothing. Land fraction coarsens to decade steps, which
# locates a crossover to within a factor of ~10 -- enough to see it move, not enough to quote.
# Both keep values the fine grids already use, so no run is orphaned between designs.
COARSE_INSTELLATION = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
COARSE_LAND_FRACTIONS = [0.3, 0.03, 0.003, 0.0003, 0.0]

# Three points per decade-and-a-bit on each tectonic axis, bracketing Earth. Ints where
# parameter_sweep uses ints (see the note above), so the land-free corners match the runs
# sweep_basic and sweep_basic_high_mgsi already wrote and are reused rather than recomputed.
GRID_OUTGASSING = [0.1, 1, 10]
GRID_CRUST = [0.1, 1, 10]

# Earth's mantle Mg/Si and the olivine-rich end member. 1.8 sits past the ~1.69 ceiling section
# 25.4 measured, but that predates Akermanite closing the norm (25.5); at 1.8 the assemblage sums
# to 1.0 with no mass-balance warning, Akermanite taking 11 wt%. State that if these are published.
GRID_MG_SI = [MG_SI_EARTH, 1.8]

# The reactive-area scaling in the seafloor weathering law. This is the axis worth resolving:
# in the kinetic limit the seafloor flux is LINEAR in alpha and continental weathering does not
# see alpha at all, so f* (the crossover land fraction) should go as alpha^1. Habitable planets
# are measured to be kinetic (median Da 0.007 over the steady states on disk), so that linearity
# should hold across the habitable population rather than only in a corner of it.
#
# It matters because alpha is NOT identifiable from Earth (development history 28.2) and spans
# 1.1-50 here -- a 45x range on a quantity nothing observable pins down. If f* really is linear
# in alpha then alpha, a model parameter, is a larger control on the crossover than any planetary
# property in the grid, which is worth knowing explicitly rather than by inference.
#
# Same three values as parameter_sweep.alpha, so these runs sit in the same family as the
# land-free alpha arm; all three stay in the kinetic limit (Da <= 0.13).
GRID_ALPHA = [ALPHA_CALIB, 10, 50]

# Which sweep __main__ runs -- edit this rather than passing a flag.
#   'baseline' : the two-arm instellation line (land 0.3 and 0), Earth on every other axis
#   'land'     : the land-fraction series at Earth outgassing and crust production
#   'grid'     : the coarse instellation x land x outgassing x crust x Mg/Si factorial
#   'alpha'    : instellation x land x outgassing x alpha, at Earth crust production and Mg/Si
SWEEP = 'all'

# These are ints on purpose. `_run_name` interpolates them with plain str(), so 1 and 1.0 give
# 'crust_1' and 'crust_1.0' -- two names for one config, and the ocean arm would stop matching
# the runs `sweep_basic` already wrote. Match the types parameter_sweep uses.
OUTGASSING = 1                # x Earth
CRUST_PRODUCTION = 1          # x Earth
OCEAN_DEPTH = 3000            # m

REVERSE_WEATHERING = True
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
        (s, OUTGASSING, CRUST_PRODUCTION, OCEAN_DEPTH, REVERSE_WEATHERING, MG_SI_EARTH, DELTA_IW,
         ALPHA_CALIB, KD_MG_CALIB, K_NA_CALIB, pe, land)
        for s, land, pe in itertools.product(instellation or INSTELLATION,
                                             land_arms or LAND_ARMS,
                                             pe_states or PE_STATES)
    ]
    combos.sort(key=lambda c: (ps._cost_rank(c), c[11] == 0.0))
    return combos


def _grid_combos(instellation=None, lands=None, outgassing=None, crust=None,
                 mg_si=None, alpha=None, pe_states=None):
    """Combos for the coarse factorial: instellation x land x outgassing x crust x Mg/Si.

    Same argument order and cost ordering as `_combos`; only the axes differ. The land-free
    corners reproduce the names `sweep_basic` and `sweep_basic_high_mgsi` already wrote, so they
    are reused off disk and only the land-bearing runs actually cost anything.
    """
    combos = [
        (s, o, c, OCEAN_DEPTH, REVERSE_WEATHERING, mg, DELTA_IW,
         a, KD_MG_CALIB, K_NA_CALIB, pe, land)
        for s, land, o, c, mg, a, pe in itertools.product(
            instellation or COARSE_INSTELLATION,
            lands or COARSE_LAND_FRACTIONS,
            outgassing or GRID_OUTGASSING,
            crust or GRID_CRUST,
            mg_si or GRID_MG_SI,
            alpha or [ALPHA_CALIB],
            pe_states or PE_STATES)
    ]
    combos.sort(key=lambda x: (ps._cost_rank(x), x[11] == 0.0))
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
    # Every axis is read off the COMBOS, never off the module constants. The constants describe
    # the baseline sweep only, so printing them made the grid sweep's log claim it had run at
    # Earth outgassing, Earth crust production and Earth Mg/Si while it was doing nothing of the
    # kind. The log is the record of what was run, so it has to be derived from what was run.
    def _axis(i):
        return ', '.join(f'{v:g}' for v in sorted({c[i] for c in combos}))

    print(f"  instellation:   {_axis(0)}")
    print(f"  land_fraction:  {_axis(11)}")
    print(f"  outgassing:     {_axis(1)}")
    print(f"  crust prod.:    {_axis(2)}")
    print(f"  mantle Mg/Si:   {_axis(5)}")
    print(f"  ocean depth:    {_axis(3)} m      dIW: {_axis(6)}")
    print(f"  pe: {[f'{v:g} ({_pe_label(v)})' for v in sorted({c[10] for c in combos})]}")
    print(f"  alpha:          {_axis(7)}")
    print(f"  kd_mg_ht={KD_MG_CALIB:g}  k_na={K_NA_CALIB:g}")
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


def _plot_results():
    import plot_results
    return plot_results


def _drop_clamped(group, pr):
    """Drop out-of-domain rows whose stored T is a box limit rather than a computed value."""
    clamped = (group['termination'].isin(pr.OUT_OF_DOMAIN) &
               ((group['T'] >= T_CLAMP_HOT) | (group['T'] <= T_CLAMP_COLD)))
    return group[~clamped]


def _draw_arm(axes, group, colour, cols, pr):
    """Draw one arm, and report whether any of it is habitable.

    `plot_results._plot_group_on_axes` returns without drawing when a group contains no habitable
    run at all, which is the right call for a facet of a larger figure but wrong here: at Earth's
    outgassing rate EVERY land-free run leaves the domain, and silently omitting the line would
    hide the very comparison this figure exists to make.

    Such an arm is drawn faint, with its clamp sentinels removed and plot_results' hollow
    per-termination markers on every point, so it reads as "measured, but not habitable". The
    line stays SOLID deliberately: all four line styles are already spoken for by DA_LEGEND, so a
    dashed fallback would read as "Da >= 1" rather than "not habitable". Colour separates the
    arms; the hollow markers say these are not habitable states.
    """
    if group['termination'].isin(pr.HABITABLE).any():
        pr._plot_group_on_axes(axes, group, colour, show_markers=False, cols=cols)
        return True
    shown = _drop_clamped(group, pr)
    for ax, col in zip(axes, cols):
        ax.plot(shown['instellation'], shown[col], color=colour, linewidth=1.2, alpha=0.5,
                zorder=2)
        for _, row in shown.iterrows():
            if np.isfinite(row[col]):
                ax.scatter(row['instellation'], row[col],
                           marker=pr.FAILED_MARKERS.get(row['termination'], 'x'), s=22,
                           facecolors='none', edgecolors=colour, linewidths=1.0, zorder=4)
    return False


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

    # A temperature window alone is not a habitability test in this model. Two states pass
    # `T < T_RUNAWAY` without being habitable at all: one past the runaway greenhouse, whose T is
    # read off the hot branch (see `_past_runaway`), and one above the OLR fit's 350 K ceiling,
    # where the climate model is extrapolating. Both are excluded here rather than by moving
    # T_RUNAWAY, which is a plot_results convention shared with every other figure.
    S_ok = np.array([not _past_runaway(s, p) for s, p in
                     zip(S, g['P_CO2'].to_numpy(dtype=float))])
    hab = trusted & (T > pr.T_SNOWBALL) & (T < pr.T_RUNAWAY) & (T <= OLR_FIT_T_MAX) & S_ok
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
        elif wall[i1 + 1] == 'hot' or not S_ok[i1 + 1]:
            # `not S_ok` is the runaway greenhouse: the neighbour has no cool-branch solution, so
            # the edge lies in this interval. Not interpolated -- the neighbour's T is on the hot
            # branch, so a line drawn through it has no meaning.
            inner, inner_kind = 0.5 * (S[i1] + S[i1 + 1]), 'bracketed'

    return float(outer), float(inner), outer_kind, inner_kind


def plot_baseline_vs_ocean(arms, output_path, pr):
    """T, pCO2, pH and salinity against instellation, continental arm against ocean arm.

    Both lines carry plot_results' Damkohler styling, so the comparison also shows whether the
    two arms sit in the same weathering regime -- which they do not: continental weathering is
    transport-limited on Earth (Da >> 1) while the land-free worlds are kinetically limited.
    """
    habitable = {land: group['termination'].isin(pr.HABITABLE).any()
                 for land, group in arms.items()}
    handles = [pr.Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6,
                         alpha=1.0 if habitable[l] else 0.5,
                         marker='' if habitable[l] else 's', markerfacecolor='none',
                         label=ARM_LABELS[l] + ('' if habitable[l]
                                                else ' — never habitable'))
               for l in arms] + list(pr.DA_LEGEND)

    for cols, sfx in pr._panel_groups(True):
        fig, axes = pr.plt.subplots(len(cols), 1, sharex=True,
                                    figsize=pr.figure_size('single', n_rows=len(cols),
                                                           row_height=2.0))
        for land, group in arms.items():
            _draw_arm(axes, group, ARM_COLOURS[land], cols, pr)
        pr._style_axes(axes, cols)
        for ax, col in zip(axes, cols):
            ax.scatter(EARTH['S'], EARTH[col], marker='*', s=180, color='gold',
                       edgecolors='k', linewidths=0.7, zorder=6)
        pr._add_figure_legend(fig, axes, handles)
        pr._save_fig(fig, pr.figure_path(output_path, f'continental_vs_ocean{sfx}.png'))


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

    habitable = {}
    for land, group in arms.items():
        habitable[land] = _draw_arm([ax], group, ARM_COLOURS[land], ['T'], pr)
    pr._style_axes([ax], ['T'])
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

    handles = [pr.Line2D([0], [0], color=ARM_COLOURS[l], linewidth=1.6,
                         alpha=1.0 if habitable[l] else 0.5,
                         marker='' if habitable[l] else 's', markerfacecolor='none',
                         label=ARM_LABELS[l] + ('' if habitable[l]
                                                else ' — never habitable'))
               for l in arms]
    if any(k == 'open' for e in edges.values() for k in e[2:]):
        handles.append(pr.Line2D([0], [0], color='k', linestyle='none', marker='>', markersize=5,
                                 label='Edge is a bound (sweep limit)'))
    pr._add_figure_legend(fig, [ax, ax_z], handles)
    pr._save_fig(fig, pr.figure_path(output_path, 'continental_habitable_zone.png'))
    return edges


def _report(arms, edges, pr):
    """Print the habitable-zone edges, with how each was located. See `hz_edges`."""
    print(f"\nHabitable zone (T between {pr.T_SNOWBALL:.0f} K and {pr.T_RUNAWAY:.0f} K), "
          f"{_pe_label(pr.REF_PE)} ocean, {OCEAN_DEPTH/1000:g} km, outgassing "
          f"{OUTGASSING:g}x, crust {CRUST_PRODUCTION:g}x Earth:")
    print(f"  {'arm':>32s} {'outer S':>8s} {'inner S':>8s} {'width':>7s}   how located")
    for land in arms:
        if land not in edges:
            walls = arms[land]['domain_wall'].dropna().value_counts()
            why = ', '.join(f'{n} {pr.WALL_LABELS.get(w, w)}' for w, n in walls.items())
            print(f"  {ARM_LABELS[land]:>32s} {'--':>8} {'--':>8} {'none':>7}"
                  f"   no habitable run ({why})")
            continue
        lo, hi, lo_kind, hi_kind = edges[land]
        print(f"  {ARM_LABELS[land]:>32s} {lo:8.3f} {hi:8.3f} {hi - lo:7.3f}"
              f"   outer {lo_kind}, inner {hi_kind}")
    if len(edges) == 2:
        (a_lo, a_hi, _, _), (b_lo, b_hi, _, _) = edges[LAND_FRACTION], edges[0.0]
        print(f"  continental zone is {(a_hi - a_lo) - (b_hi - b_lo):+.3f} S wide relative to "
              f"the ocean world ({(a_hi - a_lo) / (b_hi - b_lo):.2f}x)")

    # plot_results draws these edges as vertical lines on every other instellation figure, from
    # its own hardcoded copy. Say so loudly when the two disagree rather than let every figure
    # in the paper quote a stale zone -- the same reason parameter_sweep._warn_constant_drift
    # exists for the chemistry constants.
    if LAND_FRACTION in edges:
        lo, hi = edges[LAND_FRACTION][:2]
        for label, here, there in (('CONTINENTAL_HZ_OUTER', lo, pr.CONTINENTAL_HZ_OUTER),
                                   ('CONTINENTAL_HZ_INNER', hi, pr.CONTINENTAL_HZ_INNER)):
            if abs(here - there) > 5e-4:
                print(f"  NOTE plot_results.{label} = {there:g}, but this sweep measures "
                      f"{here:.3f}. Update it, or the HZ lines on every other figure are stale.")


def _land_series(df, output_path, pr):
    """{land_fraction: instellation-sorted rows with diagnostics}, for the land fractions present."""
    series = {}
    for land in LAND_FRACTIONS:
        sub = _arm(df, land, pr)
        if not sub.empty:
            series[land] = pr._add_diag_columns(sub, output_path).sort_values('instellation')
    return series


def _land_colours(lands, pr):
    """Colour per land fraction, log-scaled over the positive ones.

    0 cannot sit on a log scale and is not just 'a bit less land' -- it is the land-free ocean
    world every other sweep runs. It keeps the baseline figures' blue and its own legend entry.
    """
    positive = sorted(l for l in lands if l > 0)
    cmap = pr.LAND_FRACTION_CMAP
    # A single positive value gives LogNorm(vmin == vmax), which cannot be normalised. That
    # happens whenever only the baseline's own land fraction is on disk -- i.e. before the
    # land-fraction sweep has been run -- so it is the ordinary case, not an error.
    norm = (pr.mcolors.LogNorm(vmin=min(positive), vmax=max(positive))
            if len(positive) > 1 else None)
    if norm is not None:
        colours = {l: cmap(norm(l)) for l in positive}
    else:
        colours = {l: ARM_COLOURS[LAND_FRACTION] for l in positive}
    if 0.0 in lands:
        colours[0.0] = ARM_COLOURS[0.0]
    return colours, cmap, norm


def plot_land_fraction_series(series, output_path, pr):
    """T, pCO2, pH and salinity against instellation, one line per land fraction."""
    if len(series) < 2:
        print("Fewer than two land fractions on disk -- skipping the land-fraction series.")
        return
    colours, cmap, norm = _land_colours(series, pr)

    for cols, sfx in pr._panel_groups(True):
        fig, axes = pr.plt.subplots(len(cols), 1, sharex=True,
                                    figsize=pr.figure_size('single', n_rows=len(cols),
                                                           row_height=2.0))
        for land, group in sorted(series.items(), reverse=True):
            _draw_arm(axes, group, colours[land], cols, pr)
        pr._style_axes(axes, cols)
        if norm is not None:
            positive = sorted(l for l in series if l > 0)
            pr._add_colorbar(fig, list(axes), cmap, norm, 'Land fraction',
                             ticks=positive, ticklabels=[f'{v:g}' for v in positive],
                             aspect=len(cols) * 7.5)
        handles = []
        if 0.0 in series:
            handles.append(pr.Line2D([0], [0], color=ARM_COLOURS[0.0], linewidth=1.6,
                                     label='Land free (0)'))
        handles += list(pr.DA_LEGEND)
        pr._add_figure_legend(fig, axes, handles)
        pr._save_fig(fig, pr.figure_path(output_path, f'land_fraction_series{sfx}.png'))


def plot_weathering_crossover(series, output_path, pr):
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
        steady = group[group['termination'].isin(pr.HABITABLE)]
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
    fig, (ax, ax_c) = pr.plt.subplots(1, 2, figsize=pr.figure_size('double', height=2.8))

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
    pr._save_fig(fig, pr.figure_path(output_path, 'weathering_crossover.png'))

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
              f"{min(ratios_all):.3g} to {max(ratios_all):.3g}; extend LAND_FRACTIONS to bracket it.")
    return crossings


def plot_weathering_ratio_map(series, output_path, pr, levels=13):
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
        steady = group[group['termination'].isin(pr.HABITABLE)]
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
    norm = pr.mcolors.TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=max(hi, step))
    cmap = pr.cmr.fusion_r          # diverging, neutral (white) midpoint at ratio = 1

    fig, ax = pr.plt.subplots(1, 1, figsize=pr.figure_size('single', height=2.9))
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

    pr._save_fig(fig, pr.figure_path(output_path, 'weathering_ratio_map.png'))

    print("\nSeafloor / continental alkalinity flux (steady states only):")
    print(f"  {'land':>8} " + ' '.join(f'{s:>8.2f}' for s in s_vals))
    for i, land in enumerate(reversed(lands)):
        row = Z[len(lands) - 1 - i]
        cells = [(f"{10 ** v:8.2g}" if np.isfinite(v) else f"{'--':>8}") for v in row]
        print(f"  {land:8g} " + ' '.join(cells))
    if dropped_pts:
        print(f"  ({len(dropped_pts)} cell(s) blank: no steady state)")
    return ratio


def _grid_slice(df, pr, outgassing, crust, mg_si):
    """{land_fraction: rows} for one (outgassing, crust production, Mg/Si) cell of the grid.

    Restricted to the COARSE axes even where finer runs exist. The Earth-reference cell
    (out 1x, crust 1x, Mg/Si 1.25) is also where the land-fraction series ran, so without this it
    would carry ~3x the instellation samples and two extra land fractions. Contour interpolation
    depends on sampling density, so that one panel would be smoother and reach further in
    instellation than its neighbours -- a difference in the sampling, read as a difference in the
    physics. The fine runs keep their own figure (plot_weathering_ratio_map).
    """
    sub = df[
        df['instellation'].isin(COARSE_INSTELLATION) &
        df['land_fraction'].apply(
            lambda v: any(np.isclose(v, l) for l in COARSE_LAND_FRACTIONS)) &
        pr._ref_redox(df) &
        pr._ref_chem(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == OCEAN_DEPTH) &
        (df['outgassing'] == outgassing) &
        (df['crust_production'] == crust) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], DELTA_IW)
    ]
    return {float(l): sub[np.isclose(sub['land_fraction'], l)].sort_values('instellation')
            for l in sorted(sub['land_fraction'].unique())}


def _ratio_cells(series, pr, output_path):
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
        group = pr._add_diag_columns(group, output_path)
        steady = group[group['termination'].isin(pr.HABITABLE)]
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


def plot_weathering_ratio_grid(df, output_path, pr, step=0.5):
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
    mg_vals = [m for m in GRID_MG_SI if np.isclose(df['mg_si'], m).any()]
    if not mg_vals:
        print("No grid runs on disk -- skipping the faceted ratio map.")
        return None

    cells = {}
    for mg in mg_vals:
        for c in GRID_CRUST:
            for o in GRID_OUTGASSING:
                series = _grid_slice(df, pr, o, c, mg)
                if series:
                    cells[(mg, c, o)] = _ratio_cells(series, pr, output_path)

    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state grid runs with both fluxes positive -- skipping.")
        return None
    lo = np.floor(min(allv) / step) * step
    hi = np.ceil(max(allv) / step) * step
    bands = np.arange(lo, hi + 0.5 * step, step)
    norm = pr.mcolors.TwoSlopeNorm(vmin=lo, vcenter=0.0, vmax=max(hi, step))
    cmap = pr.cmr.fusion_r
    ticks = [t for t in range(-9, 10) if lo <= t <= hi]

    for mg in mg_vals:
        fig, axes = pr.plt.subplots(len(GRID_CRUST), len(GRID_OUTGASSING), sharex=True,
                                    sharey=True, squeeze=False,
                                    figsize=pr.figure_size('double', height=5.0))
        cf = None
        for i, c in enumerate(reversed(GRID_CRUST)):
            for j, o in enumerate(GRID_OUTGASSING):
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
                if i == len(GRID_CRUST) - 1:
                    ax.set_xlabel('Instellation (S/S0)')

        if cf is not None:
            cbar = fig.colorbar(cf, ax=list(axes.ravel()), pad=0.02, aspect=30, ticks=ticks)
            cbar.set_label('Seafloor / continental alkalinity flux')
            cbar.set_ticklabels([('1' if t == 0 else f'$10^{{{t}}}$') for t in ticks])
            cbar.ax.axhline(0, color='k', linewidth=1.2)
        fig.suptitle(f'Mantle Mg/Si = {mg:g}', fontsize=9)
        pr._save_fig(fig, pr.figure_path(output_path, f'weathering_ratio_grid_mgsi{mg:g}.png'))

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


def _alpha_slice(df, pr, outgassing, alpha, crust=None, mg_si=None):
    """{land_fraction: rows} for one (outgassing, alpha) cell.

    Deliberately does NOT use pr._ref_chem: that helper pins alpha to the most-run value, which
    would silently discard the alpha = 10 and 50 arms and leave a "sweep" of one point. kd_mg and
    k_na are still pinned, explicitly, to the calibrated values.
    """
    crust = CRUST_PRODUCTION if crust is None else crust
    mg_si = MG_SI_EARTH if mg_si is None else mg_si
    sub = df[
        df['instellation'].isin(COARSE_INSTELLATION) &
        df['land_fraction'].apply(
            lambda v: any(np.isclose(v, l) for l in COARSE_LAND_FRACTIONS)) &
        pr._ref_redox(df) &
        df['reverse_weathering'] &
        (df['ocean_depth'] == OCEAN_DEPTH) &
        (df['outgassing'] == outgassing) &
        (df['crust_production'] == crust) &
        (df['f_HT'] == 0.0) &
        np.isclose(df['mg_si'], mg_si) &
        np.isclose(df['delta_iw'], DELTA_IW) &
        np.isclose(df['alpha'], alpha) &
        np.isclose(df['kd_mg'], KD_MG_CALIB) &
        np.isclose(df['k_na'], K_NA_CALIB)
    ]
    return {float(l): sub[np.isclose(sub['land_fraction'], l)].sort_values('instellation')
            for l in sorted(sub['land_fraction'].unique())}


def plot_alpha_scaling(df, output_path, pr):
    """Does the crossover land fraction really go as alpha^1?

    In the kinetic limit the seafloor flux is linear in alpha while continental weathering does
    not see alpha at all, so f* should scale as alpha^1. The climate feedback should DAMP that:
    raising alpha strengthens the sink, which cools the planet and draws CO2 down, weakening both
    fluxes again. An exponent below 1 is therefore the expected outcome, and its size is the
    result -- it says how much of alpha's nominal leverage survives the feedback.
    """
    rows = []
    for o in GRID_OUTGASSING:
        for a in GRID_ALPHA:
            series = _alpha_slice(df, pr, o, a)
            if not series:
                continue
            ratio, _, _ = _ratio_cells(series, pr, output_path)
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

    fig, ax = pr.plt.subplots(1, 1, figsize=pr.figure_size('single', height=3.0))
    cmap = pr.OUTGASSING_CMAP
    norm = pr.mcolors.LogNorm(vmin=min(GRID_OUTGASSING), vmax=max(GRID_OUTGASSING))

    print("\nCrossover land fraction f* against alpha (geometric mean over instellation):")
    print(f"  {'outgassing':>10} " + ' '.join(f'{a:>10.4g}' for a in GRID_ALPHA)
          + f" {'exponent':>9}")
    exponents, anchor = {}, None
    for o in GRID_OUTGASSING:
        g = tab[tab.outgassing == o]
        if g.empty:
            continue
        # Geometric mean over instellation: f* spans decades, so an arithmetic mean would be
        # dominated by whichever instellation sits nearest the runaway.
        means = {a: float(np.exp(np.log(g[g.alpha == a].f_star).mean()))
                 for a in GRID_ALPHA if (g.alpha == a).any()}
        cells = [(f'{means[a]:10.4g}' if a in means else f'{"--":>10}') for a in GRID_ALPHA]
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
        xs = np.array(GRID_ALPHA, dtype=float)
        ax.plot(xs, f0 * xs / a0, color='0.4', linestyle=(0, (6, 3)), linewidth=1.2,
                label=r'$\propto \alpha$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Reactive area scaling $\alpha$')
    ax.set_ylabel(r'Crossover land fraction $f^*$')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=7, frameon=False, title='Outgassing', title_fontsize=7)
    pr._save_fig(fig, pr.figure_path(output_path, 'alpha_scaling.png'))

    if exponents:
        v = list(exponents.values())
        print(f"  exponent d log f* / d log alpha: {min(v):.2f} to {max(v):.2f} "
              f"(alpha^1 would be 1.00)")
    tab.to_csv(os.path.join(output_path, 'alpha_crossover.csv'), index=False)
    return tab


def plot_alpha_ratio_grid(df, output_path, pr, step=0.5):
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
    ROWS, COLS = GRID_ALPHA, GRID_OUTGASSING          # rows: alpha, columns: outgassing
    row_label, col_label = r'$\alpha$', 'outgassing'

    cells = {}
    for a in ROWS:
        for o in COLS:
            series = _alpha_slice(df, pr, o, a)
            if series:
                cells[(a, o)] = _ratio_cells(series, pr, output_path)

    allv = [np.log10(v) for r, *_ in cells.values() for v in r.values()]
    if not allv:
        print("No steady-state alpha runs with both fluxes positive -- skipping.")
        return None
    lo = np.floor(min(allv) / step) * step
    hi = np.ceil(max(allv) / step) * step
    bands = np.arange(lo, hi + 0.5 * step, step)
    norm = pr.mcolors.TwoSlopeNorm(vmin=min(lo, -step), vcenter=0.0, vmax=max(hi, step))
    cmap = pr.cmr.fusion_r

    fig, axes = pr.plt.subplots(len(ROWS), len(COLS), sharex=True, sharey=True, squeeze=False,
                                figsize=pr.figure_size('double', height=5.0))
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
    pr._save_fig(fig, pr.figure_path(output_path, 'weathering_ratio_alpha_grid.png'))
    return cells


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
    _report(arms, edges or {}, pr)

    # The land-fraction series: only draws once intermediate land fractions are on disk, so this
    # is a no-op until the sweep has been run with RUN_LAND_FRACTION_SWEEP.
    series = _land_series(df, output_path, pr)
    if len(series) > len(LAND_ARMS):
        print(f"\n  land fractions on disk: {sorted(series, reverse=True)}")
    plot_land_fraction_series(series, output_path, pr)
    plot_weathering_crossover(series, output_path, pr)
    plot_weathering_ratio_map(series, output_path, pr)
    plot_weathering_ratio_grid(df, output_path, pr)
    plot_alpha_scaling(df, output_path, pr)
    plot_alpha_ratio_grid(df, output_path, pr)
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
    parser.add_argument('--both-redox', action='store_true', default=False,
                        help=f'Also run the oxidising arm (pe = {PE_OXIDISING:g}). The figures '
                             f'are drawn at one pe either way -- see --pe.')
    parser.add_argument('--pe', type=float, default=None,
                        help='Ocean pe the figures are drawn at (default: the model reference, '
                             f'{PE_REDUCING:g}, reducing).')
    args = parser.parse_args()

    if not args.plot_only:
        pe_states = [PE_REDUCING, PE_OXIDISING] if args.both_redox else PE_STATES
        if SWEEP == 'grid':
            combos = _grid_combos(pe_states=pe_states)
        elif SWEEP == 'alpha':
            # Crust production and Mg/Si pinned to Earth: the grid sweep already showed crust
            # production is a secondary control and composition a minor one, and holding them
            # fixed keeps this factorial to a size worth running.
            combos = _grid_combos(crust=[CRUST_PRODUCTION], mg_si=[MG_SI_EARTH],
                                  alpha=GRID_ALPHA, pe_states=pe_states)
        elif SWEEP == 'land':
            combos = _combos(land_arms=LAND_FRACTIONS, pe_states=pe_states)
        elif SWEEP == 'baseline':
            combos = _combos(land_arms=LAND_ARMS, pe_states=pe_states)
        elif SWEEP == 'all':
            combos = _grid_combos(crust=[CRUST_PRODUCTION], mg_si=[MG_SI_EARTH], alpha=GRID_ALPHA, pe_states=pe_states)
        else:
            raise SystemExit(f"SWEEP must be 'baseline', 'land' or 'grid', not {SWEEP!r}")
        print(f"sweep: {SWEEP}")

        if SWEEP == 'all':
            combos = _grid_combos(crust=[CRUST_PRODUCTION], mg_si=[MG_SI_EARTH], alpha=GRID_ALPHA, pe_states=pe_states)
            run(combos, output_path=args.path)
            combos = _combos(land_arms=LAND_FRACTIONS, pe_states=pe_states)
            run(combos, output_path=args.path)
            combos = _combos(land_arms=LAND_ARMS, pe_states=pe_states)
            run(combos, output_path=args.path)
        else:
            run(combos, output_path=args.path)

    if not args.no_plots:
        make_plots(output_path=args.path, pe=args.pe)
    print("Done.")
