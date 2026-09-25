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

The ocean starts blank, as in `parameter_sweep` (set KAMINO_SEED_OCEAN=1 for the seawater seed);
runs on disk from the other initial condition are re-run, not reused.

    # run the sweep, then draw the figures (resumes; runs already on disk are reused)
    /data/pt426/big-venv/bin/python experiments/continental_baseline.py

    # re-draw the figures from runs already on disk
    /data/pt426/big-venv/bin/python experiments/continental_baseline.py --plot-only

The 'earth' sweep is the exception: it reproduces the calibration setup (experiments/calibrate_earth.py),
a seawater-seeded ocean with Earth's Cl outgassing ratio run for 4 Gyr, so its S = 1 run is the
calibrated Earth. The continental_baseline_{tp,chem,ions} figures are drawn from it.
"""

import argparse
import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import parameter_sweep as ps
from parameter_sweep import (ALPHA_CALIB, KD_MG_CALIB, K_NA_CALIB, OUTPUT_PATH, WORKERS,
                             PE_REDUCING, PE_OXIDISING, _pe_label, _run_name, run_simulation)
from kamino.constants import EARTH_MANTLE_MG_SI, EARTH_DELTA_IW, EARTH_CL_OUTGASSING_RATIO

# ── The baseline planet: Earth, on every axis ─────────────────────────────────────────────────
# Earth's land fraction. `Planet` scales the continental weathering and the aeolian dust flux
# (Jickells et al. 2005) with land area relative to Earth's, so 0.3 is Earth's reference point.
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
GRID_ALPHA = [ALPHA_CALIB * 0.1, ALPHA_CALIB, ALPHA_CALIB * 10]

# Which sweep __main__ runs -- edit this rather than passing a flag.
#   'baseline' : the two-arm instellation line (land 0.3 and 0), Earth on every other axis
#   'land'     : the land-fraction series at Earth outgassing and crust production
#   'grid'     : the coarse instellation x land x outgassing x crust x Mg/Si factorial
#   'alpha'    : instellation x land x outgassing x alpha, at Earth crust production and Mg/Si
#   'earth'    : the instellation line at land 0.3 in the calibration setup (seeded, Earth Cl, 4 Gyr)
SWEEP = 'all'

# The 'earth' sweep: the conditions calibrate_earth.py fits the constants in, so S = 1 is the calibrated Earth.
EARTH_CL_RATIO = EARTH_CL_OUTGASSING_RATIO   # Cl/C outgassing ratio (the sweeps run with none)
EARTH_T_END_GYR = 4.0                        # calibrate_earth.T_END
EARTH_SEED = True                            # seawater initial ocean: Cl 546, SO4 28.2, K 10.2 mM

# These are ints on purpose. `_run_name` interpolates them with plain str(), so 1 and 1.0 give
# 'crust_1' and 'crust_1.0' -- two names for one config, and the ocean arm would stop matching
# the runs `sweep_basic` already wrote. Match the types parameter_sweep uses.
OUTGASSING = 1                # x Earth
CRUST_PRODUCTION = 1          # x Earth
OCEAN_DEPTH = 3700            # m

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


def run(combos, output_path=OUTPUT_PATH, cl=ps.CL_OUTGASSING_RATIO, t_end_gyr=ps.T_END_GYR, seed=None):
    """Execute a combo list. Mirrors `parameter_sweep.run_combos`, carrying land fraction through.

    `cl`, `t_end_gyr` and `seed` apply to every combo; the defaults are the sweep's (no Cl, 2 Gyr,
    ps.SEED_OCEAN).
    """
    if not output_path.endswith('/'):
        output_path += '/'
    ps.p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    # Distinct configs must map to distinct filenames or one silently overwrites the other, and
    # the resume path then hands back the survivor's result for both (the fast_13 resume trap).
    names = [_run_name(*combo, cl=cl, t_end_gyr=t_end_gyr) for combo in combos]
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
    print(f"  Cl ratio: {cl:g}   t_end: {t_end_gyr:g} Gyr   initial ocean: "
          f"{'seawater seed' if (ps.SEED_OCEAN if seed is None else seed) else 'blank'}")
    ps._warn_constant_drift()

    completed = aborted = 0
    with ProcessPoolExecutor(max_workers=WORKERS, mp_context=mp.get_context('spawn')) as executor:
        futures = [executor.submit(run_simulation, *combo[:11], output_path, combo[11],
                                   cl=cl, t_end_gyr=t_end_gyr, seed=seed)
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


def make_plots(output_path=OUTPUT_PATH, pe=None):
    # Imported lazily: spawned workers re-import this module and never draw anything.
    import plot_results as pr
    if pe is not None:
        pr.REF_PE = pe
    df = pr.load_data(output_path)
    if df.empty:
        print(f"No runs found in {output_path}.")
        return
    pr.plot_continental(df, output_path)


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
        elif SWEEP == 'earth':
            combos = _combos(land_arms=[LAND_FRACTION], pe_states=pe_states)
        elif SWEEP == 'all':
            combos = _grid_combos(crust=[CRUST_PRODUCTION], mg_si=[MG_SI_EARTH], alpha=GRID_ALPHA, pe_states=pe_states)
        else:
            raise SystemExit(f"SWEEP must be 'baseline', 'land', 'grid', 'alpha', 'earth' or 'all', "
                             f"not {SWEEP!r}")
        print(f"sweep: {SWEEP}")

        if SWEEP == 'all':
            combos = _grid_combos(crust=[CRUST_PRODUCTION], mg_si=[MG_SI_EARTH], alpha=GRID_ALPHA, pe_states=pe_states)
            run(combos, output_path=args.path)
            combos = _combos(land_arms=LAND_FRACTIONS, pe_states=pe_states)
            run(combos, output_path=args.path)
            combos = _combos(land_arms=LAND_ARMS, pe_states=pe_states)
            run(combos, output_path=args.path)
            combos = _combos(land_arms=[LAND_FRACTION], pe_states=pe_states)
            run(combos, output_path=args.path, cl=EARTH_CL_RATIO, t_end_gyr=EARTH_T_END_GYR,
                seed=EARTH_SEED)
        elif SWEEP == 'earth':
            run(combos, output_path=args.path, cl=EARTH_CL_RATIO, t_end_gyr=EARTH_T_END_GYR,
                seed=EARTH_SEED)
        else:
            run(combos, output_path=args.path)

    if not args.no_plots:
        make_plots(output_path=args.path, pe=args.pe)
    print("Done.")
