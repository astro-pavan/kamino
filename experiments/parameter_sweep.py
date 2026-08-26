import itertools
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet as p2
from kamino.planet import Planet, KD_MG_HT, K_NA_CONT_REMOVAL
from kamino.weathering import ALPHA_REF
from kamino.crust_composition import mineral_composition
from kamino.constants import M_EARTH, R_EARTH, EARTH_MANTLE_MG_SI, EARTH_DELTA_IW
from kamino.mineral_info import *

RERUN = False
MAX_CHEMISTRY_FALLBACKS = 5000

# ── Calibrated constants ──────────────────────────────────────────────────────────────────────
# From experiments/calibrate_earth.py, 2026-08-26, re-run after the hedenbergite adoption (§27)
# invalidated the §22 fit. 25 evaluations, converged on its own tolerance, best cost 0.1077:
#   Earth: converged, T = 294.4 K, pH 7.76, pCO2 694 ppm
#          Na -8.3%, Ca +3.0%, Mg +37.0%, Alk +33.1%, C +32.1%
# The Mg residual is a consequence of §27, not a solver failure: the deleted Fe->fayalite
# exchange had been manufacturing diopside (a Ca source) out of forsterite, so removing it cut
# the rate-weighted Ca supply 9.5% and raised Mg supply 12.3%.
KD_MG_CALIB = 1.394362e-02
K_NA_CALIB  = 4.272026e-03

# ALPHA is NOT from the fit. The Earth calibration cannot identify it -- measured, ocean
# concentrations move <6% across a 41x change -- because Earth is transport-limited (Da >> 1)
# while the land-free worlds these sweeps target are kinetically limited (measured Da ~ 0.005
# over 19 pilot states, 0/19 with Da > 1), where F is linear in alpha.
#
# calibrate_earth.py's own diagnostic puts the ~1 Tmol/yr primary seafloor dissolution anchor
# that weathering.py documents at alpha = 10.4. That is worth CITING as where alpha would sit,
# but it is NOT adopted as the production value, because it censors the data:
#
#     S     Mg/Si   alpha=2   alpha=10
#     0.8   1.25     296.93     279.99   (-16.9 K)
#     1.0   1.25     308.46     298.62   ( -9.8 K)
#     1.0   0.50     339.87     316.75   (-23.1 K)
#     1.2   1.25   out_of_domain  out_of_domain
#
# alpha = 10 does not rescue the hot end (S = 1.2 still leaves the domain at 389 K) and pushes
# the cold end under. CROSS_INSTELLATION starts at 0.50, and a uniform -17 K puts S <= 0.7 at or
# below freezing, removing the cold half of the feedback curve -- the figure this sweep exists
# to produce.
#
# So production runs at alpha = 2 and the alpha ARM carries the argument instead. That is the
# stronger claim anyway: the FEEDBACK STRENGTH is alpha-invariant -- d ln F/dT measured
# 0.078 / 0.081 / 0.075 at alpha = 0.49 / 2 / 20, flat to 7% over a 40x range (degrading only
# near 200, where the system starts leaving the kinetic limit). Showing the Mg/Si and dIW
# orderings unchanged across the arm answers the "alpha is unconstrained" objection directly,
# where picking a single anchored value would not.
ALPHA_CALIB = 2.0
alpha = [2, 10, 50]   # the sensitivity arm; all three stay in the kinetic limit (Da <= 0.13)

def _warn_constant_drift():
    """The sweep used alpha=2, kd=0.02, k_na=0.004 for a year while planet.py shipped kd=0.07 and
    k_na=2.19e-3 -- the two drifted apart silently and every sweep result was attributable only by
    reading this file. Say so loudly rather than let it happen again."""
    for label, here, there in (('alpha', ALPHA_CALIB, ALPHA_REF),
                               ('kd_mg_ht', KD_MG_CALIB, KD_MG_HT),
                               ('k_na_cont_removal', K_NA_CALIB, K_NA_CONT_REMOVAL)):
        if abs(here / there - 1) > 1e-6:
            print(f"  NOTE {label}: sweep uses {here:g}, module default is {there:g} "
                  f"-- runs are tagged accordingly and are NOT comparable to untagged ones.")


# Wall-clock cap per run, seconds. The fallback cap alone does not bound the sweep: a run whose
# chemistry keeps CONVERGING but whose integrator crawls spends no fallback budget and never
# trips chemistry_void either, so it grinds indefinitely. In fast_18 the top 5% of runs consumed
# 57% of all CPU time and the slowest single run took 3.2 hours.
#
# Depth-dependent, because deep oceans are legitimately slower rather than pathological: a 20 km
# ocean carries ~7x the water mass, so its carbon cycle equilibrates over 1.0-1.4 Gyr against
# ~200 Myr for 3 km (measured from the saved trajectories). Capping them at the shallow budget
# would truncate real physics, not waste.
#
# NOTE this is a COST control, not a physics control. Do NOT shorten t_end instead: at 500 Myr a
# 20 km ocean at S = 1.0 reports 328 K against its converged 344 K -- a 16 K error, in exactly
# the regime the water-world results need.
#
# WALL_SECONDS_DEEP raised 1800 -> 2700 after the 20-run pilot: one deep run (S = 1.0, 20 km,
# Mg/Si 1.6, dIW = -2.0) hit the 1800 s cap at t = 1.72 of 2.0 Gyr -- 86% through, and still
# converging. That is a truncation of real physics, which is what these budgets exist to avoid.
# The headroom is cheap now that tau_prec scales with ocean depth (planet.TAU_PREC_REF): the deep
# runs that motivated the cap take 3-6x fewer steps than the pilot's did.
WALL_SECONDS_SHALLOW = 900
WALL_SECONDS_DEEP = 2700
DEEP_OCEAN_M = 10000


def wall_budget(ocean_depth):
    """Wall-clock budget for a run, in seconds. See WALL_SECONDS_* above."""
    return WALL_SECONDS_DEEP if ocean_depth >= DEEP_OCEAN_M else WALL_SECONDS_SHALLOW


# Output directory. The fast_18 path (/data/pt426/...) was on a 26-worker machine that is no
# longer available; default to a sibling of the repo so a sweep runs anywhere.
OUTPUT_PATH = os.environ.get(
    'KAMINO_SWEEP_OUTPUT',
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'sweep_output'))

# Worker processes. fast_18 used 26. Set from the machine unless overridden: PHREEQC is
# CPU-bound and hyperthreads help little, so use physical cores, leaving one free so the box
# stays usable.
WORKERS = int(os.environ.get('KAMINO_SWEEP_WORKERS',
                             max(1, (os.cpu_count() or 2) // 2 - 1)))


def _tag(value, reference, prefix):
    """Name fragment for a swept parameter; empty at its reference so old names are reproduced."""
    return '' if value == reference else f'_{prefix}{value:g}'


def _run_name(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na):
    """Run name. Every parameter that differs from the Planet default MUST appear, or two configs
    would share a filename and RERUN=False would silently return the first one's result."""
    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}'

    if rw:
        run_name += '_rw'

    # The crust axes replace the old `_mt{T_p}` tag. Both are always tagged, so a name from
    # this script can never collide with a legacy `_mt` name from an older sweep.
    run_name += f'_mgsi{mgsi:g}_diw{diw:g}'
    run_name += _tag(alpha, ALPHA_REF, 'a')
    run_name += _tag(kd_mg, KD_MG_HT, 'kmg')
    run_name += _tag(k_na, K_NA_CONT_REMOVAL, 'kna')

    return run_name


def run_simulation(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    run_name = _run_name(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na)

    if not RERUN:
        json_path = os.path.join(output_path, f'{run_name}.json')
        if os.path.exists(json_path):
            try:
                with open(json_path) as fh:
                    existing = json.load(fh)
                if 'termination' in existing:
                    return run_name, None, existing.get('T'), existing['termination']
            except Exception:
                pass  # corrupt/incomplete file — fall through and re-run

    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=c,
            outgassing=o,
            ocean_depth=d,
            reverse_weathering=rw,
            mantle_mg_si=mgsi,
            delta_iw=diw,
            alpha=alpha,
            kd_mg_ht=kd_mg,
            k_na_cont_removal=k_na,
            name=run_name
        )
        p.time_evolve(max_chemistry_fallbacks=MAX_CHEMISTRY_FALLBACKS,
                      max_wall_seconds=wall_budget(d))
        with open(p._output_filename) as fh:  # time_evolve records T and termination here
            result = json.load(fh)
        return run_name, None, result.get('T'), result.get('termination')
    except Exception as e:
        return run_name, str(e), None, None


def cross_combos(instellation, mantle_mg_si, delta_iw, ocean_depth,
                 mg_si_ref=EARTH_MANTLE_MG_SI, diw_ref=EARTH_DELTA_IW,
                 outgassing=0.1, crust=1.0, rw=True,
                 alpha=None, kd_mg=None, k_na=None):
    """Combos for a CROSS design: vary one composition axis at a time through the Earth centre.

    A full factorial over instellation x Mg/Si x dIW is the obvious design and it is mostly
    waste: the three figures the crust-composition results need are the feedback curve, the Mg/Si
    effect and the dIW effect, and each is a one-axis cut through the centre. For 14 instellations
    x 5 Mg/Si x 5 dIW the factorial is 350 runs against this design's 126 -- and the runs the
    factorial adds are the off-centre corners, which no planned figure uses.

    Returns a deduplicated list of combos in run_simulation's argument order, cheapest first
    (see `_cost_rank`) so a sweep that is interrupted or short of time has already covered the
    most ground.
    """
    # Default to the CALIBRATED constants, not the module ones. Defaulting to ALPHA_REF etc.
    # would silently run the cross sweep at the shipped values while every other sweep here used
    # the calibration -- the exact drift `_warn_constant_drift` exists to catch.
    alpha = ALPHA_CALIB if alpha is None else alpha
    kd_mg = KD_MG_CALIB if kd_mg is None else kd_mg
    k_na = K_NA_CALIB if k_na is None else k_na

    seen, combos = set(), []
    for depth in ocean_depth:
        for S in instellation:
            axes = ([(mg, diw_ref) for mg in mantle_mg_si]
                    + [(mg_si_ref, dw) for dw in delta_iw])
            for mg, dw in axes:
                key = (S, outgassing, crust, depth, rw, mg, dw, alpha, kd_mg, k_na)
                if key in seen:
                    continue
                seen.add(key)
                combos.append(key)
    combos.sort(key=_cost_rank)
    return combos


def _cost_rank(combo):
    """Cheap-first ordering. From fast_18's timings, cost rises steeply with instellation up to
    S ~ 1.10 and collapses beyond it (those runs leave the domain in seconds), and deep oceans
    are uniformly slower. Ordering this way means an interrupted sweep loses its most expensive
    runs, not a random slice of the grid."""
    S, _out, _crust, depth = combo[0], combo[1], combo[2], combo[3]
    s_cost = 0.0 if S > 1.12 else S
    return (depth >= DEEP_OCEAN_M, s_cost)


def run_sweep(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering,
              mantle_mg_si, delta_iw, alpha=(ALPHA_REF,), kd_mg=(KD_MG_HT,),
              k_na=(K_NA_CONT_REMOVAL,), output_path=OUTPUT_PATH):

    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    workers = WORKERS

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering, mantle_mg_si, delta_iw, alpha, kd_mg, k_na))
    return run_combos(combos, output_path=output_path)


def run_combos(combos, output_path=OUTPUT_PATH):
    """Execute an explicit list of combos (run_simulation argument order, minus output_path).

    Shared by run_sweep (full factorial) and the cross design, so both get the same collision
    check, resume behaviour and reporting.
    """
    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)
    workers = WORKERS
    alpha = sorted({c[7] for c in combos})
    kd_mg = sorted({c[8] for c in combos})
    k_na = sorted({c[9] for c in combos})
    total = len(combos)

    # Distinct configs must map to distinct filenames, or one silently overwrites the other and
    # RERUN=False then returns the survivor's result for both (the fast_13 resume trap).
    names = [_run_name(*combo) for combo in combos]
    if len(set(names)) != len(names):
        duplicated = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(
            f"{len(names) - len(set(names))} run name collision(s) in this grid, e.g. "
            f"{duplicated[:3]}. Two configs would share an output file."
        )

    cap_str = MAX_CHEMISTRY_FALLBACKS if MAX_CHEMISTRY_FALLBACKS is not None else 'disabled'
    print(f"Running {total} simulations with {workers} worker processes "
          f"(fallback cap: {cap_str})...")
    print(f"Output: {output_path}")
    for label, values in (('alpha', alpha), ('kd_mg_ht', kd_mg), ('k_na_cont_removal', k_na)):
        print(f"  {label}: {list(values)}")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(run_simulation, s, o, c, d, rw, mgsi, diw, a, kmg, kna, output_path): (s, o, c, d, rw, mgsi, diw, a, kmg, kna)
            for s, o, c, d, rw, mgsi, diw, a, kmg, kna in combos
        }

        completed = 0
        aborted = 0
        for future in as_completed(futures):
            completed += 1
            run_name, error, T, termination = future.result()
            if error:
                print(f"[{completed}/{total}] FAILED {run_name}: {error}")
            else:
                if termination == 'fallback_limit':
                    aborted += 1
                T_str = f"{T:.1f} K" if T is not None else "T unknown"
                print(f"[{completed}/{total}] Done: {run_name} ({T_str}, {termination or 'unknown'})")

    print("All simulations complete.")
    if aborted:
        # Not an error -- these are the runs the cap was added to stop. A large fraction here
        # means the chemistry is unhealthy over a wide part of the grid, not that the cap is wrong.
        print(f"{aborted}/{total} run(s) hit the fallback cap ({MAX_CHEMISTRY_FALLBACKS}) "
              f"and were recorded as 'fallback_limit'.")

fast = False

if fast:
    instellation = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    outgassing = [0.03, 0.1, 1, 3]
    crust_production_rate = [0.01, 0.1, 1, 10]
    ocean_depth = [300, 3000, 30000]
    mantle_mg_si = [0.8, 1.25, 1.6]
    delta_iw = [-3.0, -2.0, -1.5]
else:
    instellation = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2]
    outgassing = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    crust_production_rate = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    ocean_depth = [300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000]
    # Crust composition axes. Bounded by the mineralogical limits of Guimond et al. (2024):
    # mantles go olivine-free below Mg/Si ~0.8 and orthopyroxene-free above ~1.6, and dIW > -1
    # puts mantle FeO past the ~25 wt% ceiling where the thermodynamic models fail.
    mantle_mg_si = [0.5, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0]
    delta_iw = [-4.0, -3.0, -2.5, -2.0, -1.5, -1.0]


reverse_weathering = [False, True]
k_mg = [KD_MG_CALIB, 0]   # 0 disables the Mg->Ca exchange entirely
k_na = [K_NA_CALIB, 0]    # 0 disables the Na sink entirely

crust_production_rate_default = [1]
outgassing_default = [0.1]
ocean_depth_default = [3000]
# The water-world counterpart of ocean_depth_default, used by the *_deep sweeps. 20 km matches
# CROSS_DEPTHS so deep results from either design are directly comparable.
ocean_depth_deep_default = [20000]
reverse_weathering_default = [True]
alpha_default = [ALPHA_CALIB]
k_mg_default = [KD_MG_CALIB]
k_na_default = [K_NA_CALIB]
mantle_mg_si_default = [EARTH_MANTLE_MG_SI]
delta_iw_default = [EARTH_DELTA_IW]


# ── The crust-composition cross sweep ──────────────────────────────────────────
# Axes for the three figures the crust-composition results need: the feedback curve, the Mg/Si
# effect and the dIW effect.
#
# Instellation stops at 1.15: fast_18 put every run at S >= 1.15 out_of_domain within seconds,
# so they cost nothing but say nothing either. It starts at 0.50 because below that the planet
# is frozen and also leaves the domain.
CROSS_INSTELLATION = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05,
                      1.10, 1.15]
# Mg/Si stops at 1.6: above it every composition is mass-violating in the norm and the melts are
# ultracalcic (development history 25.4). 0.5 is included deliberately -- it is the quartz-crust
# corner, and whether an inert crust breaks the feedback is a result either way.
CROSS_MG_SI = [0.5, 0.8, 1.0, EARTH_MANTLE_MG_SI, 1.6]
# dIW spans mantle FeO 0.26 wt% (enstatite-chondrite-like) to 24 wt% (the ~25 wt% model ceiling).
CROSS_DELTA_IW = [-5.0, -4.0, -3.0, EARTH_DELTA_IW, -1.0]
# 3 km is the Earth-like reference; 20 km is the water world. Deep runs get the longer wall
# budget because they equilibrate over 1.0-1.4 Gyr rather than ~200 Myr.
CROSS_DEPTHS = [3000, 20000]


def run_cross(depths=CROSS_DEPTHS, output_path=OUTPUT_PATH):
    """The cross sweep. Cheapest runs first, so an interrupted run still covers the most ground."""
    combos = cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW, depths)
    return run_combos(combos, output_path=output_path)


# ── Named sweeps ──────────────────────────────────────────────────────────────────────────────
# Each entry is (description, callable). Select with KAMINO_SWEEPS, comma-separated:
#
#     KAMINO_SWEEPS=composition,cross python experiments/parameter_sweep.py
#
# This replaces the previous comment-out-to-choose block. That pattern meant the file had to be
# edited to run anything, only one sweep was ever active, and which one had run was recoverable
# only from the git history of this file.

def sweep_basic(output_path=OUTPUT_PATH):
    """Instellation x outgassing x crust production, at the Earth-reference crust."""
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, output_path=output_path)


def sweep_basic_deep(output_path=OUTPUT_PATH):
    """As sweep_basic, on the 20 km water world.

    This is the expensive one. Deep runs cost ~5.7x shallow (pilot: 154.6 min for 10 deep against
    27.2 min for 10 shallow), partly recovered by the depth-scaled tau_prec (§26, 3-6x fewer
    steps). At 931 combos it is still the largest sweep here by CPU time -- cost it with the
    preamble's estimate before launching, and consider trimming `outgassing` or
    `crust_production_rate` rather than running the full factorial.
    """
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_deep_default,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, output_path=output_path)


def sweep_depth(output_path=OUTPUT_PATH):
    """Instellation x ocean depth."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, output_path=output_path)


def sweep_composition(output_path=OUTPUT_PATH):
    """Instellation x Mg/Si x dIW -- the full composition factorial at 3 km.

    This is the factorial the CROSS design deliberately avoids (see cross_combos). Use it when
    the off-centre corners matter -- i.e. when the question is whether Mg/Si and dIW interact,
    rather than what each does on its own through the Earth centre.
    """
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si, delta_iw,
                     alpha_default, k_mg_default, k_na_default, output_path=output_path)


def sweep_composition_deep(output_path=OUTPUT_PATH):
    """As sweep_composition, on the 20 km water world."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_deep_default, reverse_weathering_default, mantle_mg_si, delta_iw,
                     alpha_default, k_mg_default, k_na_default, output_path=output_path)


def sweep_alpha(output_path=OUTPUT_PATH):
    """The alpha sensitivity arm at the Earth-reference crust. alpha is a CHOICE (see
    ALPHA_CALIB), so any result sensitive to the absolute CO2 level needs this reported."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha, k_mg_default, k_na_default, output_path=output_path)


def sweep_alpha_composition(output_path=OUTPUT_PATH):
    """alpha x Mg/Si and alpha x dIW: does the composition signal survive the alpha choice?

    This is the sweep that answers the referee question directly -- if the Mg/Si and dIW
    orderings are the same at alpha = 2, 10 and 50, the conclusions do not rest on alpha.
    """
    combos = []
    for a in alpha:
        combos += cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                               ocean_depth_default, alpha=a,
                               kd_mg=k_mg_default[0], k_na=k_na_default[0])
    return run_combos(combos, output_path=output_path)


def sweep_chemistry(output_path=OUTPUT_PATH):
    """kd_mg_ht and k_na on/off, to isolate what each sink contributes."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha_default, k_mg, k_na, output_path=output_path)


def sweep_cross(output_path=OUTPUT_PATH):
    """The cross design at 3 km: one composition axis at a time through the Earth centre.

    Split from the deep half deliberately. The two halves cost very differently (126 runs at
    ~2.7 min against 126 at ~15 min) and answer different questions -- 3 km is the Earth-like
    reference the Mg/Si and dIW figures are built on, 20 km is the water-world case. Running the
    shallow half alone secures those figures in a few hours; the deep half can follow.
    """
    return run_cross(depths=ocean_depth_default, output_path=output_path)


def sweep_cross_deep(output_path=OUTPUT_PATH):
    """The cross design on the 20 km water world. See sweep_cross."""
    return run_cross(depths=ocean_depth_deep_default, output_path=output_path)


SWEEPS = {
    'basic':             ('instellation x outgassing x crust production, 3 km', sweep_basic),
    'basic_deep':        ('instellation x outgassing x crust production, 20 km', sweep_basic_deep),
    'depth':             ('instellation x ocean depth', sweep_depth),
    'composition':       ('instellation x Mg/Si x dIW factorial, 3 km', sweep_composition),
    'composition_deep':  ('instellation x Mg/Si x dIW factorial, 20 km', sweep_composition_deep),
    'cross':             ('cross design: Mg/Si and dIW cuts, 3 km', sweep_cross),
    'cross_deep':        ('cross design: Mg/Si and dIW cuts, 20 km', sweep_cross_deep),
    'alpha':             ('alpha sensitivity arm', sweep_alpha),
    'alpha_composition': ('alpha x composition -- does the signal survive alpha?',
                          sweep_alpha_composition),
    'chemistry':         ('kd_mg_ht / k_na on-off', sweep_chemistry),
}

DEFAULT_SWEEPS = 'cross'


# Measured per-run wall cost, from the 20-run pilot (2026-08-25): 27.2 min for 10 shallow runs,
# 154.6 min for 10 deep. The deep figure PREDATES the depth-scaled tau_prec (§26), which cut deep
# step counts 3-6x, so treat the deep estimate as an upper bound.
MINUTES_PER_RUN_SHALLOW = 2.72
MINUTES_PER_RUN_DEEP = 15.46


def _sweep_shape(name):
    """(n_shallow, n_deep) for a sweep, so cost can be estimated before launching."""
    n = _sweep_size(name)
    if name in ('basic_deep', 'composition_deep', 'cross_deep'):
        return 0, n
    if name == 'depth':
        deep_frac = sum(1 for d in ocean_depth if d >= DEEP_OCEAN_M)
        return len(instellation) * (len(ocean_depth) - deep_frac), len(instellation) * deep_frac
    return n, 0


def _sweep_cost_hours(name):
    """Estimated CPU-hours. Divide by WORKERS for wall time."""
    ns, nd = _sweep_shape(name)
    return (ns * MINUTES_PER_RUN_SHALLOW + nd * MINUTES_PER_RUN_DEEP) / 60.0


def _sweep_size(name):
    """Run count without executing anything, so a sweep can be costed before it is launched."""
    sizers = {
        'basic':            len(instellation)*len(outgassing)*len(crust_production_rate),
        'basic_deep':       len(instellation)*len(outgassing)*len(crust_production_rate),
        'depth':            len(instellation)*len(ocean_depth),
        'composition':      len(instellation)*len(mantle_mg_si)*len(delta_iw),
        'composition_deep': len(instellation)*len(mantle_mg_si)*len(delta_iw),
        'cross':            len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                                             ocean_depth_default)),
        'cross_deep':       len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                                             ocean_depth_deep_default)),
        'alpha':            len(instellation)*len(alpha),
        'alpha_composition': sum(len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI,
                                                  CROSS_DELTA_IW, ocean_depth_default, alpha=a))
                                 for a in alpha),
        'chemistry':        len(instellation)*len(k_mg)*len(k_na),
    }
    return sizers.get(name, 0)


if __name__ == "__main__":
    requested = [n.strip() for n in
                 os.environ.get('KAMINO_SWEEPS', DEFAULT_SWEEPS).split(',') if n.strip()]
    unknown = [n for n in requested if n not in SWEEPS]
    if unknown:
        raise SystemExit(f"Unknown sweep(s) {unknown}. Available: {sorted(SWEEPS)}")

    print(f"Sweeps requested: {requested}   (set KAMINO_SWEEPS to change)")
    print(f"  alpha={ALPHA_CALIB:g}  kd_mg_ht={KD_MG_CALIB:g}  k_na={K_NA_CALIB:g}")
    _warn_constant_drift()
    total = sum(_sweep_size(n) for n in requested)
    total_h = sum(_sweep_cost_hours(n) for n in requested)
    print(f"  {'sweep':19s} {'runs':>5s} {'CPU-h':>7s} {'wall-h':>7s}   description")
    for n in requested:
        h = _sweep_cost_hours(n)
        print(f"  {n:19s} {_sweep_size(n):5d} {h:7.1f} {h/max(WORKERS,1):7.1f}   {SWEEPS[n][0]}")
    print(f"  {'TOTAL':19s} {total:5d} {total_h:7.1f} {total_h/max(WORKERS,1):7.1f}"
          f"   at {WORKERS} workers")
    print("  (deep estimate predates the depth-scaled tau_prec and is an upper bound; see §26)")

    for n in requested:
        print(f"\n{'='*78}\n  {n}: {SWEEPS[n][0]}\n{'='*78}")
        SWEEPS[n][1]()
