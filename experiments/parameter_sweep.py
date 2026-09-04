import itertools
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet as p2
from kamino.planet import Planet, KD_MG_HT, K_NA_CONT_REMOVAL, PE_DEFAULT
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

# alpha is STILL NOT identified by the Earth fit, even though it is now the same number the
# fit reports. Measured, ocean concentrations move <6% across a 41x change in alpha -- because
# Earth is transport-limited (Da >> 1) while the land-free worlds these sweeps target are
# kinetically limited (measured Da ~ 0.005 over 19 pilot states, 0/19 with Da > 1), where F is
# linear in alpha. So ALPHA_REF is a value the joint least-squares solver happened to land on
# while fitting Na/Ca/Mg (which alpha barely affects), not a measurement of alpha itself -- see
# development_history.md section 28.2.
#
# 2026-09-01: production now runs at ALPHA_REF rather than a separately pinned round number, so
# the main sweeps and the module default can never drift apart again (the drift check in
# _warn_constant_drift becomes a tautology for alpha specifically, which is the point).
#
# Re-ran the domain-coverage check (section 28.2) with ALPHA_REF (~1.10015) added as a column:
#
#     S     Mg/Si   alpha=2    alpha=1.10   alpha=10
#     0.8   1.25    298.33 K   304.49 K     280.01 K
#     1.0   1.25    321.31 K   327.22 K     310.35 K
#     1.0   0.50    346.11 K   348.34 K     335.99 K
#     1.2   1.25    out_of_domain (all three, alpha-independent)
#
# ALPHA_REF is WARMER than alpha=2 at every point (lower alpha -> weaker weathering -> less
# cooling), so nothing that stayed in-domain at alpha=2 falls out at ALPHA_REF; the S=1.2 wall is
# unrelated to alpha entirely. This move relaxes the cold-end constraint, it does not tighten it.
ALPHA_CALIB = ALPHA_REF
alpha = [ALPHA_REF, 10, 50]   # the sensitivity arm; all three stay in the kinetic limit (Da <= 0.13)

# ── Ocean redox ───────────────────────────────────────────────────────────────────────────────
# Every sweep is run under BOTH redox states, because the model has no basis for preferring one
# and the difference is large (+11.5 K at S = 1.0, 3 km, Mg/Si 1.25 -- see development history
# section 28.4).
#
#   PE_REDUCING  -3.0  an abiotic planet, which is what this model actually is: no oxygenic
#                      photosynthesis, so no free O2. Ferrous iron stays dissolved and Siderite
#                      (FeCO3) is the iron sink. Matches anoxic marine pore water (-3 to -5) and
#                      Archean ocean reconstructions (-3 to 0).
#   PE_OXIDISING  +4.0 an oxygenated ocean. Ferric Goethite is supersaturated and strips every
#                      mole of dissolved iron. This is ALSO the value PHREEQC silently defaulted
#                      to before pe was a parameter, so it reproduces every result on disk from
#                      before 2026-08-27 -- which is why it is the comparison arm rather than
#                      just an alternative.
#
# Modern oxic seawater is nearer pe = +12.5, but +4 is used for the oxidising arm because it is
# the value the pre-2026-08-27 sweeps implicitly ran at; the iron system is saturated by then
# anyway (Goethite SI +7.64 at pe 4 against +7.65 at pe 12).
#
# The effect SATURATES below pe ~ 0 (pe of 0, -3 and -6 all give 321.3 K), so this is a binary
# oxic/anoxic distinction rather than a continuum that needs resolving. pe = -3 is also the
# numerically cleanest of the reducing values (0 chemistry fallbacks against 9 at pe = -6).
PE_REDUCING = -3.0
PE_OXIDISING = 4.0
PE_DEFAULT_SWEEP = PE_REDUCING          # the physically-motivated arm for an abiotic planet
PE_STATES = [PE_REDUCING, PE_OXIDISING]


def _pe_label(pe):
    """Short name for a redox state, for logs and figure captions."""
    if pe is None:
        return 'PHREEQC default (+4, oxidising)'
    return 'reducing' if pe < 2.0 else 'oxidising'

# The pe sensitivity arm, for the one sweep that resolves the axis rather than bracketing it.
# Spans oxic seawater to below the Goethite saturation boundary (~ -5.8).
pe_arm = [12.0, 4.0, 0.0, -3.0, -6.0]

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


def _run_name(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na, pe=PE_DEFAULT_SWEEP):
    """Run name. Every parameter that differs from the Planet default MUST appear, or two configs
    would share a filename and RERUN=False would silently return the first one's result.

    `pe` is tagged against planet.PE_DEFAULT rather than against this module's sweep default, so
    a run at the model's own default keeps an untagged name and every other redox state is
    distinguishable. Without this the oxidising and reducing arms of the SAME sweep would share a
    filename and one would silently overwrite the other.
    """
    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}'

    if rw:
        run_name += '_rw'

    # The crust axes replace the old `_mt{T_p}` tag. Both are always tagged, so a name from
    # this script can never collide with a legacy `_mt` name from an older sweep.
    run_name += f'_mgsi{mgsi:g}_diw{diw:g}'
    run_name += _tag(alpha, ALPHA_REF, 'a')
    run_name += _tag(kd_mg, KD_MG_HT, 'kmg')
    run_name += _tag(k_na, K_NA_CONT_REMOVAL, 'kna')
    run_name += _tag(pe, PE_DEFAULT, 'pe')

    return run_name


def run_simulation(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na, pe, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    run_name = _run_name(s, o, c, d, rw, mgsi, diw, alpha, kd_mg, k_na, pe)

    if not RERUN:
        json_path = os.path.join(output_path, f'{run_name}.json')
        if os.path.exists(json_path):
            try:
                with open(json_path) as fh:
                    existing = json.load(fh)
                # Resume guard. A run at the model default is deliberately UNTAGGED, so a file
                # written before `pe` existed has the same name as a reducing run -- but it was
                # produced at PHREEQC's implicit pe = 4, i.e. OXIDISING. Reusing it would label
                # 2000 oxidising runs as reducing. Any file whose recorded pe does not match what
                # was asked for is re-run rather than trusted (the fast_13 resume trap, again).
                stored_pe = existing.get('pe', 'ABSENT')
                if stored_pe == 'ABSENT':
                    stale = pe is not None      # pre-pe output: only valid if pe was unset
                else:
                    stale = not (stored_pe is None and pe is None) and stored_pe != pe
                if stale:
                    print(f"  re-running {run_name}: stored pe={stored_pe!r} != requested "
                          f"pe={pe!r} (output predates the pe parameter?)", flush=True)
                elif 'termination' in existing:
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
            pe=pe,
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
                 alpha=None, kd_mg=None, k_na=None, pe=PE_DEFAULT_SWEEP):
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
                key = (S, outgassing, crust, depth, rw, mg, dw, alpha, kd_mg, k_na, pe)
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
              k_na=(K_NA_CONT_REMOVAL,), pe=(PE_DEFAULT_SWEEP,), output_path=OUTPUT_PATH):

    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    workers = WORKERS

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth,
                                   reverse_weathering, mantle_mg_si, delta_iw, alpha, kd_mg,
                                   k_na, pe))
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
    pe = sorted({c[10] for c in combos})
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
    for label, values in (('alpha', alpha), ('kd_mg_ht', kd_mg), ('k_na_cont_removal', k_na),
                          ('pe', [f'{v:g} ({_pe_label(v)})' for v in pe])):
        print(f"  {label}: {list(values)}")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(run_simulation, s, o, c, d, rw, mgsi, diw, a, kmg, kna, pe_,
                            output_path): (s, o, c, d, rw, mgsi, diw, a, kmg, kna, pe_)
            for s, o, c, d, rw, mgsi, diw, a, kmg, kna, pe_ in combos
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

# The two Mg/Si end-members the basic sweep is repeated at, bracketing Earth's 1.25.
# 0.8 is the low limit of the composition axes: below it the mantle is olivine-free and the crust
# quartz-normative (Guimond et al. 2024; development history 25.4).
# 1.8 sits at mantle (Mg+Fe)/Si = 1.98, PAST the ~1.69 ceiling 25.4 measured -- but that was
# measured before Akermanite closed the norm (25.5); at 1.8 the assemblage now sums to 1.0 with
# no mass-balance warning, Akermanite taking 11 wt%. State this if these runs are published.
mantle_mg_si_low = [0.8]
mantle_mg_si_high = [1.8]


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


def run_cross(depths=CROSS_DEPTHS, output_path=OUTPUT_PATH, pe=PE_STATES):
    """The cross sweep. Cheapest runs first, so an interrupted run still covers the most ground."""
    combos = []
    for pe_ in ([pe] if isinstance(pe, (int, float)) else pe):
        combos += cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW, depths, pe=pe_)
    combos.sort(key=_cost_rank)
    return run_combos(combos, output_path=output_path)


# ── Named sweeps ──────────────────────────────────────────────────────────────────────────────
# Each entry is (description, callable). Select with KAMINO_SWEEPS, comma-separated:
#
#     KAMINO_SWEEPS=composition,cross python experiments/parameter_sweep.py
#
# This replaces the previous comment-out-to-choose block. That pattern meant the file had to be
# edited to run anything, only one sweep was ever active, and which one had run was recoverable
# only from the git history of this file.

def sweep_basic(output_path=OUTPUT_PATH, pe=PE_STATES):
    """Instellation x outgassing x crust production, at the Earth-reference crust."""
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_basic_deep(output_path=OUTPUT_PATH, pe=PE_STATES):
    """As sweep_basic, on the 20 km water world.

    This is the expensive one. Deep runs cost ~5.7x shallow (pilot: 154.6 min for 10 deep against
    27.2 min for 10 shallow), partly recovered by the depth-scaled tau_prec (§26, 3-6x fewer
    steps). At 931 combos it is still the largest sweep here by CPU time -- cost it with the
    preamble's estimate before launching, and consider trimming `outgassing` or
    `crust_production_rate` rather than running the full factorial.
    """
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_deep_default,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_basic_low_mgsi(output_path=OUTPUT_PATH, pe=PE_STATES):
    """As sweep_basic, on the low-Mg/Si (0.8) crust -- the olivine-free, silica-rich end-member."""
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default,
                     reverse_weathering_default, mantle_mg_si_low, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_basic_high_mgsi(output_path=OUTPUT_PATH, pe=PE_STATES):
    """As sweep_basic, on the high-Mg/Si (1.8) crust -- the orthopyroxene-free, olivine-rich end.

    Paired with sweep_basic_low_mgsi, this is the outgassing x crust-production plane repeated at
    both ends of the composition axis, so the feedback strength can be read as a function of
    crust chemistry rather than only along the one-axis cut the cross design takes.
    """
    return run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default,
                     reverse_weathering_default, mantle_mg_si_high, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_depth(output_path=OUTPUT_PATH, pe=PE_STATES):
    """Instellation x ocean depth."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth,
                     reverse_weathering_default, mantle_mg_si_default, delta_iw_default,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_composition(output_path=OUTPUT_PATH, pe=PE_STATES):
    """Instellation x Mg/Si x dIW -- the full composition factorial at 3 km.

    This is the factorial the CROSS design deliberately avoids (see cross_combos). Use it when
    the off-centre corners matter -- i.e. when the question is whether Mg/Si and dIW interact,
    rather than what each does on its own through the Earth centre.
    """
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si, delta_iw,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_composition_deep(output_path=OUTPUT_PATH, pe=PE_STATES):
    """As sweep_composition, on the 20 km water world."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_deep_default, reverse_weathering_default, mantle_mg_si, delta_iw,
                     alpha_default, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_basic_oxidised(output_path=OUTPUT_PATH):
    """As sweep_basic, oxidising arm (pe = +4) only.

    basic/composition/depth already cover both redox states by default (§31), so this is not a
    new axis -- it lets the oxidising arm be launched on its own, e.g. split across a machine from
    the reducing arm, rather than paying for both every time the sweep is (re)started.
    """
    return sweep_basic(output_path=output_path, pe=(PE_OXIDISING,))


def sweep_composition_oxidised(output_path=OUTPUT_PATH):
    """As sweep_composition, oxidising arm (pe = +4) only. See sweep_basic_oxidised."""
    return sweep_composition(output_path=output_path, pe=(PE_OXIDISING,))


def sweep_depth_oxidised(output_path=OUTPUT_PATH):
    """As sweep_depth, oxidising arm (pe = +4) only. See sweep_basic_oxidised."""
    return sweep_depth(output_path=output_path, pe=(PE_OXIDISING,))


def sweep_alpha(output_path=OUTPUT_PATH, pe=PE_STATES):
    """The alpha sensitivity arm at the Earth-reference crust. alpha is a CHOICE (see
    ALPHA_CALIB), so any result sensitive to the absolute CO2 level needs this reported."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha, k_mg_default, k_na_default, pe=pe, output_path=output_path)


def sweep_alpha_composition(output_path=OUTPUT_PATH, pe=PE_STATES):
    """alpha x Mg/Si and alpha x dIW: does the composition signal survive the alpha choice?

    This is the sweep that answers the referee question directly -- if the Mg/Si and dIW
    orderings are the same at alpha = ALPHA_REF, 10 and 50, the conclusions do not rest on alpha.
    """
    combos = []
    for pe_ in pe:
        for a in alpha:
            combos += cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                                   ocean_depth_default, alpha=a,
                                   kd_mg=k_mg_default[0], k_na=k_na_default[0], pe=pe_)
    return run_combos(combos, output_path=output_path)


def sweep_chemistry(output_path=OUTPUT_PATH, pe=PE_STATES):
    """kd_mg_ht and k_na on/off, to isolate what each sink contributes."""
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha_default, k_mg, k_na, pe=pe, output_path=output_path)


def sweep_pe(output_path=OUTPUT_PATH, pe=None):
    """The redox sensitivity arm: instellation x pe at the Earth-reference crust.

    Resolves the axis that every other sweep only brackets. Spans modern oxic seawater (+12) to
    below the Goethite saturation boundary (~ -5.8), so it shows both the plateau and where the
    iron sink actually switches from ferric Goethite to ferrous Siderite (around pe = 0).

    The effect is expected to SATURATE below pe ~ 0 -- measured at one point, pe of 0, -3 and -6
    all give 321.3 K -- so the interesting structure is between +12 and 0, and the reducing tail
    is there to confirm the plateau rather than to resolve it.
    """
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha_default, k_mg_default, k_na_default,
                     pe=(pe_arm if pe is None else pe), output_path=output_path)


def sweep_pe_deep(output_path=OUTPUT_PATH, pe=None):
    """As sweep_pe, on the 20 km water world.

    Worth running separately: the deep ocean is where the carbonate sink is pressure-suppressed
    (section 29.3), and Siderite -- the reducing-arm iron sink -- is a CARBONATE, so the redox
    switch and the depth effect are not obviously independent.
    """
    return run_sweep(instellation, outgassing_default, crust_production_rate_default,
                     ocean_depth_deep_default, reverse_weathering_default, mantle_mg_si_default,
                     delta_iw_default, alpha_default, k_mg_default, k_na_default,
                     pe=(pe_arm if pe is None else pe), output_path=output_path)


def sweep_pe_composition(output_path=OUTPUT_PATH, pe=None):
    """pe x Mg/Si and pe x dIW: does the composition signal survive the redox choice?

    The dIW half is the one that matters. Section 29.2 measured dIW as a weak, phase-boundary-gated
    control, but that was at pe = 4 where dissolved iron is stripped by Goethite and so cannot
    affect the carbon cycle at all. Under anoxia iron stays in solution, and dIW is the parameter
    that sets how much of it there is -- so the dIW result may be substantially redox-dependent.
    """
    combos = []
    for pe_ in (pe_arm if pe is None else pe):
        combos += cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                               ocean_depth_default, pe=pe_)
    combos.sort(key=_cost_rank)
    return run_combos(combos, output_path=output_path)


def sweep_cross(output_path=OUTPUT_PATH, pe=PE_STATES):
    """The cross design at 3 km: one composition axis at a time through the Earth centre.

    Split from the deep half deliberately. The two halves cost very differently (126 runs at
    ~2.7 min against 126 at ~15 min) and answer different questions -- 3 km is the Earth-like
    reference the Mg/Si and dIW figures are built on, 20 km is the water-world case. Running the
    shallow half alone secures those figures in a few hours; the deep half can follow.
    """
    return run_cross(depths=ocean_depth_default, output_path=output_path, pe=pe)


def sweep_cross_deep(output_path=OUTPUT_PATH, pe=PE_STATES):
    """The cross design on the 20 km water world. See sweep_cross."""
    return run_cross(depths=ocean_depth_deep_default, output_path=output_path, pe=pe)


SWEEPS = {
    'basic':             ('instellation x outgassing x crust production, 3 km', sweep_basic),
    'basic_deep':        ('instellation x outgassing x crust production, 20 km', sweep_basic_deep),
    'basic_low_mgsi':    ('basic at Mg/Si = 0.8, 3 km', sweep_basic_low_mgsi),
    'basic_high_mgsi':   ('basic at Mg/Si = 1.8, 3 km', sweep_basic_high_mgsi),
    'depth':             ('instellation x ocean depth', sweep_depth),
    'composition':       ('instellation x Mg/Si x dIW factorial, 3 km', sweep_composition),
    'composition_deep':  ('instellation x Mg/Si x dIW factorial, 20 km', sweep_composition_deep),
    'basic_oxidised':    ('basic, oxidising arm only', sweep_basic_oxidised),
    'composition_oxidised': ('composition, oxidising arm only', sweep_composition_oxidised),
    'depth_oxidised':    ('depth, oxidising arm only', sweep_depth_oxidised),
    'cross':             ('cross design: Mg/Si and dIW cuts, 3 km', sweep_cross),
    'cross_deep':        ('cross design: Mg/Si and dIW cuts, 20 km', sweep_cross_deep),
    'alpha':             ('alpha sensitivity arm', sweep_alpha),
    'alpha_composition': ('alpha x composition -- does the signal survive alpha?',
                          sweep_alpha_composition),
    'pe':                ('redox sensitivity arm, 3 km', sweep_pe),
    'pe_deep':           ('redox sensitivity arm, 20 km', sweep_pe_deep),
    'pe_composition':    ('pe x composition -- does the signal survive redox?',
                          sweep_pe_composition),
    'chemistry':         ('kd_mg_ht / k_na on-off', sweep_chemistry),
}

DEFAULT_SWEEPS = 'basic_high_mgsi,basic_low_mgsi'


# Measured per-run wall cost, from the 20-run pilot (2026-08-25): 27.2 min for 10 shallow runs,
# 154.6 min for 10 deep. The deep figure PREDATES the depth-scaled tau_prec (§26), which cut deep
# step counts 3-6x, so treat the deep estimate as an upper bound.
MINUTES_PER_RUN_SHALLOW = 2.72
MINUTES_PER_RUN_DEEP = 15.46


def _sweep_shape(name):
    """(n_shallow, n_deep) for a sweep, so cost can be estimated before launching."""
    n = _sweep_size(name)
    if name in ('basic_deep', 'composition_deep', 'cross_deep', 'pe_deep'):
        return 0, n
    if name in ('depth', 'depth_oxidised'):
        deep = sum(1 for d in ocean_depth if d >= DEEP_OCEAN_M)
        per_state = len(instellation)
        states = 1 if name == 'depth_oxidised' else len(PE_STATES)
        return (per_state * (len(ocean_depth) - deep) * states,
                per_state * deep * states)
    return n, 0


def _sweep_cost_hours(name):
    """Estimated CPU-hours. Divide by WORKERS for wall time."""
    ns, nd = _sweep_shape(name)
    return (ns * MINUTES_PER_RUN_SHALLOW + nd * MINUTES_PER_RUN_DEEP) / 60.0


def _sweep_size(name):
    """Run count without executing anything, so a sweep can be costed before it is launched."""
    n_redox, n_pe = len(PE_STATES), len(pe_arm)
    n_cross = len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                               ocean_depth_default))
    # Every sweep except the pe arms runs under BOTH redox states, so it doubles. The pe sweeps
    # resolve that axis themselves and are NOT doubled.
    sizers = {
        'basic':            len(instellation)*len(outgassing)*len(crust_production_rate)*n_redox,
        'basic_deep':       len(instellation)*len(outgassing)*len(crust_production_rate)*n_redox,
        'basic_low_mgsi':   len(instellation)*len(outgassing)*len(crust_production_rate)*n_redox,
        'basic_high_mgsi':  len(instellation)*len(outgassing)*len(crust_production_rate)*n_redox,
        'depth':            len(instellation)*len(ocean_depth)*n_redox,
        'composition':      len(instellation)*len(mantle_mg_si)*len(delta_iw)*n_redox,
        'composition_deep': len(instellation)*len(mantle_mg_si)*len(delta_iw)*n_redox,
        'basic_oxidised':   len(instellation)*len(outgassing)*len(crust_production_rate),
        'composition_oxidised': len(instellation)*len(mantle_mg_si)*len(delta_iw),
        'depth_oxidised':   len(instellation)*len(ocean_depth),
        'cross':            n_cross*n_redox,
        'cross_deep':       len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI, CROSS_DELTA_IW,
                                             ocean_depth_deep_default))*n_redox,
        'alpha':            len(instellation)*len(alpha)*n_redox,
        'alpha_composition': sum(len(cross_combos(CROSS_INSTELLATION, CROSS_MG_SI,
                                                  CROSS_DELTA_IW, ocean_depth_default, alpha=a))
                                 for a in alpha)*n_redox,
        'chemistry':        len(instellation)*len(k_mg)*len(k_na)*n_redox,
        'pe':               len(instellation)*n_pe,
        'pe_deep':          len(instellation)*n_pe,
        'pe_composition':   n_cross*n_pe,
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
