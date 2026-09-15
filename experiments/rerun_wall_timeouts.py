"""Re-run every `wall_timeout` run in a sweep directory to completion.

A `wall_timeout` is the only termination in this model that is purely an artefact of the cost
cap: the run's chemistry was converging and its integrator was still stepping when
`max_wall_seconds` fired (see `Planet.time_evolve` / `WallClockLimitExceeded`). It is not a
statement about physics, and unlike `fallback_limit` or `chemistry_void` there is nothing wrong
with the state -- it just needs more wall clock. `timeout` (reaching t_end) is a settled run and
is NOT re-run here.

Usage
-----
    python experiments/rerun_wall_timeouts.py --path /data/pt426/sweep_output
    python experiments/rerun_wall_timeouts.py --path ... --dry-run     # list, run nothing
    python experiments/rerun_wall_timeouts.py --path ... --wall 43200  # 12 h per run

History: this file used to be a one-off carrying a hardcoded filter for the runs feeding the
ratio figures (development_history.md section 33.11). It is now general -- every `wall_timeout`
in the directory -- with the subset behaviour available through `--filter`.
"""
import argparse
import glob
import json
import os
import re
import sys

# ── Environment MUST be set before importing parameter_sweep ──────────────────────────────────
# `RERUN`, `WALL_SECONDS_SHALLOW` and `WALL_SECONDS_DEEP` are read at MODULE IMPORT time, and
# ProcessPoolExecutor spawns workers that re-import the module. A value assigned to
# `ps.RERUN` after import would therefore never reach a worker: the parent would think it was
# re-running while every child silently reused the very output it was told to replace, and the
# script would report success having recomputed nothing. The environment is inherited across
# spawn, so it is the only channel that works -- hence this block sitting above the imports,
# which is not an accident and should not be tidied downwards.
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument('--wall', type=int, default=3600*6)
_pre.add_argument('--wall-deep', type=int, default=None)
_known, _ = _pre.parse_known_args()
_WALL = _known.wall
_WALL_DEEP = _known.wall_deep if _known.wall_deep is not None else _WALL

os.environ['KAMINO_RERUN'] = '1'
os.environ.setdefault('KAMINO_WALL_SHALLOW', str(_WALL))
os.environ.setdefault('KAMINO_WALL_DEEP', str(_WALL_DEEP))

import numpy as np  # noqa: E402

import parameter_sweep as ps  # noqa: E402
import continental_baseline as cb  # noqa: E402


def _int(v):
    """Ints where the run name needs them, so a rebuilt combo reproduces its original filename.

    `_run_name` interpolates outgassing / crust production / depth with plain str(), so 1 and 1.0
    give different filenames. Getting this wrong does not error -- it writes a NEW file alongside
    the original and leaves the wall_timeout in place, which looks like the re-run silently
    failing to fix anything.
    """
    return int(v) if float(v).is_integer() else float(v)


def wall_timeout_combos(path, filter_expr=None):
    """Every wall_timeout run in `path`, as combos in `run_simulation` argument order.

    Reads the JSONs directly rather than going through `plot_results.load_data`, which parses
    every trajectory to derive salinity -- unnecessary here and slow over thousands of files.
    The termination is grepped first so only the matching files are actually parsed.
    """
    combos, skipped = [], []
    for f in sorted(glob.glob(os.path.join(path, '*.json'))):
        # Cheap pre-filter: avoid json.load on files that cannot match.
        with open(f) as fh:
            head = fh.read()
        if '"termination"' not in head or 'wall_timeout' not in head:
            continue
        try:
            d = json.loads(head)
        except json.JSONDecodeError:
            skipped.append((os.path.basename(f), 'corrupt JSON'))
            continue
        if d.get('termination') != 'wall_timeout':
            continue          # 'wall_timeout' appeared in termination_raw or a message, not here
        try:
            combo = (
                float(d['instellation']),
                _int(d['outgassing']),
                _int(d['crust_production_rate']),
                _int(d['ocean_depth']),
                bool(d.get('reverse_weathering', False)),
                float(d['mantle_mg_si']),
                float(d['delta_iw']),
                float(d['alpha']),
                float(d['kd_mg_ht']),
                float(d['k_na_cont_removal']),
                float(d['pe']) if d.get('pe') is not None else ps.PE_DEFAULT,
                float(d.get('land_fraction', 0.0)),
            )
        except KeyError as e:
            skipped.append((os.path.basename(f), f'missing {e}'))
            continue
        if filter_expr and not eval(filter_expr, {'np': np}, dict(  # noqa: S307 - operator-supplied
                s=combo[0], out=combo[1], crust=combo[2], depth=combo[3], rw=combo[4],
                mgsi=combo[5], diw=combo[6], alpha=combo[7], kd_mg=combo[8], k_na=combo[9],
                pe=combo[10], land=combo[11])):
            continue
        combos.append(combo)
    return combos, skipped


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--path', default=ps.OUTPUT_PATH, help='Sweep directory to scan and rewrite.')
    p.add_argument('--wall', type=int, default=86400,
                   help='Per-run wall budget in seconds (default 86400 = 24 h).')
    p.add_argument('--wall-deep', type=int, default=None,
                   help='Separate budget for ocean_depth >= 10 km (default: same as --wall).')
    p.add_argument('--filter', default=None,
                   help="Python expression over s, out, crust, depth, rw, mgsi, diw, alpha, "
                        "kd_mg, k_na, pe, land. e.g. \"depth == 3000 and land > 0\"")
    p.add_argument('--dry-run', action='store_true', help='List what would be re-run, then stop.')
    p.add_argument('--limit', type=int, default=None, help='Re-run at most this many (cheapest first).')
    args = p.parse_args()

    path = args.path.rstrip('/')
    if not os.path.isdir(path):
        sys.exit(f"no such directory: {path}")

    assert ps.RERUN, "KAMINO_RERUN did not reach parameter_sweep -- existing output would be reused"
    print(f"scanning {path}")
    combos, skipped = wall_timeout_combos(path, args.filter)
    for name, why in skipped:
        print(f"  SKIPPED {name}: {why}")
    if not combos:
        print("no wall_timeout runs found -- nothing to do.")
        return

    # Names must round-trip, or the re-run writes new files beside the originals and every
    # wall_timeout stays exactly where it was.
    names = [cb._run_name(*c) for c in combos]
    missing = [n for n in names if not os.path.exists(os.path.join(path, f'{n}.json'))]
    if missing:
        print(f"\n{len(missing)} of {len(names)} rebuilt names do not match a file on disk, e.g.:")
        for n in missing[:5]:
            print("   ", n)
        sys.exit("aborting: run names do not round-trip, fix _int() first")

    if len(set(names)) != len(names):
        dup = sorted({n for n in names if names.count(n) > 1})
        sys.exit(f"aborting: {len(names) - len(set(names))} name collision(s), e.g. {dup[:3]}")

    combos.sort(key=lambda c: (ps._cost_rank(c), c[11] == 0.0))
    if args.limit:
        combos = combos[:args.limit]

    deep = sum(1 for c in combos if c[3] >= ps.DEEP_OCEAN_M)
    print(f"\n{len(combos)} wall_timeout run(s) to redo ({deep} deep, {len(combos) - deep} shallow)")
    print(f"wall budget: {ps.WALL_SECONDS_SHALLOW} s shallow / {ps.WALL_SECONDS_DEEP} s deep")
    worst = ps.WALL_SECONDS_DEEP * deep + ps.WALL_SECONDS_SHALLOW * (len(combos) - deep)
    print(f"worst case if every run uses its full budget: "
          f"{worst / 3600 / max(1, cb.WORKERS):.1f} h on {cb.WORKERS} workers")

    if args.dry_run:
        print("\n--dry-run: listing only\n")
        for c, n in zip(combos, [cb._run_name(*c) for c in combos]):
            print(f"  S={c[0]:<5g} depth={c[3]:<6g} land={c[11]:<6g} alpha={c[7]:<8g} {n}")
        return

    cb.run(combos, output_path=path)


if __name__ == '__main__':
    # The __main__ guard is load-bearing: ProcessPoolExecutor uses spawn, so every worker
    # re-imports this module. Without the guard each worker would re-execute the launch and try
    # to start its own pool, which multiprocessing turns into a BrokenProcessPool.
    main()
