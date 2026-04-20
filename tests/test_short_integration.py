"""
Tier 3 short-integration stress test.

Runs time_evolve_to_steady_state capped at MAX_STEPS across the full parameter
grid in parallel.  The goal is not convergence but survival: does each
combination make it through early dynamics without crashes, non-finite values,
or accumulating PHREEQC errors?

Most failure modes (charge balance blowup, climate solver divergence, Ca
runaway) appear within the first few hundred steps.  Full convergence runs are
unnecessary to detect them.

Usage:
    /data/pt426/big-venv/bin/python tests/test_short_integration.py
"""

import sys
import os
import io
import itertools
import contextlib
import tempfile
import shutil
import traceback
import concurrent.futures
import multiprocessing

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
os.environ['JAX_PLATFORMS'] = 'cpu'

import numpy as np
from tqdm import tqdm

from kamino.constants import M_EARTH, R_EARTH, YR
from kamino.planet import Planet

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MAX_STEPS   = 300    # ~6 Myr at dt=20 kyr; gives 3 convergence checks
MAX_WORKERS = 20

instellation_range = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
tectonics_range    = [0.01, 0.1, 1, 10, 100]
ocean_depth_range  = [100, 500, 1000, 3000, 5000, 10000]

BACKGROUND_PRESSURE = 1e5
FIXED_DEPTH         = 3000
FIXED_TECTONICS     = 1.0

# ---------------------------------------------------------------------------
# Worker (must be top-level for multiprocessing spawn)
# ---------------------------------------------------------------------------

def _run_one(args):
    inst, tect, depth, label = args

    # Redirect stdout so PHREEQC noise and convergence prints don't mangle
    # the progress bar in the parent process.
    captured = io.StringIO()

    tmpdir = tempfile.mkdtemp(prefix='kamino_tier3_')
    status = 'ok'
    convergence_status = None
    T_final = float('nan')
    P_CO2_final = float('nan')
    phreeqc_error_count = 0
    exception_msg = None

    try:
        p = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, inst, tect, depth,
                   name='tier3')

        with contextlib.redirect_stdout(captured):
            p.time_evolve_to_steady_state(
                output_dir=tmpdir,
                max_steps=MAX_STEPS,
            )

        # Parse summary JSON
        import json, glob
        jsons = glob.glob(os.path.join(tmpdir, '*.json'))
        if jsons:
            with open(jsons[0]) as f:
                summary = json.load(f)
            convergence_status = summary['convergence']['status']
            T_final     = summary['final_state']['T_K']
            P_CO2_final = summary['final_state']['P_CO2_Pa']

        # Count PHREEQC error lines in captured stdout
        for line in captured.getvalue().splitlines():
            if line.startswith('ERROR'):
                phreeqc_error_count += 1

        if convergence_status == 'nonfinite':
            status = 'nonfinite'
        elif phreeqc_error_count > 0:
            status = 'phreeqc_warn'
        else:
            status = 'ok'

    except Exception:
        status = 'exception'
        exception_msg = traceback.format_exc(limit=4)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    return {
        'label':              label,
        'inst':               inst,
        'tect':               tect,
        'depth':              depth,
        'status':             status,
        'convergence_status': convergence_status,
        'T_final':            T_final,
        'P_CO2_Pa':           P_CO2_final,
        'phreeqc_errors':     phreeqc_error_count,
        'exception':          exception_msg,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all():
    tasks = []

    for inst, tect in itertools.product(instellation_range, tectonics_range):
        label = f'inst={inst:.2f} tect={tect:6.2f} depth={FIXED_DEPTH}'
        tasks.append((inst, tect, FIXED_DEPTH, label))

    for inst, depth in itertools.product(instellation_range, ocean_depth_range):
        label = f'inst={inst:.2f} tect={FIXED_TECTONICS:6.2f} depth={depth}'
        tasks.append((inst, FIXED_TECTONICS, depth, label))

    total = len(tasks)
    print(f"Running {total} short integrations ({MAX_STEPS} steps each, "
          f"{MAX_WORKERS} workers)...\n")

    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = {ex.submit(_run_one, t): t for t in tasks}
        with tqdm(total=total, unit='run') as pbar:
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
                pbar.update(1)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_ok      = sum(1 for r in results if r['status'] == 'ok')
    n_warn    = sum(1 for r in results if r['status'] == 'phreeqc_warn')
    n_nonfin  = sum(1 for r in results if r['status'] == 'nonfinite')
    n_exc     = sum(1 for r in results if r['status'] == 'exception')

    print(f'\n{"="*72}')
    print(f'  {n_ok} ok   {n_warn} phreeqc_warn   {n_nonfin} nonfinite   '
          f'{n_exc} exception   / {total} total')
    print(f'{"="*72}\n')

    problems = [r for r in results if r['status'] != 'ok']
    if problems:
        for stype in ('exception', 'nonfinite', 'phreeqc_warn'):
            group = sorted([r for r in problems if r['status'] == stype],
                           key=lambda r: (r['inst'], r['tect'], r['depth']))
            if not group:
                continue
            print(f'--- {stype.upper()} ({len(group)}) ---')
            for r in group:
                print(f'  {r["label"]}  '
                      f'conv={r["convergence_status"]}  '
                      f'T={r["T_final"]:.1f} K  '
                      f'phreeqc_errors={r["phreeqc_errors"]}')
                if r['exception']:
                    for l in r['exception'].strip().splitlines()[-3:]:
                        print(f'    {l}')
            print()
    else:
        print('  All short integrations completed without errors.\n')

    # Convergence outcome breakdown
    from collections import Counter
    conv_counts = Counter(r['convergence_status'] for r in results
                          if r['convergence_status'] is not None)
    print('Convergence status breakdown:')
    for status, count in sorted(conv_counts.items(), key=lambda x: -x[1]):
        print(f'  {count:4d}  {status}')
    print()

    # PHREEQC error totals
    total_phreeqc = sum(r['phreeqc_errors'] for r in results)
    if total_phreeqc:
        print(f'Total PHREEQC error lines across all runs: {total_phreeqc}')
        worst = sorted([r for r in results if r['phreeqc_errors'] > 0],
                       key=lambda r: -r['phreeqc_errors'])[:10]
        print('Worst offenders:')
        for r in worst:
            print(f'  {r["label"]}  →  {r["phreeqc_errors"]} error lines')
    else:
        print('No PHREEQC errors across any run.')

    return results


if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')
    run_all()
