"""
Tier 2 single-step stress test.

For every (instellation, tectonics, depth) combination in the parameter sweep,
constructs the planet and calls _compute_fluxes_and_derivatives once with the
standard initial state.  No time integration — just checks that the full call
chain (climate solver + PHREEQC + precipitation + seafloor weathering) can
complete without errors.

A failure here means the combination will never survive the full sweep.

Usage:
    /data/pt426/big-venv/bin/python tests/test_single_step.py
"""

import sys
import os
import io
import contextlib
import itertools
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
os.environ['JAX_PLATFORMS'] = 'cpu'

import numpy as np
from kamino.planet import Planet, T_min, T_max
from kamino.kamino_chem.ocean_chemistry import elements
from kamino.constants import M_EARTH, R_EARTH, YR

# ---------------------------------------------------------------------------
# Parameter grids — mirrors test_planet2.py exactly
# ---------------------------------------------------------------------------

instellation_range = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
tectonics_range    = [0.01, 0.1, 1, 10, 100]
ocean_depth_range  = [100, 500, 1000, 3000, 5000, 10000]

BACKGROUND_PRESSURE = 1e5   # Pa
FIXED_DEPTH         = 3000  # m  (used when sweeping tectonics)
FIXED_TECTONICS     = 1.0   #    (used when sweeping depth)

DT = 20000 * YR   # same default as time_evolve_to_steady_state

# ---------------------------------------------------------------------------
# Helper: build the standard initial Y vector for a planet
# ---------------------------------------------------------------------------

def make_initial_Y(planet: Planet) -> np.ndarray:
    Y = np.zeros(2 + elements.shape[0])
    Y[0] = 300.0
    Y[1] = 280e-6 * planet.P_background
    return Y


def phreeqc_errors_from_stdout(text: str) -> list[str]:
    return [l.strip() for l in text.splitlines() if l.startswith('ERROR')]


# ---------------------------------------------------------------------------
# Single-planet test
# ---------------------------------------------------------------------------

def test_one(label: str, planet: Planet) -> dict:
    Y = make_initial_Y(planet)
    captured = io.StringIO()
    status = 'ok'
    phreeqc_errors = []
    exception_msg = None
    dY = None

    try:
        with contextlib.redirect_stdout(captured):
            dY, fb = planet._compute_fluxes_and_derivatives(Y, DT)

        phreeqc_errors = phreeqc_errors_from_stdout(captured.getvalue())

        if not np.all(np.isfinite(dY)):
            status = 'nonfinite_dY'
        elif phreeqc_errors:
            status = 'phreeqc_warn'
        else:
            status = 'ok'

    except Exception:
        status = 'exception'
        exception_msg = traceback.format_exc(limit=3)
        phreeqc_errors = phreeqc_errors_from_stdout(captured.getvalue())

    T_solved = float(fb['T_new']) if status not in ('exception',) else float('nan')
    P_CO2_solved = float(fb['P_CO2_new']) if status not in ('exception',) else float('nan')

    return {
        'label':      label,
        'status':     status,
        'T_new':      T_solved,
        'P_CO2_ppm':  P_CO2_solved / (BACKGROUND_PRESSURE * 1e-6) if np.isfinite(P_CO2_solved) else float('nan'),
        'phreeqc':    phreeqc_errors[:1],
        'exception':  exception_msg,
    }


# ---------------------------------------------------------------------------
# Run both sweeps
# ---------------------------------------------------------------------------

def run_all():
    tasks = []

    # Tectonics sweep
    for inst, tect in itertools.product(instellation_range, tectonics_range):
        label = f'inst={inst:.2f}  tect={tect:6.2f}  depth={FIXED_DEPTH}'
        p = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, inst, tect, FIXED_DEPTH,
                   name='stress')
        tasks.append((label, p))

    # Depth sweep
    for inst, depth in itertools.product(instellation_range, ocean_depth_range):
        label = f'inst={inst:.2f}  tect={FIXED_TECTONICS:6.2f}  depth={depth}'
        p = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, inst, FIXED_TECTONICS, depth,
                   name='stress')
        tasks.append((label, p))

    total = len(tasks)
    print(f"Running {total} single-step checks ...\n")

    results = []
    for i, (label, planet) in enumerate(tasks):
        r = test_one(label, planet)
        results.append(r)
        symbol = {'ok': '.', 'phreeqc_warn': 'W', 'nonfinite_dY': 'X',
                  'exception': 'E'}.get(r['status'], '?')
        print(symbol, end='', flush=True)
        if (i + 1) % 50 == 0:
            print(f'  {i+1}/{total}')

    print(f'\n\n{"="*70}')

    n_ok      = sum(1 for r in results if r['status'] == 'ok')
    n_warn    = sum(1 for r in results if r['status'] == 'phreeqc_warn')
    n_nonfin  = sum(1 for r in results if r['status'] == 'nonfinite_dY')
    n_exc     = sum(1 for r in results if r['status'] == 'exception')

    print(f'  {n_ok} ok   {n_warn} phreeqc_warn   {n_nonfin} nonfinite   {n_exc} exception   / {total} total')
    print(f'{"="*70}\n')

    problems = [r for r in results if r['status'] != 'ok']

    if not problems:
        print('  All combinations passed the single-step check cleanly.\n')
    else:
        # Group by status
        for status_type in ('exception', 'nonfinite_dY', 'phreeqc_warn'):
            group = [r for r in problems if r['status'] == status_type]
            if not group:
                continue
            print(f'--- {status_type.upper()} ({len(group)}) ---')
            for r in group:
                print(f'  {r["label"]}')
                if r['phreeqc']:
                    print(f'    PHREEQC: {r["phreeqc"][0]}')
                if r['exception']:
                    # Print last 3 lines of traceback only
                    tb_lines = r['exception'].strip().splitlines()
                    for l in tb_lines[-3:]:
                        print(f'    {l}')
            print()

    # Snowball / runaway summary from successful runs
    snowball  = [r for r in results if r['status'] != 'exception' and r['T_new'] <= T_min + 1]
    runaway   = [r for r in results if r['status'] != 'exception' and r['T_new'] >= T_max - 1]
    habitable = [r for r in results if r['status'] != 'exception'
                 and T_min + 1 < r['T_new'] < T_max - 1]

    print(f'Climate outcome at step 0:')
    print(f'  Habitable  : {len(habitable):3d}')
    print(f'  Snowball   : {len(snowball):3d}  (T ≤ {T_min+1:.0f} K)')
    print(f'  Runaway    : {len(runaway):3d}  (T ≥ {T_max-1:.0f} K)')
    print()

    return results


if __name__ == '__main__':
    run_all()
