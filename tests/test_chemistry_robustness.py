"""
Tier 1 chemistry robustness stress test.

Tests PHREEQC-facing functions directly over a grid of ocean states — no time
integration.  Completes in seconds and maps the valid operating envelope before
running any full parameter sweeps.

Usage:
    /data/pt426/big-venv/bin/python tests/test_chemistry_robustness.py
"""

import sys
import os
import io
import itertools
import contextlib
import traceback

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
os.environ['JAX_PLATFORMS'] = 'cpu'

import numpy as np
from kamino.constants import EARTH_ATM
from kamino.kamino_chem.ocean_chemistry import (
    elements, get_precipitation_flux, get_P_CO2, get_weathering_flux,
    carbonate_minerals, secondary_sink_minerals,
)

# ---------------------------------------------------------------------------
# Grid definition
# ---------------------------------------------------------------------------

T_values       = [200.0, 240.0, 274.0, 288.0, 310.0, 350.0]
P_values       = [1e5, 3e5, 1e7, 3e7, 1e8]            # Pa  (1 bar seafloor, 3 bar, 100 bar)
P_CO2_values   = [1.0, 100.0, 1e4, 1e5]    # Pa  (~10 ppb, 1 ppm, 100 ppm, 1 bar)

# Alkalinity spans freshwater to evolved ocean
Alk_values     = [1e-3, 5e-3, 0.05, 0.1, 0.2]   # mol eq / kgw

# Ca/Alk ratios: approach the charge-balance limit (2Ca = Alk → ratio 0.5)
Ca_Alk_ratios  = [0.05, 0.2, 0.4, 0.45, 0.49, 0.499]

# C is set to match Alk (realistic DIC)
C_fraction     = 0.95    # C = C_fraction * Alk

# Fixed trace concentrations
Si_fixed = 1e-4
Mg_fixed = 1e-5
Al_fixed = 1e-9
Fe_fixed = 1e-9
Na_fixed = 0.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
_C_idx   = int(np.where(elements == 'C')[0][0])
_Si_idx  = int(np.where(elements == 'Si')[0][0])
_Al_idx  = int(np.where(elements == 'Al')[0][0])
_Fe_idx  = int(np.where(elements == 'Fe')[0][0])
_Ca_idx  = int(np.where(elements == 'Ca')[0][0])
_Mg_idx  = int(np.where(elements == 'Mg')[0][0])
_Na_idx  = int(np.where(elements == 'Na')[0][0])


def make_b(alk, ca_alk_ratio):
    b = np.zeros(len(elements))
    b[_alk_idx] = alk
    b[_C_idx]   = alk * C_fraction
    b[_Ca_idx]  = alk * ca_alk_ratio
    b[_Mg_idx]  = Mg_fixed
    b[_Si_idx]  = Si_fixed
    b[_Al_idx]  = Al_fixed
    b[_Fe_idx]  = Fe_fixed
    b[_Na_idx]  = Na_fixed
    return b


@contextlib.contextmanager
def capture_stderr():
    """Capture text written to stderr (PHREEQC prints errors there via the C library)."""
    old = sys.stderr
    sys.stderr = buf = io.StringIO()
    try:
        yield buf
    finally:
        sys.stderr = old


def check_phreeqc_output(captured_text):
    """Return list of distinct PHREEQC error lines from captured output."""
    errors = []
    for line in captured_text.splitlines():
        if line.startswith('ERROR'):
            errors.append(line.strip())
    return errors


# ---------------------------------------------------------------------------
# Test runners
# ---------------------------------------------------------------------------

def run_get_precipitation(T, P, b):
    minerals = carbonate_minerals + secondary_sink_minerals
    result, _ = get_precipitation_flux(P, T, b, precipitating_minerals=minerals)
    return result


def run_get_P_CO2(P, T, b):
    return get_P_CO2(P, T, b)


def run_get_weathering(P, T, P_CO2, b):
    flux, _ = get_weathering_flux(P, T, P_CO2, b)
    return flux


TESTS = [
    ('get_precipitation_flux', run_get_precipitation,
     lambda T, P, P_CO2, b: (T, P, b)),
    ('get_P_CO2',              run_get_P_CO2,
     lambda T, P, P_CO2, b: (P, T, b)),
    ('get_weathering_flux',    run_get_weathering,
     lambda T, P, P_CO2, b: (P, T, P_CO2, b)),
]


def run_all():
    grid = list(itertools.product(T_values, P_values, P_CO2_values,
                                  Alk_values, Ca_Alk_ratios))
    total = len(grid) * len(TESTS)

    results = []   # (fn_name, T, P, P_CO2, alk, ca_alk, status, errors)

    passed = failed = phreeqc_warned = 0

    print(f"Running {total} checks ({len(grid)} states × {len(TESTS)} functions)...\n")

    for T, P, P_CO2, alk, ca_alk in grid:
        b = make_b(alk, ca_alk)

        for fn_name, fn, arg_builder in TESTS:
            args = arg_builder(T, P, P_CO2, b)
            status = 'ok'
            phreeqc_errors = []
            exception_msg = None

            # PHREEQC writes errors to stdout via print() in ocean_chemistry.py
            captured_stdout = io.StringIO()
            try:
                with contextlib.redirect_stdout(captured_stdout):
                    result = fn(*args)
                if not np.all(np.isfinite(result)):
                    status = 'nonfinite'
                    failed += 1
                else:
                    passed += 1
            except Exception as e:
                status = 'exception'
                exception_msg = f'{type(e).__name__}: {e}'
                failed += 1

            phreeqc_out = captured_stdout.getvalue()
            phreeqc_errors = check_phreeqc_output(phreeqc_out)
            if phreeqc_errors:
                phreeqc_warned += 1
                if status == 'ok':
                    status = 'phreeqc_warn'

            results.append({
                'fn':        fn_name,
                'T':         T,
                'P':         P,
                'P_CO2':     P_CO2,
                'Alk':       alk,
                'Ca_Alk':    ca_alk,
                '2Ca+2Mg':   2*b[_Ca_idx] + 2*b[_Mg_idx],
                'status':    status,
                'phreeqc':   phreeqc_errors[:1],   # first error only for brevity
                'exception': exception_msg,
            })

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"{'='*70}")
    print(f"  RESULTS:  {passed} passed   {failed} failed   {phreeqc_warned} with PHREEQC warnings")
    print(f"{'='*70}\n")

    failures    = [r for r in results if r['status'] in ('nonfinite', 'exception')]
    phreeqc_w   = [r for r in results if r['status'] == 'phreeqc_warn']

    if failures:
        print(f"--- FAILURES ({len(failures)}) ---")
        for r in failures:
            print(f"  {r['fn']:30s}  T={r['T']:5.0f}K  P={r['P']:.0e}Pa"
                  f"  P_CO2={r['P_CO2']:.0e}Pa  Alk={r['Alk']:.3f}"
                  f"  Ca/Alk={r['Ca_Alk']:.3f}  2Ca+2Mg={r['2Ca+2Mg']:.3f}")
            if r['exception']:
                print(f"    → {r['exception']}")
        print()

    if phreeqc_w:
        # Group by function and first error message
        from collections import Counter
        key = lambda r: (r['fn'], r['phreeqc'][0] if r['phreeqc'] else '')
        counts = Counter(key(r) for r in phreeqc_w)
        print(f"--- PHREEQC WARNINGS ({len(phreeqc_w)} cases, grouped) ---")
        for (fn, msg), n in sorted(counts.items(), key=lambda x: -x[1]):
            print(f"  {n:4d}×  {fn:30s}  {msg}")

        print(f"\n  Worst offenders by (T, Ca/Alk):")
        from collections import defaultdict
        by_state = defaultdict(int)
        for r in phreeqc_w:
            by_state[(r['T'], r['Ca_Alk'])] += 1
        for (T, ca_alk), n in sorted(by_state.items(), key=lambda x: -x[1])[:10]:
            print(f"    T={T:5.0f}K  Ca/Alk={ca_alk:.3f}  →  {n} warnings")
        print()

    if not failures and not phreeqc_w:
        print("  All checks passed cleanly — no PHREEQC errors, no exceptions.\n")

    return results


if __name__ == '__main__':
    run_all()
