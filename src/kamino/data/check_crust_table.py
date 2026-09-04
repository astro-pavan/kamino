"""Validate a generated crust_compositions.csv, end to end.

Runs the checks that a merged table has to pass before the model is pointed at it:

  1. structural  -- complete rectangular grid, monotonic axes, no NaNs
  2. anchor      -- Mg/Si 1.25, dIW -2 reproduces mantle FeO 8.05 and a basaltic melt
  3. norm        -- every grid point passes cipw_norm with no warnings, and mass closes
  4. petrology   -- CaO/Al2O3 flagged where it exceeds the ultracalcic threshold
  5. chemistry   -- PHREEQC equilibrates the crust at the anchor and both Mg/Si extremes

Step 5 is the one that matters most, and section 23.9 of the development history explains why:
the Python-side stoichiometry parser and PHREEQC read the same database differently, so a phase
can look present to Python and be invisible to the solver. Validate with a solve, never with
the parser.

    python check_crust_table.py [--table PATH] [--skip-chemistry]
"""

import argparse
import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))

from kamino.constants import EARTH_DELTA_IW, EARTH_MANTLE_MG_SI
from kamino.crust_composition import (CRUST_TABLE, cipw_norm, feo_from_delta_iw,
                                      load_crust_table, mineral_composition, oxide_composition)
from kamino.mineral_info import MINERAL_MOLAR_MASS

ULTRACALCIC_WARN = 1.2          # must match make_crust_compositions.jl
EARTH_MANTLE_FEO = 8.05

# Oxide components per formula unit, for the mass-closure check. NOTE these are per FORMULA
# UNIT, not per cation: Forsterite is Mg2SiO4, so 1 SiO2 and 2 MgO. cipw_norm counts Si per
# CATION internally (0.5 Si per Mg in olivine), and conflating the two silently breaks the check
# rather than the norm.
FORMULA = {
    'Quartz':       {'SiO2': 1},
    'Albite':       {'Na2O': 0.5, 'Al2O3': 0.5, 'SiO2': 3},
    'Nepheline':    {'Na2O': 0.5, 'Al2O3': 0.5, 'SiO2': 1},
    'Anorthite':    {'CaO': 1, 'Al2O3': 1, 'SiO2': 2},
    'Diopside':     {'CaO': 1, 'MgO': 1, 'SiO2': 2},
    # Fe endmembers of the two pyroxenes, emitted directly by the norm rather than exchanged into
    # fayalite, so mass closure has to account for them.
    'Hedenbergite': {'CaO': 1, 'FeOt': 1, 'SiO2': 2},
    'Ferrosilite':  {'FeOt': 1, 'SiO2': 1},
    'Akermanite':   {'CaO': 2, 'MgO': 1, 'SiO2': 2},
    'Wollastonite': {'CaO': 1, 'SiO2': 1},
    'Enstatite':    {'MgO': 1, 'SiO2': 1},
    'Forsterite':   {'MgO': 2, 'SiO2': 1},
    'Fayalite':     {'FeOt': 2, 'SiO2': 1},
    'K-Feldspar':   {'K2O': 0.5, 'Al2O3': 0.5, 'SiO2': 3},
}
OXIDE_MASS = {'SiO2': 60.084, 'Al2O3': 101.961, 'FeOt': 71.844, 'MgO': 40.304,
              'CaO': 56.077, 'Na2O': 61.979, 'K2O': 94.196}
# Oxides the norm cannot express and therefore drops: no rutile, no chromite, and K is not a
# tracked ion so K-Feldspar is off.
NORM_DROPS = ('TiO2', 'Cr2O3', 'K2O')


def reconstruct(assemblage):
    """Oxide analysis (wt%) implied by a weight-fraction mineral assemblage."""
    out = {}
    for m, f in assemblage.items():
        n = f / MINERAL_MOLAR_MASS[m]
        for ox, k in FORMULA[m].items():
            out[ox] = out.get(ox, 0.0) + n * k * OXIDE_MASS[ox]
    return {k: v * 100 for k, v in out.items()}

_fails: list[str] = []


def check(ok: bool, label: str, detail: str = '') -> bool:
    print(f'  [{"PASS" if ok else "FAIL"}] {label}{"  -- " + detail if detail else ""}')
    if not ok:
        _fails.append(label)
    return ok


def main(table: str, skip_chemistry: bool) -> int:
    print(f'Checking {table}\n')
    interp, axes = load_crust_table(table)
    df, mg_si, d_iw = axes['table'], axes['mg_si'], axes['delta_iw']

    print('1. Structure')
    check(len(df) == len(mg_si) * len(d_iw), 'complete rectangular grid',
          f'{len(df)} rows = {len(mg_si)} Mg/Si x {len(d_iw)} dIW')
    check(bool(np.all(np.diff(mg_si) > 0)) and bool(np.all(np.diff(d_iw) > 0)), 'axes strictly increasing')
    oxcols = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O']
    check(not df[oxcols].isna().to_numpy().any(), 'no NaN oxides')
    check(bool(np.allclose(df[oxcols + ['Cr2O3']].sum(axis=1), 100.0, atol=0.5)),
          'oxides sum to 100 wt%')
    print(f'      Mg/Si {mg_si.min():g} .. {mg_si.max():g}   dIW {d_iw.min():+g} .. {d_iw.max():+g}')

    print('\n2. Earth anchor')
    row = df[np.isclose(df['mg_si'], EARTH_MANTLE_MG_SI) & np.isclose(df['delta_iw'], EARTH_DELTA_IW)]
    if check(len(row) == 1, 'Earth grid point present'):
        r = row.iloc[0]
        check(abs(r['mantle_feo'] - EARTH_MANTLE_FEO) < 1e-2, 'mantle FeO = 8.05 wt%',
              f"{r['mantle_feo']:.3f}")
        check(abs(feo_from_delta_iw(EARTH_DELTA_IW) - EARTH_MANTLE_FEO) < 1e-9,
              'python/julia FeO mappings agree')
        check(44 <= r['SiO2'] <= 52, 'Earth melt is basaltic', f"SiO2 {r['SiO2']:.2f} wt%")
        # T_melt in isobaric tables, T_p in the older isentropic ones.
        t_col = 'T_melt' if 'T_melt' in r else 'T_p'
        print(f"      {t_col} {r[t_col]:.0f} C   F {r['melt_fraction']:.3f}   "
              f"MgO {r['MgO']:.2f}   Al2O3 {r['Al2O3']:.2f}   CaO {r['CaO']:.2f}   "
              f"CaO/Al2O3 {r['CaO_Al2O3']:.2f}")
        print(f"      PRIMELT reference:  SiO2 48.76  MgO 11.27  Al2O3 17.05  CaO 11.91")
        misfit = (abs(r['SiO2'] - 48.76) + abs(r['MgO'] - 11.27)
                  + abs(r['Al2O3'] - 17.05) + abs(r['CaO'] - 11.91))
        print(f"      misfit to PRIMELT = {misfit:.2f} wt% "
              f"(section 24.2 got 2.22 at F = 0.117; F = 0.20 is expected to be worse)")

    print('\n3. CIPW norm over every grid point')
    bad_norm, bad_mass, quartz_rows, neph_rows = [], [], 0, 0
    for _, r in df.iterrows():
        ox = {'SiO2': r['SiO2'], 'TiO2': r['TiO2'], 'Al2O3': r['Al2O3'], 'FeOt': r['FeO'],
              'MgO': r['MgO'], 'CaO': r['CaO'], 'Na2O': r['Na2O'], 'K2O': r['K2O']}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            c = cipw_norm(ox, emit_quartz=True)
            msgs = [str(x.message) for x in w]
        if msgs:
            bad_norm.append((r['mg_si'], r['delta_iw'], msgs))
        # Mass closure on EVERY oxide, not just silica: the norm can be silica-balanced and
        # still have lost Al to the peraluminous drop or Ca to the wollastonite overflow.
        # Compare against the melt renormalised over the oxides the norm can express.
        kept = {k: v for k, v in ox.items() if k not in NORM_DROPS}
        tot = sum(kept.values())
        expect = {k: v / tot * 100 for k, v in kept.items()}
        got = reconstruct(c)
        worst = max((abs(got.get(k, 0.0) - expect.get(k, 0.0)), k) for k in expect)
        if worst[0] > 0.05:
            bad_mass.append((r['mg_si'], r['delta_iw'], worst[1], worst[0]))
        quartz_rows += 'Quartz' in c
        neph_rows += 'Nepheline' in c
    check(not bad_norm, 'no cipw_norm warnings on any grid point',
          f'{len(bad_norm)} rows warned' if bad_norm else '')
    for m, d, msgs in bad_norm[:5]:
        print(f'        Mg/Si={m:.2f} dIW={d:+.1f}: {"; ".join(msgs)}')
    check(not bad_mass, 'assemblage mass-balances the melt on every oxide (<0.05 wt%)',
          f'{len(bad_mass)} rows off' if bad_mass else '')
    for m, d, ox_worst, err in bad_mass[:5]:
        print(f'        Mg/Si={m:.2f} dIW={d:+.1f}: worst oxide {ox_worst} off by {err:.2f} wt%')
    print(f'      Quartz-bearing rows: {quartz_rows}/{len(df)}   '
          f'Nepheline-bearing rows: {neph_rows}/{len(df)}')

    print('\n4. Petrology flags (reported, not failed -- F = 0.20 melts past cpx-out at high Mg/Si)')
    ultra = df[df['CaO_Al2O3'] > ULTRACALCIC_WARN]
    print(f'      ultracalcic (CaO/Al2O3 > {ULTRACALCIC_WARN}): {len(ultra)}/{len(df)} rows')
    if len(ultra):
        print(f'      Mg/Si range affected: {ultra["mg_si"].min():g} .. {ultra["mg_si"].max():g}, '
              f'max CaO/Al2O3 {ultra["CaO_Al2O3"].max():.2f}')
    silicic = df[df['SiO2'] > 60]
    print(f'      silicic melts (SiO2 > 60 wt%): {len(silicic)}/{len(df)} rows'
          + (f', Mg/Si <= {silicic["mg_si"].max():g}' if len(silicic) else ''))

    print('\n5. End-to-end chemistry (PHREEQC must equilibrate, not just parse)')
    if skip_chemistry:
        print('      skipped (--skip-chemistry)')
    else:
        from kamino.chemistry import get_b_eq
        for mgsi in (float(mg_si.min()), EARTH_MANTLE_MG_SI, float(mg_si.max())):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    crust = mineral_composition(mgsi, EARTH_DELTA_IW)
                    _, pH = get_b_eq(300e5, 300.0, 280 * 101325, crust,
                                     dissolve_only=True, water_rock_ratio=3.0)
                check(True, f'Mg/Si={mgsi:g} equilibrates', f'pH {float(pH):.2f}, '
                      f'{len(crust)} phases: {", ".join(sorted(crust))}')
            except Exception as e:
                check(False, f'Mg/Si={mgsi:g} equilibrates', f'{type(e).__name__}: {e}')

    print(f'\n{"ALL CHECKS PASSED" if not _fails else str(len(_fails)) + " CHECK(S) FAILED: " + "; ".join(_fails)}')
    return 1 if _fails else 0


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=CRUST_TABLE)
    ap.add_argument('--skip-chemistry', action='store_true')
    a = ap.parse_args()
    sys.exit(main(a.table, a.skip_chemistry))
