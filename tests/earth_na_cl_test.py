"""
Test Earth-like configuration with Na and Cl chemistry.
Compares final steady-state ocean chemistry against modern seawater targets.

Two initial condition strategies run side-by-side for each config:
  ramp  — blank ocean, instellation rises 0.7→1.0 over 4.5 Gyr (physically motivated cold start)
  seeded — ocean pre-loaded at modern Earth values, fixed instellation=1.0
If both converge to the same state, the attractor is unique.
If they diverge, the system has multiple stable states.
"""

import sys
import numpy as np
import json

sys.path.insert(0, '/home/pt426/Code/kamino/src')

from kamino.constants import M_EARTH, R_EARTH, YR
from kamino.chemistry import elements, alk_idx, c_idx, ca_idx, mg_idx, na_idx, cl_idx
from kamino.planet import Planet

TARGETS = {
    'T':      288.0,    # K
    'P_CO2':  28.0,     # Pa  (280 ppm × 101325 Pa/atm)
    'pH':     8.1,
    'Alk':    2.3e-3,   # mol/kgw
    'DIC':    2.0e-3,   # mol/kgw (C)
    'Ca':     10.3e-3,  # mol/kgw
    'Mg':     52.8e-3,  # mol/kgw
    'Na':     480e-3,   # mol/kgw
    'Cl':     550e-3,   # mol/kgw
}

BACKGROUND_PRESSURE = 1e5   # Pa (1 bar N2)
OCEAN_DEPTH         = 3700  # m (mean)

T_RAMP  = 4.5e9 * YR
T_END   = 5e9   * YR

def instellation_ramp(t):
    return 0.7 + 0.3 * min(t / T_RAMP, 1.0)

b0_blank  = None

b0_seeded = np.zeros(len(elements))
b0_seeded[alk_idx] = TARGETS['Alk']
b0_seeded[c_idx]   = TARGETS['DIC']
b0_seeded[ca_idx]  = TARGETS['Ca']
b0_seeded[mg_idx]  = TARGETS['Mg']
b0_seeded[na_idx]  = TARGETS['Na']
b0_seeded[cl_idx]  = TARGETS['Cl']

base_configs = [
    dict(name='bio_fc0.0',  land_fraction=0.3, outgassing=1.0, reverse_weathering=True,
         f_HT=0.0, cl_outgassing_ratio=0.0133, cl_subduction=1.0, f_carb=0.0,  f_bio=1.0, tau_prec_yr=3e5),
    dict(name='bio_fc0.38', land_fraction=0.3, outgassing=1.0, reverse_weathering=True,
         f_HT=0.0, cl_outgassing_ratio=0.0133, cl_subduction=1.0, f_carb=0.38, f_bio=1.0, tau_prec_yr=3e5),
    dict(name='abiotic',    land_fraction=0.3, outgassing=1.0, reverse_weathering=True,
         f_HT=0.0, cl_outgassing_ratio=0.0133, cl_subduction=1.0, f_carb=0.38, f_bio=0.0, tau_prec_yr=3e5),
]

runs = []
for cfg in base_configs:
    # ramp: start high-pCO2 (0.3 bar) so T>260K is maintained at S=0.7 before outgassing builds CO2
    runs.append(dict(**cfg, name=cfg['name'] + '_ramp',   b0=b0_blank,  instellation=instellation_ramp, initial_pco2=3e4))
    runs.append(dict(**cfg, name=cfg['name'] + '_seeded', b0=b0_seeded, instellation=1.0,               initial_pco2=1000))

header = f"{'Config':<28} {'T':>6} {'pCO2':>8} {'pH':>5} {'Alk':>8} {'DIC':>8} {'Ca':>8} {'Mg':>8} {'Na':>8} {'Cl':>8}  {'term':<12}"
units  = f"{'':28} {'K':>6} {'ppm':>8} {'':>5} {'meq/kg':>8} {'mmol/kg':>8} {'mmol/kg':>8} {'mmol/kg':>8} {'mmol/kg':>8} {'mmol/kg':>8}"
sep    = '-' * len(header)
print(header)
print(units)
print(sep)

t = TARGETS
print(f"{'TARGET':<28} {t['T']:>6.1f} {t['P_CO2']/101325*1e6:>8.0f} {t['pH']:>5.2f} "
      f"{t['Alk']*1e3:>8.2f} {t['DIC']*1e3:>8.2f} {t['Ca']*1e3:>8.2f} "
      f"{t['Mg']*1e3:>8.2f} {t['Na']*1e3:>8.1f} {t['Cl']*1e3:>8.1f}  {'':12}")
print(sep)

prev_base = None
for run in runs:
    name        = run['name']
    b0          = run['b0']
    instellation = run['instellation']

    base = name.rsplit('_', 1)[0]
    if prev_base and base != prev_base:
        print()
    prev_base = base

    print(f"Running {name}...", flush=True)

    p = Planet(
        M_EARTH, R_EARTH,
        BACKGROUND_PRESSURE,
        instellation=instellation,
        crust_production_rate=1.0,
        outgassing=run['outgassing'],
        ocean_depth=OCEAN_DEPTH,
        land_fraction=run['land_fraction'],
        reverse_weathering=run['reverse_weathering'],
        f_HT=run['f_HT'],
        cl_outgassing_ratio=run['cl_outgassing_ratio'],
        cl_subduction=run['cl_subduction'],
        f_carb=run['f_carb'],
        tau_prec=run['tau_prec_yr'],
        f_bio=run['f_bio'],
        verbose=True,
        name=name,
    )

    p.time_evolve(t_end=T_END, b0=b0, initial_pco2=run['initial_pco2'])

    with open(f'/home/pt426/Code/kamino/output/{name}.json') as f:
        out = json.load(f)

    y_end = np.array(out['data']['y'])[:, -1]
    b_end = y_end[2:-1]

    T    = out['T']
    pco2 = out['P_CO2'] * 1e5 / 101325 * 1e6
    pH   = out['pH']
    alk  = b_end[alk_idx] * 1e3
    dic  = b_end[c_idx]   * 1e3
    ca   = b_end[ca_idx]  * 1e3
    mg   = b_end[mg_idx]  * 1e3
    na   = b_end[na_idx]  * 1e3
    cl   = b_end[cl_idx]  * 1e3
    term = out['termination']

    print(f"{name:<28} {T:>6.1f} {pco2:>8.0f} {pH:>5.2f} {alk:>8.2f} {dic:>8.2f} "
          f"{ca:>8.2f} {mg:>8.2f} {na:>8.1f} {cl:>8.1f}  {term:<12}")
