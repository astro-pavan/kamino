"""
Run bio_fc0.0 and bio_fc0.38 with a 100x tighter convergence criterion
to find the true steady state — testing if the model naturally approaches
Earth values when given more time to converge.
"""
import sys, numpy as np, json

sys.path.insert(0, '/home/pt426/Code/kamino/src')

from kamino.constants import M_EARTH, R_EARTH, YR
from kamino.chemistry import elements, alk_idx, c_idx, ca_idx, mg_idx, na_idx, cl_idx
from kamino.planet import Planet

TARGETS = {'T': 288.0, 'P_CO2': 28.0, 'pH': 8.1, 'Alk': 2.3e-3, 'Ca': 10.3e-3, 'Mg': 52.8e-3}
BACKGROUND_PRESSURE = 1e5
OCEAN_DEPTH = 3700

b0 = np.zeros(len(elements))
b0[na_idx] = 480e-3
b0[cl_idx] = 550e-3

configs = [
    dict(name='bio_fc0.0_tight',  land_fraction=0.3, outgassing=1.0, reverse_weathering=True,
         f_HT=0.0, cl_outgassing_ratio=0.02, cl_subduction=1.0, f_carb=0.0,  f_bio=1.0, tau_prec_yr=3e6),
    dict(name='bio_fc0.38_tight', land_fraction=0.3, outgassing=1.0, reverse_weathering=True,
         f_HT=0.0, cl_outgassing_ratio=0.02, cl_subduction=1.0, f_carb=0.38, f_bio=1.0, tau_prec_yr=3e6),
]

header = f"{'Config':<26} {'T':>6} {'pCO2':>8} {'pH':>5} {'Alk':>8} {'Ca':>8} {'Mg':>8}  {'term':<12}"
print(header)
t = TARGETS
print(f"{'TARGET':<26} {t['T']:>6.1f} {t['P_CO2']:>8.0f} {t['pH']:>5.2f} {t['Alk']*1e3:>8.2f} {t['Ca']*1e3:>8.2f} {t['Mg']*1e3:>8.2f}")
print('-' * len(header))

for cfg in configs:
    name = cfg.pop('name')
    print(f"Running {name}...", flush=True)

    p = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE,
               instellation=1.0, crust_production_rate=1.0,
               ocean_depth=OCEAN_DEPTH,
               land_fraction=cfg.pop('land_fraction'),
               outgassing=cfg.pop('outgassing'),
               reverse_weathering=cfg.pop('reverse_weathering'),
               f_HT=cfg.pop('f_HT'),
               cl_outgassing_ratio=cfg.pop('cl_outgassing_ratio'),
               cl_subduction=cfg.pop('cl_subduction'),
               f_carb=cfg.pop('f_carb'),
               tau_prec=cfg.pop('tau_prec_yr'),
               f_bio=cfg.pop('f_bio'),
               verbose=True, name=name)

    # 100x tighter convergence: 0.001/Gyr instead of 0.1/Gyr
    p.time_evolve(t_end=5e9 * YR, b0=b0, convergence_threshold=0.001)

    with open(f'/home/pt426/Code/kamino/output/{name}.json') as f:
        out = json.load(f)

    y_end = np.array(out['data']['y'])[:, -1]
    b_end = y_end[2:-1]

    T    = out['T']
    pco2 = out['P_CO2'] * 1e5 / 101325 * 1e6
    pH   = out['pH']
    alk  = b_end[alk_idx] * 1e3
    ca   = b_end[ca_idx]  * 1e3
    mg   = b_end[mg_idx]  * 1e3
    term = out['termination']

    print(f"\n{name:<26} {T:>6.1f} {pco2:>8.0f} {pH:>5.2f} {alk:>8.2f} {ca:>8.2f} {mg:>8.2f}  {term:<12}\n")
