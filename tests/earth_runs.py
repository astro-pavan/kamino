"""
Earth calibration runs — 5 configurations to explore the parameter space.
All use S=1.0 (solar), basalt_49 crust, 3800 m ocean depth.

Run with: /data/pt426/big-venv/bin/python tests/earth_runs.py
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import kamino.planet as p_mod
p_mod.output_path = os.path.join(os.path.dirname(__file__), '../output/')

from kamino.planet import Planet
from kamino.mineral_info import basalt_49
from kamino.constants import M_EARTH, R_EARTH

CONFIGS = [
    # (name_tag, land, f_HT, f_carb, outgassing, rw, desc)
    ('ow_rw',         0.0, 0.0, 0.0, 1.0, True,
     'Ocean world: seafloor only, rw=True'),
    ('cont_out1',     0.3, 0.0, 0.0, 1.0, True,
     'Continental (30%), full outgassing'),
    ('cont_out05',    0.3, 0.0, 0.0, 0.5, True,
     'Continental (30%), half outgassing (organic-C burial proxy)'),
    ('cont_ht_out05', 0.3, 0.3, 0.0, 0.5, True,
     'Continental + HT (30%), half outgassing'),
    ('cont_full',     0.3, 0.3, 0.5, 0.5, True,
     'Full: continental + HT + carbonate weathering, half outgassing'),
]

results = {}

for tag, land, fht, fcarb, out, rw, desc in CONFIGS:
    name = f'earth_{tag}'
    print(f'\n{"="*60}')
    print(f'  {name}: {desc}')
    print(f'{"="*60}')
    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=1.0,
            crust_production_rate=1.0,
            outgassing=out,
            ocean_depth=3800,
            land_fraction=land,
            crust_composition=basalt_49,
            reverse_weathering=rw,
            f_HT=fht,
            f_carb=fcarb,
            name=name,
            verbose=True,
        )
        p.time_evolve(t_end=2e9 * 3.16e7)
    except Exception as e:
        print(f'  ERROR: {e}')

# Print summary table
import json
print(f'\n\n{"="*80}')
print(f'  EARTH CALIBRATION SUMMARY')
print(f'{"="*80}')
print(f'  {"Config":<22}  {"Term":>10}  {"T(K)":>6}  {"pCO2(ppm)":>10}  {"pH":>5}  {"Alk(meq)":>8}  {"C(mmol)":>8}  {"Ca(uM)":>7}  {"Mg(uM)":>7}')
print(f'  {"-"*22}  {"-"*10}  {"-"*6}  {"-"*10}  {"-"*5}  {"-"*8}  {"-"*8}  {"-"*7}  {"-"*7}')

targets = {'T': 288, 'pCO2': 280, 'pH': 8.1, 'Alk': 2.3, 'C': 2.0}

for tag, *_ in CONFIGS:
    name = f'earth_{tag}'
    fpath = os.path.join(p_mod.output_path, f'{name}.json')
    try:
        with open(fpath) as f:
            d = json.load(f)
        term = d.get('termination', '?')
        T = d.get('T', float('nan'))
        P_CO2_bar = d.get('P_CO2', float('nan'))
        P_CO2_ppm = P_CO2_bar * 1e5 / 101325 * 1e6 if P_CO2_bar else float('nan')
        pH = d.get('pH', float('nan'))
        y = d.get('data', {}).get('y', [])
        if y and len(y) >= 9:
            alk = float(y[2][-1]) * 1e3      # meq/kg
            c   = float(y[3][-1]) * 1e3      # mmol/kg
            ca  = float(y[7][-1]) * 1e6      # umol/kg
            mg  = float(y[8][-1]) * 1e6      # umol/kg
        else:
            alk = c = ca = mg = float('nan')
        print(f'  {name:<22}  {term:>10}  {T:>6.1f}  {P_CO2_ppm:>10.1f}  {pH:>5.2f}  {alk:>8.2f}  {c:>8.2f}  {ca:>7.1f}  {mg:>7.1f}')
    except Exception as e:
        print(f'  {name:<22}  ERROR: {e}')

print(f'\n  {"TARGET (Earth)":<22}  {"converged":>10}  {288:>6.1f}  {280:>10.1f}  {8.1:>5.1f}  {2.3:>8.1f}  {2.0:>8.1f}  {"10000":>7}  {"53000":>7}')
print(f'\n  Note: climate model gives T=280K at 280ppm (8K cold bias), Ca/Mg in model')
print(f'        are low due to absent Na/K/Cl/SO4 in charge balance.')
