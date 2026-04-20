import sys
sys.path.insert(0, 'src')

from kamino.planet import Planet
from kamino.constants import *

# In this model there are no Cl/SO4, so all cation charge must be balanced
# by alkalinity: 2*Ca + 2*Mg + Na < Alk is required for PHREEQC to work.
# Real seawater values (Ca~10 mM, Mg~54 mM) are incompatible.
#
# For P_CO2 ~ 280 ppm at pH ~ 8.2:
#   [H2CO3*] = KH * P_CO2 ~ 9.5e-6 mol/kg  =>  DIC ~ 2 mmol, Alk ~ 2.1 meq
# Ca set just at calcite saturation: [Ca] = Ksp / [CO3^2-] ~ 7e-5 mol/kg
# This ensures carbonate precipitation is active from the start.
# All 2*Ca + 2*Mg = 1.4e-4 + trace << Alk = 2.1e-3  =>  PHREEQC-consistent.

b_init = {
    'Alkalinity': 2.1e-3,   # eq/kgw  — Alk > 2*Ca to satisfy PHREEQC charge balance
    'C':          2.0e-3,   # mol/kgw — DIC for ~280 ppm with above Alk
    'Ca':         7.0e-5,   # mol/kgw — at calcite saturation (2*Ca << Alk)
    'Mg':         1.0e-5,   # mol/kgw — trace
    'Si':         1.0e-4,   # mol/kgw
    'Al':         1.0e-9,
    'Fe':         1.0e-9,
}

p = Planet(
    mass=M_EARTH,
    radius=R_EARTH,
    background_pressure=0.8e5,
    instellation=1.0,
    tectonics=1.0,
    ocean_depth=3800,
    land_fraction=0.3,
    name='earth_validation',
)
p.time_evolve_to_steady_state(output_dir='output', b_ocean_init=b_init)

import json
with open('output/earth_validation_summary.json') as f:
    s = json.load(f)

print("\n=== Earth Validation Results ===")
print(f"  Status        : {s['convergence']['status']}")
print(f"  Time (Myr)    : {s['convergence']['time_to_converge_Myr']:.1f}")
print(f"  T (K)         : {s['final_state']['T_K']:.2f}   (target ~288)")
ppm = s['final_state']['P_CO2_Pa'] / 101325 * 1e6
print(f"  P_CO2 (ppm)   : {ppm:.1f}   (target ~280)")
chem = s['final_state']['ocean_chemistry_mol_kgw']
print(f"  Alk (meq/kg)  : {chem['Alkalinity']*1e3:.3f}")
print(f"  C   (mmol/kg) : {chem['C']*1e3:.3f}")
print(f"  Ca  (mmol/kg) : {chem['Ca']*1e3:.4f}  (model equilibrium ~0.07)")
print(f"  Mg  (mmol/kg) : {chem['Mg']*1e3:.4f}")
print(f"  Da            : {s['final_diagnostics']['Da']:.3e}")
print(f"  supply_eff    : {s['final_diagnostics']['supply_efficiency']:.1%}")
