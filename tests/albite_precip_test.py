"""
Test whether albite precipitates as a secondary mineral in seafloor pore fluids.

Key question: after other basalt minerals dissolve into seawater, is the pore
fluid supersaturated with respect to Albite? And if we allow Albite to
precipitate, does it consume Na?

Note: Albite as an *equilibriating* phase (unlimited source) causes PHREEQC
non-convergence at low T — seawater is too far from albite equilibrium.
So we test Albite only as a *precipitating* (secondary) phase here.
"""
import sys
sys.path.insert(0, '/home/pt426/Code/kamino/src')

import numpy as np
from kamino.chemistry import solve_solution, elements, p_LT
from kamino.mineral_info import clay_minerals
from kamino.constants import EARTH_ATM

# Earth seawater at steady state (mol/kgw)
b_ocean = np.zeros(len(elements))
b_ocean[0] = 2.3e-3    # Alk
b_ocean[1] = 2.1e-3    # C
b_ocean[2] = 0.1e-3    # Si
b_ocean[3] = 1e-9      # Al (trace)
b_ocean[4] = 1e-9      # Fe (trace)
b_ocean[5] = 10.3e-3   # Ca
b_ocean[6] = 52.8e-3   # Mg
b_ocean[7] = 480e-3    # Na
b_ocean[8] = 550e-3    # Cl
b_ocean[9] = 28e-3     # S

T_pore = 280.0
P_pore = 3e7
P_CO2  = 280e-6 * EARTH_ATM

na_idx = int(np.where(elements == 'Na')[0][0])
al_idx = int(np.where(elements == 'Al')[0][0])

def si_of(output, mineral):
    key = f'si_{mineral}'
    v = float(output.get(key, [-999.999])[-1])
    return float('nan') if abs(v + 999.999) < 1e-3 else v

def run(label, eq_minerals, prec_minerals, b_input=b_ocean):
    try:
        out = solve_solution(P_pore, T_pore, b_input, P_CO2=P_CO2,
                             equilibriating_minerals=eq_minerals,
                             precipitating_minerals=prec_minerals)
        na  = float(out['Na(mol/kgw)'][-1])
        al  = float(out['Al(mol/kgw)'][-1])
        ph  = float(out['pH'][-1])
        return dict(
            na=na, al=al, ph=ph,
            si_albite  = si_of(out, 'Albite'),
            si_kaolin  = si_of(out, 'Kaolinite'),
            si_smectna = si_of(out, 'Smectite-Na'),
            ok=True)
    except Exception:
        err = p_LT.GetErrorString()
        print(f"  [{label}] PHREEQC failed: {err[:200]}")
        return dict(ok=False)

# ─── Scenario 0: plain seawater, no minerals ───────────────────────────────
print("="*60)
print("Scenario 0: plain seawater — SI of Albite before any dissolution")
r0 = run("0", [], [], b_ocean)
if r0['ok']:
    print(f"  pH = {r0['ph']:.2f}   Al = {r0['al']*1e9:.2f} nmol/kg")
    print(f"  SI(Albite)    = {r0['si_albite']:.3f}")
    print(f"  SI(Kaolinite) = {r0['si_kaolin']:.3f}")

# ─── Non-Albite crust minerals from basalt_49 at low T ─────────────────────
# Exclude: Anorthite, Forsterite, Enstatite (model excludes these at low T)
# Exclude: Albite (to isolate the secondary precipitation question)
# Keep:    Wollastonite, Fayalite
other_crust = ['Wollastonite', 'Fayalite']

print()
print("="*60)
print("Scenario 1: Wollastonite + Fayalite dissolve, no secondary Albite")
print(f"  (Equilibriating: {other_crust}, Precipitating: {clay_minerals})")
r1 = run("1", other_crust, clay_minerals)
if r1['ok']:
    print(f"  pH = {r1['ph']:.2f}   Al = {r1['al']*1e6:.3f} µmol/kg")
    print(f"  Na = {r1['na']*1e3:.2f} mmol/kg  (ΔNa = {(r1['na']-b_ocean[na_idx])*1e3:+.3f})")
    print(f"  SI(Albite)      = {r1['si_albite']:.3f}  ({'SUPERSATURATED' if r1['si_albite']>0 else 'undersaturated'})")
    print(f"  SI(Kaolinite)   = {r1['si_kaolin']:.3f}")
    print(f"  SI(Smectite-Na) = {r1['si_smectna']:.3f}")

print()
print("="*60)
print("Scenario 2: same dissolution, WITH secondary Albite allowed to precipitate")
r2 = run("2", other_crust, clay_minerals + ['Albite'])
if r2['ok']:
    print(f"  pH = {r2['ph']:.2f}   Al = {r2['al']*1e6:.3f} µmol/kg")
    print(f"  Na = {r2['na']*1e3:.2f} mmol/kg  (ΔNa = {(r2['na']-b_ocean[na_idx])*1e3:+.3f})")
    print(f"  SI(Albite)      = {r2['si_albite']:.3f}")
    print(f"  SI(Kaolinite)   = {r2['si_kaolin']:.3f}")
    if r1['ok']:
        print(f"  Net Na from secondary Albite: {(r2['na']-r1['na'])*1e3:+.3f} mmol/kg")

print()
print("="*60)
print("Scenario 3: artificially elevated Al (simulate Al from anorthite dissolution)")
# At high T, Anorthite dissolves and releases Al. Simulate that pore fluid.
b_high_al = b_ocean.copy()
b_high_al[al_idx] = 1e-4   # 0.1 mmol/kg Al (elevated from crust dissolution)
b_high_al[2] = 1e-3        # 1 mmol/kg Si (also elevated)

r3a = run("3a-noAlbite",  [], clay_minerals,            b_high_al)
r3b = run("3b-withAlbite", [], clay_minerals + ['Albite'], b_high_al)
if r3a['ok']:
    print(f"  High-Al pore fluid (no secondary Albite):")
    print(f"  pH = {r3a['ph']:.2f}   SI(Albite) = {r3a['si_albite']:.3f}  SI(Kaolinite) = {r3a['si_kaolin']:.3f}")
if r3b['ok'] and r3a['ok']:
    print(f"  High-Al pore fluid (with secondary Albite):")
    print(f"  pH = {r3b['ph']:.2f}   Na = {r3b['na']*1e3:.2f} mmol/kg")
    print(f"  Net Na consumed by Albite: {(r3b['na']-r3a['na'])*1e3:+.3f} mmol/kg")
    print(f"  Net Al consumed by Albite: {(r3b['al']-r3a['al'])*1e6:+.3f} µmol/kg")
