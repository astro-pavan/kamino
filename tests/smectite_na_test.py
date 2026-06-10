"""
Two tests:
1. Does Smectite-Na precipitate at Earth ocean conditions? What Na flux does it provide?
2. How does that compare to the parameterised K_NA reverse weathering sink?
   i.e., if K_NA were set to zero, would Smectite-Na alone handle the Na budget?
"""
import sys
sys.path.insert(0, '/home/pt426/Code/kamino/src')

import numpy as np
from kamino.chemistry import elements, solve_solution
from kamino.precipitation import get_precipitation
from kamino.weathering import J_ref_normalised, rate_ref
from kamino.mineral_info import reverse_weathering_minerals, carbonate_minerals, clay_minerals, silica_minerals, evaporite_minerals
from kamino.constants import EARTH_ATM, YR, R_EARTH

na_idx  = int(np.where(elements == 'Na')[0][0])
mg_idx  = int(np.where(elements == 'Mg')[0][0])
ca_idx  = int(np.where(elements == 'Ca')[0][0])
alk_idx = int(np.where(elements == 'Alkalinity')[0][0])

# ── Earth ocean parameters ──────────────────────────────────────────────────
b_ocean = np.zeros(len(elements))
b_ocean[0] = 2.3e-3    # Alk
b_ocean[1] = 2.1e-3    # C
b_ocean[2] = 0.1e-3    # Si
b_ocean[3] = 1e-9      # Al
b_ocean[4] = 1e-9      # Fe
b_ocean[5] = 10.3e-3   # Ca
b_ocean[6] = 52.8e-3   # Mg
b_ocean[7] = 480e-3    # Na
b_ocean[8] = 550e-3    # Cl
b_ocean[9] = 28e-3     # S

T_surface  = 287.0                    # K
T_seafloor = max(1.02 * T_surface - 16.7, 274.0)   # planet.py formula
P_CO2      = 280e-6 * EARTH_ATM
P_surface  = EARTH_ATM
ocean_depth = 3800.0                  # m
gravity     = 9.81                    # m/s²
P_pore      = P_surface + 1000 * gravity * ocean_depth   # ~3.8e7 Pa

tau_prec    = 100e3 * YR             # 100 kyr

surface_area   = 4 * np.pi * R_EARTH**2
ocean_water_mass = ocean_depth * surface_area * 1000   # kg

J_total = J_ref_normalised           # Earth crust production rate=1, rate=rate_ref

K_NA_REVERSE_WEATHERING = 3.4e-3    # from planet.py

print(f"Earth conditions: T_seafloor={T_seafloor:.1f} K, P_pore={P_pore/1e6:.2f} MPa")
print(f"tau_prec = {tau_prec/YR/1e3:.0f} kyr")

# ── Parameterised Na RW flux (planet.py) ───────────────────────────────────
F_na_rw_param = (K_NA_REVERSE_WEATHERING
                 * b_ocean[na_idx]
                 * J_total
                 * surface_area
                 / ocean_water_mass)   # mol/kgw/s

print()
print("="*60)
print("PARAMETERISED Na RW flux (K_NA × [Na] × J × A / M)")
print("="*60)
print(f"  F_na_rw = {F_na_rw_param*1e12:.3f}  pmol/kgw/s")
print(f"  Global  = {F_na_rw_param * ocean_water_mass / 1e12 * YR:.3f}  Tmol/yr")

# ── TEST 1: SI of Smectite-Na in Earth seawater ────────────────────────────
print()
print("="*60)
print("TEST 1: SI of Smectite-Na and Saponite-Mg in Earth seawater")
print("="*60)

out_sw = solve_solution(P_pore, T_seafloor, b_ocean, P_CO2=P_CO2)
for mineral in ['Smectite-Na', 'Saponite-Mg', 'Kaolinite', 'Calcite']:
    key = f'si_{mineral}'
    v = float(out_sw.get(key, [-999.999])[-1])
    if abs(v + 999.999) > 1e-3:
        flag = 'SUPERSATURATED' if v > 0 else 'undersaturated'
        print(f"  SI({mineral:16s}) = {v:+.3f}  {flag}")
    else:
        print(f"  SI({mineral:16s}) = N/A")

# ── TEST 2: PHREEQC Smectite-Na precipitation flux ─────────────────────────
print()
print("="*60)
print("TEST 2: get_precipitation with reverse_weathering_minerals")
print(f"  Minerals: {reverse_weathering_minerals}")
print("="*60)

try:
    F_prec, pH, SI = get_precipitation(
        P_pore, T_seafloor, b_ocean,
        precipitating_minerals=reverse_weathering_minerals,
        precipitation_timescale=tau_prec,
    )
    print(f"  pH after precipitation = {pH:.2f}")
    print(f"  ΔNa  = {F_prec[na_idx]*1e15:.3f}  fmol/kgw/s  ({F_prec[na_idx]*ocean_water_mass/1e12*YR:.4f} Tmol/yr)")
    print(f"  ΔAlk = {F_prec[alk_idx]*1e15:.3f}  fmol/kgw/s")
    print(f"  ΔMg  = {F_prec[mg_idx]*1e15:.3f}  fmol/kgw/s")
    print(f"  ΔCa  = {F_prec[ca_idx]*1e15:.3f}  fmol/kgw/s")
    for mineral, si_val in SI.items():
        if abs(si_val + 999.999) > 1e-3:
            print(f"  SI_initial({mineral}) = {si_val:+.3f}")

    print()
    ratio = abs(F_prec[na_idx]) / abs(F_na_rw_param) if F_na_rw_param != 0 else float('inf')
    print(f"  Smectite-Na flux / K_NA param flux = {ratio:.4f}")
    print(f"  → PHREEQC Smectite-Na provides {ratio*100:.2f}% of the parameterised Na sink")

except Exception as e:
    print(f"  FAILED: {e}")

# ── TEST 3: Full ocean precipitating minerals (as in the model) ─────────────
print()
print("="*60)
print("TEST 3: Full ocean_precipitating_minerals (all minerals incl. RW)")
ocean_prec = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals + reverse_weathering_minerals
print(f"  Minerals: {ocean_prec}")
print("="*60)

try:
    F_prec_full, pH_full, SI_full = get_precipitation(
        P_pore, T_seafloor, b_ocean,
        precipitating_minerals=ocean_prec,
        precipitation_timescale=tau_prec,
    )
    print(f"  pH = {pH_full:.2f}")
    print(f"  ΔNa  = {F_prec_full[na_idx]*1e15:.3f}  fmol/kgw/s  ({F_prec_full[na_idx]*ocean_water_mass/1e12*YR:.4f} Tmol/yr)")
    print(f"  ΔMg  = {F_prec_full[mg_idx]*1e15:.3f}  fmol/kgw/s")
    for mineral in ['Smectite-Na', 'Saponite-Mg', 'Kaolinite', 'Calcite']:
        si_val = SI_full.get(mineral, -999.999)
        if abs(si_val + 999.999) > 1e-3:
            print(f"  SI_initial({mineral:16s}) = {si_val:+.3f}")
    ratio_full = abs(F_prec_full[na_idx]) / abs(F_na_rw_param) if F_na_rw_param != 0 else float('inf')
    print(f"  Smectite-Na vs K_NA ratio = {ratio_full:.4f}  ({ratio_full*100:.2f}%)")
except Exception as e:
    print(f"  FAILED: {e}")
