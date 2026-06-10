"""
Characterise Al sources in the model:
1. What is b_eq[Al] in the pore fluid from low-T seafloor weathering?
   (This is the waterworld Al source.)
2. What continental Al/Alk ratio would give ~20 nM steady-state ocean Al?
3. With that Al flux, how much Smectite-Na precipitates from the ocean?
"""
import sys
sys.path.insert(0, '/home/pt426/Code/kamino/src')

import numpy as np
from kamino.chemistry import elements, get_b_eq, solve_solution
from kamino.precipitation import get_precipitation
from kamino.mineral_info import basalt_49, clay_minerals, reverse_weathering_minerals
from kamino.constants import EARTH_ATM, YR, R_EARTH
from kamino.weathering import J_ref_normalised, rate_ref

na_idx  = int(np.where(elements == 'Na')[0][0])
al_idx  = int(np.where(elements == 'Al')[0][0])
mg_idx  = int(np.where(elements == 'Mg')[0][0])
alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
si_idx  = int(np.where(elements == 'Si')[0][0])

b_ocean = np.zeros(len(elements))
b_ocean[0] = 2.3e-3; b_ocean[1] = 2.1e-3; b_ocean[2] = 0.1e-3
b_ocean[3] = 1e-9;   b_ocean[4] = 1e-9;   b_ocean[5] = 10.3e-3
b_ocean[6] = 52.8e-3; b_ocean[7] = 480e-3; b_ocean[8] = 550e-3; b_ocean[9] = 28e-3

T_surface   = 287.0
T_seafloor  = max(1.02 * T_surface - 16.7, 274.0)
T_pore      = T_seafloor + 9
P_CO2       = 280e-6 * EARTH_ATM
gravity     = 9.81
ocean_depth = 3800.0
P_pore      = EARTH_ATM + 1000 * gravity * ocean_depth

surface_area     = 4 * np.pi * R_EARTH**2
ocean_water_mass = ocean_depth * surface_area * 1000
A_seafloor       = 0.7 * surface_area
tau_prec         = 100e3 * YR

K_NA_REVERSE_WEATHERING = 3.4e-3
J_total = J_ref_normalised

F_na_rw_param = (K_NA_REVERSE_WEATHERING * b_ocean[na_idx]
                 * J_total * surface_area / ocean_water_mass)   # mol/kgw/s

# ── 1. Seafloor b_eq[Al] — waterworld Al source ────────────────────────────
print("="*60)
print("1. Pore fluid equilibrium Al (waterworld source)")
print("   Excludes Anorthite, Forsterite, Enstatite at low T")
print("="*60)

try:
    b_eq, pH_eq = get_b_eq(
        P_pore, T_pore, P_CO2,
        composition={'Wollastonite': basalt_49['Wollastonite'],
                     'Fayalite':     basalt_49['Fayalite']},
        b_input=b_ocean,
        precipitating_minerals=clay_minerals,
    )
    al_eq   = b_eq[al_idx]
    F_al_sf = J_total * (al_eq - b_ocean[al_idx]) * A_seafloor  # mol/s
    print(f"  b_eq[Al] = {al_eq*1e9:.3f} nmol/kg")
    print(f"  b_ocean[Al] = {b_ocean[al_idx]*1e9:.3f} nmol/kg")
    print(f"  Seafloor Al source = {F_al_sf / 1e12 * YR:.4f} Tmol/yr")
except Exception as e:
    print(f"  get_b_eq failed: {e}")
    # Fallback: just check Al solubility at kaolinite equilibrium
    try:
        out = solve_solution(P_pore, T_pore, b_ocean, P_CO2=P_CO2,
                             precipitating_minerals=['Kaolinite'])
        al_eq_sol = float(out['Al(mol/kgw)'][-1])
        print(f"  Kaolinite-buffered Al = {al_eq_sol*1e9:.3f} nmol/kg")
    except Exception as e2:
        print(f"  Fallback also failed: {e2}")

# ── 2. Continental Al flux needed for ~20 nM steady-state Al ───────────────
print()
print("="*60)
print("2. Continental Al flux — calibration to 20 nM steady-state")
print("="*60)

# Target: ~20 nM = 2e-8 mol/kg average ocean Al
# At steady state: F_Al_in = F_Al_out (precipitation)
# Residence time of Al ~ 200 yr (dominated by scavenging/reverse weathering)
# F_Al_in = M_ocean × [Al]_ss / tau_Al
tau_Al = 200 * YR   # ~200 yr residence time
al_target = 20e-9   # mol/kg

F_al_needed = ocean_water_mass * al_target / tau_Al   # mol/s
F_alk_earth = 8e12 / YR                                # mol/s (Earth continental flux)

alpha_Al = F_al_needed / F_alk_earth
print(f"  Target [Al]_ss  = {al_target*1e9:.0f} nmol/kg")
print(f"  tau_Al          = {tau_Al/YR:.0f} yr")
print(f"  F_Al needed     = {F_al_needed/1e12*YR:.3f} Tmol/yr")
print(f"  Earth Alk flux  = {F_alk_earth/1e12*YR:.1f} Tmol/yr")
print(f"  alpha_Al (Al/Alk ratio) = {alpha_Al:.4f}")

# ── 3. Smectite-Na precipitation with realistic Al ─────────────────────────
print()
print("="*60)
print("3. Smectite-Na flux with realistic ocean [Al] = 20 nM")
print("="*60)

b_al_realistic = b_ocean.copy()
b_al_realistic[al_idx] = al_target   # 20 nM

try:
    F_prec, pH_p, SI = get_precipitation(
        P_pore, T_seafloor, b_al_realistic,
        precipitating_minerals=reverse_weathering_minerals,
        precipitation_timescale=tau_prec,
    )
    F_na_smect = F_prec[na_idx]   # mol/kgw/s (negative = removal)
    F_na_smect_global = F_na_smect * ocean_water_mass   # mol/s
    print(f"  pH = {pH_p:.2f}")
    print(f"  ΔNa   = {F_na_smect*1e15:.3f} fmol/kgw/s  ({F_na_smect_global/1e12*YR:.4f} Tmol/yr)")
    print(f"  ΔMg   = {F_prec[mg_idx]*1e15:.3f} fmol/kgw/s")
    print(f"  ΔAlk  = {F_prec[alk_idx]*1e15:.3f} fmol/kgw/s")
    for mineral, si_val in SI.items():
        if abs(si_val + 999.999) > 1e-3:
            print(f"  SI_initial({mineral}) = {si_val:+.3f}")
    ratio = abs(F_na_smect_global) / (abs(F_na_rw_param) * ocean_water_mass)
    print(f"  Smectite-Na / K_NA param = {ratio:.4f}  ({ratio*100:.2f}%)")
except Exception as e:
    print(f"  FAILED: {e}")

# ── 4. Scan: how does Smectite-Na Na flux scale with [Al]? ─────────────────
print()
print("="*60)
print("4. Na flux from Smectite-Na vs [Al] (pH dependence proxy)")
print("="*60)
print(f"  {'[Al] (nM)':>12}  {'ΔNa (Tmol/yr)':>15}  {'% of K_NA':>12}  pH")

for al_conc in [1e-9, 5e-9, 20e-9, 100e-9, 500e-9, 1e-6]:
    b_test = b_ocean.copy()
    b_test[al_idx] = al_conc
    try:
        Fp, pHp, _ = get_precipitation(
            P_pore, T_seafloor, b_test,
            precipitating_minerals=reverse_weathering_minerals,
            precipitation_timescale=tau_prec,
        )
        val = Fp[na_idx] * ocean_water_mass / 1e12 * YR
        pct = abs(val) / (abs(F_na_rw_param)*ocean_water_mass/1e12*YR) * 100
        print(f"  {al_conc*1e9:>12.1f}  {val:>15.4f}  {pct:>11.2f}%  {pHp:.2f}")
    except Exception:
        print(f"  {al_conc*1e9:>12.1f}  {'FAILED':>15}")
