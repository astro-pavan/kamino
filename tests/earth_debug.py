"""
Earth validation diagnostic script.

Runs five offline checks to identify why the Earth validation fails
without running the full time evolution. Prints a structured report.
"""

import sys
sys.path.insert(0, 'src')

import numpy as np
from kamino.kamino_chem.ocean_chemistry import (
    elements, get_precipitation_flux, get_weathering_flux, get_P_CO2,
    get_continental_weathering_flux, carbonate_minerals, secondary_sink_minerals,
    EARTH_CONTINENTAL_WEATHERING_REF, solve_solution,
)
from kamino.speedy_climate.clima_interpolator import get_T_surface
from kamino.utils import august_roche_magnus_formula
from kamino.constants import *

YR = 3.16e7
R_EARTH = 6371000.0
SURFACE_AREA = 4 * np.pi * R_EARTH**2
OCEAN_DEPTH = 3800.0
OCEAN_MASS = OCEAN_DEPTH * SURFACE_AREA * 1000.0  # kg

LAND_FRACTION = 0.3
LAND_AREA = LAND_FRACTION * SURFACE_AREA
SEA_AREA = (1 - LAND_FRACTION) * SURFACE_AREA
P_BACKGROUND = 0.8e5

# Reference ocean composition (from earth_validation.py b_init)
b_ref = np.zeros(len(elements))
for e, v in [
    ('Alkalinity', 2.1e-3), ('C', 2.0e-3), ('Ca', 7.0e-5),
    ('Mg', 1.0e-5),         ('Si', 1.0e-4), ('Al', 1.0e-9),
    ('Fe', 1.0e-9),
]:
    b_ref[int(np.where(elements == e)[0][0])] = v

_alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
_C_idx   = int(np.where(elements == 'C')[0][0])
_Ca_idx  = int(np.where(elements == 'Ca')[0][0])
_Mg_idx  = int(np.where(elements == 'Mg')[0][0])
_Si_idx  = int(np.where(elements == 'Si')[0][0])

def tmol_yr(flux_mol_kgw_s):
    """Convert mol/kgw/s flux to Tmol/yr total."""
    return flux_mol_kgw_s * OCEAN_MASS * YR / 1e12

def sep(title):
    print()
    print('=' * 70)
    print(f'  {title}')
    print('=' * 70)

# ─────────────────────────────────────────────────────────────────────────────
sep('STEP 1 — STEADY-STATE FLUX AUDIT AT REFERENCE CONDITIONS')
# ─────────────────────────────────────────────────────────────────────────────
print('Reference state: T=288 K, P_CO2=280 ppm, b = b_init')
print()

T_ref    = 288.0
P_CO2_ref = EARTH_ATM * 280e-6
P_H2O_ref = august_roche_magnus_formula(T_ref) * 0.5
P_atm_ref = P_BACKGROUND + P_CO2_ref + P_H2O_ref

# --- Outgassing ---
F_out = np.zeros(len(elements))
F_out[_C_idx] = (EARTH_OUTGASSING / YR) * SURFACE_AREA * 1.0  # tectonics=1
F_out_norm = F_out / OCEAN_MASS  # mol/kgw/s

# --- Continental weathering ---
F_cont_per_area = get_continental_weathering_flux(T_ref, P_CO2_ref)
F_cont = F_cont_per_area * LAND_AREA  # mol/s total
F_cont_norm = F_cont / OCEAN_MASS

# --- Seafloor weathering ---
T_seafloor_ref = max(1.02 * T_ref - 16.7, 274.0)
T_pore_ref     = T_seafloor_ref + 9
P_pore_ref     = P_atm_ref + 1000 * (G * 5.972e24 / R_EARTH**2) * OCEAN_DEPTH
crust_rate     = EARTH_CRUST_PRODUCTION_RATE_PER_AREA  # tectonics=1
hydro_flux     = EARTH_HYDROTHERMAL_FLUX_PER_AREA

F_diss_per_area, w_diag = get_weathering_flux(
    P_pore_ref, T_pore_ref, P_CO2_ref, b_ref,
    alpha=2, rate=crust_rate, clog=False,
)
F_diss = F_diss_per_area * SURFACE_AREA
F_diss_norm = F_diss / OCEAN_MASS

# --- Precipitation (surface) ---
P_surf_ref = P_atm_ref
delta_b_prec_surf = get_precipitation_flux(
    P_surf_ref, T_ref, b_ref,
    precipitating_minerals=carbonate_minerals,
)
# tau_prec = 1e4 yr, use limit dt→∞: F_prec = delta_b / tau_prec
tau_prec = 1e4 * YR
F_prec_surf_norm = delta_b_prec_surf / tau_prec
F_prec_surf = F_prec_surf_norm * OCEAN_MASS

# --- Precipitation (seafloor) ---
delta_b_prec_deep = get_precipitation_flux(
    P_pore_ref, T_seafloor_ref, b_ref,
    precipitating_minerals=carbonate_minerals + secondary_sink_minerals,
)
F_prec_deep_norm = delta_b_prec_deep / tau_prec
F_prec_deep = F_prec_deep_norm * OCEAN_MASS

F_prec = F_prec_surf + LAND_FRACTION * F_prec_surf  # NEW: shallow uses land_fraction
# Recalculate properly
F_prec_total     = F_prec_deep + LAND_FRACTION * F_prec_surf
F_prec_total_norm = F_prec_total / OCEAN_MASS

# --- Net ---
F_net = F_out + F_cont + F_diss + F_prec_total

print(f'  {"Flux":20s}  {"Alk":>10s}  {"C":>10s}  {"Ca":>10s}  {"Mg":>10s}  unit')
print(f'  {"-"*20}  {"-"*10}  {"-"*10}  {"-"*10}  {"-"*10}  ----')

def frow(name, arr):
    scale = YR / 1e12  # mol/s → Tmol/yr
    a = arr[_alk_idx] * scale
    c = arr[_C_idx]   * scale
    ca= arr[_Ca_idx]  * scale
    mg= arr[_Mg_idx]  * scale
    print(f'  {name:20s}  {a:>10.3f}  {c:>10.3f}  {ca:>10.3f}  {mg:>10.3f}  Tmol/yr')

frow('F_out (volcanic)',   F_out)
frow('F_cont (WHK ref)',   F_cont)
frow('F_diss (seafloor)',  F_diss)
frow('F_prec_deep',        F_prec_deep)
frow('F_prec_surf×fland',  LAND_FRACTION * F_prec_surf)
frow('F_NET',              F_net)

print()
print('  Needed for C steady state:')
needed_prec_C = -(F_out[_C_idx] + F_cont[_C_idx] + F_diss[_C_idx])
print(f'    F_prec_C must be {needed_prec_C/1e12*YR:.2f} Tmol/yr')
actual_prec_C = F_prec_total[_C_idx]
print(f'    Actual F_prec_C = {actual_prec_C/1e12*YR:.4f} Tmol/yr')
print(f'    Deficit factor  = {abs(needed_prec_C/actual_prec_C) if actual_prec_C != 0 else float("inf"):.1f}×')

print()
print('  Urey cycle check (F_out should ≈ net silicate weathering CO2 sink):')
print(f'    F_out_C          = {F_out[_C_idx]/1e12*YR:.2f} Tmol C/yr')
print(f'    F_cont net CO2 sink = F_cont_Ca × 1  = {F_cont[_Ca_idx]/1e12*YR:.2f} Tmol/yr')
print(f'    (Each mol CaSiO3 weathers → net -1 mol CO2 from atm, not -2)')
print(f'    Imbalance (out - net_sink) = {(F_out[_C_idx] - F_cont[_Ca_idx])/1e12*YR:.2f} Tmol/yr')
print(f'    Note: ~half of Earth outgassing is balanced by organic C burial')
print(f'    (absent in this model) → structural ~3.5 Tmol/yr excess')

# ─────────────────────────────────────────────────────────────────────────────
sep('STEP 2 — PRECIPITATION SENSITIVITY SWEEP')
# ─────────────────────────────────────────────────────────────────────────────
print('Sweeping Alk/DIC ratio and Ca at T=288K, P=1 atm')
print('(How supersaturated does the ocean need to be for adequate precipitation?)')
print()

# Target: remove 15 Tmol C/yr → in mol/kgw/s:
target_Cflux_norm = 15e12 / (OCEAN_MASS * YR)  # mol/kgw/s
# With tau_prec = 1e4 yr: delta_b_C needed = target_Cflux × tau_prec
target_delta_b_C = target_Cflux_norm * tau_prec
print(f'  Target delta_b_C needed: {target_delta_b_C:.4e} mol/kgw')
print(f'  (to remove 15 Tmol C/yr via precipitation)')
print()

print(f'  {"Alk(meq)":>9} {"DIC(mmol)":>9} {"Ca(mmol)":>8} {"SI_Cal":>7} {"db_C(mol/kgw)":>14} {"F_prec_C(Tmol/yr)":>18}')
cases = [
    # Alk,   DIC,    Ca
    (2.1e-3, 2.0e-3, 7e-5),    # initial b_init
    (2.3e-3, 2.1e-3, 1.0e-2),  # modern seawater Ca
    (5.0e-3, 2.0e-3, 7e-5),    # high Alk
    (5.0e-3, 4.0e-3, 5e-4),    # high Alk, high DIC, high Ca
    (10e-3,  2.0e-3, 1e-3),    # very high Alk, low DIC
    (10e-3,  9.0e-3, 1e-3),    # very high Alk, high DIC, high Ca
    (20e-3,  10e-3,  5e-3),    # extreme
]
for alk, dic, ca in cases:
    b_test = b_ref.copy()
    b_test[_alk_idx] = alk
    b_test[_C_idx]   = dic
    b_test[_Ca_idx]  = ca
    try:
        out = solve_solution(P_surf_ref, T_ref, b_test, verbose=False)
        import pandas as pd
        df = pd.DataFrame(out)
        si_cal = float(df['si_Calcite'].iloc[-1]) if 'si_Calcite' in df.columns else float('nan')
        db = get_precipitation_flux(P_surf_ref, T_ref, b_test,
                                     precipitating_minerals=carbonate_minerals)
        db_C = db[_C_idx]
        fp_C = db_C * OCEAN_MASS * YR / (tau_prec * 1e12)
        print(f'  {alk*1e3:>9.2f} {dic*1e3:>9.2f} {ca*1e3:>8.3f} {si_cal:>7.3f} {db_C:>14.4e} {fp_C:>18.3f}')
    except Exception as ex:
        print(f'  {alk*1e3:>9.2f} {dic*1e3:>9.2f} {ca*1e3:>8.3f}  FAILED: {ex}')

# ─────────────────────────────────────────────────────────────────────────────
sep('STEP 3 — Ca/Mg CHARGE BALANCE AUDIT')
# ─────────────────────────────────────────────────────────────────────────────
print('Which flux drives 2Ca+2Mg toward Alk (causing PHREEQC failure)?')
print()

scale = YR / 1e12
print(f'  {"Flux":25s}  {"Δ(2Ca+2Mg) (Tmol/yr)":>22s}  {"ΔAlk (Tmol/yr)":>16s}  {"net_ratio_change":>16s}')
def cation_charge(arr):
    return 2*arr[_Ca_idx] + 2*arr[_Mg_idx]

fluxes = [
    ('F_out',              F_out),
    ('F_cont (WHK ref)',   F_cont),
    ('F_diss (seafloor)',  F_diss),
    ('F_prec_deep',        F_prec_deep),
    ('F_prec_surf×fland',  LAND_FRACTION * F_prec_surf),
]
for name, arr in fluxes:
    d_cat = cation_charge(arr) * scale
    d_alk = arr[_alk_idx] * scale
    # How much does Alk/cation ratio change? (positive = safer, negative = worse)
    # If Alk/(2Ca+2Mg) at ref = alk_ref/cation_ref, then d(ratio) ≈ (d_alk*cation - alk*d_cat)/cation²
    alk_ref  = b_ref[_alk_idx]
    cat_ref  = 2*b_ref[_Ca_idx] + 2*b_ref[_Mg_idx]
    d_ratio  = (d_alk * cat_ref - alk_ref * d_cat) / cat_ref**2 if cat_ref > 0 else float('nan')
    print(f'  {name:25s}  {d_cat:>22.4f}  {d_alk:>16.4f}  {d_ratio:>16.3f}')

print()
print('  Current Alk/(2Ca+2Mg) ratio at b_init:',
      f'{b_ref[_alk_idx] / (2*b_ref[_Ca_idx]+2*b_ref[_Mg_idx]):.1f}')
print('  PHREEQC fails when ratio < 1.0')

# ─────────────────────────────────────────────────────────────────────────────
sep('STEP 4 — UREY CYCLE / OUTGASSING CALIBRATION CHECK')
# ─────────────────────────────────────────────────────────────────────────────
print('Checking consistency between outgassing, weathering, and burial fluxes')
print()

F_out_C_tmol  = F_out[_C_idx] * YR / 1e12
F_cont_Ca_tmol = F_cont[_Ca_idx] * YR / 1e12
F_cont_C_tmol  = F_cont[_C_idx] * YR / 1e12
F_cont_alk_tmol = F_cont[_alk_idx] * YR / 1e12

print(f'  Volcanic outgassing F_out_C = {F_out_C_tmol:.2f} Tmol C/yr')
print(f'  WHK reference: F_cont_Alk   = {F_cont_alk_tmol:.2f} Tmol eq/yr')
print(f'               F_cont_C      = {F_cont_C_tmol:.2f} Tmol C/yr')
print(f'               F_cont_Ca     = {F_cont_Ca_tmol:.2f} Tmol Ca/yr')
print()
print(f'  CaCO3 stoichiometry (per mol): -1 Ca, -1 C (net), -2 Alk')
print(f'  For C balance:   F_burial_C  = F_out_C + F_cont_C = {F_out_C_tmol + F_cont_C_tmol:.2f} Tmol/yr')
print(f'  For Ca balance:  F_burial_Ca = F_cont_Ca = {F_cont_Ca_tmol:.2f} Tmol/yr')
print(f'  C/Ca mismatch in burial: {(F_out_C_tmol+F_cont_C_tmol)/F_cont_Ca_tmol:.2f}×  (should be 1.0 for CaCO3)')
print()
print(f'  For Alk balance: F_burial_Alk = 2 × F_burial_C = {2*(F_out_C_tmol+F_cont_C_tmol):.2f} Tmol/yr')
print(f'  Available Alk:   F_cont_Alk + F_diss_Alk = {F_cont_alk_tmol + F_diss[_alk_idx]*YR/1e12:.2f} Tmol/yr')
print(f'  Alk shortfall: {2*(F_out_C_tmol+F_cont_C_tmol) - (F_cont_alk_tmol + F_diss[_alk_idx]*YR/1e12):.2f} Tmol/yr')
print()
print('  Diagnosis: F_cont adds equal C and Alk, but burial removes 2 Alk per 1 C.')
print('  The only way to close the Alk budget is:')
alk_needed_from_seafloor = 2*(F_out_C_tmol + F_cont_C_tmol) - F_cont_alk_tmol
print(f'    (a) Seafloor Alk flux >= {alk_needed_from_seafloor:.1f} Tmol/yr (currently {F_diss[_alk_idx]*YR/1e12:.2f})')
print(f'    (b) F_cont_C should be 0, not = F_cont_Alk (rivers deliver Alk not C directly)')
print(f'        → then Alk shortfall = {2*F_out_C_tmol - F_cont_alk_tmol:.2f} Tmol/yr')
print(f'    (c) Reduce F_out_C to {F_cont_Ca_tmol:.2f} Tmol/yr (= Ca supply from weathering)')
print(f'        → represents only the silicate-weathering-paired outgassing')

# ─────────────────────────────────────────────────────────────────────────────
sep('STEP 5 — FIX OPTIONS ASSESSMENT')
# ─────────────────────────────────────────────────────────────────────────────
print('Testing specific fix hypotheses to see which closes the budget')
print()

# Fix A: Don't add C from continental weathering (rivers deliver Alk+Ca but
#        C was drawn from atmosphere — the ocean-atm equilibrium handles it)
print('--- Fix A: F_cont_C = 0 (rivers add Alk+Ca only, C tracked via equilibrium) ---')
F_net_fixA_C = F_out[_C_idx] + 0 + F_diss[_C_idx]
F_prec_C_needed_A = -F_net_fixA_C
print(f'  Needed F_prec_C = {F_prec_C_needed_A*YR/1e12:.2f} Tmol/yr')
F_alk_avail_A = F_cont[_alk_idx] + F_diss[_alk_idx]
F_alk_from_burial_A = 2 * F_prec_C_needed_A * (-1)  # negative because removal
print(f'  Alk removed by burial = {abs(F_alk_from_burial_A)*YR/1e12:.2f} Tmol/yr')
print(f'  Alk supplied = {F_alk_avail_A*YR/1e12:.2f} Tmol/yr')
print(f'  Alk balance: {"CLOSED" if abs(F_alk_avail_A*YR/1e12 - abs(F_alk_from_burial_A)*YR/1e12) < 1 else "STILL OPEN by %.2f Tmol/yr" % abs((F_alk_avail_A + F_alk_from_burial_A)*YR/1e12)}')

print()
print('--- Fix B: Halve F_out_C (organic burial accounts for other half) ---')
F_out_C_half = F_out[_C_idx] * 0.5
F_net_fixB_C = F_out_C_half + F_cont[_C_idx] + F_diss[_C_idx]
F_prec_C_needed_B = -F_net_fixB_C
F_alk_from_burial_B = 2 * F_prec_C_needed_B * (-1)
print(f'  Needed F_prec_C = {F_prec_C_needed_B*YR/1e12:.2f} Tmol/yr')
print(f'  Alk removed by burial = {abs(F_alk_from_burial_B)*YR/1e12:.2f} Tmol/yr')
print(f'  Alk supplied = {F_alk_avail_A*YR/1e12:.2f} Tmol/yr')
print(f'  Alk balance: {"CLOSED" if abs((F_alk_avail_A - abs(F_alk_from_burial_B)/(-1))*YR/1e12) < 1 else "STILL OPEN by %.2f Tmol/yr" % abs((F_alk_avail_A + F_alk_from_burial_B)*YR/1e12)}')

print()
print('--- Fix C: Fix A + Fix B combined ---')
F_net_fixC_C = F_out_C_half + 0 + F_diss[_C_idx]
F_prec_C_needed_C = -F_net_fixC_C
F_alk_from_burial_C = 2 * abs(F_prec_C_needed_C)
print(f'  Needed F_prec_C = {F_prec_C_needed_C*YR/1e12:.2f} Tmol/yr')
print(f'  Alk removed by burial = {F_alk_from_burial_C*YR/1e12:.2f} Tmol/yr')
print(f'  Alk supplied = {F_alk_avail_A*YR/1e12:.2f} Tmol/yr')
alk_balance_C = F_alk_avail_A*YR/1e12 - F_alk_from_burial_C*YR/1e12
print(f'  Alk balance: {"CLOSED" if abs(alk_balance_C) < 1 else "OPEN by %.2f Tmol/yr" % alk_balance_C}')

print()
print('--- Precipitation magnitude check for Fix C ---')
print(f'  Target F_prec_C = {abs(F_prec_C_needed_C)*YR/1e12:.2f} Tmol/yr')
print(f'  → delta_b_C needed = {abs(F_prec_C_needed_C)*tau_prec/OCEAN_MASS:.4e} mol/kgw')
print(f'  Actual delta_b_C at b_init (surface) = {delta_b_prec_surf[_C_idx]:.4e} mol/kgw')
ratio = abs(F_prec_C_needed_C)*tau_prec/OCEAN_MASS / abs(delta_b_prec_surf[_C_idx]) if delta_b_prec_surf[_C_idx] != 0 else float('inf')
print(f'  Still need {ratio:.0f}× more precipitation signal than currently available.')
print(f'  → SI(Calcite) = 0.04 is too low. Ocean must be much more supersaturated.')
print()
print('  What Ca concentration gives adequate precipitation?')
print(f'  {"Ca(mmol)":>9} {"SI_Cal":>7} {"db_C":>14} {"F_prec_C(Tmol/yr)":>18} {"vs target":>10}')
target_db_C = abs(F_prec_C_needed_C) * tau_prec / OCEAN_MASS
for ca_test in [7e-5, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1]:
    b_test = b_ref.copy()
    b_test[_Ca_idx] = ca_test
    try:
        out2 = solve_solution(P_surf_ref, T_ref, b_test, verbose=False)
        import pandas as pd
        df2 = pd.DataFrame(out2)
        si_cal2 = float(df2['si_Calcite'].iloc[-1]) if 'si_Calcite' in df2.columns else float('nan')
        db2 = get_precipitation_flux(P_surf_ref, T_ref, b_test,
                                      precipitating_minerals=carbonate_minerals)
        db_C2 = db2[_C_idx]
        fp_C2 = abs(db_C2) * OCEAN_MASS * YR / (tau_prec * 1e12)
        print(f'  {ca_test*1e3:>9.2f} {si_cal2:>7.3f} {db_C2:>14.4e} {fp_C2:>18.3f} {fp_C2/abs(F_prec_C_needed_C*YR/1e12)*100:>9.1f}%')
    except Exception as ex:
        print(f'  {ca_test*1e3:>9.2f}  FAILED: {ex}')

sep('SUMMARY')
print("""
  ROOT CAUSES IDENTIFIED:

  1. STOICHIOMETRIC IMBALANCE (structural):
     F_cont adds equal Alk and C (WHK: CaSiO3+2CO2 → Ca+2HCO3-).
     CaCO3 burial removes 2 Alk per 1 C.
     With F_out=7.5 and F_cont_C=8 Tmol/yr, burial must remove 15.5 Tmol C/yr,
     which requires 31 Tmol/yr Alk removal — but only 8+1=9 Tmol/yr is supplied.
     Alk shortfall: ~22 Tmol/yr. Model CANNOT reach steady state as calibrated.

  2. MISSING ORGANIC CARBON (structural):
     Earth's outgassing (7.5 Tmol/yr) is half balanced by organic burial.
     Without it, inorganic silicate weathering alone cannot close the C budget.
     Fix: Either halve F_out_C, or add organic burial scaling with land_fraction.

  3. F_cont_C DOUBLE-COUNTS ATMOSPHERIC CO2 (conceptual):
     Rivers deliver HCO3- (DIC) to ocean, but that C came FROM the atmosphere.
     Adding it to ocean DIC is correct ONLY IF the atmospheric CO2 loss is also
     tracked. In this model, the ocean IS the atmosphere (P_CO2 from PHREEQC),
     so adding F_cont_C raises P_CO2 instead of lowering it — backward!
     Fix: Set F_cont_C = 0 and let ocean-atmosphere equilibrium handle it.

  4. PRECIPITATION SIGNAL TOO SMALL:
     Even with a balanced budget, delta_b_C at SI=0.04 is negligible.
     Need SI >> 1 (achieved with high Ca, ~10-50 mmol/kgw modern seawater Ca)
     for adequate precipitation to act as the thermostat.

  RECOMMENDED FIXES (in priority order):
  A. Set F_cont_C = 0 in get_continental_weathering_flux
     (rivers add Alk and Ca, not DIC — atmosphere-ocean equilibrium handles C)
  B. Scale outgassing to inorganic-only: multiply EARTH_OUTGASSING by ~0.5,
     OR add an organic burial term proportional to land_fraction.
  C. Increase initial Ca to modern seawater (~10 mmol/kgw) for adequate
     precipitation signal — but check PHREEQC charge balance.
""")
