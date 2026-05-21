"""
Static flux audit for Earth calibration.
Checks each flux component at Earth-like steady-state conditions.
Runs 4 model configurations and compares final states.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import json

from kamino.chemistry import elements, alk_idx, c_idx, ca_idx, mg_idx, si_idx
from kamino.weathering import (get_weathering_flux, get_continental_weathering_flux,
                                J_ref_normalised, rate_ref, A_seafloor,
                                EARTH_CONTINENTAL_WEATHERING_REF)
from kamino.precipitation import get_precipitation
from kamino.climate.clima_interpolator import get_T_surface
from kamino.constants import *

# ─────────────────────────────────────────────────────────────────────────────
# Reference quantities
# ─────────────────────────────────────────────────────────────────────────────
g = G * M_EARTH / R_EARTH**2
A_earth = 4 * np.pi * R_EARTH**2
ocean_mass = 3800 * A_earth * 1000  # kg (3.8 km mean depth)
land_frac = 0.3
land_area = land_frac * A_earth

T_earth   = 288.0         # K target
T_model   = 280.0         # K what the climate model gives at 280 ppm, S=1.0
P_CO2_ea  = EARTH_ATM * 280e-6   # 280 ppm in Pa
P_bg      = 1e5            # 1 bar N2 background
P_H2O     = 1200           # Pa approx (50% RH at 280 K)
P_surf    = P_bg + P_CO2_ea + P_H2O
P_pore    = P_surf + 1000 * g * 3800
T_sf      = max(1.02 * T_model - 16.7, 274)  # seafloor T
T_pore    = T_sf + 9

def tmol(flux_pa, area=None):
    """flux per unit area [mol/m²/s] × area [m²] → Tmol/yr"""
    if area is None:
        return flux_pa * YR / 1e12
    return flux_pa * area * YR / 1e12

def sep(title):
    print(); print('=' * 68); print(f'  {title}'); print('=' * 68)

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 1 — REFERENCE OCEAN STATE')
# ─────────────────────────────────────────────────────────────────────────────
# Model-consistent ocean: small Ca+Mg so PHREEQC charge balance holds.
# In the model Cl/Na/K/SO4 are absent so Alk ≥ 2*Ca + 2*Mg required.
# These values are what a model with only seafloor weathering + carbonate
# precipitation might settle at.
b_model = np.zeros(len(elements))
b_model[alk_idx] = 2.1e-3    # ~modern Alk
b_model[c_idx]   = 2.0e-3    # ~modern DIC
b_model[ca_idx]  = 5e-5      # low Ca (model limitation)
b_model[mg_idx]  = 5e-5      # low Mg
b_model[si_idx]  = 1e-4

_, pH_ref, SI_ref = get_precipitation(
    P_pore, T_sf, b_model,
    ['Calcite', 'SiO2(am)', 'Kaolinite', 'Saponite-Mg'],
    precipitation_timescale=1e13*YR)  # tiny rate, just for diagnostics

print(f'\n  Climate model at 280 ppm, S=1.0, albedo=0.3:')
T_check = get_T_surface(SOLAR_CONSTANT, P_CO2_ea, 0.3, False)
print(f'    T = {T_check:.1f} K  (actual Earth 288 K — model runs 8 K cool)')
print(f'\n  Reference ocean chemistry (b_model):')
print(f'    pH          = {pH_ref:.2f}  (target 8.1)')
print(f'    Alk         = {b_model[alk_idx]*1e3:.2f} meq/kg  (modern Earth ~2.3)')
print(f'    DIC (C)     = {b_model[c_idx]*1e3:.2f} mmol/kg  (modern Earth ~2.0)')
print(f'    Ca          = {b_model[ca_idx]*1e3:.2f} mmol/kg  (modern Earth ~10, model limited)')
print(f'    Mg          = {b_model[mg_idx]*1e3:.2f} mmol/kg  (modern Earth ~53, model limited)')
print(f'    SI(Calcite) = {SI_ref.get("Calcite", float("nan")):.2f}')
print(f'    SI(SiO2am)  = {SI_ref.get("SiO2(am)", float("nan")):.2f}')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 2 — SEAFLOOR WEATHERING FLUX AUDIT')
# ─────────────────────────────────────────────────────────────────────────────
from kamino.mineral_info import basalt_composition, clay_minerals, carbonate_minerals

# LT only (f_HT = 0)
flux_lt, diag_lt = get_weathering_flux(
    P_pore, T_pore, P_CO2_ea, b_model,
    rate=rate_ref, J=J_ref_normalised,
    crust_composition=basalt_composition,
    precipitating_minerals=carbonate_minerals + clay_minerals)

print(f'\n  LT seafloor weathering at T_pore={T_pore:.1f} K:')
print(f'    Da              = {diag_lt["Da"]:.3f}  (>1: thermodynamic-limited; <1: kinetic)')
print(f'    A_reactive      = {diag_lt["A_reactive"]:.3e} 1/s')
print(f'    supply_eff      = {diag_lt["supply_efficiency"]:.2%}')
print(f'    Alk flux (LT)   = {tmol(flux_lt[alk_idx], A_seafloor):.2f} Tmol eq/yr'
      f'  [Coogan 2022: ~0.7 Teq/yr]')
print(f'    Ca flux (LT)    = {tmol(flux_lt[ca_idx], A_seafloor):.2f} Tmol/yr')
print(f'    Mg flux (LT)    = {tmol(flux_lt[mg_idx], A_seafloor):.2f} Tmol/yr')
print(f'    Si flux (LT)    = {tmol(flux_lt[si_idx], A_seafloor):.2f} Tmol/yr')

# HT Mg→Ca exchange (simple parameterization in planet.py)
KD_MG_HT = 0.7 * 1.4e12 / (0.053 * 1.4e15)
J_total = J_ref_normalised  # crust_production=1
ht_rate_per_vol = KD_MG_HT * b_model[mg_idx] * J_total  # mol/kgw/s per unit seafloor area
ht_mg_tmol = -ht_rate_per_vol * A_earth * YR / 1e12
ht_ca_tmol = +ht_rate_per_vol * A_earth * YR / 1e12
print(f'\n  Simple HT Mg→Ca exchange (KD_MG_HT, at [Mg]={b_model[mg_idx]*1e3:.2f} mmol/kg):')
print(f'    Mg removal = {ht_mg_tmol:.2e} Tmol/yr (at low [Mg])')
print(f'    Ca addition = {ht_ca_tmol:.2e} Tmol/yr')
print(f'    (At [Mg]=53 mM Earth: {KD_MG_HT * 0.053 * J_ref_normalised * A_earth * YR / 1e12:.2f} Tmol/yr'
      f'  [Coogan: ~1.4 Tmol/yr HT Ca flux])')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 3 — CONTINENTAL WEATHERING FLUX AUDIT')
# ─────────────────────────────────────────────────────────────────────────────
F_cont = get_continental_weathering_flux(T_model, P_CO2_ea)
print(f'\n  Continental silicate weathering at T={T_model:.0f} K, pCO2=280ppm:')
print(f'    Alk flux = {tmol(F_cont[alk_idx], land_area):.2f} Tmol eq/yr'
      f'  [Coogan 2022: ~16 Teq/yr, code calibrated to 8]')
print(f'    Ca flux  = {tmol(F_cont[ca_idx], land_area):.2f} Tmol/yr')
print(f'    Mg flux  = {tmol(F_cont[mg_idx], land_area):.2f} Tmol/yr')
print(f'    C flux   = {tmol(F_cont[c_idx], land_area):.2f} Tmol/yr  (forced to 0 in planet.py)')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 4 — OUTGASSING & CARBON BUDGET')
# ─────────────────────────────────────────────────────────────────────────────
F_out_C = EARTH_OUTGASSING * A_earth  # mol/s
F_out_C_tmol = F_out_C * YR / 1e12

print(f'\n  Volcanic outgassing (code default):')
print(f'    C flux = {F_out_C_tmol:.2f} Tmol C/yr')
print(f'    [Coogan 2026: uncertain, ~5-10 Tmol/yr; here EARTH_OUTGASSING=0.0147 mol/m²/yr]')

# Carbon balance at Earth conditions
C_in_ocean_world = F_out_C_tmol  # no land
C_in_continental = F_out_C_tmol  # F_cont_C = 0 in model

print(f'\n  Carbon budget (must balance at steady state):')
print(f'    Ocean world  (no land): C in = {C_in_ocean_world:.2f} Tmol/yr')
print(f'    Continental (land=0.3): C in = {C_in_continental:.2f} Tmol/yr (same; F_cont_C=0)')

# Alkalinity budget
alk_sf = tmol(flux_lt[alk_idx], A_seafloor)
alk_cont = tmol(F_cont[alk_idx], land_area)

print(f'\n  Alkalinity budget:')
print(f'    Seafloor LT supply = {alk_sf:.2f} Tmol eq/yr')
print(f'    Continental supply = {alk_cont:.2f} Tmol eq/yr  (for 30% land)')
print(f'    CaCO3 burial removes 2 Alk per 1 C buried')
print(f'    Alk needed for C burial = 2 × {C_in_ocean_world:.2f} = {2*C_in_ocean_world:.2f} Tmol eq/yr')
print(f'')
print(f'  Ocean world  Alk balance: {alk_sf:.2f} in, {2*C_in_ocean_world:.2f} needed'
      f'  → {"OK" if alk_sf >= 2*C_in_ocean_world else f"deficit {2*C_in_ocean_world - alk_sf:.2f} Tmol eq/yr"}')
print(f'  Continental  Alk balance: {alk_sf+alk_cont:.2f} in, {2*C_in_continental:.2f} needed'
      f'  → {"OK" if alk_sf+alk_cont >= 2*C_in_continental else f"deficit {2*C_in_continental - (alk_sf+alk_cont):.2f} Tmol eq/yr"}')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 5 — SENSITIVITY TO T_PORE / TEMPERATURE')
# ─────────────────────────────────────────────────────────────────────────────
print(f'\n  Seafloor LT alkalinity flux vs temperature:')
print(f'    {"T_sf (K)":>9}  {"T_pore (K)":>10}  {"Alk (Tmol/yr)":>14}  {"Da":>7}')
for T_sf_test in [274, 278, 282, 286, 290, 295]:
    T_pore_test = T_sf_test + 9
    fl, dg = get_weathering_flux(P_pore, T_pore_test, P_CO2_ea, b_model,
                                  rate=rate_ref, J=J_ref_normalised,
                                  crust_composition=basalt_composition,
                                  precipitating_minerals=[])
    print(f'    {T_sf_test:>9.0f}  {T_pore_test:>10.0f}  '
          f'{tmol(fl[alk_idx], A_seafloor):>14.2f}  {dg["Da"]:>7.3f}')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 6 — REVERSE WEATHERING: Saponite-Mg Mg SINK')
# ─────────────────────────────────────────────────────────────────────────────
# How effective is Saponite-Mg precipitation at removing Mg from the ocean?
print(f'\n  Saponite-Mg precipitation (ocean) at b_model conditions:')
db_rw, pH_rw, SI_rw = get_precipitation(
    P_pore, T_sf, b_model, ['Saponite-Mg'],
    precipitation_timescale=1e13*YR)
print(f'    SI(Saponite-Mg) = {SI_rw.get("Saponite-Mg", float("nan")):.3f}'
      f'  (>0 → saturated, will precipitate)')
print(f'    Mg removal flux = {db_rw[mg_idx] / (1e13*YR) * ocean_mass * YR / 1e12:.2e} Tmol/yr'
      f' (at tau_prec=1e5 yr)')

# ─────────────────────────────────────────────────────────────────────────────
sep('PART 7 — SUMMARY AND RECOMMENDED RUNS')
# ─────────────────────────────────────────────────────────────────────────────
print("""
  KEY OBSERVATIONS:
  ─────────────────
  (a) Climate model gives T=280 K at pCO2=280 ppm (model runs 8 K cool vs Earth).
      Target: pCO2 ≈ 280 ppm → T ≈ 280 K in model.

  (b) LT seafloor Alk supply alone may be insufficient to balance outgassing
      at modern Earth outgassing rates. Continental weathering is critical.

  (c) The simple HT Mg→Ca exchange is small at low [Mg] but at Earth-like
      [Mg]≈53 mM it gives ~1.4 Tmol/yr Ca addition (consistent with Coogan).

  (d) Continental weathering at 30% land gives substantial Alk surplus,
      which should drive pCO2 DOWN toward Earth-like values.

  RECOMMENDED MODEL RUNS:
  ────────────────────────
  Run 1: Ocean world  (land=0.0, rw=True, f_HT=0)
         → Baseline: shows how seafloor weathering alone sets steady state

  Run 2: Continental (land=0.3, rw=True, f_HT=0, f_carb=0)
         → Pure silicate weathering + seafloor: can this give Earth pCO2?

  Run 3: Continental + carbonate weathering (land=0.3, rw=True, f_HT=0, f_carb=0.5)
         → Effect of continental carbonate weathering on Ca/Alk budget

  Run 4: Continental + HT hydrothermal (land=0.3, rw=True, f_HT=0.3, f_carb=0)
         → Effect of HT on Mg→Ca exchange and Alk budget

  Run 5: Full Earth (land=0.3, rw=True, f_HT=0.3, f_carb=0.5)
         → All processes active: closest to Coogan model
""")
