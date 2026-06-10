"""
Check seafloor weathering Alk flux at two states:
  1. The bio_fc0.0 converged state (Alk=17, Ca=5)
  2. Earth target conditions (Alk=2.3, Ca=10.3)
to see whether the seafloor is giving anomalously high flux away from calibration.
"""
import sys, numpy as np
sys.path.insert(0, '/home/pt426/Code/kamino/src')

from kamino.constants import R_EARTH, YR, G, M_EARTH, EARTH_CRUST_PRODUCTION_RATE_PER_AREA
from kamino.chemistry import elements, alk_idx, ca_idx, mg_idx, si_idx, na_idx, cl_idx, c_idx
from kamino.weathering import get_weathering_flux, J_ref_normalised, rate_ref, ALPHA_REF
from kamino.mineral_info import *

OCEAN_DEPTH = 3700
BACKGROUND_PRESSURE = 1e5

gravity = G * M_EARTH / R_EARTH**2
surface_area = 4 * np.pi * R_EARTH**2
ocean_water_mass = OCEAN_DEPTH * surface_area * 1000

crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * 1.0
J_norm = J_ref_normalised * (crust_production_rate / rate_ref)

def seafloor_alk(P_CO2_pa, T_surface_K, b_ocean_labels, b_ocean_vals):
    T_sf = max(1.02 * T_surface_K - 16.7, 274.0)
    T_pore = T_sf + 9
    P_pore = BACKGROUND_PRESSURE + P_CO2_pa + 1000 * gravity * OCEAN_DEPTH

    b_ocean = np.zeros(len(elements))
    for k, v in zip(b_ocean_labels, b_ocean_vals):
        b_ocean[k] = v

    pore_prec = [m for m in carbonate_minerals if m != 'Calcite'] + clay_minerals
    flux_lt, diag = get_weathering_flux(P_pore, T_pore, P_CO2_pa, b_ocean,
                                        rate=crust_production_rate, J=J_norm,
                                        precipitating_minerals=pore_prec)
    F_diss = (flux_lt * surface_area) / ocean_water_mass
    F_diss[alk_idx] -= F_diss[na_idx]  # balanced Na cycle correction

    alk_tmol = F_diss[alk_idx] * ocean_water_mass * YR / 1e12
    ca_tmol  = F_diss[ca_idx]  * ocean_water_mass * YR / 1e12
    mg_tmol  = F_diss[mg_idx]  * ocean_water_mass * YR / 1e12
    si_tmol  = F_diss[si_idx]  * ocean_water_mass * YR / 1e12

    print(f"  T_surface={T_surface_K:.1f}K  T_pore={T_pore:.1f}K  P_CO2={P_CO2_pa:.1f}Pa ({P_CO2_pa/101.325:.0f}ppm)")
    print(f"  Seafloor Alk = {alk_tmol:+.3f} Tmol/yr")
    print(f"  Seafloor Ca  = {ca_tmol:+.3f} Tmol/yr")
    print(f"  Seafloor Mg  = {mg_tmol:+.3f} Tmol/yr")
    print(f"  Seafloor Si  = {si_tmol:+.3f} Tmol/yr")
    print(f"  Da (Mg proxy)= {diag['Da']:.3f}")
    return alk_tmol

print("=== 1. bio_fc0.0 converged state (Alk=17.34 meq/kg, Ca=5.03 mmol/kg) ===")
# pCO2=823ppm → 823e-6 atm × 101325 Pa = 83.4 Pa
seafloor_alk(83.4, 287.5,
    [alk_idx, ca_idx, mg_idx, si_idx, na_idx, cl_idx, c_idx],
    [17.34e-3, 5.03e-3, 56.8e-3, 0.965e-3, 525.6e-3, 559.1e-3, 14.43e-3])

print()
print("=== 2. Earth target conditions (Alk=2.3 meq/kg, Ca=10.3 mmol/kg) ===")
# pCO2=280ppm → 280e-6 × 101325 = 28.4 Pa
seafloor_alk(28.4, 288.0,
    [alk_idx, ca_idx, mg_idx, si_idx, na_idx, cl_idx, c_idx],
    [2.3e-3, 10.3e-3, 52.8e-3, 0.1e-3, 480e-3, 550e-3, 2.0e-3])

print()
print("=== 3. Calibration conditions (empty ocean, T=280K, pCO2=280ppm) ===")
seafloor_alk(28.4, 273.0,   # T_surface that gives T_sf=282, T_pore≈291? Let's use direct T_surface so T_pore=280+9=289
    [], [])
# Actually the calibration uses T_pore=280 directly. Let's set T_surface such that T_sf=271, T_pore=280
# T_sf = 1.02*T_s - 16.7 = 280 → T_s = (280+16.7)/1.02 = 290.9K
seafloor_alk(28.4, 290.9, [], [])
