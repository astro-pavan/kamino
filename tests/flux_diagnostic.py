"""
Flux budget diagnostic at the bio_fc0.38 converged state.
Computes each term in the Alk and Ca balance to diagnose the true steady state.
"""
import sys, json
import numpy as np

sys.path.insert(0, '/home/pt426/Code/kamino/src')

from kamino.constants import M_EARTH, R_EARTH, YR, EARTH_OUTGASSING, EARTH_CRUST_PRODUCTION_RATE_PER_AREA, G
from kamino.chemistry import elements, alk_idx, c_idx, ca_idx, mg_idx, si_idx, na_idx, cl_idx
from kamino.weathering import get_weathering_flux, get_continental_weathering_flux, J_ref_normalised, rate_ref
from kamino.precipitation import get_precipitation
from kamino.mineral_info import *
from kamino.climate.clima_interpolator import get_T_surface
from kamino.utils import august_roche_magnus_formula
from kamino.planet import KD_MG_HT, K_BIO_CA, K_BIO_SI, K_NA_REVERSE_WEATHERING, K_CL_SUBDUCTION

BACKGROUND_PRESSURE = 1e5
OCEAN_DEPTH = 3700

# Load saved results
configs_to_check = ['bio_fc0.0', 'bio_fc0.38', 'abiotic']

for name in configs_to_check:
    with open(f'/home/pt426/Code/kamino/output/{name}.json') as f:
        out = json.load(f)

    # Parameters from the saved planet config (stored at top of json)
    cfg = out
    f_bio = cfg['f_bio']
    f_carb = cfg['f_carb']
    land_fraction = cfg['land_fraction']
    outgassing = cfg['outgassing']
    cl_outgassing_ratio = cfg['cl_outgassing_ratio']
    cl_subduction = cfg['cl_subduction']
    tau_prec = cfg['tau_prec_yr'] * YR

    # Final state
    y_end = np.array(out['data']['y'])[:, -1]
    P_CO2 = y_end[0]
    P_H2O = y_end[1]
    b_ocean = y_end[2:-1]

    # Planet geometry
    mass = M_EARTH; radius = R_EARTH
    gravity = G * mass / radius**2
    surface_area = 4 * np.pi * radius**2
    ocean_water_mass = OCEAN_DEPTH * surface_area * 1000
    land_area = land_fraction * surface_area

    P_surface = BACKGROUND_PRESSURE + P_CO2 + P_H2O
    T_surface = get_T_surface(1.0 * 1361, max(P_CO2, 1.0), 0.3, False)
    T_seafloor = max(1.02 * T_surface - 16.7, 274)
    T_pore = T_seafloor + 9
    P_pore = P_surface + 1000 * gravity * OCEAN_DEPTH

    crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * 1.0
    J_total_norm = J_ref_normalised * (crust_production_rate / rate_ref)
    J_total = J_total_norm  # per unit area; ×surface_area to get global

    # Precipitation / SI
    rw_minerals = carbonate_minerals + clay_minerals + silica_minerals + reverse_weathering_minerals + evaporite_minerals
    pore_prec_minerals = [m for m in carbonate_minerals if m != 'Calcite'] + clay_minerals
    shelf_prec_minerals = carbonate_minerals

    F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean,
                                       precipitating_minerals=rw_minerals,
                                       precipitation_timescale=tau_prec)

    # Continental weathering
    F_sil_cont = get_continental_weathering_flux(T_surface, P_CO2).copy()
    F_carb_w = np.zeros(elements.shape)
    if f_carb > 0:
        F_carb_w[ca_idx]  = f_carb * F_sil_cont[ca_idx]
        F_carb_w[alk_idx] = 2.0 * f_carb * F_sil_cont[ca_idx]
        F_carb_w[c_idx]   = f_carb * F_sil_cont[ca_idx]
    F_cont = (F_sil_cont + F_carb_w) * land_area / ocean_water_mass

    # Shelf precipitation
    P_shelf = P_surface + 1000 * gravity * 1000.0
    F_shelf_prec, _, _ = get_precipitation(P_shelf, T_seafloor, b_ocean,
                                           precipitating_minerals=shelf_prec_minerals,
                                           precipitation_timescale=tau_prec)

    # Seafloor weathering
    F_carb_sed = max(0.0, -F_prec[c_idx])
    F_sil_sed  = max(0.0, -F_prec[si_idx])
    ocean_water_per_area = OCEAN_DEPTH * 1000.0
    S_sed = (F_carb_sed * 0.100/2710.0 + F_sil_sed * 0.060/2650.0) * ocean_water_per_area

    flux_lt, _ = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean,
                                     rate=crust_production_rate, J=J_total_norm,
                                     sedimentation_rate=S_sed,
                                     precipitating_minerals=pore_prec_minerals)
    F_diss = (flux_lt * surface_area) / ocean_water_mass
    F_diss[alk_idx] -= F_diss[na_idx]  # balanced Na cycle

    # HT Mg→Ca exchange
    _ht_rate = KD_MG_HT * b_ocean[mg_idx] * J_total * surface_area / ocean_water_mass
    F_ht = np.zeros(elements.shape)
    F_ht[mg_idx] = -_ht_rate
    F_ht[ca_idx] = +_ht_rate

    # Outgassing
    F_vol = np.zeros(elements.shape)
    F_vol[c_idx]   = (EARTH_OUTGASSING / YR) * surface_area * outgassing / ocean_water_mass
    F_vol[cl_idx]  = F_vol[c_idx] * cl_outgassing_ratio
    F_vol[alk_idx] = -F_vol[cl_idx]

    # Cl subduction
    F_cl_subd = np.zeros(elements.shape)
    F_cl_subd[cl_idx] = -cl_subduction * K_CL_SUBDUCTION * b_ocean[cl_idx] * J_total * surface_area / ocean_water_mass

    # Na reverse weathering
    F_na_rw = np.zeros(elements.shape)
    F_na_rw[na_idx] = -K_NA_REVERSE_WEATHERING * b_ocean[na_idx] * J_total * surface_area / ocean_water_mass

    # Biogenic burial
    F_bio = np.zeros(elements.shape)
    if f_bio > 0 and SI.get('Calcite', -1.0) > 0.0:
        _bio_ca = f_bio * K_BIO_CA * b_ocean[ca_idx]
        _bio_si = f_bio * K_BIO_SI * b_ocean[si_idx]
        F_bio[ca_idx]  = -_bio_ca
        F_bio[alk_idx] = -2.0 * _bio_ca
        F_bio[c_idx]   = -_bio_ca
        F_bio[si_idx]  = -_bio_si

    # Convert to Tmol/yr globally (multiply mol/kgw/s × ocean_mass × YR / 1e12)
    scale = ocean_water_mass * YR / 1e12

    def tmol(F, idx):
        return F[idx] * scale

    print(f"\n{'='*70}")
    print(f"Config: {name}  (f_bio={f_bio}, f_carb={f_carb})")
    print(f"  State: T={T_surface:.1f}K  P_CO2={P_CO2:.1f}Pa ({P_CO2/101325*1e6:.0f}ppm)  pH={pH:.2f}")
    print(f"  b_ocean: Ca={b_ocean[ca_idx]*1e3:.2f}  Mg={b_ocean[mg_idx]*1e3:.1f}  Alk={b_ocean[alk_idx]*1e3:.2f}  DIC={b_ocean[c_idx]*1e3:.2f}  Si={b_ocean[si_idx]*1e3:.3f}  mmol/kg")
    print(f"  Calcite SI = {SI.get('Calcite', float('nan')):.3f}")
    print()
    print(f"  Alk budget (Tmol/yr):  [positive = input, negative = sink]")
    print(f"    Continental sil:    {tmol(F_cont, alk_idx):+.3f}")
    print(f"    Carb weathering:    {tmol(F_carb_w * land_area / ocean_water_mass, alk_idx):+.3f}")
    print(f"    Seafloor weathering:{tmol(F_diss, alk_idx):+.3f}")
    print(f"    Outgassing (HCl):   {tmol(F_vol, alk_idx):+.3f}")
    print(f"    Ocean precipitation:{tmol(F_prec, alk_idx):+.3f}")
    print(f"    Shelf precipitation:{tmol(F_shelf_prec, alk_idx):+.3f}")
    print(f"    Biogenic burial:    {tmol(F_bio, alk_idx):+.3f}")
    alk_net = tmol(F_cont + F_carb_w*land_area/ocean_water_mass + F_diss + F_vol + F_prec + F_shelf_prec + F_bio, alk_idx)
    print(f"    NET:                {alk_net:+.3f}  {'<-- Alk accumulating' if alk_net > 0 else '<-- Alk declining'}")
    print()
    print(f"  Ca budget (Tmol/yr):")
    print(f"    Continental:        {tmol(F_cont + F_carb_w*land_area/ocean_water_mass, ca_idx):+.3f}")
    print(f"    Seafloor:           {tmol(F_diss, ca_idx):+.3f}")
    print(f"    HT exchange:        {tmol(F_ht, ca_idx):+.3f}")
    print(f"    Ocean precipitation:{tmol(F_prec, ca_idx):+.3f}")
    print(f"    Shelf precipitation:{tmol(F_shelf_prec, ca_idx):+.3f}")
    print(f"    Biogenic burial:    {tmol(F_bio, ca_idx):+.3f}")
    ca_net = tmol(F_cont + F_carb_w*land_area/ocean_water_mass + F_diss + F_ht + F_prec + F_shelf_prec + F_bio, ca_idx)
    print(f"    NET:                {ca_net:+.3f}  {'<-- Ca accumulating' if ca_net > 0 else '<-- Ca declining'}")
    print()
    # What Calcite SI would be at Earth conditions
    print(f"  Biogenic active: {'YES' if (f_bio>0 and SI.get('Calcite',-1)>0) else 'NO (SI<0 or f_bio=0)'}")
    print(f"  Ca residence time (biogenic): {b_ocean[ca_idx]/(max(abs(tmol(F_bio,ca_idx)),1e-9)/1e12*ocean_water_mass*YR)*YR/1e6:.1f} Myr" if f_bio>0 else "  Ca residence time (biogenic): N/A")
