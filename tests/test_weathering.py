import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

# Adjust these imports if your directory structure requires it
from kamino.constants import *
from kamino.kamino_chem.ocean_chemistry import *

def run_all_benchmarks():
    print("--- Starting Ocean Chemistry Benchmarks ---\n")
    
    # Modern Earth Seawater Approximation
    b_seawater = np.array([2.3e-3, 2.0e-3, 1e-4, 0.0, 0.0, 10.3e-3, 52.7e-3, 468e-3, 546e-3])
    
    # Standard Earth Seafloor Conditions
    P_seafloor = 1e5 + 1000 * 9.81 * 3000  # ~300 atm
    T_seafloor = 275  # ~2°C
    alpha_ref = 2 # From your least_squares comment
    
    # Indices for elements
    idx_Alk = int(np.where(elements == 'Alkalinity')[0][0])
    idx_C = int(np.where(elements == 'C')[0][0])
    idx_Ca = int(np.where(elements == 'Ca')[0][0])
    idx_Mg = int(np.where(elements == 'Mg')[0][0])
    idx_Na = int(np.where(elements == 'Na')[0][0])
    idx_Cl = int(np.where(elements == 'Cl')[0][0])

    # ---------------------------------------------------------
    # Test 1: Modern Earth Baseline
    # ---------------------------------------------------------
    print("1. Testing Modern Earth Baseline...")
    P_CO2_modern = get_P_CO2(1e5, 288, b_seawater)
    print(f"   Calculated Modern P_CO2: {P_CO2_modern / EARTH_ATM * 1e6:.1f} ppm (Expected: ~280-400 ppm)")
    
    # Calculate modern weathering flux
    flux, Da = get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_seawater, alpha_ref, rate_ref)
    
    # Total C flux in Tmol/yr
    A_seafloor = 0.7 * 4 * np.pi * R_EARTH ** 2
    flux_C_Tmol_yr = (flux[idx_C] * A_seafloor * YR) / 1e12
    print(f"   Calculated C Sink Flux: {flux_C_Tmol_yr:.4f} Tmol/yr (Expected: Negative, ~ -1.5 to -3.0)\n")

    # ---------------------------------------------------------
    # Test 2: Carbonate System Thermodynamics (Bjerrum Plot)
    # ---------------------------------------------------------
    print("2. Testing Carbonate Thermodynamics (Alkalinity Sweep)...")
    alk_multipliers = np.linspace(0.5, 2.0, 10)
    pH_vals = []
    pco2_vals = []
    
    for mult in alk_multipliers:
        b_test = b_seawater.copy()
        b_test[idx_Alk] *= mult
        
        # Get pH via get_b_eq
        _, pH = get_b_eq(1e5, 288, P_CO2_modern, basalt_composition, b_input=b_test)
        pco2 = get_P_CO2(1e5, 288, b_test)
        
        pH_vals.append(pH)
        pco2_vals.append(pco2 / EARTH_ATM * 1e6)

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(alk_multipliers, pH_vals, 'b-', label='pH')
    ax1.set_xlabel('Alkalinity Multiplier (Relative to Modern)')
    ax1.set_ylabel('pH', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    
    ax2 = ax1.twinx()
    ax2.plot(alk_multipliers, pco2_vals, 'r-', label='P_CO2')
    ax2.set_ylabel('P_CO2 (ppm)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    plt.title("Test 2: Carbonate Speciation (Constant Total C)")
    plt.savefig('test2.png')
    plt.close()
    print("   -> Plot generated. pH should rise and P_CO2 should fall as Alk increases.\n")

    # ---------------------------------------------------------
    # Test 3: Transport vs Kinetic Limits
    # ---------------------------------------------------------
    print("3. Testing Transport vs. Kinetic Limits...")
    huge_J = J_ref_normalised * 1e6
    huge_alpha = alpha_ref * 1e6
    
    flux_norm, _ = get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_seawater, alpha_ref, rate_ref)
    flux_fast_flow, _ = get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_seawater, alpha_ref, rate_ref, J=huge_J)
    flux_fast_rxn, _ = get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_seawater, huge_alpha, rate_ref)
    
    print(f"   Normal Flux (Ca):    {flux_norm[idx_Ca]:.2e}")
    print(f"   Fast Flow (Kinetic): {flux_fast_flow[idx_Ca]:.2e} (Should be > Normal)")
    print(f"   Fast Rxn (Transport):{flux_fast_rxn[idx_Ca]:.2e} (Should be > Normal)\n")

    # ---------------------------------------------------------
    # Test 4: Thermostat Dependency (Temperature Sweep)
    # ---------------------------------------------------------
    print("4. Testing Thermostat Dependency...")
    T_sweep = np.linspace(273, 350, 15)
    fluxes_Ca, fluxes_Mg, fluxes_Alk = [], [], []
    
    for T in T_sweep:
        f, _ = get_weathering_flux(P_seafloor, T, P_CO2_modern, b_seawater, alpha_ref, rate_ref)
        fluxes_Ca.append(f[idx_Ca])
        fluxes_Mg.append(f[idx_Mg])
        fluxes_Alk.append(f[idx_Alk])
        
    plt.figure(figsize=(8, 5))
    plt.plot(T_sweep, fluxes_Ca, label='Ca Flux')
    plt.plot(T_sweep, fluxes_Mg, label='Mg Flux')
    plt.plot(T_sweep, fluxes_Alk, label='Alkalinity Flux')
    plt.xlabel('Seafloor Temperature (K)')
    plt.ylabel('Weathering Flux (mol/kgw/s)')
    plt.yscale('symlog', linthresh=1e-15)
    plt.title('Test 4: Thermostat Exponential Dependency')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    print("   -> Plot generated. Fluxes should scale strongly with temperature.\n")

    # ---------------------------------------------------------
    # Test 5: Charge Balance / Alkalinity Conservation
    # ---------------------------------------------------------
    print("5. Testing Charge Balance...")
    b_pure = np.zeros(elements.shape)
    flux_pure, _ = get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_pure, alpha_ref, rate_ref)
    
    # Approx check: Alk ≈ 2*Ca + 2*Mg + Na (ignoring minor Fe/Al species for a rough check)
    cations_charge = 2 * flux_pure[idx_Ca] + 2 * flux_pure[idx_Mg] + 1 * flux_pure[idx_Na]
    alk_flux = flux_pure[idx_Alk]
    cl_flux = flux_pure[idx_Cl]

    
    print(f"   Flux Alkalinity:    {alk_flux:.4e}")
    print(f"   Charge from Anions: {cl_flux:.4e}")
    print(f"   Charge from Cations:{cations_charge:.4e}")
    print(f"   Ratio (Anions+Alk/Cations): {(alk_flux + cl_flux) / cations_charge:.2f} (Should be ~1.0)\n")

    # ---------------------------------------------------------
    # Test 6: Dead Ocean Edge Case
    # ---------------------------------------------------------
    print("6. Testing 'Dead Ocean' (Zeros) Edge Case...")
    try:
        get_P_CO2(1e5, 288, b_pure)
        get_weathering_flux(P_seafloor, T_seafloor, P_CO2_modern, b_pure, alpha_ref, rate_ref)
        print("   Success: Model handled zero-concentration arrays without crashing.\n")
    except Exception as e:
        print(f"   FAILED: Dead ocean crashed with error: {e}\n")

    print("--- Benchmarks Complete ---")

if __name__ == '__main__':
    run_all_benchmarks()