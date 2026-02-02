import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

# --- IMPORTS ---
try:
    from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate
    from kamino.ocean_chemistry.fast_chemistry import get_calcite_data_fast
except ImportError:
    print("Error: Could not import solvers. Run this script from the root directory or adjust paths.")
    exit()

def test_parameter_space():
    print("--- Starting Robustness Test (Wide Parameter Sweep) ---")

    # 1. Define Parameter Ranges
    # --------------------------
    temps = np.linspace(275, 360, 10)       # 2°C to 87°C
    cas   = np.linspace(5e-3, 50e-3, 5)     # 5mM to 50mM (Low to High Salinity Ca)
    dics  = np.linspace(1e-3, 10e-3, 5)     # 1mM to 10mM (Carbon starved to Carbon rich)
    
    # Fixed parameters
    P = 1e5
    Alk_ratio = 1.1 # Alk typically slightly > DIC in ocean

    results = []
    
    total_steps = len(temps) * len(cas) * len(dics)
    step = 0
    t_start_all = time.time()

    print(f"Testing {total_steps} unique conditions...")

    for T in temps:
        for Ca in cas:
            for DIC in dics:
                Alk = DIC * Alk_ratio
                
                # A. Run PHREEQC (The Reference)
                try:
                    t0 = time.time()
                    rate_ref, SI_ref = get_calcite_precipitation_rate(P, T, Alk, DIC, Ca)
                    t_ref = time.time() - t0
                except Exception:
                    rate_ref, SI_ref = np.nan, np.nan
                    t_ref = 0

                # B. Run Fast Solver (The Candidate)
                try:
                    t0 = time.time()
                    # scaling_factor=1.0 for direct physics comparison
                    rate_fast, SI_fast = get_calcite_data_fast(P, T, Alk, DIC, Ca, scaling_factor=1.0)
                    t_fast = time.time() - t0
                except Exception:
                    rate_fast, SI_fast = np.nan, np.nan
                    t_fast = 0

                # C. Calculate Errors
                # SI Error (Absolute difference)
                err_SI = SI_fast - SI_ref
                
                # Rate Error (Log10 difference, handle zeros)
                if rate_ref > 1e-20 and rate_fast > 1e-20:
                    err_rate = np.log10(rate_fast) - np.log10(rate_ref)
                elif rate_ref < 1e-20 and rate_fast < 1e-20:
                    err_rate = 0.0 # Both zero match
                else:
                    err_rate = np.nan # One is zero, one is not (dissolution/precip boundary mismatch)

                results.append({
                    "T_K": T,
                    "Ca_mM": Ca * 1000,
                    "DIC_mM": DIC * 1000,
                    "SI_Ref": SI_ref,
                    "SI_Fast": SI_fast,
                    "Err_SI": err_SI,
                    "Rate_Ref": rate_ref,
                    "Rate_Fast": rate_fast,
                    "Log_Err_Rate": err_rate,
                    "Time_Speedup": t_ref / t_fast if t_fast > 1e-9 else 0
                })
                
                step += 1
                if step % 20 == 0:
                    print(f"  Progress: {step}/{total_steps}...")

    df = pd.DataFrame(results)
    
    print(f"\nCompleted in {time.time() - t_start_all:.2f}s")
    print(f"Mean SI Error: {df['Err_SI'].mean():.4f}")
    print(f"Max SI Error:  {df['Err_SI'].abs().max():.4f}")
    print(f"Mean Speedup:  {df['Time_Speedup'].mean():.1f}x")

    # 2. Visualizations
    # -----------------
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot A: SI Error Heatmap (Temp vs Calcium)
    # We aggregate DIC by taking the mean error for that T/Ca cell
    pivot_si = df.pivot_table(index='Ca_mM', columns='T_K', values='Err_SI', aggfunc='mean')
    sns.heatmap(pivot_si, ax=axes[0,0], cmap="coolwarm", center=0, annot=True, fmt=".2f")
    axes[0,0].set_title("Mean SI Error (Fast - Ref)\nby Temp and Calcium")
    axes[0,0].invert_yaxis()

    # Plot B: Rate Error vs SI (Scatter)
    # Checks if rate deviates specifically at high or low saturation
    sc = axes[0,1].scatter(df['SI_Ref'], df['Log_Err_Rate'], c=df['T_K'], cmap='viridis', alpha=0.7)
    axes[0,1].set_xlabel("Reference Saturation Index")
    axes[0,1].set_ylabel("Log10 Rate Error (Fast/Ref)")
    axes[0,1].set_title("Rate Discrepancy vs Saturation")
    axes[0,1].axhline(0, color='k', linestyle='--', linewidth=1)
    plt.colorbar(sc, ax=axes[0,1], label='Temperature (K)')

    # Plot C: Rate Comparison (1:1 Line)
    axes[1,0].loglog(df['Rate_Ref'], df['Rate_Fast'], 'k.', alpha=0.3)
    min_rate = min(df['Rate_Ref'].min(), df['Rate_Fast'].min())
    max_rate = max(df['Rate_Ref'].max(), df['Rate_Fast'].max())
    axes[1,0].plot([min_rate, max_rate], [min_rate, max_rate], 'r--', label='Perfect Match')
    axes[1,0].set_xlabel("PHREEQC Rate")
    axes[1,0].set_ylabel("Fast Rate")
    axes[1,0].set_title("Precipitation Rate 1:1 Comparison")
    axes[1,0].legend()

    # Plot D: Speedup Histogram
    axes[1,1].hist(df['Time_Speedup'], bins=20, color='green', alpha=0.7)
    axes[1,1].set_xlabel("Speedup Factor (X times faster)")
    axes[1,1].set_title("Performance Gain Distribution")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    test_parameter_space()