import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate as get_phreeqc
from kamino.ocean_chemistry.precipitation_interpolator import get_calcite_data_interpolated as get_interpolator
from kamino.ocean_chemistry.precipitation_interpolator import PrecipitationInterpolator

def run_accuracy_test(n_samples=500):
    print(f"--- Starting Accuracy Test (N={n_samples} random points) ---")
    print("Comparing: Interpolator vs. PHREEQC (Ground Truth)")

    # 1. Generate Random Test Set
    # ---------------------------
    np.random.seed(42) # Reproducible results
    
    # Randomly sample within the Interpolator's valid bounds
    # Note: Using Log-uniform sampling for concentrations to cover orders of magnitude
    Ps   = np.geomspace(1e5, 500e5, n_samples)      # 1 to 500 bar
    Ts   = np.random.uniform(275, 370, n_samples)   # 275 to 370 K
    DICs = np.geomspace(1e-4, 0.1, n_samples)       # 0.1 mM to 100 mM
    Cas  = np.geomspace(1e-4, 0.1, n_samples)       # 0.1 mM to 100 mM
    
    # Alkalinity Ratio (0.5 to 2.5 typical for carbonates)
    Ratios = np.random.uniform(0.5, 2.5, n_samples) 
    Alks   = DICs * Ratios

    results = []
    start_time = time.time()

    # 2. Run Comparison Loop
    # ----------------------
    for i in range(n_samples):
        P, T, Alk, DIC, Ca = Ps[i], Ts[i], Alks[i], DICs[i], Cas[i]

        # A. Run PHREEQC (The Reference)
        # Note: This is slow/expensive
        try:
            rate_ref, SI_ref = get_phreeqc(P, T, Alk, DIC, Ca)
        except Exception:
            rate_ref, SI_ref = np.nan, np.nan

        # B. Run Interpolator (The Approximation)
        # Fast look-up
        try:
            # Note: The interpolator returns 'Kinetic Potential' (k=1). 
            # If your PHREEQC code uses a specific rate constant, the magnitudes might differ,
            # but the Saturation Index (SI) must match.
            rate_interp, SI_interp = get_interpolator(P, T, Alk, DIC, Ca)
        except Exception as e:
            print(f"Interpolator Error at {T:.1f}K, {P:.1e}Pa: {e}")
            rate_interp, SI_interp = np.nan, np.nan

        # C. (Optional) Run Fast Physics directly to check interpolation error specifically
        # rate_phys, SI_phys = get_fast_physics(P, T, Alk, DIC, Ca, rate_constant=1.0)

        results.append({
            "P_bar": P/1e5,
            "T_K": T,
            "SI_Ref": SI_ref,
            "SI_Interp": SI_interp,
            "Rate_Ref": rate_ref,
            "Rate_Interp": rate_interp,
            "Error_SI": SI_interp - SI_ref
        })

        if i % 50 == 0:
            print(f"  Progress: {i}/{n_samples}...")

    df = pd.DataFrame(results).dropna()
    
    print(f"\nCompleted in {time.time() - start_time:.2f}s")
    
    # 3. Analyze Results
    # ------------------
    # Calculate SI Error Statistics
    mae_si = df['Error_SI'].abs().mean()
    rmse_si = np.sqrt((df['Error_SI']**2).mean())
    max_si = df['Error_SI'].abs().max()

    print("\n--- Saturation Index (SI) Accuracy ---")
    print(f"  Mean Absolute Error: {mae_si:.4f}")
    print(f"  RMSE:                {rmse_si:.4f}")
    print(f"  Max Error:           {max_si:.4f}")
    
    if mae_si < 0.05:
        print("  >> RESULT: EXCELLENT MATCH")
    elif mae_si < 0.2:
        print("  >> RESULT: GOOD MATCH (Acceptable for global models)")
    else:
        print("  >> RESULT: DISCREPANCY DETECTED (Check thermodynamics)")

    # 4. Plotting
    # -----------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot A: SI Correlation
    axes[0].scatter(df['SI_Ref'], df['SI_Interp'], c=df['T_K'], cmap='coolwarm', alpha=0.6, edgecolors='none')
    axes[0].plot([df['SI_Ref'].min(), df['SI_Ref'].max()], [df['SI_Ref'].min(), df['SI_Ref'].max()], 'k--', lw=1)
    axes[0].set_xlabel('PHREEQC SI')
    axes[0].set_ylabel('Interpolator SI')
    axes[0].set_title(f'Saturation Index (SI) Agreement\nMAE={mae_si:.3f}')
    axes[0].grid(True, alpha=0.3)

    # Plot B: SI Error Distribution
    axes[1].hist(df['Error_SI'], bins=30, color='teal', edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('Error (Interpolator - PHREEQC)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('SI Error Distribution')
    axes[1].axvline(0, color='k', lw=1)
    
    # Plot C: Rate Comparison (Log-Log)
    # Filter out zeros/negatives for log plot
    valid_rates = df[(df['Rate_Ref'] > 1e-20) & (df['Rate_Interp'] > 1e-20)]
    
    if not valid_rates.empty:
        sc = axes[2].scatter(valid_rates['Rate_Ref'], valid_rates['Rate_Interp'], c=valid_rates['P_bar'], cmap='viridis', alpha=0.6)
        axes[2].set_xscale('log')
        axes[2].set_yscale('log')
        axes[2].set_xlabel('PHREEQC Rate (mol/kgw/s)')
        axes[2].set_ylabel('Interpolator Rate (Unscaled)')
        axes[2].set_title('Rate Comparison (Shape Check)')
        plt.colorbar(sc, ax=axes[2], label='Pressure (bar)')
        
        # Add note about rate scaling
        axes[2].text(0.05, 0.95, "Note: Offset expected if k=1 in table\nbut k!=1 in PHREEQC.", 
                     transform=axes[2].transAxes, fontsize=9, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        axes[2].text(0.5, 0.5, "No precipitation (SI < 0)\nin random sample", ha='center')

    # Filter for valid non-zero rates to avoid log(0) errors
    valid_mask = (df['Rate_Ref'] > 1e-25) & (df['Rate_Interp'] > 1e-25)
    df_valid = df[valid_mask].copy()

    # Calculate Log-Space Error
    df_valid['Log_Ref'] = np.log10(df_valid['Rate_Ref'])
    df_valid['Log_Interp'] = np.log10(df_valid['Rate_Interp'])
    df_valid['Log_Error'] = df_valid['Log_Interp'] - df_valid['Log_Ref']

    # Statistics
    mae_log_rate = df_valid['Log_Error'].abs().mean()
    max_log_rate = df_valid['Log_Error'].abs().max()

    print("\n--- Rate Accuracy (Log-10 Space) ---")
    print(f"  Mean Absolute Log Error: {mae_log_rate:.4f} log units")
    print(f"  Max Log Error:           {max_log_rate:.4f} log units")

    # Evaluation
    if mae_log_rate < 0.02:
        print("  >> RESULT: PERFECT (Errors < 5%)")
    elif mae_log_rate < 0.05:
        print("  >> RESULT: EXCELLENT (Errors < 12%)")
    elif mae_log_rate < 0.1:
        print("  >> RESULT: GOOD (Errors < 25%)")
    else:
        print("  >> RESULT: POOR (Systematic bias detected)")

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    _interpolator = PrecipitationInterpolator()
    _interpolator.generate_table()
    
    run_accuracy_test()