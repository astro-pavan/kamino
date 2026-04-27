import os
import glob
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Set the path to where your JSON files are saved
OUTPUT_PATH = '/home/pavan/PhD/kamino_experiments/'

def load_data(output_path):
    """Loads all config and results JSONs into a single Pandas DataFrame."""
    print("Loading data...")
    all_data = []
    
    # Find all config files
    config_files = glob.glob(os.path.join(output_path, '*_config.json'))
    
    for config_file in config_files:
        # Construct the matching results filename
        base_name = config_file.replace('_config.json', '')
        results_file = f"{base_name}_results.json"
        
        if os.path.exists(results_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
            with open(results_file, 'r') as f:
                results = json.load(f)
                
            # Combine relevant data
            row = {
                'name': config['name'],
                'instellation': config['instellation'],
                'outgassing': config['outgassing'],
                'crust_production': config['crust_production_rate'],
                'rw': config.get('reverse_weathering', False),
                'carb_frac': config.get('crust_carbonate_content', 0.0),
                'T': results.get('T', np.nan),
                'P_CO2': results.get('P_CO2', np.nan),
                'pH': results.get('pH', np.nan),
                'status': results.get('status', -1)
            }
            all_data.append(row)
        else:
            print(f"Warning: Missing results file for {config_file}")
            
    df = pd.DataFrame(all_data)
    print(f"Loaded {len(df)} successful simulation runs.")
    return df

def plot_method_1_faceted_lines(df):
    """Method 1: Stacked line plots of T, P_CO2, and pH vs. Instellation."""
    print("Generating Faceted Line Plots...")
    
    # Select a few representative crust production rates to avoid 50+ plots
    target_crust_rates = [0.1, 1.0, 10.0]
    
    # Setup the colormap for outgassing (logarithmic scale)
    outgassing_vals = df['outgassing'].unique()
    norm = mcolors.LogNorm(vmin=min(outgassing_vals), vmax=max(outgassing_vals))
    cmap = plt.get_cmap('plasma')
    
    for crust_rate in target_crust_rates:
        # Filter data for this specific crust production rate
        subset = df[(df['crust_production'] == crust_rate) & 
                    (df['rw'] == False) & 
                    (df['carb_frac'] == 0.0)].copy()
        
        if subset.empty:
            continue
            
        fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
        fig.subplots_adjust(hspace=0.05, right=0.85)
        
        # Group by outgassing so we can plot one line per outgassing value
        grouped = subset.groupby('outgassing')
        
        for outgassing, group in grouped:
            group = group.sort_values('instellation') # Ensure lines are drawn left-to-right
            color = cmap(norm(outgassing))
            
            # Panel 1: Temperature
            axes[0].plot(group['instellation'], group['T'], color=color, linewidth=2)
            # Panel 2: P_CO2
            axes[1].plot(group['instellation'], group['P_CO2'], color=color, linewidth=2)
            # Panel 3: pH
            axes[2].plot(group['instellation'], group['pH'], color=color, linewidth=2)
            
        # --- Formatting Axes ---
        
        # Temp axis
        axes[0].set_ylabel('Temperature (K)', fontsize=12)
        axes[0].axhspan(250, 273, color='blue', alpha=0.1) # Snowball shading
        axes[0].text(0.22, 260, 'Snowball', fontsize=12)
        axes[0].axhspan(340, 350, color='red', alpha=0.1)  # Runaway shading
        axes[0].text(0.22, 342, 'Runaway Greenhouse', fontsize=12)
        axes[0].set_ylim(245, 350)
        
        # P_CO2 axis
        axes[1].set_ylabel('$P_{CO2}$ (bar)', fontsize=12)
        axes[1].set_yscale('log')
        axes[1].set_ylim(1e-7, 10)
        
        # pH axis
        axes[2].set_ylabel('Ocean pH', fontsize=12)
        axes[2].set_xlabel('Instellation ($S/S_0$)', fontsize=12)
        axes[2].set_ylim(4.5, 9.5)
        
        for ax in axes:
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_xlim(0.2, 1.4)
            
        # Add a single colorbar for the whole figure
        cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.7])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label('Outgassing Factor', fontsize=12)
        
        axes[0].set_title(f'Tectonics Sweep (Crust Production = {crust_rate}x Earth)', fontsize=14)
        
        # Save the figure
        plt.savefig(f'line_plot_crust_{crust_rate}.png', dpi=300, bbox_inches='tight')
        plt.close()

def plot_method_2_heatmaps(df):
    """Method 2: Heatmaps of Outgassing vs Crust Production at fixed instellations."""
    print("Generating Heatmaps...")
    
    # Select representative instellations
    target_s_values = [0.6, 1.0, 1.2]
    
    # We will create one figure per variable (T, P_CO2, pH), with subplots for instellation
    variables = [
        ('T', 'Temperature (K)', 'viridis', None),
        ('P_CO2', '$P_{CO2}$ (bar)', 'inferno', mcolors.LogNorm(vmin=1e-6, vmax=10)),
        ('pH', 'Ocean pH', 'cividis', None)
    ]
    
    for var, label, cmap_name, norm in variables:
        fig, axes = plt.subplots(1, len(target_s_values), figsize=(15, 5), sharey=True)
        
        for i, s in enumerate(target_s_values):
            ax = axes[i]
            
            # Filter and pivot the data to create a 2D grid
            subset = df[(df['instellation'] == s) & 
                        (df['rw'] == False) & 
                        (df['carb_frac'] == 0.0)].copy()
            
            if subset.empty:
                continue
                
            # Create a 2D matrix: rows=crust_production, cols=outgassing
            grid = subset.pivot(index='crust_production', columns='outgassing', values=var)
            
            # Extract x and y arrays for pcolormesh
            X = grid.columns.values
            Y = grid.index.values
            Z = grid.values
            
            # Plot the heatmap
            mesh = ax.pcolormesh(X, Y, Z, cmap=cmap_name, norm=norm, shading='nearest')
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_title(f'Instellation S = {s}')
            ax.set_xlabel('Outgassing Factor')
            
            if i == 0:
                ax.set_ylabel('Crust Production Rate Factor')
                
        # Add colorbar
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(mesh, cax=cbar_ax, label=label)
        
        fig.suptitle(f'Phase Space Map: {label}', fontsize=16)
        plt.subplots_adjust(right=0.9)
        
        # Save figure
        clean_var_name = var.replace('_', '')
        plt.savefig(f'heatmap_{clean_var_name}.png', dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    # 1. Load the data
    df = load_data(OUTPUT_PATH)
    
    # Make sure we actually found data before trying to plot
    if not df.empty:
        # 2. Plot Method 1 (Line Plots)
        plot_method_1_faceted_lines(df)
        
        # 3. Plot Method 2 (Heatmaps)
        plot_method_2_heatmaps(df)
        
        print("Plotting complete! Check your working directory for the PNG files.")
    else:
        print("No data found. Check your OUTPUT_PATH.")