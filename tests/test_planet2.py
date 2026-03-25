import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

os.environ['JAX_PLATFORMS'] = 'cpu'

import multiprocessing
import concurrent.futures
from tqdm import tqdm

import glob
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable

from kamino.planet2 import Planet
from kamino.constants import *

import itertools

# instellation_range = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25]
instellation_range = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
tectonics_range = [0.01, 0.1, 1, 10, 100]
ocean_depth_range = [100, 500, 1000, 3000, 5000, 10000]


# 2. Fixed Earth-like constants for the un-swept parameters
BACKGROUND_PRESSURE = 1e5   # Pa (~1 bar)
OCEAN_DEPTH = 3000          # m
TECTONICS = 1.0

def run_single_planet(params, rerun=True):
    """
    Worker function to simulate a single planet.
    Must be defined at the top level for multiprocessing to work.
    """
    inst, tect, depth, base_dir = params
    run_name = f"planet_inst_{inst:.2f}_tect_{tect:.2f}_depth_{depth:.1e}"
    run_dir = os.path.join(base_dir, run_name)

    # Check if the run already exists
    if os.path.exists(run_dir) and not rerun:
        print(f'Run directory already exists: {run_dir}')
        return

    print(f'Starting {run_name}...')
    
    # Only create the directory if we are actually going to run it
    os.makedirs(run_dir, exist_ok=True)

    p = Planet(
        mass=M_EARTH,
        radius=R_EARTH,
        background_pressure=BACKGROUND_PRESSURE,
        instellation=inst,
        tectonics=tect,
        ocean_depth=depth,
        name=run_name
    )

    try:
        # Suppress standard print statements inside the worker to avoid 
        # mangling the progress bar in the terminal
        p.time_evolve_to_steady_state(output_dir=run_dir)
        return f"Success: {run_name}"
    except Exception as e:
        return f"Failed: {run_name} | Error: {e}"

def run_parameter_sweep_tectonics():
    base_output_dir = "sweep_results"
    os.makedirs(base_output_dir, exist_ok=True)

    # Build the list of task arguments
    combinations = list(itertools.product(instellation_range, tectonics_range))
    tasks = [(inst, tect, OCEAN_DEPTH, base_output_dir) for inst, tect in combinations]
    total_runs = len(tasks)

    print(f"Starting parallel parameter sweep: {total_runs} total configurations.")
    
    # Use ProcessPoolExecutor to manage the CPU cores
    # By default, it uses as many workers as you have processors on your machine
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        # Submit all tasks to the executor
        futures = {executor.submit(run_single_planet, task): task for task in tasks}
        
        # Use tqdm to track completed futures as they finish
        with tqdm(total=total_runs, unit="run", desc="Sweeping") as pbar:
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                # You can uncomment the line below if you want to see individual successes/fails
                # tqdm.write(result) 
                pbar.update(1)

    print("\nParameter sweep complete! Check the 'sweep_results' folder.")


def run_parameter_sweep_ocean_depth():
    # Use a new output directory so we don't mix with the tectonics sweep
    base_output_dir = "sweep_results_depth"
    os.makedirs(base_output_dir, exist_ok=True)

    combinations = list(itertools.product(instellation_range, ocean_depth_range))
    tasks = [(inst, TECTONICS, depth, base_output_dir) for inst, depth in combinations]
    total_runs = len(tasks)

    print(f"Starting parallel parameter sweep: {total_runs} total configurations.")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(run_single_planet, task): task for task in tasks}
        
        with tqdm(total=total_runs, unit="run", desc="Sweeping") as pbar:
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                pbar.update(1)

    print(f"\nParameter sweep complete! Check the '{base_output_dir}' folder.")


def plot_steady_states(
    results_dir="sweep_results", 
    sweep_param_key="tectonics", 
    cbar_label="Tectonics Factor",
    cmap_name="inferno",
    output_filename="steady_states_sweep.pdf"
):
    """
    Generic plotting function to handle different parameter sweeps.
    """
    # 1. Locate all summary JSON files from the specified directory
    search_pattern = os.path.join(results_dir, "**", "*_summary.json")
    json_files = glob.glob(search_pattern, recursive=True)
    
    if not json_files:
        print(f"No summary files found in '{results_dir}'. Make sure you've run the sweep first.")
        return

    data_by_param = {}
    print(f"Reading {len(json_files)} result files from {results_dir}...")
    
    for file in json_files:
        with open(file, 'r') as f:
            summary = json.load(f)
            
            inst = summary['input_parameters']['instellation_solar']
            
            # Extract the specific swept parameter dynamically
            param_value = summary['input_parameters'][sweep_param_key]
            
            T_final = summary['final_state']['T_K']
            P_CO2_Pa = summary['final_state']['P_CO2_Pa']
            
            P_CO2_bar = P_CO2_Pa / 1e5
            
            if param_value not in data_by_param:
                data_by_param[param_value] = {'inst': [], 'T': [], 'P_CO2': []}
                
            data_by_param[param_value]['inst'].append(inst)
            data_by_param[param_value]['T'].append(T_final)
            data_by_param[param_value]['P_CO2'].append(P_CO2_bar)

    # 2. Setup the plot 
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    plt.subplots_adjust(hspace=0.05) 
    
    cmap = plt.get_cmap(cmap_name)
    param_values = list(data_by_param.keys())
    
    # Use LogNorm to handle orders of magnitude for both tectonics and depth
    norm = colors.LogNorm(vmin=min(param_values), vmax=max(param_values))
    
    # 3. Plot the lines
    for val in sorted(param_values):
        sort_idx = np.argsort(data_by_param[val]['inst'])
        inst_sorted = np.array(data_by_param[val]['inst'])[sort_idx]
        T_sorted = np.array(data_by_param[val]['T'])[sort_idx]
        P_CO2_sorted = np.array(data_by_param[val]['P_CO2'])[sort_idx]
        
        color = cmap(norm(val))
        
        # zorder=2 ensures lines are drawn on top of the shaded background regions
        ax1.plot(inst_sorted, T_sorted, color=color, linewidth=1.5, zorder=2)
        ax2.plot(inst_sorted, P_CO2_sorted, color=color, linewidth=1.5)

    # 4. Formatting Top Panel (Temperature) with Shaded Regions
    ax1.set_ylabel('Temperature (K)', fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.grid(True, linestyle='--', alpha=0.5, zorder=1)
    
    # Add Snowball Region (Blue, below freezing ~273.15K) 
    ax1.axhspan(240, 273.15, color='#ccccff', alpha=0.5, zorder=0)
    ax1.text(0.25, 260, 'Snowball', fontsize=12, va='center', ha='left', zorder=3)
    
    # Add Runaway Greenhouse Region (Red, >= 340K) 
    ax1.axhspan(340, 360, color='#ffcccc', alpha=0.5, zorder=0)
    ax1.text(0.25, 345, 'Runaway Greenhouse', fontsize=12, va='center', ha='left', zorder=3)
    
    # Lock the y-axis to match the visual frame of the reference 
    ax1.set_ylim(245, 350)
    ax1.set_xlim(0.2, 1.3)
    
    # 5. Formatting Bottom Panel (P_CO2)
    ax2.set_ylabel('$P_{CO2}$ (bar)', fontsize=12)
    ax2.set_xlabel('Instellation ($S/S_0$)', fontsize=12)
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # Lock the P_CO2 limits to match the reference 
    ax2.set_ylim(1e-7, 1e1)
    
    # 6. Add the Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.03, 0.7]) 
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([]) 
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(cbar_label, fontsize=12)
    
    # Output file
    plt.savefig(output_filename, bbox_inches='tight')
    print(f"\nPlot saved successfully to: {output_filename}")
    plt.show()

if __name__ == "__main__":

    multiprocessing.set_start_method('spawn')
    run_parameter_sweep_tectonics()
    # run_parameter_sweep_ocean_depth()

    plot_steady_states()

    # plot_steady_states(
    #     results_dir="sweep_results_depth", 
    #     sweep_param_key="ocean_depth_m", 
    #     cbar_label="Ocean Depth (m)",
    #     cmap_name="viridis",
    #     output_filename="steady_states_depth.pdf"
    # )
