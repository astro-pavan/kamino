import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import concurrent.futures

import os
import importlib.resources
from pathlib import Path
import time

from kamino.speedy_climate.climate import run_HELIOS
from kamino.constants import *

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'data' / 'climate_runs'

def get_T_surface(S, x_co2, f_circ=0.25, spec_type='G2'):

    res = run_HELIOS('test', S * SOLAR_CONSTANT, spec_type, R_EARTH, M_EARTH, EARTH_ATM, x_co2, 0.01, 0.05, f_circ, clouds=0.6, verbose=False, august_roche_magnus=True)

    T_surface = res['Surface_Temp (K)']

    return T_surface

def generate_data(n: int=10):

    x_co2_range = np.logspace(-6, 0, num=n)
    S_range = np.linspace(0.1, 1.5, num=n)
    
    x_co2_arr = []
    S_arr = []
    T_arr = []

    count = 1

    for x_co2 in x_co2_range:
        for S in S_range:

            print(f'Run: {count}/{n*n}')
            x_co2_arr.append(x_co2)
            S_arr.append(S)

            T = -1
            with concurrent.futures.ThreadPoolExecutor(max_workers=6) as executor:
                future = executor.submit(get_T_surface, S, x_co2)
                try:
                    T = future.result(timeout=60)
                except concurrent.futures.TimeoutError:
                    print(f"Timeout: get_T_surface took longer than 1 minute for S={S}, x_co2={x_co2}")
                except Exception as e:
                    print(f"Error: {e} for S={S}, x_co2={x_co2}")
            
            T_arr.append(T)
            count += 1

    x_co2_arr = np.array(x_co2_arr)
    S_arr = np.array(S_arr)
    T_arr = np.array(T_arr)

    df = pd.DataFrame({'x_co2' : x_co2_arr, 'S' : S_arr, 'T' : T_arr})

    df.to_csv('dataset.csv')

def plot_dataset():

    # ---------------------------------------------------------
    # 1. Load Data
    # ---------------------------------------------------------

    # OPTION A: If loading from a real file, uncomment the line below:
    df = pd.read_csv('dataset.csv', index_col=0)
    df = df[df['T'] != -1]
    df = df[df['S'] > 0.2]

    # ---------------------------------------------------------
    # 2. Create the Plot
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define variables for clarity
    x = df['x_co2']
    y = df['S']
    z = df['T']

    # Create filled contours (tricontourf is best for column data)
    # levels=20 makes the gradient smoother
    contour_filled = ax.tricontourf(x, y, z, levels=200, cmap='plasma')
    contour_lines = ax.tricontour(x, y, z, [273, 288, 340], colors='k', linestyles='dotted')
    ax.clabel(contour_lines, fmt='%d K', colors='k', inline=False, use_clabeltext=True)
    ax.set_xscale('log')

    # ax.axvline(280e-6, label='Pre-Industrial CO2', color='k', linestyle='--')
    ax.scatter([280e-6], [1], s=10, c='black', label='Earth')
    ax.annotate('Earth', (280e-6, 1))

    # Optional: Add black contour lines on top for precision
    # contour_lines = ax.tricontour(x, y, z, levels=20, colors='k', linewidths=0.5, alpha=0.5)

    # ---------------------------------------------------------
    # 3. Formatting
    # ---------------------------------------------------------

    # Add a colorbar to show what T values the colors represent
    cbar = plt.colorbar(contour_filled, ax=ax)
    cbar.set_label('Temperature (K)')

    # Axis Labels
    ax.set_xlabel('$x_{CO2}$')
    ax.set_ylabel('S ($S_{\\oplus}$)')

    # Optional: Since x_co2 is very small, scientific notation is automatic,
    # but we can force a log scale if your x data spans many orders of magnitude.
    # ax.set_xscale('log') 

    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig('climate.png')

