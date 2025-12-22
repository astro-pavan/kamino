import pandas as pd
import numpy as np
from scipy.stats import qmc

import time

import concurrent.futures

from pathlib import Path
import os

from climate import run_HELIOS
from kamino.constants import *

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'data' / 'climate_runs'

def generate_input_parameters(n_samples, spectral_type, recirculation_factor):

    param_bounds = {
        'instellation': (0.1 * SOLAR_CONSTANT, 2.0 * SOLAR_CONSTANT),
        'log_pressure': (4, 6),
        'log_x_h2o': (-6, 0),
        'log_x_co2': (-6, 0),
        'albedo': (0, 0.5),
    }

    n_dimensions = len(param_bounds)
    l_bounds = [b[0] for b in param_bounds.values()]
    u_bounds = [b[1] for b in param_bounds.values()]

    sampler = qmc.LatinHypercube(d=n_dimensions)

    samples_unit_cube = sampler.random(n=n_samples)
    samples_scaled = qmc.scale(samples_unit_cube, l_bounds, u_bounds)

    inputs = []

    for i, sample in enumerate(samples_scaled):

        # Format: (name, instellation, spec_type, Rp, Mp, P_surf, x_CO2, x_H2O, alb, recirc)

        p = 10 ** float(sample[1])
        x_h2o = 10 ** float(sample[2])
        x_co2 = 10 ** float(sample[3])
        alb = float(sample[4])

        sample_tuple = (f'run_{i}',
                        float(sample[0]),
                        spectral_type,
                        1.0 * R_EARTH,
                        1.0 * M_EARTH,
                        p,
                        x_co2, 
                        x_h2o,
                        alb,
                        recirculation_factor)
        
        inputs.append(sample_tuple)

    return inputs


def run_batch_simulation(inputs, output_csv_name="helios_results.csv"):
    
    start_time = time.time()
    results_list = []

    print('Beginning HELIOS model sweep:')

    total_models = len(inputs)
    count = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
        
        futures = []

        for args in inputs:
            futures.append(executor.submit(run_HELIOS, *args))

        for future in concurrent.futures.as_completed(futures):
            
            result = future.result()
            results_list.append(result)

            count += 1

            if result["Status"] == "Success":
                print(f"Finished {result['Run_ID']}: {result['Surface_Temp (K)']} K [{count}/{total_models}]")
            else:
                print(f"Failed {result['Run_ID']}: {result.get('Status')} [{count}/{total_models}]")

    elapsed_time = time.time() - start_time
    print(f"Total execution time: {elapsed_time:.2f} seconds")

    if results_list:
        df = pd.DataFrame(results_list)
        
        # Ensure the directory exists before writing
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        
        save_path = OUTPUT_DIR / output_csv_name
        df.to_csv(save_path, index=False)
        print(f"\nResults successfully saved to: {save_path}")
    else:
        print("No results were generated.")
        
if __name__ == '__main__':

    inputs = generate_input_parameters(1000, 'G2', 0.25)
    run_batch_simulation(inputs, 'helios_1000_runs_earth_rapid_rotator.csv')

    inputs = generate_input_parameters(1000, 'M5', 0.6666)
    run_batch_simulation(inputs, 'helios_1000_runs_earth_tidally_locked.csv')
