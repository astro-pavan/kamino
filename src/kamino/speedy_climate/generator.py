import pandas as pd
import numpy as np
from scipy.stats import qmc

import time

import concurrent.futures

from pathlib import Path
import os
from typing import Union

from kamino.speedy_climate.climate import run_HELIOS
from kamino.constants import *

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'data' / 'climate_runs'

def generate_input_parameters(
        n_samples: int, 
        instellation_bounds: Union[tuple[float, float], float],
        log_P_background_bounds: Union[tuple[float, float], float],
        log_P_CO2_bounds: Union[tuple[float, float], float],
        log_P_H2O_bounds: Union[tuple[float, float], float],
        albedo_bounds: Union[tuple[float, float], float],
        cloud_fraction_bounds: Union[tuple[float, float], float],
        spectral_type: str, 
        recirculation_factor: float,
        M_planet: float=M_EARTH,
        R_planet: float=R_EARTH
        ):
    
    param_bounds = {}
    param_const = {}

    args = {
        'instellation': instellation_bounds,
        'log_p_background': log_P_background_bounds,
        'log_p_co2': log_P_CO2_bounds,
        'log_p_h2o': log_P_H2O_bounds,
        'albedo': albedo_bounds,
        'cloud_fraction': cloud_fraction_bounds
    }

    for key, val in args.items():
        if isinstance(val, tuple):
            param_bounds[key] = val
        else:
            param_const[key] = val

    n_dimensions = len(param_bounds)
    l_bounds = [b[0] for b in param_bounds.values()]
    u_bounds = [b[1] for b in param_bounds.values()]

    sampler = qmc.LatinHypercube(d=n_dimensions)

    samples_unit_cube = sampler.random(n=n_samples)
    samples_scaled = qmc.scale(samples_unit_cube, l_bounds, u_bounds)

    inputs = []

    for j, sample in enumerate(samples_scaled):

        # Format: (name, instellation, spec_type, Rp, Mp, P_background, P_CO2, P_H2O, alb, recirc, cloud_fraction)

        for i, key in enumerate(param_bounds.keys()):
            if key == 'instellation':
                instellation = float(sample[i])
            if key == 'log_p_background':
                log_p_background = float(sample[i])
            if key == 'log_p_co2':
                log_p_co2 = float(sample[i])
            if key == 'log_p_h2o':
                log_p_h2o = float(sample[i])
            if key == 'albedo':
                albedo = float(sample[i])
            if key == 'cloud_fraction':
                cloud_fraction = float(sample[i])

        for i, key in enumerate(param_const.keys()):
            if key == 'instellation':
                instellation = param_const[key]
            if key == 'log_p_background':
                log_p_background = param_const[key]
            if key == 'log_p_co2':
                log_p_co2 = param_const[key]
            if key == 'log_p_h2o':
                log_p_h2o = param_const[key]
            if key == 'albedo':
                albedo = param_const[key]
            if key == 'cloud_fraction':
                cloud_fraction = param_const[key]

        # instellation = float(sample[0]) if 'instellation' in param_bounds else param_const['instellation']
        # log_p_background = float(sample[1]) if 'log_p_background' in param_bounds else param_const['log_p_background']
        # log_p_co2 = float(sample[2]) if 'log_p_co2' in param_bounds else param_const['log_p_co2']
        # log_p_h2o = float(sample[3]) if 'log_p_h2o' in param_bounds else param_const['log_p_h2o']
        # albedo = float(sample[4]) if 'albedo' in param_bounds else param_const['albedo']
        # cloud_fraction = float(sample[5]) if 'cloud_fraction' in param_bounds else param_const['cloud_fraction']

        P_back = 10 ** log_p_background # type: ignore
        P_co2 = 10 ** log_p_co2 # type: ignore
        P_h2o = 10 ** log_p_h2o # type: ignore

        sample_tuple = (
            f'run_{j}',
            instellation, # type: ignore
            spectral_type,
            R_planet,
            M_planet,
            P_back,
            P_co2,
            P_h2o,
            albedo, # type: ignore
            recirculation_factor,
            cloud_fraction, # type: ignore
            True if 'log_p_h2o' in param_const else False,
            True if 'cloud_fraction' in param_const else False,
        )

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

            count += 1

            result = future.result()
            results_list.append(result)

            if result["Status"] == "Success":
                print(f"Finished {result['Run_ID']}: {result['Surface_Temp (K)']:.0f} K [{count}/{total_models}]")
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

