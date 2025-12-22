import pandas as pd
import numpy as np

import subprocess
import os
import shutil

from kamino.constants import *

WD = os.getcwd()
HELIOS_PATH = f'{WD}/src/kamino/speedy_climate/helios'

def run_HELIOS(name: str, instellation: float, spectral_type: str, R_planet: float, M_planet: float, P_surface: float, x_CO2: float, x_H2O: float, albedo: float, recirculation_factor: float, verbose: bool=False) -> dict[str, object]:

    g_surface = (G * M_planet) / (R_planet ** 2)

    R_star = SPECTRAL_TYPE_DATA[spectral_type]['Radius']
    T_star = SPECTRAL_TYPE_DATA[spectral_type]['Temperature']

    orbital_distance = np.sqrt((((R_star * R_SUN) ** 2) * STEFAN_BOLTZMANN * (T_star ** 4)) / (instellation)) / AU

    species_data = {
        'species' : ['H2O', 'CO2'],
        'absorbing' : ['yes', 'yes'],
        'scattering' : ['yes', 'yes'],
        'mixing_ratio' : [x_H2O, x_CO2]
    }

    species_df = pd.DataFrame(species_data)
    species_df.to_csv(f'{HELIOS_PATH}/input/species_{name}.dat', sep='\t', index=False)

    command = [f'{WD}/.venv/bin/python', 'helios.py']

    parameters = [
        '-name', name,
        '-boa_pressure', f'{P_surface / 10}',
        '-f_factor', f'{recirculation_factor}',
        '-surface_albedo', f'{albedo}',
        '-surface_gravity', f'{g_surface * 100}',
        '-radius_planet', f'{R_planet / R_JUPITER}',
        '-path_to_species_file', f'./input/species_{name}.dat', 
        '-radius_star', f'{R_star}',
        '-temperature_star', f'{T_star}',
        '-orbital_distance', f'{orbital_distance}',
        '-directory_with_opacity_files', f'{WD}/src/kamino/speedy_climate/data/opacity/',
        '-opacity_mixing', 'on-the-fly',
        '-stellar_spectral_model', 'blackbody',
        '-realtime_plotting', 'no',
        '-planet', 'manual',
        '-planet_type', 'rocky',
        '-number_of_layers', '25',
        '-k_coefficients_mixing_method', 'correlated-k'
    ]

    env = os.environ.copy()
    #env["PATH"] = '/data/pt426/cuda/cuda12/bin' + ":" + env["PATH"]
    #env["LD_LIBRARY_PATH"] = '/data/pt426/cuda/cuda12/lib64' + ":" + env.get("LD_LIBRARY_PATH", "")
    #env["DYLD_LIBRARY_PATH"] = '/data/pt426/cuda/cuda12/lib' + ":" + env.get("LD_LIBRARY_PATH", "")

    if verbose:
        subprocess.run(command + parameters, cwd=HELIOS_PATH, env=env) 
    else:
        subprocess.run(command + parameters, cwd=HELIOS_PATH, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    try:
        atm_df = pd.read_table(f'{HELIOS_PATH}/output/{name}/{name}_tp.dat', sep='\s+', skiprows=1)

        # P = np.array(atm_df['press.[10^-6bar]']) * 10
        T = np.array(atm_df['temp.[K]'])
        # z = np.array(atm_df['altitude[cm]']) / 100

        T_surface = float(T[0])

        result_dict = {
            "Run_ID": name,
            "Instellation (W/m^2)": instellation,
            "Spectral_Type": spectral_type,
            "R_Planet (m)": R_planet,
            "M_Planet (kg)": M_planet,
            "P_Surface (Pa)": P_surface,
            "x_CO2": x_CO2,
            "x_H2O": x_H2O,
            "Albedo": albedo,
            "Recirculation_Factor": recirculation_factor,
            "Surface_Temp (K)": T_surface,
            "Status": "Success"
        }

    except Exception as e :

        result_dict = {
            "Run_ID": name,
            # "Zenith_Angle": zenith_angle,
            "Instellation (W/m^2)": instellation,
            "Spectral_Type": spectral_type,
            "R_Planet (m)": R_planet,
            "M_Planet (kg)": M_planet,
            "P_Surface (Pa)": P_surface,
            "x_CO2": x_CO2,
            "x_H2O": x_H2O,
            "Albedo": albedo,
            "Recirculation_Factor": recirculation_factor,
            "Surface_Temp (K)": -1,
            "Status": str(e)
        }

    shutil.rmtree(f'{HELIOS_PATH}/output/{name}', ignore_errors=True)
    os.remove(f'{HELIOS_PATH}/input/species_{name}.dat')

    return result_dict

if __name__ == '__main__':

    print(run_HELIOS('test', 1.5 * SOLAR_CONSTANT, 'G2', R_EARTH, M_EARTH, EARTH_ATM, 0.9 * EARTH_ATM, 0.9 * EARTH_ATM, 0.0, 0.25, verbose=True))
