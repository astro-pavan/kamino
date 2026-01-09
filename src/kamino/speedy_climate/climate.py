import pandas as pd
import numpy as np

import subprocess
import os
import shutil

import sys
import tarfile
import importlib.resources
from pathlib import Path

from kamino.constants import *

def get_climate_resources():
    """
    Locates the HELIOS submodule and the Data directory.
    Returns: (helios_path, data_path)
    """
    # Locate 'kamino.speedy_climate' package
    try:
        pkg_path = importlib.resources.files("kamino.speedy_climate")
    except (ImportError, AttributeError):
        pkg_path = Path(__file__).parent

    helios_path = pkg_path / "HELIOS"
    data_path = pkg_path / "data"
    
    return helios_path, data_path


def ensure_opacity_data(data_path: Path):
    """
    Checks if opacity data is unpacked. If not, unpacks opacity_data.tar.xz.
    """
    opacity_dir = data_path / "opacity"
    tarball = opacity_dir / "opacity_data.tar.xz"

    # Check if directory exists and is not empty
    if opacity_dir.exists():
        files = list(opacity_dir.iterdir())
        # Check if there are files other than .tar.xz archives
        if any(f.suffix != ".xz" or not f.name.endswith(".tar.xz") for f in files):
            return

    print(f"First-time setup: Unpacking {tarball.name}...")
    
    if not tarball.exists():
        raise FileNotFoundError(f"Opacity tarball not found at {tarball}")

    # Unpack
    with tarfile.open(tarball, "r:xz") as tar:
        tar.extractall(path=opacity_dir)
        
    print("Unpacking complete.")


def run_HELIOS(name: str, instellation: float, spectral_type: str, R_planet: float, M_planet: float, P_surface: float, x_CO2: float, x_H2O: float, albedo: float, recirculation_factor: float, clouds: bool=False, verbose: bool=False) -> dict[str, object]:
    """_summary_

    Parameters
    ----------
    name : str
        Name of HELIOS run.
    instellation : float
        Instellation flux in W/m^2.
    spectral_type : str
        Spectral type of host star (can be any from SPECTRAL_TYPE_DATA dict).
    R_planet : float
        Planet radius in m.
    M_planet : float
        Planet mass in kg.
    P_surface : float
        Surface pressure in Pa.
    x_CO2 : float
        CO2 volume mixing ratio.
    x_H2O : float
        H2O volume mixing ration.
    albedo : float
        Albedo.
    recirculation_factor : float
        Recirculation factor (0.25 if rapidly rotating, 0.666 if tidally locked).
    clouds : bool
        Whether to include water clouds, by default False.
    verbose : bool, optional
        Whether to print HELIOS output to terminal, by default False.

    Returns
    -------
    dict[str, object]
        Dictionary of climate run inputs and output surface temperature.  
    """


    helios_path, data_path = get_climate_resources()
    ensure_opacity_data(data_path)

    g_surface = (G * M_planet) / (R_planet ** 2)

    R_star = SPECTRAL_TYPE_DATA[spectral_type]['Radius']
    T_star = SPECTRAL_TYPE_DATA[spectral_type]['Temperature']

    orbital_distance = np.sqrt((((R_star * R_SUN) ** 2) * STEFAN_BOLTZMANN * (T_star ** 4)) / (instellation)) / AU

    result_dict: dict[str, object] = {}

    try:

        input_dir = helios_path / "input"
        input_dir.mkdir(exist_ok=True)
        species_file = input_dir / f'species_{name}.dat'

        species_data = {
            'species' : ['H2O', 'CO2'],
            'absorbing' : ['yes', 'yes'],
            'scattering' : ['yes', 'yes'],
            'mixing_ratio' : [x_H2O, x_CO2]
        }

        pd.DataFrame(species_data).to_csv(species_file, sep='\t', index=False)

        command = [sys.executable, 'helios.py']

        # print((data_path / "opacity").as_posix() + "/water.mie")

        parameters = [
            '-name', name,
            '-boa_pressure', f'{P_surface * 10}',
            '-f_factor', f'{recirculation_factor}',
            '-surface_albedo', f'{albedo}',
            '-surface_gravity', f'{g_surface * 100}',
            '-radius_planet', f'{R_planet / R_JUPITER}',
            '-path_to_species_file', species_file.as_posix(), 
            '-radius_star', f'{R_star}',
            '-temperature_star', f'{T_star}',
            '-orbital_distance', f'{orbital_distance}',
            '-directory_with_opacity_files', (data_path / "opacity").as_posix() + "/",
            '-opacity_mixing', 'on-the-fly',
            '-stellar_spectral_model', 'blackbody',
            '-realtime_plotting', 'yes',
            '-planet', 'manual',
            '-planet_type', 'rocky',
            '-number_of_layers', '25',
            '-k_coefficients_mixing_method', 'correlated-k',
            '-number_of_cloud_decks', '1' if clouds else '0',
            '-path_to_mie_files', (data_path / "opacity").as_posix() + "/water.mie",
            '-aerosol_radius_mode', '11.0',
            # '-aerosol_radius_geometric_std_dev ', '2.0',
            '-cloud_mixing_ratio', 'file',
            '-path_to_file_with_cloud_data', (data_path / "opacity").as_posix() + "/earth_clouds.dat",
            '-aerosol_name', 'water',
            '-radiative_equilibrium_criterion', '1e-4',
            # '-adaptive_interval', '1',
            #'-tp_profile_smoothing', 'yes',
            #'-cloud_bottom_pressure', '8.5e5',
            #'-cloud_bottom_mixing_ratio', '2.1e-5',
            #'-cloud_to_gas_scale_height_ratio', '0.2'
        ]

        env = os.environ.copy()

        subprocess.run(
                command + parameters, 
                cwd=helios_path, 
                env=env,
                stdout=None if verbose else subprocess.DEVNULL, 
                stderr=None if verbose else subprocess.DEVNULL,
                check=True
        ) # type: ignore

        try:
            atm_df = pd.read_table(f'{helios_path}/output/{name}/{name}_tp.dat', sep='\s+', skiprows=1)

            # P = np.array(atm_df['press.[10^-6bar]']) * 10
            T = np.array(atm_df['temp.[K]'])
            print(T)
            z = np.array(atm_df['altitude[cm]']) / 100
            print(z)

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

    finally:

        # This block runs whether the code succeeds OR fails
        # removing the species file
        # if 'species_file' in locals() and species_file.exists():
        #     os.remove(species_file)
            
        # removing the output directory
        output_dir = helios_path / "output" / name
        # if output_dir.exists():
        #     shutil.rmtree(output_dir, ignore_errors=True)
        print(output_dir)

    return result_dict
