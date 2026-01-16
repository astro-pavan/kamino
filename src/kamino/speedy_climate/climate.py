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
        if any(f.suffix == ".h5" for f in files):
            return

    print(f"First-time setup: Unpacking {tarball.name}...")
    
    if not tarball.exists():
        raise FileNotFoundError(f"Opacity tarball not found at {tarball}")

    # Unpack
    with tarfile.open(tarball, "r:xz") as tar:
        tar.extractall(path=opacity_dir)
        
    print("Unpacking complete.")


def run_HELIOS(
        name: str, 
        instellation: float, 
        spectral_type: str, 
        R_planet: float, 
        M_planet: float, 
        P_background: float, 
        P_CO2: float, 
        P_H2O: float, 
        albedo: float, 
        recirculation_factor: float,
        clouds: float=0,
        august_roche_magnus: bool=False,
        cloud_destruction: bool=True,
        rainout: bool=True,
        relative_humidity: float=0.77,
        moist_convection: bool=True,
        verbose: bool=False
        ) -> dict[str, object]:
    """
    Runs the HELIOS radiative convective climate code.

    Parameters
    ----------
    name : str
        Name of HELIOS run.
    instellation : float
        Instellation flux in W/m^2.
    spectral_type : str
        Spectral type of host star (can be any from SPECTRAL_TYPE_DATA dict in constants.py).
    R_planet : float
        Planet radius in m.
    M_planet : float
        Planet mass in kg.
    P_surface : float
        Pressure of background transparent gas (most likely N2) in Pa.
    P_CO2 : float
        CO2 partial pressure in Pa.
    P_H2O : float
        H2O partial pressure in Pa.
    albedo : float
        Albedo.
    recirculation_factor : float
        Recirculation factor (0.25 if rapidly rotating, 0.666 if tidally locked).
    clouds : float
        Fraction of planet covered in cloud, by default 0.
    august_roche_magnus : bool
        Whether to set surface x_H2O with the August-Roches-Magnus formula (overwrites provided x_H2O), by default False. 
    cloud_destruction : bool
        Whether to apply temperature dependant cloud destruction, by default True.
    rainout : bool
        Whether to apply H2O rainout, by default True.
    relative_humidity : float
        Value of relative humidity used in rainout claculations, by default 0.77.
    moist_convection : bool
        Whether to apply moist convection, by default True.
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

    P_surface = P_background + P_CO2 + P_H2O
    x_CO2 = P_CO2 / P_surface
    x_H2O = P_H2O / P_surface

    result_dict = {
        "Run_ID": name,
        "Instellation (W/m^2)": instellation,
        "Spectral_Type": spectral_type,
        "R_Planet (m)": R_planet,
        "M_Planet (kg)": M_planet,
        "P_Surface (Pa)": P_surface,
        "P_CO2 (Pa)": P_CO2,
        "P_H2O (Pa)": P_H2O,
        "Albedo": albedo,
        "Recirculation_Factor": recirculation_factor,
        "Rainout": rainout,
        "Relative Humidity": relative_humidity,
        "August Roche Magnus": august_roche_magnus,
        "Cloud Cover": clouds,
        "Cloud Destruction": cloud_destruction,
        "Moist Convection": moist_convection
    }

    try:
        
        input_dir = helios_path / "input"
        input_dir.mkdir(exist_ok=True)
        species_file = input_dir / f'species_{name}.dat'
        cloud_flie = input_dir / f'clouds_{name}.dat'

        species_data = {
            'species' : ['H2O', 'CO2'],
            'absorbing' : ['yes', 'yes'],
            'scattering' : ['yes', 'yes'],
            'mixing_ratio' : [x_H2O, x_CO2]
        }

        pd.DataFrame(species_data).to_csv(species_file, sep='\t', index=False)

        P_surf_bar = P_surface / 1e5
        
        # Cloud Limits (0.7 to 0.85 P_surf)
        cloud_base_P = 0.85 * P_surf_bar
        cloud_top_P  = 0.70 * P_surf_bar
        
        # Cloud Density (Mixing Ratio in g/g)
        max_mixing_ratio = 5.0e-5 
        
        pressure_grid = np.linspace(P_surf_bar, 1e-6, 20)
        cloud_profile = np.zeros_like(pressure_grid)
        
        # Fill in the cloud deck (Block / Slab profile)
        # select indices where Pressure is between Top and Base
        in_cloud = (pressure_grid <= cloud_base_P) & (pressure_grid >= cloud_top_P)
        cloud_profile[in_cloud] = max_mixing_ratio
        
        with open(cloud_flie.as_posix(), 'w') as f:
            # Header (Required for names=True in HELIOS)
            # Column 1: Pressure, Column 2: Cloud_Deck_1
            f.write("pressure\tlow_cloud\n")

            for p, mix in zip(pressure_grid, cloud_profile):
                # Format: Scientific notation, tab separated
                f.write(f"{p:.4f}\t{mix:.6e}\n")

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
            '-realtime_plotting', 'no',
            '-planet', 'manual',
            '-planet_type', 'rocky',
            '-number_of_layers', '25',
            '-k_coefficients_mixing_method', 'correlated-k',
            '-number_of_cloud_decks', '1' if clouds > 0 else '0',
            '-path_to_mie_files', (data_path / "opacity").as_posix() + "/water.mie",
            '-aerosol_radius_mode', '11.0',
            # '-aerosol_radius_geometric_std_dev ', '2.0', # this is set in the parameter file
            '-cloud_mixing_ratio', 'file',
            '-path_to_file_with_cloud_data', cloud_flie.as_posix(),
            '-aerosol_name', 'water',
            '-radiative_equilibrium_criterion', '1e-4',
            # parameters for added physics
            '-cloud_cover', f'{np.minimum(1.0, clouds)}',
            '-relative_humidity', f'{np.minimum(1.0, relative_humidity)}',
            '-rainout', 'yes' if rainout else 'no',
            '-moist_convection', 'yes' if moist_convection else 'no',
            '-convective_damping_parameter', '20',
            '-set_h2o_with_august_roche_magnus', 'yes' if august_roche_magnus else 'no',
            '-cloud_destruction', 'yes' if cloud_destruction else 'no'
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
            z = np.array(atm_df['altitude[cm]']) / 100

            if verbose:
                print(T)
                print(z)

            T_surface = float(T[0])

            result_dict['Surface_Temp (K)'] = T_surface
            result_dict['Status'] = "Success"

        except Exception as e :

            result_dict['Surface_Temp (K)'] = -1
            result_dict['Status'] = str(e)

    finally:

        # This block runs whether the code succeeds OR fails
        # removing the species file
        if 'species_file' in locals() and species_file.exists():
            os.remove(species_file)

        if 'cloud_file' in locals() and cloud_flie.exists():
            os.remove(cloud_flie)
            
        # removing the output directory
        output_dir = helios_path / "output" / name
        if output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)

    return result_dict
