import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel, Matern

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

from scipy.optimize import least_squares
from scipy.interpolate import CubicSpline

import joblib
import os
import importlib.resources
from pathlib import Path

from kamino.constants import *

def get_data_paths():
    """
    Returns the paths to the climate_runs and emulators directories.
    """
    # 1. Locate the package data directory
    try:
        # For Python 3.9+ standard package structure
        pkg_root = importlib.resources.files("kamino.speedy_climate")
    except (ImportError, AttributeError):
        # Fallback if not installed as a package
        pkg_root = Path(__file__).parent

    data_root = pkg_root / "data"
    
    # 2. Define subdirectories
    runs_dir = data_root / "climate_runs"
    emulators_dir = data_root / "emulators"
    
    # Ensure they exist (prevents FileNotFoundError when saving)
    emulators_dir.mkdir(parents=True, exist_ok=True)
    
    return runs_dir, emulators_dir

def august_roche_magnus_formula(T: float) -> float:
    """
    Implementation of August-Roche-Magnus formula to calculate H2O partial pressure.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        H2O partial pressure in Pa
    """

    T_celsius = T + ABSOLUTE_ZERO
    return 610.94 * np.exp((17.625 * T_celsius)/(T_celsius + 243.04))

class climate_emulator:

    def __init__(self, emulator_name: str, climate_data_file: str="", make_accuracy_plot: bool=False, force_retraining: bool=False):
        """
        Climate emulator class. 

        Parameters
        ----------
        emulator_name : str
            Name of emulator. If the emulator doesn't exist, it will train one from the climate data file provided.
        climate_data_file : str, optional
            Name of climate data file for the model to be trained on, by default "".
        make_accuracy_plot : bool, optional
            Whether to make an accuracy plot, by default False.
        force_retraining : bool, optional
            Whether to force retraining the emulator

        Raises
        ------
        FileNotFoundError
            Training data or emulator data file not found.
        """

        # Get the standard directories
        runs_dir, emulators_dir = get_data_paths()
        
        # Define paths for the emulator files
        gp_path = emulators_dir / f'{emulator_name}_gp_climate_emulator.pkl'
        x_scaler_path = emulators_dir / f'{emulator_name}_inputs_scaler.pkl'
        y_scaler_path = emulators_dir / f'{emulator_name}_output_scaler.pkl'

        emulator_files_exist: bool = (gp_path.exists() and x_scaler_path.exists() and y_scaler_path.exists()) # type: ignore
        climate_data_provided: bool = climate_data_file != ""

        if (not emulator_files_exist or force_retraining) and climate_data_provided: # type: ignore

            print(f"Training new emulator from: {climate_data_file}")
            
            # Robustly find the CSV
            csv_path = runs_dir / climate_data_file
            if not csv_path.exists():
                raise FileNotFoundError(f"Could not find training data at: {csv_path}")
            data = pd.read_csv(csv_path)
            data = data[data['Surface_Temp (K)'] != -1]

            # Identify input features: columns with non-constant values, excluding specified columns
            excluded_cols = {'Surface_Temp (K)', 'Status', 'Run_ID'}
            bool_cols = set(data.select_dtypes(include=['bool']).columns)
            candidate_cols = [col for col in data.columns if col not in excluded_cols and col not in bool_cols]

            input_features = []
            for col in candidate_cols:
                if data[col].nunique() > 1:
                    input_features.append(col)

            # Get indexes of pressure columns in input_features
            pressure_columns = ['P_Surface (Pa)', 'x_CO2', 'x_H2O']
            pressure_indexes = [input_features.index(col) for col in pressure_columns if col in input_features]

            # input_features = ['Instellation (W/m^2)', 'P_Surface (Pa)', 'x_CO2', 'x_H2O']
            output_targets = ['Surface_Temp (K)']

            X = data[input_features].values
            y = data[output_targets].values

            y = np.log10(y) # log scale y

            # log scale pressures
            for i in pressure_indexes:
                X[:, i] = np.log10(X[:, i])

            X_train_full, X_test, y_train_full, y_test = train_test_split(X, y, test_size=0.1, random_state=42)
            X_train, X_val, y_train, y_val = train_test_split(X_train_full, y_train_full, test_size=0.1, random_state=42)

            x_scaler = StandardScaler()
            y_scaler = StandardScaler()

            # FIT *ONLY* ON THE TRAINING DATA
            x_scaler.fit(X_train)
            y_scaler.fit(y_train)

            # TRANSFORM all three datasets
            X_train_scaled = x_scaler.transform(X_train)
            X_val_scaled   = x_scaler.transform(X_val)
            X_test_scaled  = x_scaler.transform(X_test)

            y_train_scaled = y_scaler.transform(y_train)
            y_val_scaled   = y_scaler.transform(y_val)
            y_test_scaled  = y_scaler.transform(y_test)

            # Set length_scale size to match number of input features
            n_features = X_train.shape[1]
            length_scale = [1.0] * n_features

            kernel = (
                ConstantKernel(1.0, (1e-3, 1e5)) * RBF(length_scale=length_scale, length_scale_bounds=(1e-2, 1e2)) 
                # + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-10, 1e-1))
            )

            gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, alpha=1e-3)

            print('Fitting Gaussian process...')
            gaussian_process.fit(X_train_scaled, y_train_scaled)
            print('Fitted Gaussian process')

            print(f"Learned kernel: {gaussian_process.kernel_}")
            print(f"Log-marginal-likelihood: {gaussian_process.log_marginal_likelihood(gaussian_process.kernel_.theta):.3f}")

            # 1. Score the model on the Validation set
            print(f"Training R2: {gaussian_process.score(X_train_scaled, y_train_scaled):.4f}")
            print(f"Validation R2: {gaussian_process.score(X_val_scaled, y_val_scaled):.4f}")

            # 2. Visual Check (Optional but recommended)

            if make_accuracy_plot:

                # 1. Get Predictions on Validation Data
                y_val_pred_scaled, y_val_std_scaled = gaussian_process.predict(X_val_scaled, return_std=True)
                
                # 2. Inverse Transform to get Kelvin
                y_val_pred_k = y_scaler.inverse_transform(y_val_pred_scaled.reshape(-1, 1))
                y_val_k = y_val # y_val is already unscaled in the split
                
                # 3. Calculate Residuals (Error)
                residuals = y_val_pred_k - y_val_k
                
                # 4. Set up the plotting grid (2x2 for 4 features)
                fig, axes = plt.subplots(2, 2, figsize=(12, 10))
                feature_names = input_features
                
                # X_val contains: [Instellation, log_P, log_CO2, log_H2O] 
                # (logs applied in lines 80-82 of emulator.py)
                
                axes = axes.flatten()
                
                for i, ax in enumerate(axes):
                    # Scatter plot: Input Feature vs Residual
                    ax.scatter(X_val[:, i], residuals, alpha=0.6, c='blue', edgecolors='k', s=20)
                    
                    # Reference lines
                    ax.axhline(0, color='black', lw=1, linestyle='-') # Perfect prediction
                    ax.axhline(1, color='red', lw=1, linestyle='--')  # +1 K limit
                    ax.axhline(-1, color='red', lw=1, linestyle='--') # -1 K limit
                    
                    ax.set_xlabel(feature_names[i])
                    ax.set_ylabel("Error (Predicted - Actual) [K]")
                    ax.set_title(f"Residuals vs {feature_names[i]}")
                    ax.grid(True, alpha=0.3)

                plt.tight_layout()
                plt.savefig("emulator_diagnostics.png")
                print("Saved emulator_diagnostics.png")
                
            self.gaussian_process = gaussian_process
            self.x_scaler = x_scaler
            self.y_scaler = y_scaler

            # Save to the package directory for future use
            print(f"Saving emulator to {emulators_dir}...")
            joblib.dump(self.gaussian_process, gp_path)
            joblib.dump(self.x_scaler, x_scaler_path)
            joblib.dump(self.y_scaler, y_scaler_path)
        
        else:
    
            print("Loading emulator...")

            if not gp_path.exists():
                raise FileNotFoundError(
                    f"Emulator '{emulator_name}' not found at {gp_path}.\n"
                    "Please provide a 'climate_data_file' to train it first."
                )

            self.gaussian_process = joblib.load(gp_path)
            self.x_scaler = joblib.load(x_scaler_path)
            self.y_scaler = joblib.load(y_scaler_path)

            print("Emulator loaded.")

    def get_temperature_from_emulator(self, instellation: float, P_surface: float, x_CO2: float, x_H2O: float) -> tuple[float, float]:
        """
        Calculates the surface temperature from the emulator with the direct input parameters. 

        Parameters
        ----------
        instellation : float
            Instellation flux in W/m^2.
        P_surface : float
            Surface pressure in Pa.
        x_CO2 : float
            CO2 volume mixing ratio.
        x_H2O : float
            H2O volume mixing ratio.
        albedo : float
            Albedo.

        Returns
        -------
        tuple[float, float]
            Surface temperature in K, Uncertainty
        """

        log_p = float(np.log10(P_surface))
        log_x_co2 = float(np.log10(x_CO2))
        log_x_h2o = float(np.log10(x_H2O))

        features_raw = np.array([[instellation, log_p, log_x_co2, log_x_h2o]])
        features_scaled = self.x_scaler.transform(features_raw)

        log_temp_scaled, log_std_scaled = self.gaussian_process.predict(features_scaled, return_std=True)
        
        # Unscale to get "Real" Log10(Temperature)
        # Note: We unscale the mean using the scaler's mean/scale logic
        log_temp = self.y_scaler.inverse_transform(log_temp_scaled.reshape(-1, 1))[0][0]
        
        # The standard deviation must be unscaled by multiplying by the scale factor
        # (Standard deviation is a width, not a point, so we don't add the mean)
        log_std = log_std_scaled[0] * self.y_scaler.scale_[0]

        # Convert to Kelvin
        temperature_kelvin = 10 ** log_temp
        
        # Calculate Uncertainty in Kelvin
        # We calculate the temperature at (Log_Mean + Log_Std) and subtract the mean
        temp_upper_bound = 10 ** (log_temp + log_std)
        uncertainty_kelvin = temp_upper_bound - temperature_kelvin

        # temp_scaled, std_scaled = self.gaussian_process.predict(features_scaled, return_std=True)

        # temp_kelvin = 10 ** self.y_scaler.inverse_transform(temp_scaled.reshape(-1, 1))

        # uncertainty_kelvin = std_scaled[0] * self.y_scaler.scale_[0]

        return float(temperature_kelvin), float(uncertainty_kelvin)
    
    def get_temperature(self, instellation: float, P_background: float, P_CO2: float, P_H2O: float) -> float:
        """
        Calculates the surface temperature from the emulator with partial pressures as input parameters.

        Parameters
        ----------
        instellation : float
            Instellation flux in W/m^2.
        P_background : float
            Partial pressure of background gas with no opacity (typically N2) in Pa.
        P_CO2 : float
            Partial pressure of CO2 in Pa.
        P_H2O : float
            Partial pressure of H2O in Pa.
        albedo : float
            Albedo.

        Returns
        -------
        float
            Surface temperature in K.
        """

        P_surface = P_background + P_CO2 + P_H2O
        x_CO2 = P_CO2 / P_surface
        x_H2O = P_H2O / P_surface

        return self.get_temperature_from_emulator(instellation, P_surface, x_CO2, x_H2O)[0]
    
    def make_temperature_pco2_interpolator(self, instellation: float, P_background: float):
        """
        Makes an interpolator for the surface temperature as function of P_CO2 only, calculating P_H2O with August-Roche-Magnus formula and all other parameters kept constant.

        Parameters
        ----------
        instellation : float
            Instellation flux in W/m^2.
        P_background : float
            Partial pressure of background gas with no opacity (typically N2) in Pa.
        albedo : float
            Albedo.

        Raises
        ------
        ValueError
            One of the data points for the intepolator didn't converge.
        """

        log_P_CO2_range = np.linspace(-6, 0)
        T_results = np.zeros_like(log_P_CO2_range)

        def residual(log_P_CO2: float, T: float) -> float:

            P_H2O = august_roche_magnus_formula(T)
            P_CO2 = 10 ** log_P_CO2
            T_from_climate = self.get_temperature(instellation, P_background, P_CO2, P_H2O)

            return T_from_climate - T
        
        for i, log_P_CO2 in enumerate(log_P_CO2_range):
            
            target_func = lambda T: residual(log_P_CO2, T) # type: ignore

            sol = least_squares(target_func, 300)

            if sol.success:
                T_results[i] = sol.x
            else:
                raise ValueError

        self.T_pco2_interpolator = CubicSpline(log_P_CO2_range, T_results)
    
    def get_temperature_from_pco2(self, P_CO2: float) -> float:
        """
        Calculates the surface temperature from the interpolator as a function of P_CO2 only. 

        Parameters
        ----------
        P_CO2 : float
            Partial pressure of CO2 in Pa.

        Returns
        -------
        float
            Surface temperature in K.
        """
        P_CO2 = np.maximum(P_CO2, 1e-6) # makes sure P_CO2 doesn't go below the minimum value of the interpolator
        return self.T_pco2_interpolator(np.log10(P_CO2)) # type: ignore
