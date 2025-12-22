import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

import joblib
import os

class climate_emulator:

    def __init__(self, emulator_name, climate_data_file=None, make_accuracy_plot=False):
        
        if climate_data_file is not None:

            try:
                data = pd.read_csv(f'speedy_climate_data/climate_runs/{climate_data_file}')
            except FileNotFoundError:
                data = pd.read_csv(f'{os.getcwd()}/src/kamino/speedy_climate/data/climate_runs/{climate_data_file}')

            input_features = ['Instellation (W/m^2)', 'P_Surface (Pa)', 'x_CO2', 'x_H2O', 'Albedo']
            output_targets = ['Surface_Temp (K)']

            X = data[input_features].values
            y = data[output_targets].values

            # log scale P_surface, x_CO2 and x_H2O
            X[:, 1] = np.log10(X[:, 1])
            X[:, 2] = np.log10(X[:, 2])
            X[:, 3] = np.log10(X[:, 3])

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

            kernel = (
                ConstantKernel(1.0, (1e-3, 1e5)) * RBF(length_scale=[1.0, 1.0, 1.0, 1.0, 1.0], length_scale_bounds=(1e-2, 1e2)) 
                + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-10, 1e-1))
            )

            gaussian_process = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=20, alpha=1e-10)

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

                # Predict on validation data
                y_val_pred, y_val_std = gaussian_process.predict(X_val_scaled, return_std=True) # type: ignore

                # Plot Actual vs Predicted
                plt.figure(figsize=(6, 6))
                plt.scatter(y_val, y_scaler.inverse_transform(y_val_pred.reshape(-1, 1)), alpha=0.5)
                plt.plot([y_val.min(), y_val.max()], [y_val.min(), y_val.max()], 'r--', lw=2) # Perfect fit line
                plt.xlabel("Actual Surface Temp (K)")
                plt.ylabel("Predicted Surface Temp (K)")
                plt.title("Model Accuracy Check")
                plt.savefig("accuracy_check.png")
                print("Saved accuracy_check.png")

            joblib.dump(gaussian_process, f'speedy_climate_data/emulators/{emulator_name}_gp_climate_emulator.pkl')
            joblib.dump(x_scaler, f'speedy_climate_data/emulators/{emulator_name}_inputs_scaler.pkl')
            joblib.dump(y_scaler, f'speedy_climate_data/emulators/{emulator_name}_output_scaler.pkl')
        
        print("Loading emulator...")

        try:
            self.gaussian_process = joblib.load(f'speedy_climate_data/emulators/{emulator_name}_gp_climate_emulator.pkl')
            self.x_scaler = joblib.load(f'speedy_climate_data/emulators/{emulator_name}_inputs_scaler.pkl')
            self.y_scaler = joblib.load(f'speedy_climate_data/emulators/{emulator_name}_output_scaler.pkl')
        except FileNotFoundError:
            self.gaussian_process = joblib.load(f'{os.getcwd()}/src/kamino/speedy_climate/data/emulators/{emulator_name}_gp_climate_emulator.pkl')
            self.x_scaler = joblib.load(f'{os.getcwd()}/src/kamino/speedy_climate/data/emulators/{emulator_name}_inputs_scaler.pkl')
            self.y_scaler = joblib.load(f'{os.getcwd()}/src/kamino/speedy_climate/data/emulators/{emulator_name}_output_scaler.pkl')

        print("Emulator loaded.")

    def get_temperature(self, instellation, P_surface, x_CO2, x_H2O, albedo):

        log_p = np.log10(P_surface)
        log_x_co2 = np.log10(x_CO2)
        log_x_h2o = np.log10(x_H2O)

        features_raw = np.array([[instellation, log_p, log_x_co2, log_x_h2o, albedo]])
        features_scaled = self.x_scaler.transform(features_raw)

        temp_scaled, std_scaled = self.gaussian_process.predict(features_scaled, return_std=True)

        temp_kelvin = self.y_scaler.inverse_transform(temp_scaled.reshape(-1, 1))

        uncertainty_kelvin = std_scaled[0] * self.y_scaler.scale_[0]

        return temp_kelvin[0][0], uncertainty_kelvin


if __name__ == '__main__':

    cem = climate_emulator("earth_rapid_rotator")

    print(cem.get_temperature(1300, 1e5, 0.003, 0.01, 0.3))