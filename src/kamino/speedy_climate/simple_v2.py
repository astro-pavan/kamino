import pandas as pd
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RBFInterpolator, RegularGridInterpolator

from kamino.utils import *

data_path_original = '/home/pt426/Code/kamino/src/kamino/speedy_climate/data/climate_runs/dataset_4d.csv'
data_path = '/home/pt426/Code/kamino/src/kamino/speedy_climate/data/climate_runs/dataset_4d_v2.csv'

def fill_in_gaps():

    # Load data
    df = pd.read_csv(data_path_original)

    # 1. Filter Data
    df_success = df[df['Status'] == 'Success']
    df_failed = df[df['Status'] != 'Success']

    # 2. Prepare Inputs/Outputs
    input_cols = ['Instellation (W/m^2)', 'P_CO2 (Pa)', 'P_H2O (Pa)', 'Albedo']
    output_col = 'Surface_Temp (K)'

    X_train = df_success[input_cols].values
    y_train = df_success[output_col].values
    X_target = df_failed[input_cols].values

    # 3. Normalize Data (Crucial for multi-scale inputs)
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1 # Prevent division by zero

    X_train_norm = (X_train - mean) / std
    X_target_norm = (X_target - mean) / std

    # 4. Interpolate
    # LinearND is best for filling internal holes in a grid
    interp = LinearNDInterpolator(X_train_norm, y_train)
    y_pred = interp(X_target_norm)

    # Handle potential extrapolation (NaNs) with RBF if LinearND fails
    if np.isnan(y_pred).any():
        print("LinearND produced NaNs (extrapolation needed). Using RBF for remaining points.")
        rbf = RBFInterpolator(X_train_norm, y_train, kernel='linear')
        y_rbf = rbf(X_target_norm)
        
        # Fill NaNs in y_pred with RBF results
        mask = np.isnan(y_pred)
        y_pred[mask] = y_rbf[mask]

    # 5. Update DataFrame
    df.loc[df['Status'] != 'Success', output_col] = y_pred
    df.loc[df['Status'] != 'Success', 'Status'] = 'Interpolated'

    df.to_csv(data_path, index=False)

# Load the filled dataset (or use the code below to fill it on the fly)
df = pd.read_csv(data_path)

# 1. Define your grid axes
dims = ['Instellation (W/m^2)', 'P_CO2 (Pa)', 'P_H2O (Pa)', 'Albedo']
target = 'Surface_Temp (K)'

# Extract unique, sorted values for each dimension
unique_vals = [np.sort(df[col].unique()) for col in dims]

# 2. Reshape the data into the grid structure
# It is CRITICAL to sort the dataframe so the values line up with the grid axes
df_sorted = df.sort_values(by=dims)
grid_values = df_sorted[target].values.reshape([len(u) for u in unique_vals])

# 3. Create the Interpolator
# bounds_error=False allows extrapolation (returns NaN or nearest), 
# fill_value=None with bounds_error=False extrapolates using the grid edge.
rgi = RegularGridInterpolator(unique_vals, grid_values, bounds_error=True, method='linear')

def get_T_surface(S, P_CO2, P_H2O, albedo):
    P_CO2 = smooth_max(1e-1, P_CO2)
    P_H2O = smooth_max(1e-1, P_H2O)
    point = [S, float(P_CO2), float(P_H2O), albedo]
    temp = rgi(point)
    assert ~np.isnan(temp)
    return float(temp) - 15
