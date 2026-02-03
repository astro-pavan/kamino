import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from kamino.utils import *
from kamino.constants import*

data_path_rapid = 'src/kamino/speedy_climate/data/climate_runs/clima_data_grid_rapid.csv'
data_path_tidal = 'src/kamino/speedy_climate/data/climate_runs/clima_data_grid_tidal.csv'

try:
    df_rapid = pd.read_csv(data_path_rapid)
    df_tidal = pd.read_csv(data_path_tidal)
except FileNotFoundError:
    df_rapid = pd.read_csv('../' + data_path_rapid)
    df_tidal = pd.read_csv('../' + data_path_tidal)

# 1. Define your grid axes
dims = ['Instellation (S0)', 'P_CO2 (bar)', 'Albedo']
target = 'Surface_Temp (K)'

# Extract unique, sorted values for each dimension
unique_vals = [np.sort(df_rapid[col].unique()) for col in dims]

# 2. Reshape the data into the grid structure
# It is CRITICAL to sort the dataframe so the values line up with the grid axes
df_sorted = df_rapid.sort_values(by=dims)
grid_values = df_sorted[target].values.reshape([len(u) for u in unique_vals])

# 3. Create the Interpolator
# bounds_error=False allows extrapolation (returns NaN or nearest), 
# fill_value=None with bounds_error=False extrapolates using the grid edge.
rgi_rapid = RegularGridInterpolator(unique_vals, grid_values, bounds_error=False, method='linear')

# Extract unique, sorted values for each dimension
unique_vals = [np.sort(df_tidal[col].unique()) for col in dims]

# 2. Reshape the data into the grid structure
# It is CRITICAL to sort the dataframe so the values line up with the grid axes
df_sorted = df_tidal.sort_values(by=dims)
grid_values = df_sorted[target].values.reshape([len(u) for u in unique_vals])

# 3. Create the Interpolator
# bounds_error=False allows extrapolation (returns NaN or nearest), 
# fill_value=None with bounds_error=False extrapolates using the grid edge.
rgi_tidal = RegularGridInterpolator(unique_vals, grid_values, bounds_error=False, method='linear')

def get_T_surface(S, P_CO2, albedo, tidally_locked=False):
    
    P_CO2 = smooth_max(1e-2, P_CO2) # 0.01 Pa minimum
    albedo = np.clip(albedo, 0.0, 0.5) # the interpolator only goes to 0.5
    point = [S / SOLAR_CONSTANT, float(P_CO2) * 1e-5, albedo]

    if tidally_locked:
        temp = rgi_tidal(point)[0]
    else:
        temp = rgi_rapid(point)[0]

    assert ~np.isnan(temp) # NEDD TO REPLACE WITH BETTER ERROR HANDLING
    
    return float(temp)