import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from kamino.utils import *
from kamino.constants import*

data_path = 'src/kamino/speedy_climate/data/climate_runs/clima_data_grid_rapid.csv'

try:
    df = pd.read_csv(data_path)
except FileNotFoundError:
    df = pd.read_csv('../' + data_path)

# 1. Define your grid axes
dims = ['Instellation (S0)', 'P_CO2 (bar)', 'Albedo']
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

def get_T_surface(S, P_CO2, albedo):
    P_CO2 = smooth_max(1e-1, P_CO2)
    point = [S / SOLAR_CONSTANT, float(P_CO2) * 1e-5, albedo]
    temp = rgi(point)[0]
    assert ~np.isnan(temp)
    return float(temp)