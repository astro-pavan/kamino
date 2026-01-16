from scipy.interpolate import CloughTocher2DInterpolator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from kamino.constants import *
from kamino.utils import *

data_path = '/home/pt426/Code/kamino/src/kamino/speedy_climate/data/climate_runs/dataset_2d_v2.csv'

# Load data
df = pd.read_csv(data_path)
df = df[df['Surface_Temp (K)'] != -1]

# Prepare data for interpolation
p_co2_log = np.log10(df['P_CO2 (Pa)'].values)
S = df['Instellation (W/m^2)'].values
T = df['Surface_Temp (K)'].values

fig, ax = plt.subplots(figsize=(10, 6))
contour_filled = ax.tricontourf(p_co2_log, S, T, levels=200, cmap='plasma')
ax.scatter(p_co2_log, S, c=T, cmap='plasma')
cbar = plt.colorbar(contour_filled, ax=ax)
plt.savefig('Interpolator.png')

# Create grid for interpolation
points = np.array([p_co2_log, S]).T

interp = CloughTocher2DInterpolator(points, T, fill_value=np.nan)

def get_T_surface(instellation, p_co2):
    p_co2 = smooth_max(1.1, p_co2)
    log_p_co2 = np.log10(p_co2)
    # log_x_co2 = smooth_max(-5.5, log_p_co2)
    log_p_co2 = smooth_min(6, log_p_co2)
    T_val = interp(log_p_co2, instellation)
    return float(T_val)