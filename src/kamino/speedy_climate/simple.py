from scipy.interpolate import CloughTocher2DInterpolator
import numpy as np
import pandas as pd

from kamino.constants import *
from kamino.utils import *

data_path = '/home/pt426/Code/kamino/src/kamino/speedy_climate/data/climate_runs/dataset_2d.csv'

# Load data
df = pd.read_csv(data_path)
df = df[df['T'] != -1]

# Prepare data for interpolation
x_co2_log = np.log10(df['x_co2'].values)
S = df['S'].values
T = df['T'].values

# Create grid for interpolation
points = np.array([x_co2_log, S]).T

interp = CloughTocher2DInterpolator(points, T)

def get_T_surface(instellation, x_co2):
    log_x_co2 = np.log10(x_co2)
    log_x_co2 = smooth_max(-5.5, log_x_co2)
    log_x_co2 = smooth_min(-0.01, log_x_co2)
    T_val = interp(log_x_co2, instellation / SOLAR_CONSTANT)
    return float(T_val)