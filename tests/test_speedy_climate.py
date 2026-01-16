import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.emulator import climate_emulator, august_roche_magnus_formula
from kamino.constants import *

import numpy as np
import pandas as pd

cem = climate_emulator("earth_rapid_rotator", "helios_3000_runs_earth_rapid_rotator.csv", make_accuracy_plot=True, force_retraining=False)

print(cem.get_temperature_from_emulator(SOLAR_CONSTANT, EARTH_ATM, 400e-6, 0.01))

cem.make_temperature_pco2_interpolator(SOLAR_CONSTANT, EARTH_ATM)
T = cem.get_temperature_from_pco2(40)

print(T)

df = pd.read_csv('/home/pt426/Code/kamino/src/kamino/speedy_climate/data/climate_runs/helios_3000_runs_earth_rapid_rotator.csv')

# 2. Define your Target Point (Earth-like)
target_instellation = SOLAR_CONSTANT
target_log_p = np.log10(1e5)
target_log_co2 = np.log10(400e-6)
target_log_h2o = np.log10(0.01)

# 3. Calculate "Distance" to every point in the dataset
# We normalize by the rough range of each parameter to make distance meaningful
df['dist'] = np.sqrt(
    ((df['Instellation (W/m^2)'] - target_instellation) / 1360)**2 +
    ((np.log10(df['P_Surface (Pa)']) - target_log_p) / 2)**2 + 
    ((np.log10(df['x_CO2']) - target_log_co2) / 6)**2 +
    ((np.log10(df['x_H2O']) - target_log_h2o) / 6)**2
)

# 4. Get the 5 Closest Neighbors
nearest = df.sort_values('dist').head(5)

print("--- Nearest Training Points to Earth-like Test Case ---")
print(nearest[['Instellation (W/m^2)', 'P_Surface (Pa)', 'x_CO2', 'x_H2O', 'Surface_Temp (K)', 'dist']])