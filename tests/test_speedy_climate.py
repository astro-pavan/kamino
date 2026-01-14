import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.emulator import climate_emulator, august_roche_magnus_formula
from kamino.constants import *

cem = climate_emulator("earth_rapid_rotator", "helios_1000_runs_earth_rapid_rotator.csv", make_accuracy_plot=True, force_retraining=True)

print(cem.get_temperature_from_emulator(1300, 1e5, 0.0001, 0.01))

cem.make_temperature_pco2_interpolator(SOLAR_CONSTANT, 1e5)

T = cem.get_temperature_from_pco2(30)
P_H2O = august_roche_magnus_formula(T)

print(f'T from interpolator : {T:.4f}')
print(f'P_H2O : {P_H2O:.2e}')

T_test = cem.get_temperature(SOLAR_CONSTANT, 1e5, 0.003 * 1e5, P_H2O)

print(f'T_test : {T_test:.4f}')