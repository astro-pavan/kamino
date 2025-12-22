import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.emulator import climate_emulator

cem = climate_emulator("earth_rapid_rotator", "helios_1000_runs_earth_rapid_rotator.csv")

print(cem.get_temperature(1300, 1e5, 0.003, 0.01, 0.3))