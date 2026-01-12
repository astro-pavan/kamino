import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.generator import generate_input_parameters, run_batch_simulation

inputs = generate_input_parameters(1000, 'G2', 0.25, 0.05)
run_batch_simulation(inputs, 'helios_1000_runs_earth_rapid_rotator.csv')

inputs = generate_input_parameters(1000, 'M5', 0.6666, 0.05)
run_batch_simulation(inputs, 'helios_1000_runs_earth_tidally_locked.csv')