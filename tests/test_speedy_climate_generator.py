import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.generator import generate_input_parameters, run_batch_simulation
from kamino.constants import *

inputs = generate_input_parameters(
    200,
    (0.1 * SOLAR_CONSTANT, 1.5 * SOLAR_CONSTANT),
    5,
    (0, 6),
    3,
    0.05,
    0.6,
    'G2',
    0.25
)

run_batch_simulation(inputs, 'helios_n200_2d_earth_rapid_rotator.csv')