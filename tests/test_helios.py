import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.climate import run_HELIOS
from kamino.constants import *

print(run_HELIOS('test', 1 * SOLAR_CONSTANT, 'G2', R_EARTH, M_EARTH, EARTH_ATM, 0.0004, 0.004, 0.47, 0.25, verbose=True))

