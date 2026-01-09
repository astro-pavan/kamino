import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.climate import run_HELIOS
from kamino.constants import *

print(run_HELIOS('test', 1.3 * SOLAR_CONSTANT, 'G2', R_EARTH, M_EARTH, EARTH_ATM, 0.0004, 0.00004, 0.0, 0.25, clouds=False, verbose=True))

