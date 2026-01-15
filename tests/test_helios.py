import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.climate import run_HELIOS
from kamino.constants import *

print(run_HELIOS('test', SOLAR_CONSTANT, 'G2', R_EARTH, M_EARTH, EARTH_ATM, 400 * 1e6, 0.01, 0.05, 0.25, clouds=0.55, verbose=False, august_roche_magnus=True))

