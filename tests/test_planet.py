import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

import matplotlib.pyplot as plt

pl = planet(1e5, 10000, 0.9, 1, 0.05, 50e6, 100, use_WHAK_weathering=False)
T_s, P_CO2, Co, Ao, Cao, T_weather = pl.find_steady_state_no_evolution(solve_chemistry=True)
pl.stability_analysis(T_s, P_CO2, Co, Ao, Cao)

