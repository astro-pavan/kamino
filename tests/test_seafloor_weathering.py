import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

# import kamino.seafloor_weathering.chili.plot_example

from kamino.seafloor_weathering.weathering import *
from kamino.constants import YR

w_test = get_weathering_rate(1e5, 277, 0.02, 0.05 / YR, 100, 50e6 * YR)

print(w_test * YR)

