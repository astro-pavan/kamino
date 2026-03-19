import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np

from kamino.kamino_chem.ocean_chemistry import *

b_seawater = np.array([2.3e-3, 2.0e-3, 1e-4, 0.0, 0.0, 10.3e-3, 52.7e-3, 468e-3])

get_P_CO2(1e5, 300, b_seawater)