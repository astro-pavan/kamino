import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

# Adjust these imports if your directory structure requires it
from kamino.constants import *
from kamino.kamino_chem.ocean_chemistry import *

seafloor_flux = 1e12 / YR # 1 Tmol / yr
seafloor_flux_normalised = seafloor_flux / A_seafloor

print(f'Earth reference weathering flux: {seafloor_flux_normalised:.2e} mol/m^2/s')

T_ref = 280
P_ref = 1000 * 10 * 3000
P_CO2_ref = EARTH_ATM * 280e-6

print(f'Normalised hydrothermal flux: {J_ref_normalised:.2e} kg/s/m^2')

residual = lambda a: (get_weathering_flux(P_ref, T_ref, P_CO2_ref, np.array([]), a, rate_ref, J_ref_normalised)[0][0] - seafloor_flux_normalised) / seafloor_flux_normalised

from scipy.optimize import least_squares

root = least_squares(residual, 100)

print(root)

alpha_ref = float(root.x[0])

print(f'Alpha required for reference flux: {alpha_ref:.2e}')