import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

# import kamino.seafloor_weathering.chili.plot_example

from kamino.seafloor_weathering.weathering import *
from kamino.constants import YR

T_range = np.linspace(274, 340, num=20)
pco2_range = np.logspace(-1, 4, num=20)
W = np.zeros((20,20))
W2 = np.zeros((20,20))
W3 = np.zeros((20,20))

for i, T in enumerate(T_range):
    for j, pco2 in enumerate(pco2_range):
        W[i, j] = get_weathering_rate(1e6, T, pco2 / 1e5, 0.05, 100, 50e6)
        W2[i, j] = get_weathering_rate_old(1e6, T, pco2 / 1e5)
        W3[i, j] = get_land_weathering_rate_old(1e6, T, pco2 / 1e5)

# plt.contourf(T_range, pco2_range, np.log10(W2.T), 200)
# plt.colorbar()
# plt.xlabel('T (K)')
# plt.ylabel('P_CO2 (Pa)')
# plt.yscale('log')
# plt.show()

plt.contourf(T_range, pco2_range, np.log10(W.T), 200, cmap='turbo')
plt.colorbar(label='log[Weathering rate (mol/m^2/yr)]')
plt.xlabel('T (K)')
plt.ylabel('P_CO2 (Pa)')
plt.yscale('log')
plt.show()

plt.contourf(T_range, pco2_range, ((W-W2)/W2).T, 200, cmap='turbo')
plt.colorbar(label='Weathering rate difference between H21 and KT18 (mol/m^2/yr)')
plt.contour(T_range, pco2_range, ((W-W2)/W2).T, [1], cmap='turbo')
plt.xlabel('T (K)')
plt.ylabel('P_CO2 (Pa)')
plt.yscale('log')
plt.show()

plt.contourf(T_range, pco2_range, ((W-W3)/W3).T, 200, cmap='turbo')
plt.colorbar(label='Weathering rate difference between H21 and WHAK (mol/m^2/yr)')
plt.contour(T_range, pco2_range, ((W-W3)/W3).T, [1], cmap='turbo')
plt.xlabel('T (K)')
plt.ylabel('P_CO2 (Pa)')
plt.yscale('log')
plt.show()