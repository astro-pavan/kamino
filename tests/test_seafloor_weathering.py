import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

# import kamino.seafloor_weathering.chili.plot_example

from kamino.seafloor_weathering.weathering import *
from kamino.constants import YR

T_range = np.linspace(274, 340, num=20)
pco2_range = np.logspace(-1, 5, num=20)
W = np.zeros((20,20))
W2 = np.zeros((20,20))

T_range = np.linspace(274, 340, num=20)
pco2_range = np.logspace(-1, 5, num=20)

k_eff = np.zeros((20, 20))
C_eq = np.zeros((20, 20))

for i, T in enumerate(T_range):
    for j, pco2 in enumerate(pco2_range):
        
        P = 1e5
        x_CO2 = pco2 / P
        P = float(P / 1e5) # convert to bar
        P = np.clip(P, pr.P.min(), pr.P.max())
        x_CO2 = np.clip(x_CO2, pr.xCO2.min(), pr.xCO2.max())
        T = np.clip(T, pr.T.min(), pr.T.max())

        C_eq[j, i] = get_C_eq(1e6, T, x_CO2)
        k_eff[j, i] = get_k_eff(1e6, T, x_CO2)

plt.contourf(T_range, pco2_range, np.log10(k_eff), 200, cmap='inferno')
plt.colorbar(label='log[K_eff]')
plt.xlabel('T (K)')
plt.ylabel('P_CO2 (Pa)')
plt.yscale('log')
plt.show()

plt.contourf(T_range, pco2_range, np.log10(C_eq), 200, cmap='inferno')
plt.colorbar(label='log[C_eq]')
plt.xlabel('T (K)')
plt.ylabel('P_CO2 (Pa)')
plt.yscale('log')
plt.show()

# for i, T in enumerate(T_range):
#     for j, pco2 in enumerate(pco2_range):
#         W[i, j] = get_weathering_rate(1e6, T, pco2 / 1e5, 0.05, 100, 50e6)
#         W2[i, j] = get_weathering_rate_KT18(1e6, T, pco2 / 1e5)

# plt.contourf(T_range, pco2_range, np.log10(W.T), 200, cmap='turbo')
# plt.colorbar(label='log[Weathering rate (mol/m^2/yr)]')
# plt.xlabel('T (K)')
# plt.ylabel('P_CO2 (Pa)')
# plt.yscale('log')
# plt.show()

# plt.contourf(T_range, pco2_range, ((W-W2)/W2).T, 200, cmap='turbo')
# plt.colorbar(label='Difference between H21 and KT18 (mol/m^2/yr)')
# plt.contour(T_range, pco2_range, ((W-W2)/W2).T, [1], cmap='turbo')
# plt.xlabel('T (K)')
# plt.ylabel('P_CO2 (Pa)')
# plt.yscale('log')
# plt.show()

# t_rock_range = np.logspace(1, 8, num=20)
# flow_rate_range = np.logspace(-3, 0, num=20)
# W3 = np.zeros((20,20))

# for i, t_rock in enumerate(t_rock_range):
#     for j, flow_rate in enumerate(flow_rate_range):
#         W3[i, j] = get_weathering_rate(1e6, 288, 280 * 1e-6, flow_rate, 100, t_rock)

# plt.contourf(t_rock_range, flow_rate_range, np.log10(W3.T), 200, cmap='turbo')
# plt.colorbar(label='log[Weathering rate (mol/m^2/yr)]')
# plt.xlabel('t_rock (yr)')
# plt.ylabel('Flow rate (m/yr)')
# plt.yscale('log')
# plt.xscale('log')
# plt.show()

# T_range = np.linspace(274, 340, num=20)
# pco2_range = np.logspace(-1, 4, num=20)
# W = np.zeros((20,20))
# W2 = np.zeros((20,20))

# for i, T in enumerate(T_range):
#     for j, pco2 in enumerate(pco2_range):
#         W[i, j] = get_weathering_rate(1e6, T, pco2 / 1e5, 1, 10, 0)
#         W2[i, j] = get_weathering_rate_KT18(1e6, T, pco2 / 1e5)

# plt.contourf(T_range, pco2_range, np.log10(W.T), 200, cmap='turbo')
# plt.colorbar(label='log[Weathering rate (mol/m^2/yr)]')
# plt.xlabel('T (K)')
# plt.ylabel('P_CO2 (Pa)')
# plt.yscale('log')
# plt.show()

# plt.contourf(T_range, pco2_range, ((W-W2)/W2).T, 200, cmap='turbo')
# plt.colorbar(label='Difference between H21 and KT18 (mol/m^2/yr)')
# plt.contour(T_range, pco2_range, ((W-W2)/W2).T, [1], cmap='turbo')
# plt.xlabel('T (K)')
# plt.ylabel('P_CO2 (Pa)')
# plt.yscale('log')
# plt.show()
