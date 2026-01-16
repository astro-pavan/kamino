import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.speedy_climate.simple import get_T_surface
from kamino.constants import *

import numpy as np
import matplotlib.pyplot as plt

x_co2_range = np.logspace(-5.5, 0, num=30)
S_range = np.linspace(0.3, 1.5, num=30) * SOLAR_CONSTANT

x_co2_arr = []
S_arr = []
T_arr = []

for x_co2 in x_co2_range:
    for S in S_range:

        x_co2_arr.append(x_co2)
        S_arr.append(S)
        T_arr.append(get_T_surface(S, x_co2))

fig, ax = plt.subplots(figsize=(10, 6))

# Define variables for clarity
x = np.array(x_co2_arr)
y = np.array(S_arr) / SOLAR_CONSTANT
z = np.array(T_arr)

# Create filled contours (tricontourf is best for column data)
contour_filled = ax.tricontourf(x, y, z, levels=200, cmap='plasma')
contour_lines = ax.tricontour(x, y, z, [273, 288, 340], colors='k', linestyles='dotted')
ax.clabel(contour_lines, fmt='%d K', colors='k', inline=False, use_clabeltext=True)
ax.set_xscale('log')

# ax.axvline(280e-6, label='Pre-Industrial CO2', color='k', linestyle='--')
ax.scatter([280e-6], [1], s=10, c='black', label='Earth')
ax.annotate('Earth', (280e-6, 1))

# Add a colorbar to show what T values the colors represent
cbar = plt.colorbar(contour_filled, ax=ax)
cbar.set_label('Temperature (K)')

# Axis Labels
ax.set_xlabel('$x_{CO2}$')
ax.set_ylabel('S ($S_{\\oplus}$)')

# Optional: Since x_co2 is very small, scientific notation is automatic,
# but we can force a log scale if your x data spans many orders of magnitude.
# ax.set_xscale('log') 

plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig('climate.png')