import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

# import kamino.seafloor_weathering.chili.plot_example

from kamino.seafloor_weathering.weathering import *
from kamino.constants import YR

T_range = np.linspace(274, 340, num=20)
W = []
W2 = []

for T in T_range:
    W.append(get_weathering_rate(1e5, T, 280 * 1e-6, 0.05, 100, 50e6))
    W2.append(get_weathering_rate_old(1e5, T, 280 * 1e-6))

plt.plot(T_range, W)
plt.plot(T_range, W2)
plt.yscale('log')
plt.show()