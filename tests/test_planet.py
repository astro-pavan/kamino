import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

import matplotlib.pyplot as plt

# p1 = planet(1e5, 3000, 1.0)
# T_s, st = p1.find_steady_state(1000000)

# p1 = planet(1e5, 300, 1.0)
# T_s, st = p1.find_steady_state(5e6)
# T_s, pco2, Co, Cp, Ao, Ap, Cao, Cap = p1.find_steady_state_simple()

S = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
T = []
co2 = []
DIC = []
stable = []

plt.figure(figsize=(20, 20))

for s in S:

    p1 = planet(1e5, 300, s, 0.05, 50e6)
    # T_s, st = p1.find_steady_state(7e4)
    T_s, pco2 = p1.find_steady_state_no_evolution()
    T.append(T_s)
    co2.append(pco2)
    # DIC.append(Co)
    # stable.append(st)

# plt.savefig('Rates2.png')

fig, ax1 = plt.subplots(figsize=(20, 20))

ax1.plot(S, T, label='Temperature (T)', marker='o', color='tab:blue')
ax1.set_xlabel('S')
ax1.set_ylabel('Temperature (T)', color='tab:blue')
ax1.axhline(273, color='tab:blue')
ax1.axhline(340, color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax2 = ax1.twinx()
ax2.plot(S, co2, label='P_CO2', marker='s', color='tab:orange')
ax2.set_ylabel('P_CO2', color='tab:orange')
ax2.set_yscale('log')
ax2.tick_params(axis='y', labelcolor='tab:orange')

fig.tight_layout()
plt.show()
