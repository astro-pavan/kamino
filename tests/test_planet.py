import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

import matplotlib.pyplot as plt

plot_width = 15 # in
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.titlesize': 18})
plt.rcParams.update({'lines.linewidth': 2})

pl_H21 = planet(1e5, 3000, 0.95, 1.0, 0.05, 50e6, 100, tidally_locked=False, use_KT18_weathering=False)
T_s, pco2, C, A, Ca, weathering = pl_H21.find_steady_state_no_evolution(diagnostic_plots=False, solve_chemistry=True)

res_H21_1 = pl_H21.stability_analysis(0.1, T_s, pco2, C, A, Ca)
res_H21_2 = pl_H21.stability_analysis(-0.1, T_s, pco2, C, A, Ca)

# res_H21_1 = pl_H21.stability_analysis(0.0.05, T_s, pco2, C, A, Ca)
# res_H21_2 = pl_H21.stability_analysis(-0.1, T_s, pco2, C, A, Ca)

pl_WHAK = planet(1e5, 3000, 0.78, 1.0, 0.05, 50e6, 100, tidally_locked=False, use_KT18_weathering=True)
T_s, pco2, C, A, Ca, weathering = pl_WHAK.find_steady_state_no_evolution(diagnostic_plots=False, solve_chemistry=True)


res_WHAK_1 = pl_WHAK.stability_analysis(0.1, T_s, pco2, C, A, Ca)
res_WHAK_2 = pl_WHAK.stability_analysis(-0.1, T_s, pco2, C, A, Ca)



plt.subplots(figsize=(plot_width, 0.3 * plot_width))
plt.plot(res_WHAK_1['time'] / 1e6 - 0.5, res_WHAK_1['T_surface'], color='red', linestyle='--', label='Increased instellation (KT18)')
plt.plot(res_WHAK_2['time'] / 1e6 - 0.5, res_WHAK_2['T_surface'], color='darkblue', linestyle='--', label='Decreased instellation (KT18)')
plt.plot(res_H21_1['time'] / 1e6 - 0.5, res_H21_1['T_surface'], color='red', label='Increased instellation (H21)')
plt.plot(res_H21_2['time'] / 1e6 - 0.5, res_H21_2['T_surface'], color='darkblue', label='Decreased instellation (H21)')
plt.legend()
plt.xlim([0, 2.5])
plt.ylabel('T (K)')
plt.xlabel('t (Myr)')
plt.savefig('evol.pdf', bbox_inches='tight')
plt.show()

