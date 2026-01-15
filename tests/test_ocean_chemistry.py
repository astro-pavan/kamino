import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.ocean_chemistry.aqueous_geochemistry import *
from kamino.ocean_chemistry.precipitation import *
from kamino.ocean_chemistry.co2 import *

debug_mode = True

sal = 1

comp: dict[str, float] = {
    'Cl': 0.546 * sal,
    'Na': 0.469 * sal,
    'Mg': 0.0528 * sal,
    'S': 0.0282 * sal,
    'Ca': 0.0103 * sal,
    'K': 0.0102 * sal,
    'Si': 0.0000001 * sal,
    'Al': 0.0000001 * sal,
    'C': 0.002 * sal
}

# output = kinetics(10 * EARTH_ATM, 280, comp, None, {'Calcite' : 0.0}, 1e6)

# print(output)

k, SI = get_calcite_precipitation_rate(20 * EARTH_ATM, 274, 0.001, 0.0002, 0.0002)
print(k)
print(SI)

# p_CO2 = get_P_CO2(1e5, 280, 0.002, 0.002, 0)
# print(f'{(p_CO2 / 1e5) * 1e6} ppm')