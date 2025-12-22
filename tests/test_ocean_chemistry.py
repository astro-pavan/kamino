import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

# from kamino.ocean_chemistry.build_phreeqc import build

# build()

from kamino.ocean_chemistry.aqueous_geochemistry import *

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

output = kinetics(10 * EARTH_ATM, 280, comp, None, {'Calcite' : 0.0}, 1e6)

print(output)