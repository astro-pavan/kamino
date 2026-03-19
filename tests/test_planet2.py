import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet2 import Planet
from kamino.constants import *

p1 = Planet(M_EARTH, R_EARTH, 1e5, 1, 1, 3000)
res = p1.solve_steady_state()

print(res)