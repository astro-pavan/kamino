import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet2 import Planet
from kamino.constants import *

p1 = Planet(M_EARTH, R_EARTH, 1e5, 1.0, 0.1, 3000)
p1.time_evolve_to_steady_state()