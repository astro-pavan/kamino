import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

p1 = planet(1e5, 3000, 0.9)

p1.run_simulation(1e7, True, [0.2, 0.2, 0.2, 0.2, 0.2, 0.2])

# p1.find_steady_state(1e6)

# df = p1.run_simulation(100)
# df.to_csv('tests/test.csv')