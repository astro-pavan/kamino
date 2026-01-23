import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

p1 = planet(1e5, 3000, 1.1)

# p1.run_simulation(1e3, True, [0.00005, 0.00005, 0.00005, 0.00005, 0.00005, 0.00005])

p1.find_steady_state(2000000)

# df = p1.run_simulation(100)
# df.to_csv('tests/test.csv')