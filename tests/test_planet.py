import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.planet import planet
from kamino.constants import *

p1 = planet(1e5, 3000, 0.9)

p1.find_steady_state()

# df = p1.run_simulation(100)
# df.to_csv('tests/test.csv')