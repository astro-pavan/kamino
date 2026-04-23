from kamino.planet import Planet
from kamino.constants import *

BACKGROUND_PRESSURE = 1e5   # Pa (~1 bar)
OCEAN_DEPTH = 3000          # m
TECTONICS = 1.0

# instellation_range = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25]
instellation_range = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
# tectonics_range = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
ocean_depth_range = [100, 500, 1000, 3000, 5000, 10000]

for s in instellation_range:
    print(f's = {s}')
    p1 = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, s, 1.0, 3000, name=f'test_s_{s}')
    p1.time_evolve()
    print('')