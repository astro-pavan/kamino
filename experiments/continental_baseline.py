"""Run the Earth-like continental baseline sweep.

Fixes: basalt_49, rw=True, out=1×, crust=1×, depth=3000 m, land_fraction=0.3.
Sweeps instellation from 0.3 to 1.4 in the same steps as the main parameter sweep.
"""

import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import kamino.planet as p_mod
from kamino.planet import Planet
from kamino.constants import M_EARTH, R_EARTH
from kamino.mineral_info import basalt_49

OUTPUT_PATH = '/data/pt426/kamino_experiments_fast_3/'

INSTELLATION = [
    0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
    0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4,
]


def run_simulation(s, output_path):
    p_mod.output_path = output_path
    run_name = f'planet_s_{s}_out_1_crust_1_depth_3000_comp_basalt_49_rw_land'
    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=1.0,
            outgassing=1.0,
            ocean_depth=3000,
            land_fraction=0.3,
            crust_composition=basalt_49,
            reverse_weathering=True,
            name=run_name,
            f_HT=0.01
        )
        p.time_evolve()
        return run_name, None
    except Exception as e:
        return run_name, str(e)


if __name__ == '__main__':
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    p_mod.output_path = OUTPUT_PATH

    total = len(INSTELLATION)
    print(f"Running {total} continental baseline simulations...")

    with ProcessPoolExecutor(max_workers=total, mp_context=mp.get_context('spawn')) as executor:
        futures = {executor.submit(run_simulation, s, OUTPUT_PATH): s for s in INSTELLATION}
        completed = 0
        for future in as_completed(futures):
            completed += 1
            run_name, error = future.result()
            if error:
                print(f"[{completed}/{total}] FAILED {run_name}: {error}")
            else:
                print(f"[{completed}/{total}] Done: {run_name}")

    print("All simulations complete.")
