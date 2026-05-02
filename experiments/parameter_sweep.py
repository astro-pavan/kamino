import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet as p2
from kamino.planet import Planet
from kamino.constants import M_EARTH, R_EARTH


def run_simulation(s, o, c, d, rw, cc, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}'
    if rw:
        run_name += '_rw'
    if cc > 0.0:
        run_name += f'_carb_{cc}'

    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=c,
            outgassing=o,
            ocean_depth=d,
            crust_carbonate_content=cc,
            reverse_weathering=rw,
            name=run_name
        )
        p.time_evolve()
        return run_name, None
    except Exception as e:
        return run_name, str(e)


def main():

    output_path = '/data/pt426/kamino_experiments'
    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    workers = 26

    instellation = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4]
    outgassing = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    crust_production_rate = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    ocean_depth = [100, 300, 1000, 3000, 10000]
    reverse_weathering = [False]
    crust_carbonate = [0]

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering, crust_carbonate))
    total = len(combos)

    print(f"Running {total} simulations with {workers} worker processes...")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(run_simulation, s, o, c, d, rw, cc, output_path): (s, o, c, d, rw, cc)
            for s, o, c, d, rw, cc in combos
        }

        completed = 0
        for future in as_completed(futures):
            completed += 1
            run_name, error = future.result()
            if error:
                print(f"[{completed}/{total}] FAILED {run_name}: {error}")
            else:
                print(f"[{completed}/{total}] Done: {run_name}")

    print("All simulations complete.")


if __name__ == "__main__":
    main()
