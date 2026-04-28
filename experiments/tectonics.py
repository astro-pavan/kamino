import argparse
import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet2 as p2
from kamino.planet2 import Planet
from kamino.constants import M_EARTH, R_EARTH


def run_simulation(s, o, c, d, reverse_weathering, crust_carbonate_fraction, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}'
    if reverse_weathering:
        run_name += '_rw'
    if crust_carbonate_fraction > 0.0:
        run_name += f'_carb_{crust_carbonate_fraction}'

    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=c,
            outgassing=o,
            ocean_depth=d,
            crust_carbonate_content=crust_carbonate_fraction,
            reverse_weathering=reverse_weathering,
            name=run_name
        )
        p.time_evolve()
        return run_name, None
    except Exception as e:
        return run_name, str(e)


def main():
    parser = argparse.ArgumentParser(description="Sweep parameters for planet tectonics.")
    parser.add_argument('--reverse-weathering', action='store_true',
                        help="Enable reverse weathering.")
    parser.add_argument('--crust-carbonate-fraction', type=float, default=0.0,
                        help="Set the carbonate crust fraction (default: 0.0).")
    parser.add_argument('--workers', type=int, default=os.cpu_count(),
                        help="Number of parallel worker processes (default: CPU count).")
    args = parser.parse_args()

    output_path = '/home/pt426/kamino_experiments'
    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    instellation = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]
    outgassing = [0.01, 0.1, 1, 10]
    crust_production_rate = [0.01, 0.1, 1, 10]
    ocean_depth = [3000]

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth))
    total = len(combos)

    print(f"Running {total} simulations with {args.workers} worker processes...")

    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(
                run_simulation, s, o, c, d,
                args.reverse_weathering, args.crust_carbonate_fraction, output_path
            ): (s, o, c, d)
            for s, o, c, d in combos
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
