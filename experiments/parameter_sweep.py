import itertools
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet as p2
from kamino.planet import Planet
from kamino.constants import M_EARTH, R_EARTH
from kamino.mineral_info import *


def run_simulation(s, o, c, d, rw, cc, fht, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    cc_name, cc_dict = cc
    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}_comp_{cc_name}'
    if rw:
        run_name += '_rw'
    if fht > 0.0:
        run_name += f'_fht_{fht}'

    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=c,
            outgassing=o,
            ocean_depth=d,
            crust_composition=cc_dict,
            f_HT=fht,
            reverse_weathering=rw,
            name=run_name
        )
        p.time_evolve()
        return run_name, None
    except Exception as e:
        return run_name, str(e)


def run_sweep(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering, crust_composition, high_T_flux):

    output_path = '/data/pt426/kamino_experiments'
    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    workers = 26

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering, crust_composition, high_T_flux))
    total = len(combos)

    print(f"Running {total} simulations with {workers} worker processes...")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(run_simulation, s, o, c, d, rw, cc, fht, output_path): (s, o, c, d, rw, cc, fht)
            for s, o, c, d, rw, cc, fht in combos
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


instellation = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4]
outgassing = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
crust_production_rate = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
ocean_depth = [100, 300, 1000, 3000, 10000]
crust_composition = [
    ('komatiite_42', komatiite_42),
    ('komatiite_44', komatiite_44),
    ('basalt_47', basalt_47),
    ('basalt_49', basalt_49),
    ('basalt_51', basalt_51),
]
reverse_weathering = [False, True]
f_HT = [0, 0.001, 0.005, 0.01]

ocean_depth_default = [3000]
crust_composition_default = [('basalt_49', basalt_49)]
reverse_weathering_default = [True]
f_HT_default = [0]

calcium_only_crust = {'Wollastonite': 1}

basalt_49_no_iron = {
    'Anorthite': 0.1,
    'Wollastonite': 0.333 * 0.8,
    'Enstatite': 0.666 * 0.8,
    'Forsterite':  0.1,
} # 10% Olivine 80% Pyroxene 10% Plagioclase

magnesium_test_crust = [
    ('no_mg', calcium_only_crust),
    ('mg', basalt_49_no_iron)
]


if __name__ == "__main__":

    # Sweep 1: basic

    run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, reverse_weathering_default, crust_composition_default, f_HT_default)

    # Sweep 2: ocean depth

    run_sweep(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering_default, crust_composition_default, f_HT_default)

    # Sweep 3: crust composition

    run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, reverse_weathering_default, crust_composition, f_HT_default)

    # Sweep 4: reverse weathering (magnesium chemistry)

    run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, [True], [('mg', basalt_49_no_iron)], f_HT_default)
    run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, [False], [('no_mg', calcium_only_crust)], f_HT_default)

    # Sweep 5: high temperature hydrothermal

    # run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, reverse_weathering, crust_composition_default, f_HT)