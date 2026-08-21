import itertools
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault('JAX_PLATFORMS', 'cpu')

import kamino.planet as p2
from kamino.planet import Planet, KD_MG_HT, K_NA_CONT_REMOVAL
from kamino.weathering import ALPHA_REF
from kamino.crust_composition import mineral_composition, MG_SI_REF
from kamino.constants import M_EARTH, R_EARTH
from kamino.mineral_info import *

RERUN = False
MAX_CHEMISTRY_FALLBACKS = 5000

OUTPUT_PATH = '/data/pt426/kamino_experiments_fast_18'


def _tag(value, reference, prefix):
    """Name fragment for a swept parameter; empty at its reference so old names are reproduced."""
    return '' if value == reference else f'_{prefix}{value:g}'


def _run_name(s, o, c, d, rw, mt, alpha, kd_mg, k_na):
    """Run name. Every parameter that differs from the Planet default MUST appear, or two configs
    would share a filename and RERUN=False would silently return the first one's result."""
    run_name = f'planet_s_{s}_out_{o}_crust_{c}_depth_{d}'

    if rw:
        run_name += '_rw'

    run_name += f'_mt{mt}'
    run_name += _tag(alpha, ALPHA_REF, 'a')
    run_name += _tag(kd_mg, KD_MG_HT, 'kmg')
    run_name += _tag(k_na, K_NA_CONT_REMOVAL, 'kna')

    return run_name


def run_simulation(s, o, c, d, rw, mt, alpha, kd_mg, k_na, output_path):
    p2.output_path = output_path  # each subprocess imports a fresh module; set path here

    run_name = _run_name(s, o, c, d, rw, mt, alpha, kd_mg, k_na)

    if not RERUN:
        json_path = os.path.join(output_path, f'{run_name}.json')
        if os.path.exists(json_path):
            try:
                with open(json_path) as fh:
                    existing = json.load(fh)
                if 'termination' in existing:
                    return run_name, None, existing.get('T'), existing['termination']
            except Exception:
                pass  # corrupt/incomplete file — fall through and re-run

    try:
        p = Planet(
            mass=M_EARTH,
            radius=R_EARTH,
            background_pressure=1e5,
            instellation=s,
            crust_production_rate=c,
            outgassing=o,
            ocean_depth=d,
            reverse_weathering=rw,
            mantle_potential_temperature=mt,
            alpha=alpha,
            kd_mg_ht=kd_mg,
            k_na_cont_removal=k_na,
            name=run_name
        )
        p.time_evolve(max_chemistry_fallbacks=MAX_CHEMISTRY_FALLBACKS)
        with open(p._output_filename) as fh:  # time_evolve records T and termination here
            result = json.load(fh)
        return run_name, None, result.get('T'), result.get('termination')
    except Exception as e:
        return run_name, str(e), None, None


def run_sweep(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering,
              mantle_temperature, alpha=(ALPHA_REF,), kd_mg=(KD_MG_HT,),
              k_na=(K_NA_CONT_REMOVAL,), output_path=OUTPUT_PATH):

    if not output_path.endswith('/'):
        output_path += '/'
    p2.output_path = output_path
    os.makedirs(output_path, exist_ok=True)

    workers = 26

    combos = list(itertools.product(instellation, outgassing, crust_production_rate, ocean_depth, reverse_weathering, mantle_temperature, alpha, kd_mg, k_na))
    total = len(combos)

    # Distinct configs must map to distinct filenames, or one silently overwrites the other and
    # RERUN=False then returns the survivor's result for both (the fast_13 resume trap).
    names = [_run_name(*combo) for combo in combos]
    if len(set(names)) != len(names):
        duplicated = sorted({n for n in names if names.count(n) > 1})
        raise ValueError(
            f"{len(names) - len(set(names))} run name collision(s) in this grid, e.g. "
            f"{duplicated[:3]}. Two configs would share an output file."
        )

    cap_str = MAX_CHEMISTRY_FALLBACKS if MAX_CHEMISTRY_FALLBACKS is not None else 'disabled'
    print(f"Running {total} simulations with {workers} worker processes "
          f"(fallback cap: {cap_str})...")
    print(f"Output: {output_path}")
    for label, values in (('alpha', alpha), ('kd_mg_ht', kd_mg), ('k_na_cont_removal', k_na)):
        print(f"  {label}: {list(values)}")

    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as executor:
        futures = {
            executor.submit(run_simulation, s, o, c, d, rw, mt, a, kmg, kna, output_path): (s, o, c, d, rw, mt, a, kmg, kna)
            for s, o, c, d, rw, mt, a, kmg, kna in combos
        }

        completed = 0
        aborted = 0
        for future in as_completed(futures):
            completed += 1
            run_name, error, T, termination = future.result()
            if error:
                print(f"[{completed}/{total}] FAILED {run_name}: {error}")
            else:
                if termination == 'fallback_limit':
                    aborted += 1
                T_str = f"{T:.1f} K" if T is not None else "T unknown"
                print(f"[{completed}/{total}] Done: {run_name} ({T_str}, {termination or 'unknown'})")

    print("All simulations complete.")
    if aborted:
        # Not an error -- these are the runs the cap was added to stop. A large fraction here
        # means the chemistry is unhealthy over a wide part of the grid, not that the cap is wrong.
        print(f"{aborted}/{total} run(s) hit the fallback cap ({MAX_CHEMISTRY_FALLBACKS}) "
              f"and were recorded as 'fallback_limit'.")

fast = False

if fast:
    instellation = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    outgassing = [0.03, 0.1, 1, 3]
    crust_production_rate = [0.01, 0.1, 1, 10]
    ocean_depth = [300, 3000, 30000]
    mantle_temp = [1350, 1475, 1600]
else:
    instellation = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2]
    outgassing = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    crust_production_rate = [0.01, 0.03, 0.1, 0.3, 1, 3, 10]
    ocean_depth = [300, 500, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000]
    mantle_temp = [1350, 1412, 1475, 1538, 1600]


reverse_weathering = [False, True]
alpha = [2, 20, 200]
k_mg = [0.02, 0]   # 0 disables the Mg->Ca exchange entirely
k_na = [0.004, 0]  # 0 disables the Na sink entirely

crust_production_rate_default = [1]
outgassing_default = [0.1]
ocean_depth_default = [3000]
reverse_weathering_default = [True]
# Earth calibration at alpha = 19.92 gave kd_mg_ht = 1.894210e-02, k_na = 3.893608e-03.
alpha_default = [2]
k_mg_default = [0.02]
k_na_default = [0.004]
mantle_temp_default = [1350]


if __name__ == "__main__":

    # Sweep 1: basic

    # print('Running basic sweep')
    # run_sweep(instellation, outgassing, crust_production_rate, ocean_depth_default, reverse_weathering_default, mantle_temp_default, alpha_default, k_mg_default, k_na_default)

    # Sweep 2: ocean depth

    # print('Running ocean depth sweep')
    # run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth, reverse_weathering_default, mantle_temp_default, alpha_default, k_mg_default, k_na_default)

    # Sweep 3: mantle temp

    # print('Running mantle temp sweep')
    # run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth_default, reverse_weathering_default, mantle_temp, alpha_default, k_mg_default, k_na_default)

    # Sweep 4: alpha
    # print('Running alpha sweep')
    # run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth_default, reverse_weathering_default, mantle_temp_default, alpha, k_mg_default, k_na_default)

    # Sweep 5: chemistry constants
    print('Running chemistry constant sweep')
    run_sweep(instellation, outgassing_default, crust_production_rate_default, ocean_depth_default, reverse_weathering_default, mantle_temp_default, alpha_default, k_mg, k_na)

