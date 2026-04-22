import numpy as _np
if not hasattr(_np, 'in1d'):
    _np.in1d = lambda ar1, ar2, assume_unique=False, invert=False: _np.isin(ar1, ar2, assume_unique=assume_unique, invert=invert).ravel()

import exo_k as xk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

from kamino.utils import august_roche_magnus_formula
from kamino.constants import *

import csv
import itertools
import multiprocessing
import os

k_data_path = os.path.join(os.path.dirname(__file__), 'data/data_tutorial_exo-k_interior/')
runs_data_path = os.path.join(os.path.dirname(__file__), 'data/climate_runs/exo_k_sweep.csv')

xk.Settings().set_mks(True)

xk.Settings().set_search_path(k_data_path+'corrk', path_type='ktable')
xk.Settings().set_search_path(k_data_path+'cia', path_type='cia')

molecules = ['H2O', 'CO2']
k_db = xk.Kdatabase(molecules, search_path=k_data_path+'corrk/R70_homogeneous_from_R500', remove_zeros=True)

cia_db = xk.CIAdatabase(molecule_pairs=[['CO2','CO2']], search_path=k_data_path+'cia')
cia_db.sample(k_db)
cia_db.convert_to_mks()
print(cia_db)

def run_exo_k_atmosphere(instellation, P_background, P_H2O, P_CO2, surface_albedo, bond_albedo=0.3, recirculation_factor=0.25):

    grav_earth = 9.81
    mean_insolation = recirculation_factor * instellation * SOLAR_CONSTANT

    # convert to Pa (keep bar values for output)
    P_background_Pa = P_background * 1e5
    P_H2O_Pa = P_H2O * 1e5
    P_CO2_Pa = P_CO2 * 1e5

    P_surface = P_background_Pa + P_H2O_Pa + P_CO2_Pa

    H2O_vmr = P_H2O_Pa / P_surface
    CO2_vmr = P_CO2_Pa / P_surface

    composition = {'N2': 'background', 'H2O': H2O_vmr, 'CO2': CO2_vmr}

    atm = xk.Atm_evolution(
        Nlay=50,                  # Number of vertical layers
        Tsurf=288.0,              # Initial surface temperature guess
        Tstrat=200.0,             # Initial stratosphere temperature guess
        psurf=P_surface,          # Surface pressure in Pa
        ptop=10.0,                # Top of atmosphere pressure (0.1 mbar)
        grav=grav_earth,          # Surface gravity
        albedo_surf=surface_albedo,
        rcp=0.28,                 # R/Cp for the gas mixture
        bg_vmr=composition,
        k_database=k_db,
        cia_database=cia_db,        # Can be replaced with xk.CIAdatabase(...) if you have CIA tables
        internal_flux=0.0,        # No internal heat for an Earth-like rocky planet
        flux_top_dw=mean_insolation * (1.0 - bond_albedo), 
        verbose=False,
        rayleigh=True,            # Enable Rayleigh scattering
        convection=True,          # Enable convective adjustment
        acceleration_mode=4,
        dTmax_use_kernel=50.
    )

    atm.equilibrate(
        Fnet_tolerance=0.1,       # Target net flux tolerance in W/m^2
        N_timestep_ini=100,       # Initial batch of timesteps before kernel evaluation
        N_iter_max=20,            # Max kernel iterations
        N_timestep_max=40000,     # Max total timesteps
        verbose=False
    )

    T_surface = atm.tlay[-1]

    result = {
        'Instellation': instellation,
        'P_background (bar)': P_background,
        'P_H2O (bar)': P_H2O,
        'P_CO2 (bar)': P_CO2,
        'Surface Albedo': surface_albedo,
        'Bond Albedo': bond_albedo,
        'Recirculation Factor': recirculation_factor,
        'Surface Temperature (K)': T_surface
    }

    return result

def _run_single(args):
    instellation, P_background, P_H2O, P_CO2, surface_albedo = args
    try:
        return run_exo_k_atmosphere(*args)
    except Exception as e:
        print(f'\nFailed ({instellation:.2f}, {P_background:.2e}, {P_H2O:.2e}, {P_CO2:.2e}): {e}')
        return {
            'Instellation': instellation,
            'P_background (bar)': P_background,
            'P_H2O (bar)': P_H2O,
            'P_CO2 (bar)': P_CO2,
            'Surface Albedo': surface_albedo,
            'Bond Albedo': float('nan'),
            'Recirculation Factor': float('nan'),
            'Surface Temperature (K)': float('nan'),
        }

def parameter_sweep(csv_path=runs_data_path, n_workers=None):
    instellation_range = np.linspace(0.1, 1.4, num=14)
    P_H2O_range = np.logspace(-6, 1, num=8)
    P_CO2_range = np.logspace(-6, 1, num=8)
    surface_albedo = 0.05
    P_background = 1.0  # bar

    args = [
        (inst, P_background, p_h2o, p_co2, surface_albedo)
        for inst, p_h2o, p_co2 in itertools.product(instellation_range, P_H2O_range, P_CO2_range)
    ]
    total = len(args)
    print(f'Running {total} simulations with {n_workers or os.cpu_count()} workers...')

    results = []
    with multiprocessing.Pool(processes=n_workers) as pool:
        for i, result in enumerate(pool.imap_unordered(_run_single, args), 1):
            results.append(result)
            print(f'  {i}/{total} complete', end='\r', flush=True)
    print()

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f'Saved {len(results)} results to {csv_path}')
    return results

def make_basic_interpolator(csv_path=runs_data_path):

    df = pd.read_csv(csv_path)

    instellation_vals = np.sort(df['Instellation'].unique())
    P_H2O_vals       = np.sort(df['P_H2O (bar)'].unique())
    P_CO2_vals        = np.sort(df['P_CO2 (bar)'].unique())

    # sort rows to match grid order, then reshape into (inst, H2O, CO2)
    df_sorted = df.sort_values(['Instellation', 'P_H2O (bar)', 'P_CO2 (bar)'])
    T_grid = df_sorted['Surface Temperature (K)'].values.reshape(len(instellation_vals), len(P_H2O_vals), len(P_CO2_vals))

    # fill the handful of NaN cells (failed runs) with nearest-neighbour values
    # so they don't poison bilinear interpolation in neighbouring regions
    from scipy.interpolate import NearestNDInterpolator
    inst_g, h2o_g, co2_g = np.meshgrid(
        instellation_vals, np.log10(P_H2O_vals), np.log10(P_CO2_vals), indexing='ij'
    )
    valid = ~np.isnan(T_grid)
    if not valid.all():
        nn = NearestNDInterpolator(
            np.column_stack([inst_g[valid], h2o_g[valid], co2_g[valid]]),
            T_grid[valid],
        )
        T_grid[~valid] = nn(
            np.column_stack([inst_g[~valid], h2o_g[~valid], co2_g[~valid]])
        )

    _interp = RegularGridInterpolator(
        (instellation_vals, np.log10(P_H2O_vals), np.log10(P_CO2_vals)),
        T_grid,
        method='linear',
        bounds_error=False,
        fill_value=np.nan,
    )

    def interpolator(instellation, P_H2O, P_CO2):
        inst = np.asarray(instellation, dtype=float)
        h2o  = np.log10(np.asarray(P_H2O, dtype=float))
        co2  = np.log10(np.asarray(P_CO2, dtype=float))
        inst, h2o, co2 = np.broadcast_arrays(inst, h2o, co2)
        pts = np.stack([inst.ravel(), h2o.ravel(), co2.ravel()], axis=-1)
        result = _interp(pts).reshape(inst.shape)
        return float(result) if result.ndim == 0 else result

    return interpolator

if not os.path.exists(runs_data_path):
    parameter_sweep(n_workers=16)

basic_interpolator = make_basic_interpolator()

def solve_P_H2O(relative_humidity=0.8):
    from scipy.optimize import brentq

    instellation_range = np.linspace(0.1, 1.4, num=14)
    P_CO2_range = np.logspace(-6, 1, num=8)

    def target_function(instellation, P_CO2, P_H2O):
        T = basic_interpolator(instellation, P_H2O, P_CO2)
        # august_roche_magnus_formula returns Pa; convert to bar
        P_H2O_new = august_roche_magnus_formula(T) * relative_humidity / 1e5
        return P_H2O - P_H2O_new

    scan_pts = np.logspace(-6, 1, 60)

    P_H2O_eq = np.full((len(instellation_range), len(P_CO2_range)), np.nan)

    for i, inst in enumerate(instellation_range):
        for j, p_co2 in enumerate(P_CO2_range):
            fvals = np.array([target_function(inst, p_co2, p) for p in scan_pts])
            valid = ~np.isnan(fvals)
            if not np.any(valid):
                continue
            vpts, vf = scan_pts[valid], fvals[valid]
            sign_changes = np.where(np.diff(np.sign(vf)))[0]
            if len(sign_changes) == 0:
                continue
            k = sign_changes[0]  # lowest-P_H2O root
            try:
                P_H2O_eq[i, j] = brentq(
                    lambda p: target_function(inst, p_co2, p),
                    vpts[k], vpts[k + 1],
                )
            except ValueError:
                pass

    from scipy.interpolate import LinearNDInterpolator

    valid = ~np.isnan(P_H2O_eq)
    inst_grid, co2_grid = np.meshgrid(instellation_range, np.log10(P_CO2_range), indexing='ij')
    _interp = LinearNDInterpolator(
        np.column_stack([inst_grid[valid], co2_grid[valid]]),
        P_H2O_eq[valid],
        fill_value=np.nan,
    )

    def interpolator(instellation, P_CO2):
        inst = np.asarray(instellation, dtype=float)
        co2  = np.log10(np.asarray(P_CO2, dtype=float))
        inst, co2 = np.broadcast_arrays(inst, co2)
        pts = np.stack([inst.ravel(), co2.ravel()], axis=-1)
        result = _interp(pts).reshape(inst.shape)
        return float(result) if result.ndim == 0 else result

    # contour plot of T_surface(instellation, P_CO2) at equilibrium P_H2O
    inst_plot = np.linspace(instellation_range[0], instellation_range[-1], 200)
    co2_plot  = np.logspace(np.log10(P_CO2_range[0]), np.log10(P_CO2_range[-1]), 200)
    inst_mesh, co2_mesh = np.meshgrid(inst_plot, co2_plot)

    p_h2o_mesh = interpolator(inst_mesh, co2_mesh)
    T_mesh = basic_interpolator(inst_mesh, p_h2o_mesh, co2_mesh)

    fig, ax = plt.subplots(figsize=(8, 5))
    cf = ax.contourf(inst_mesh, co2_mesh, T_mesh, levels=20, cmap='RdYlBu_r')
    cs = ax.contour(inst_mesh, co2_mesh, T_mesh, levels=20, colors='k', linewidths=0.5, alpha=0.4)
    ax.clabel(cs, fmt='%d K', fontsize=7)
    plt.colorbar(cf, ax=ax, label='Surface Temperature (K)')
    ax.set_yscale('log')
    ax.set_xlabel('Instellation (S/S☉)')
    ax.set_ylabel('P$_{CO_2}$ (bar)')
    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(runs_data_path), 'T_surface_contour.png'), dpi=150)
    plt.show()

    return interpolator


climate_interpolator = solve_P_H2O()

