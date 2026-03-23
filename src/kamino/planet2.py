import numpy as np
import numpy.typing as npt
from scipy.optimize import newton, bisect
from tqdm import tqdm
import matplotlib.pyplot as plt
import csv
import os

from kamino.constants import *
from kamino.kamino_chem.ocean_chemistry import *
from kamino.speedy_climate.clima_interpolator import get_T_surface
from kamino.utils import *

T_min = 275
T_max = 350

tau_prec = 1e4 * YR
tau_atm = 1000 * YR

class Planet:

    def __init__(self, mass, radius, surface_pressure, instellation, tectonics, ocean_depth):

        self.mass = mass
        self.radius = radius
        self.gravity = (G * self.mass) / (self.radius ** 2)
        self.surface_area = 4 * np.pi * self.radius ** 2
        self.P_surface = surface_pressure

        self.ocean_depth = ocean_depth
        self.ocean_water_mass = self.ocean_depth * self.surface_area * 1000

        self.crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * tectonics
        self.hydrothermal_flux = EARTH_HYDROTHERMAL_FLUX_PER_AREA * tectonics
        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[1] = (EARTH_OUTGASSING / YR) * self.surface_area * tectonics

        self.alpha = 2

        self.tidally_locked = False

        land_fraction = 0
        ocean_albedo = 0.3
        land_albedo = 0.3

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

    def solve_climate_from_chemistry(self, b_ocean: npt.NDArray[np.float64], T_init: float=288) -> tuple[float, float]:

        def T_s_residual(T_guess):
            T_guess = np.clip(T_guess, T_min, 500)
            P_CO2 = get_P_CO2(self.P_surface, T_guess, b_ocean)
            T_calc = get_T_surface(self.instellation, P_CO2, self.albedo, tidally_locked=self.tidally_locked)
            return T_guess - T_calc
        
        try:
            T_s = newton(T_s_residual, T_init, maxiter=50, tol=1e-4)
        except (RuntimeError, ValueError):
            try:
                T_s = bisect(T_s_residual, T_min, T_max) 
            except ValueError:
                T_s = T_max if T_s_residual(T_max) < 0 else T_min
        
        P_CO2 = get_P_CO2(self.P_surface, T_s, b_ocean) # type: ignore

        assert ~np.isnan(T_s)

        return float(T_s), P_CO2 # type: ignore
    
    def get_seafloor_properties(self, T_s: float, P_CO2: float) -> tuple[float, float, float]:

        T_seafloor = 1.02 * T_s - 16.7

        P_H2O = august_roche_magnus_formula(T_s) * 0.5
        P_pore = (self.P_surface + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth

        T_seafloor = smooth_max(T_seafloor, 274)
        T_pore = T_seafloor + 9

        return T_seafloor, T_pore, P_pore
    
    def _compute_fluxes_and_derivatives(self, Y, dt):

        T = Y[0]
        P_CO2 = Y[1]
        b_ocean = smooth_max(Y[2:], np.zeros_like(Y[2:]))

        T_new, P_CO2_new = self.solve_climate_from_chemistry(b_ocean)
        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T_new, P_CO2_new)

        decay = np.exp(-dt / tau_atm)
        dT_dt = (T_new - T) * (1 - decay) / dt
        dP_CO2_dt = (P_CO2_new - P_CO2) * (1 - decay) / dt

        F_out = self.outgassing_flux / self.ocean_water_mass

        delta_b_prec = get_precipitation_flux(P_pore, T_seafloor, b_ocean, 
                                            precipitating_minerals=carbonate_minerals + secondary_sink_minerals)
        decay_prec = np.exp(-dt / tau_prec)
        F_prec = delta_b_prec * (1 - decay_prec) / dt

        _C_idx  = int(np.where(elements == 'C')[0][0])
        _Si_idx = int(np.where(elements == 'Si')[0][0])
        _ocean_water_per_area = self.ocean_depth * 1000.0  
        
        _F_carb = max(0.0, -F_prec[_C_idx])   
        _F_sil  = max(0.0, -F_prec[_Si_idx])  
        
        S_sed = (_F_carb * 0.100 / 2710.0 + _F_sil * 0.060 / 2650.0) * _ocean_water_per_area
        weathering_flux, w_diag = get_weathering_flux(P_pore, T_pore, P_CO2_new, b_ocean, 
                                                self.alpha, self.crust_production_rate, clog=False, sedimentation_rate=S_sed)
        F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

        F_net = F_out + F_diss + F_prec
        dYdt = np.zeros_like(Y)
        dYdt[0] = dT_dt
        dYdt[1] = dP_CO2_dt
        dYdt[2:] = F_net

        diagnostics = {
            'T_new': T_new, 'P_CO2_new': P_CO2_new,
            'T_seafloor': T_seafloor, 'T_pore': T_pore, 'P_pore': P_pore,
            'F_out': F_out, 'F_diss': F_diss, 'F_prec': F_prec,
            'delta_b_prec': delta_b_prec, 'F_net': F_net,
            'Da': w_diag['Da'],
            'A_reactive': w_diag['A_reactive'],
            'supply_efficiency': w_diag['supply_efficiency']
        }

        return dYdt, diagnostics
    
    def time_evolve_to_steady_state(self, dt=20000*YR, tol=1e-3, check_every=100, max_steps=100000, dt_max=200000*YR, dt_min=100*YR, debug_log='convergence_debug.csv'):

        Y = np.zeros(2 + elements.shape[0])
        Y[0] = 288
        Y[1] = 280e-6 * self.P_surface

        t_current = 0

        # Setup lists before the loop
        t = []
        P_CO2 = []
        T = []
        concentrations = [] # <--- Replaces Alk = []

        Y_prev_check = Y.copy()
        prev_rel_change = np.inf

        # Stagnation / runaway detection: track rel_change over a rolling window.
        # If it improves by less than 1% over the window, we are stuck.
        stagnation_window = 40   # number of checks
        rel_change_history: list[float] = []

        # Set up CSV debug log
        _csv_fields = (
            ['t_Myr', 'dt_kyr', 'T', 'P_CO2', 'rel_change', 'slowest', 'Da', 'A_reactive', 'supply_efficiency']
            + [f'Y_{e}' for e in elements]
            + [f'F_out_{e}' for e in elements]
            + [f'F_diss_{e}' for e in elements]
            + [f'F_prec_{e}' for e in elements]
            + [f'F_net_{e}' for e in elements]
        )
        _csv_file = open(debug_log, 'w', newline='')
        _csv_writer = csv.DictWriter(_csv_file, fieldnames=_csv_fields)
        _csv_writer.writeheader()

        pbar = tqdm(total=max_steps, unit='step')

        for i in range(max_steps):

            dY, fb = self._compute_fluxes_and_derivatives(Y, dt)

            if not np.all(np.isfinite(dY)):
                print(f"[time_evolve DEBUG] Non-finite dY at step {i}, t={t_current/YR:.2f} yr")
                print(f"  Y (before step): T={Y[0]:.4f}, P_CO2={Y[1]:.4e}, b_ocean={Y[2:]}")
                print(f"  dY: {dY}")
                break

            Y += dY * dt
            Y[2:] = np.maximum(Y[2:], 0.0)  # concentrations cannot be negative
            t_current += dt

            t.append(t_current)
            T.append(Y[0])
            P_CO2.append(Y[1])
            concentrations.append(Y[2:].copy())

            pbar.update(1)

            if i > 0 and i % check_every == 0:
                # Flux-based convergence: |F_net| / (|F_out| + |F_diss| + |F_prec|)
                # This directly measures flux balance rather than state drift, so
                # numerically oscillating-but-balanced elements (Mg, Ca) don't
                # block convergence.  T and P_CO2 use state-based check as before.
                total_flux = np.abs(fb['F_out']) + np.abs(fb['F_diss']) + np.abs(fb['F_prec'])
                safe_total = np.where(total_flux > 1e-50, total_flux, 1.0)
                flux_imbalance = np.where(total_flux > 1e-50, np.abs(fb['F_net']) / safe_total, 0.0)

                # Absolute flux floor: if |F_net| is too small to change the
                # concentration by atol over the next check interval, treat as
                # converged.  This prevents trace elements (e.g. Al at sub-pM
                # levels) from blocking convergence when their absolute flux is
                # physically negligible.
                element_atol = np.full(len(elements), 1e-9)  # mol/kgw
                flux_threshold = element_atol / (check_every * dt)   # mol/kgw/s
                flux_imbalance = np.where(np.abs(fb['F_net']) < flux_threshold, 0.0, flux_imbalance)
                
                atol = np.concatenate([[0.0, 0.0], np.full(len(elements), 1e-9)])
                Y_scale = np.maximum(np.abs(Y_prev_check), atol)
                Y_scale = np.maximum(Y_scale, 1e-30)
                state_changes = np.abs(Y[:2] - Y_prev_check[:2]) / Y_scale[:2]

                all_metrics = np.concatenate([state_changes, flux_imbalance])
                Y_labels = ['T', 'P_CO2'] + list(elements)
                rel_change = np.max(all_metrics)
                slowest_idx = int(np.argmax(all_metrics))
                slowest_label = Y_labels[slowest_idx]

                pbar.set_postfix({
                    't(Myr)': f'{t_current / YR / 1e6:.4f}',
                    'T(K)': f'{Y[0]:.0f}',
                    'rel_chg': f'{rel_change:.3%}',
                    'slowest': slowest_label,
                    'dt(kyr)': f'{dt / YR / 1e3:.2f}'
                })

                # Detailed flux breakdown for convergence diagnostics (fb already computed above)
                row = {
                    't_Myr': f'{t_current/YR/1e6:.4f}',
                    'dt_kyr': f'{dt/YR/1e3:.2f}',
                    'T': f'{Y[0]:.4f}',
                    'P_CO2': f'{Y[1]:.6e}',
                    'rel_change': f'{rel_change:.6e}',
                    'slowest': slowest_label,
                    'Da': f'{fb["Da"]:.6e}',
                    'A_reactive': f'{fb["A_reactive"]:.6e}',
                    'supply_efficiency': f'{fb["supply_efficiency"]:.4%}' # Formatted as a percentage!
                }
                
                for j, e in enumerate(elements):
                    row[f'Y_{e}']     = f'{Y[2+j]:.6e}'
                    row[f'F_out_{e}'] = f'{fb["F_out"][j]:.6e}'
                    row[f'F_diss_{e}']= f'{fb["F_diss"][j]:.6e}'
                    row[f'F_prec_{e}']= f'{fb["F_prec"][j]:.6e}'
                    row[f'F_net_{e}'] = f'{fb["F_net"][j]:.6e}'
                _csv_writer.writerow(row)
                _csv_file.flush()

                if rel_change < tol:
                    print(f"\nConverged at t = {t_current / YR / 1e6:.4f} Myr  "
                          f"(rel_change = {rel_change:.2e} < tol = {tol:.2e})")
                    break

                # Stagnation / runaway detection.
                rel_change_history.append(rel_change)
                if len(rel_change_history) > stagnation_window:
                    rel_change_history.pop(0)

                if len(rel_change_history) == stagnation_window:
                    oldest = rel_change_history[0]
                    improvement = (oldest - rel_change) / oldest if oldest > 0 else 0
                    if improvement < 0.01:
                        # Diagnose the cause before breaking.
                        element_atol = np.full(len(elements), 1e-9)
                        flux_threshold = element_atol / (check_every * dt)
                        meaningful = np.abs(fb['F_net']) >= flux_threshold

                        if abs(Y[0] - T_min) < 0.5:
                            cause = (f"T pinned at T_min = {T_min} K — planet is too cold "
                                     f"for the carbonate thermostat to operate (CO2 runaway / acid ocean).")
                        elif abs(Y[0] - T_max) < 0.5:
                            cause = (f"T pinned at T_max = {T_max} K — climate model ceiling "
                                     f"reached; consider extending the interpolation grid.")
                        else:
                            stuck = meaningful & (flux_imbalance > 0.5)
                            stuck_names = [str(elements[j]) for j in range(len(elements)) if stuck[j]]
                            cause = (f"elements {stuck_names} have persistent flux imbalance "
                                     f"> 50% with no sign of improvement — likely no accessible "
                                     f"steady state in the current mineral set.")

                        print(f"\nStagnation detected: rel_change = {rel_change:.3%} "
                              f"improved by only {improvement:.1%} over "
                              f"{stagnation_window} checks ({stagnation_window * check_every * dt / YR / 1e6:.1f} Myr).")
                        print(f"Cause: {cause}")
                        break

                # Adapt dt: grow when converging smoothly, shrink only on a
                # genuine spike (rel_change increased from previous check).
                if rel_change > prev_rel_change * 2.0 and rel_change > 0.05:
                    dt = max(dt * 0.5, dt_min)
                elif rel_change < prev_rel_change * 0.999:
                    # Monotonically improving: allow dt to grow toward dt_max.
                    dt = min(dt * 1.2, dt_max)

                prev_rel_change = rel_change

                Y_prev_check = Y.copy()

        else:
            print(f"\nMax steps ({max_steps}) reached without convergence at "
                  f"t = {t_current / YR / 1e6:.4f} Myr")

        pbar.close()
        _csv_file.close()
        print(f"\nFlux debug log written to: {os.path.abspath(debug_log)}")

        t = np.array(t)
        T = np.array(T)
        P_CO2 = np.array(P_CO2)
        concentrations = np.array(concentrations)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        for j, elem in enumerate(elements):
            ax1.plot(t / YR, concentrations[:, j], label=elem)
        
        ax1.set_yscale('log')
        ax1.set_ylabel('Concentration (mol/kgw)')
        ax1.set_title('Time Evolution of Ocean Chemistry')
        # Put the legend outside the plot so it doesn't cover the data
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, which="both", ls="--", alpha=0.5)

        # --- Bottom Plot: P_CO2 and T ---
        # Left Y-axis for P_CO2 (Log scale)
        color_co2 = 'tab:blue'
        ax2.set_xlabel('Time (years)')
        ax2.set_ylabel('P_CO2 / P_surface', color=color_co2)
        ax2.plot(t / YR, P_CO2 / self.P_surface, color=color_co2, label='P_CO2')
        ax2.tick_params(axis='y', labelcolor=color_co2)
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.5)

        # Right Y-axis for Temperature (Linear scale)
        ax3 = ax2.twinx()  
        color_t = 'tab:red'
        ax3.set_ylabel('Temperature (K)', color=color_t)
        ax3.plot(t / YR, T, color=color_t, label='Temperature')
        ax3.tick_params(axis='y', labelcolor=color_t)

        # Adjust layout so the external legend and right y-label aren't cut off
        fig.tight_layout()  
        plt.show()

        # plt.plot(t / YR, T)
        # plt.show()

        # plt.plot(t / YR, P_CO2 / self.P_surface)
        # plt.yscale('log')
        # plt.show()

        # plt.plot(t / YR, Alk)
        # plt.yscale('log')
        # plt.show()
        