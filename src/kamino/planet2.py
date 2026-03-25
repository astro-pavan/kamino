import numpy as np
import numpy.typing as npt
from scipy.optimize import newton, bisect
from tqdm import tqdm
import matplotlib.pyplot as plt
import csv
import json
import os

from kamino.constants import *
from kamino.kamino_chem.ocean_chemistry import *
from kamino.speedy_climate.clima_interpolator import get_T_surface
from kamino.utils import *

T_min = 200   # climate table lower bound ~184 K; 200 K is safely inside
T_max = 350

tau_prec = 1e4 * YR
tau_atm = 1000 * YR

class Planet:

    def __init__(self, mass: float, radius: float, background_pressure: float, instellation : float, tectonics: float, ocean_depth: float, land_fraction: float=0.0, name: str='planet', cl_chemistry: bool = False, hcl_co2_ratio: float = 0.02):        

        self.name = name
        self.mass = mass
        self.radius = radius
        self.gravity = (G * self.mass) / (self.radius ** 2)
        self.surface_area = 4 * np.pi * self.radius ** 2
        self.P_background = background_pressure

        self.ocean_depth = ocean_depth
        self.ocean_water_mass = self.ocean_depth * self.surface_area * 1000

        self.crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * tectonics
        self.hydrothermal_flux = EARTH_HYDROTHERMAL_FLUX_PER_AREA * tectonics
        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[1] = (EARTH_OUTGASSING / YR) * self.surface_area * tectonics

        self.cl_chemistry = cl_chemistry
        self.hcl_co2_ratio = hcl_co2_ratio
        if cl_chemistry:
            F_HCl = (EARTH_OUTGASSING / YR) * self.surface_area * tectonics * hcl_co2_ratio
            self.outgassing_flux[cl_idx]  += F_HCl   # HCl → Cl⁻ added to ocean
            self.outgassing_flux[alk_idx] -= F_HCl   # H⁺ released consumes alkalinity

        self.alpha = 2

        self.tidally_locked = False

        self._input_instellation = instellation  # in units of solar constant
        self._input_tectonics = tectonics

        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area

        ocean_albedo = 0.3
        land_albedo = 0.3

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

    def solve_climate_from_chemistry(self, b_ocean: npt.NDArray[np.float64], T_init: float=288, P_CO2_est: float=0.0) -> tuple[float, float]:

        def T_s_residual(T_guess):
            if not np.isfinite(T_guess):
                raise ValueError("T_guess is not finite")
            T_guess = np.clip(T_guess, T_min, 500)
            P_H2O = august_roche_magnus_formula(T_guess) * 0.5
            P_atm = self.P_background + P_CO2_est + P_H2O
            P_CO2 = get_P_CO2(P_atm, T_guess, b_ocean)
            T_calc = get_T_surface(self.instellation, P_CO2, self.albedo, tidally_locked=self.tidally_locked)
            return T_guess - T_calc

        try:
            T_s = newton(T_s_residual, T_init, maxiter=50, tol=1e-4)
        except (RuntimeError, ValueError, AssertionError):
            try:
                T_s = bisect(T_s_residual, T_min, T_max)
            except (ValueError, AssertionError):
                try:
                    T_s = T_max if T_s_residual(T_max) < 0 else T_min
                except (AssertionError, ValueError):
                    T_s = T_max  # climate table out of range even at T_max → runaway

        P_H2O = august_roche_magnus_formula(T_s) * 0.5  # type: ignore
        P_atm = self.P_background + P_CO2_est + P_H2O
        P_CO2 = get_P_CO2(P_atm, T_s, b_ocean) # type: ignore

        assert ~np.isnan(T_s)

        return float(T_s), P_CO2 # type: ignore
    
    def get_seafloor_properties(self, T_s: float, P_CO2: float) -> tuple[float, float, float]:

        T_seafloor = 1.02 * T_s - 16.7

        P_H2O = august_roche_magnus_formula(T_s) * 0.5
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth

        T_seafloor = smooth_max(T_seafloor, 274)
        T_pore = T_seafloor + 9

        return T_seafloor, T_pore, P_pore
    
    def _compute_fluxes_and_derivatives(self, Y, dt):

        T = Y[0]
        P_CO2 = Y[1]
        b_ocean = smooth_max(Y[2:], np.zeros_like(Y[2:]))

        # --- Atmospheric and Seafloor Properties ---

        T_new, P_CO2_new = self.solve_climate_from_chemistry(b_ocean, P_CO2_est=P_CO2)
        P_H2O = august_roche_magnus_formula(T_new) * 0.5
        P_surf = self.P_background + P_CO2_new + P_H2O
        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T_new, P_CO2_new)

        decay = np.exp(-dt / tau_atm)
        dT_dt = (T_new - T) * (1 - decay) / dt
        dP_CO2_dt = (P_CO2_new - P_CO2) * (1 - decay) / dt

        # --- Outgassing ---

        F_out = self.outgassing_flux / self.ocean_water_mass

        # --- Precipitation ---

        # Deep burial
        delta_b_prec, pH = get_precipitation_flux(P_pore, T_seafloor, b_ocean, precipitating_minerals=carbonate_minerals + secondary_sink_minerals)

        # Shallow burial
        if self.land_fraction > 0:
            delta_b_shallow, _ = get_precipitation_flux(P_surf, T_new, b_ocean, precipitating_minerals=carbonate_minerals)
            delta_b_prec += self.land_fraction * delta_b_shallow

        decay_prec = np.exp(-dt / tau_prec)
        F_prec = delta_b_prec * (1 - decay_prec) / dt

        # --- Weathering ---

        _ocean_water_per_area = self.ocean_depth * 1000.0  
        
        _F_carb = max(0.0, -F_prec[c_idx])   
        _F_sil  = max(0.0, -F_prec[si_idx])  
        
        S_sed = (_F_carb * 0.100 / 2710.0 + _F_sil * 0.060 / 2650.0) * _ocean_water_per_area
        weathering_flux, w_diag = get_weathering_flux(P_pore, T_pore, P_CO2_new, b_ocean, self.alpha, self.crust_production_rate, clog=False, sedimentation_rate=S_sed)
        F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

        if self.land_area > 0:
            F_cont_per_area = get_continental_weathering_flux(T_new, P_CO2_new)
            F_cont = (F_cont_per_area * self.land_area) / self.ocean_water_mass
        else:
            F_cont = np.zeros(len(elements))

        # --- Total Flux ---

        F_net = F_out + F_diss + F_prec + F_cont
        dYdt = np.zeros_like(Y)
        dYdt[0] = dT_dt
        dYdt[1] = dP_CO2_dt
        dYdt[2:] = F_net

        diagnostics = {
            'T_new': T_new, 'P_CO2_new': P_CO2_new, 'pH': pH,
            'T_seafloor': T_seafloor, 'T_pore': T_pore, 'P_pore': P_pore,
            'F_out': F_out, 'F_diss': F_diss, 'F_prec': F_prec, 'F_cont': F_cont,
            'delta_b_prec': delta_b_prec, 'F_net': F_net,
            'Da': w_diag['Da'],
            'A_reactive': w_diag['A_reactive'],
            'supply_efficiency': w_diag['supply_efficiency']
        }

        return dYdt, diagnostics
    
    def time_evolve_to_steady_state(self, dt=20000*YR, tol=5e-3, check_every=100, max_steps=30000, dt_max=200000*YR, dt_min=100*YR, output_dir='.', b_ocean_init: dict | None = None):

        output_csv = os.path.join(output_dir, f'{self.name}_timeseries.csv')
        summary_file = os.path.join(output_dir, f'{self.name}_summary.json')

        Y = np.zeros(2 + elements.shape[0])
        Y[0] = 300
        Y[1] = 280e-6 * self.P_background

        _trace_defaults = {
            'Alkalinity': 1e-3,   # mol eq/kgw  — sets a valid pH baseline
            'C':          1e-3,   # mol/kgw
            'Ca':         1e-8,
            'Mg':         1e-8,
            'Si':         1e-8,
            'Al':         1e-9,
            'Fe':         1e-9,
        }

        for elem, val in _trace_defaults.items():
            idx = int(np.where(elements == elem)[0][0])
            Y[2 + idx] = val

        if b_ocean_init is not None:
            for elem, val in b_ocean_init.items():
                idx = int(np.where(elements == elem)[0][0])
                Y[2 + idx] = val

        t_current = 0
        convergence_status = 'max_steps'
        convergence_reason = f'Max steps ({max_steps}) reached without convergence'
        last_diagnostics = None

        # Setup lists before the loop
        t = []
        P_CO2 = []
        T = []
        concentrations = []

        Y_prev_check = Y.copy()
        prev_rel_change = np.inf

        # Stagnation / runaway detection: track rel_change over a rolling window.
        # If it improves by less than 1% over the window, we are stuck.
        stagnation_window = 40   # number of checks
        rel_change_history: list[float] = []

        # Set up timeseries CSV
        _csv_fields = (
            ['t_Myr', 'dt_kyr', 'T', 'P_CO2', 'pH', 'rel_change', 'slowest', 'Da', 'A_reactive', 'supply_efficiency']
            + [f'Y_{e}' for e in elements]
            + [f'F_out_{e}' for e in elements]
            + [f'F_diss_{e}' for e in elements]
            + [f'F_prec_{e}' for e in elements]
            + [f'F_cont_{e}' for e in elements]
            + [f'F_net_{e}' for e in elements]
        )
        _csv_file = open(output_csv, 'w', newline='')
        _csv_writer = csv.DictWriter(_csv_file, fieldnames=_csv_fields)
        _csv_writer.writeheader()

        for i in range(max_steps):

            dY, diagnostics = self._compute_fluxes_and_derivatives(Y, dt)
            last_diagnostics = diagnostics

            if not np.all(np.isfinite(dY)):
                print(f"[time_evolve DEBUG] Non-finite dY at step {i}, t={t_current/YR:.2f} yr")
                print(f"  Y (before step): T={Y[0]:.4f}, P_CO2={Y[1]:.4e}, b_ocean={Y[2:]}")
                print(f"  dY: {dY}")
                convergence_status = 'nonfinite'
                convergence_reason = f'Non-finite dY at step {i}, t={t_current/YR:.2f} yr'
                break

            Y += dY * dt
            Y[2:] = np.maximum(Y[2:], 0.0)
            t_current += dt

            t.append(t_current)
            T.append(Y[0])
            P_CO2.append(Y[1])
            concentrations.append(Y[2:].copy())

            if i > 0 and i % check_every == 0:
                # Flux-based convergence: |F_net| / (|F_out| + |F_diss| + |F_prec|)
                # This directly measures flux balance rather than state drift, so
                # numerically oscillating-but-balanced elements (Mg, Ca) don't
                # block convergence.  T and P_CO2 use state-based check as before.
                total_flux = np.abs(diagnostics['F_out']) + np.abs(diagnostics['F_diss']) + np.abs(diagnostics['F_prec'])
                safe_total = np.where(total_flux > 1e-50, total_flux, 1.0)
                flux_imbalance = np.where(total_flux > 1e-50, np.abs(diagnostics['F_net']) / safe_total, 0.0)

                # Absolute flux floor: if |F_net| is too small to change the
                # concentration by atol over the next check interval, treat as
                # converged.  This prevents trace elements (e.g. Al at sub-pM
                # levels) from blocking convergence when their absolute flux is
                # physically negligible.
                element_atol = np.full(len(elements), 1e-9)  # mol/kgw
                flux_threshold = element_atol / (check_every * dt)   # mol/kgw/s
                flux_imbalance = np.where(np.abs(diagnostics['F_net']) < flux_threshold, 0.0, flux_imbalance)
                
                atol = np.concatenate([[0.0, 0.0], np.full(len(elements), 1e-9)])
                Y_scale = np.maximum(np.abs(Y_prev_check), atol)
                Y_scale = np.maximum(Y_scale, 1e-30)
                state_changes = np.abs(Y[:2] - Y_prev_check[:2]) / Y_scale[:2]

                all_metrics = np.concatenate([state_changes, flux_imbalance])
                Y_labels = ['T', 'P_CO2'] + list(elements)
                rel_change = np.max(all_metrics)
                slowest_idx = int(np.argmax(all_metrics))
                slowest_label = Y_labels[slowest_idx]

                # Detailed flux breakdown for convergence diagnostics (fb already computed above)
                row = {
                    't_Myr': f'{t_current/YR/1e6:.4f}',
                    'dt_kyr': f'{dt/YR/1e3:.2f}',
                    'T': f'{Y[0]:.4f}',
                    'P_CO2': f'{Y[1]:.6e}',
                    'pH': f"{diagnostics['pH']:.4f}",
                    'rel_change': f'{rel_change:.6e}',
                    'slowest': slowest_label,
                    'Da': f'{diagnostics["Da"]:.6e}',
                    'A_reactive': f'{diagnostics["A_reactive"]:.6e}',
                    'supply_efficiency': f'{diagnostics["supply_efficiency"]:.4%}' # Formatted as a percentage!
                }
                
                for j, e in enumerate(elements):
                    row[f'Y_{e}']      = f'{Y[2+j]:.6e}'
                    row[f'F_out_{e}']  = f'{diagnostics["F_out"][j]:.6e}'
                    row[f'F_diss_{e}'] = f'{diagnostics["F_diss"][j]:.6e}'
                    row[f'F_prec_{e}'] = f'{diagnostics["F_prec"][j]:.6e}'
                    row[f'F_cont_{e}'] = f'{diagnostics["F_cont"][j]:.6e}'
                    row[f'F_net_{e}']  = f'{diagnostics["F_net"][j]:.6e}'
                _csv_writer.writerow(row)
                _csv_file.flush()

                if rel_change < tol:
                    print(f"\nConverged at t = {t_current / YR / 1e6:.4f} Myr  "
                          f"(rel_change = {rel_change:.2e} < tol = {tol:.2e})")
                    convergence_status = 'converged'
                    convergence_reason = f'rel_change = {rel_change:.2e} < tol = {tol:.2e}'
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
                        meaningful = np.abs(diagnostics['F_net']) >= flux_threshold

                        if abs(Y[0] - T_min) < 0.5:
                            cause = (f"T pinned at T_min = {T_min} K — planet has entered a "
                                     f"snowball state (runaway glaciation).")
                            print(f"\nSnowball planet detected: T = {Y[0]:.1f} K pinned at T_min = {T_min} K. "
                                  f"rel_change = {rel_change:.3%} over "
                                  f"{stagnation_window} checks ({stagnation_window * check_every * dt / YR / 1e6:.1f} Myr).")
                            print(f"Cause: {cause}")
                            convergence_status = 'snowball'
                            convergence_reason = cause
                            break
                        elif abs(Y[0] - T_max) < 0.5:
                            cause = (f"T pinned at T_max = {T_max} K — climate model ceiling "
                                     f"reached; planet likely in runaway greenhouse state.")
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
                        convergence_status = 'stagnated'
                        convergence_reason = cause
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
            
        _csv_file.close()
        print(f"\nTimeseries written to: {os.path.abspath(output_csv)}")

        # --- Write summary JSON ---
        final_ocean_chemistry = {str(e): float(Y[2 + j]) for j, e in enumerate(elements)}
        summary = {
            'input_parameters': {
                'name': self.name,
                'mass_kg': float(self.mass),
                'radius_m': float(self.radius),
                'surface_pressure_Pa': float(self.P_background),
                'instellation_solar': float(self._input_instellation),
                'tectonics': float(self._input_tectonics),
                'land_fraction': float(self.land_fraction),
                'ocean_depth_m': float(self.ocean_depth),
                'alpha': float(self.alpha),
                'tidally_locked': bool(self.tidally_locked),
                'cl_chemistry': bool(self.cl_chemistry),
                'hcl_co2_ratio': float(self.hcl_co2_ratio),
            },
            'convergence': {
                'status': convergence_status,
                'reason': convergence_reason,
                'time_to_converge_Myr': float(t_current / YR / 1e6),
                'snowball': convergence_status == 'snowball',
                'runaway': convergence_status == 'stagnated' and abs(Y[0] - T_max) < 0.5,
                'stagnated': convergence_status == 'stagnated',
            },
            'final_state': {
                'T_K': float(Y[0]),
                'P_CO2_Pa': float(Y[1]),
                'pH': float(last_diagnostics['pH']) if last_diagnostics is not None else None,
                'ocean_chemistry_mol_kgw': final_ocean_chemistry,
            },
            'final_diagnostics': {
                'Da': float(last_diagnostics['Da']) if last_diagnostics is not None else None,
                'A_reactive_m2': float(last_diagnostics['A_reactive']) if last_diagnostics is not None else None,
                'supply_efficiency': float(last_diagnostics['supply_efficiency']) if last_diagnostics is not None else None,
            },
        }
        with open(summary_file, 'w') as _sf:
            json.dump(summary, _sf, indent=2)
        print(f"Summary written to:    {os.path.abspath(summary_file)}")

        # --- Plot results ---
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
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, which="both", ls="--", alpha=0.5)

        color_co2 = 'tab:blue'
        ax2.set_xlabel('Time (years)')
        ax2.set_ylabel('$P_{CO2}$ (bar)', color=color_co2)
        ax2.plot(t / YR, P_CO2 / 1e5, color=color_co2, label='P_CO2')
        ax2.tick_params(axis='y', labelcolor=color_co2)
        ax2.set_yscale('log')
        ax2.grid(True, alpha=0.5)

        ax3 = ax2.twinx()  
        color_t = 'tab:red'
        ax3.set_ylabel('Temperature (K)', color=color_t)
        ax3.plot(t / YR, T, color=color_t, label='Temperature')
        ax3.tick_params(axis='y', labelcolor=color_t)

        fig.tight_layout()  
        plt.savefig(os.path.join(output_dir, f'{self.name}_plot.pdf'))