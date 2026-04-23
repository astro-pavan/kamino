import numpy as np
import numpy.typing as npt
from scipy.optimize import newton, bisect, brentq
from scipy.integrate import solve_ivp
from tqdm import tqdm
import matplotlib.pyplot as plt
import csv
import json
import os

from kamino.constants import *
from kamino.chemistry.ocean_chemistry import *
from kamino.climate.clima_interpolator import get_T_surface
from kamino.utils import *

output_path = os.path.join(os.path.dirname(__file__), '../../output/')

T_min = 200   # climate table lower bound ~184 K; 200 K is safely inside
T_max = 350

tau_prec = 1e4 * YR
max_precipitation_frac = 0.0002
tau_atm = 1e4 * YR
drain_fraction = 0.01

t_max = 1e9 * YR

class Planet:

    def __init__(self, mass: float, radius: float, background_pressure: float, instellation : float, tectonics: float, ocean_depth: float, land_fraction: float=0.0, name: str='planet', use_precipitation_surrogate: bool=False, cl_chemistry: bool = False, hcl_co2_ratio: float = 0.02):        

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

        # self.cl_chemistry = cl_chemistry
        # self.hcl_co2_ratio = hcl_co2_ratio
        # if cl_chemistry:
        #     F_HCl = (EARTH_OUTGASSING / YR) * self.surface_area * tectonics * hcl_co2_ratio
        #     self.outgassing_flux[cl_idx]  += F_HCl   # HCl → Cl⁻ added to ocean
        #     self.outgassing_flux[alk_idx] -= F_HCl   # H⁺ released consumes alkalinity

        self.alpha = 2

        self.tidally_locked = False

        self._input_instellation = instellation  # in units of solar constant
        self._input_tectonics = tectonics

        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area

        ocean_albedo = 0.3
        land_albedo = 0.3

        self.relative_humidity = 0.5

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

        self.use_precipitation_surrogate = use_precipitation_surrogate

    def solve_climate_from_chemistry(self, b_ocean: npt.NDArray[np.float64], T_init: float=288, P_CO2_est: float=0.0) -> tuple[float, float]:

        # Anchor P_CO2 at T_init via PHREEQC, then build a linear P_CO2(T) model
        # using the actual PHREEQC slope.  The analytic thermo-ratio formula
        # over-estimates dP_CO2/dT by ~40%, which makes Newton overshoot and
        # land on spurious high-T roots (~350 K).
        P_atm_0 = self.P_background + P_CO2_est + august_roche_magnus_formula(T_init) * self.relative_humidity
        P_CO2_0 = get_P_CO2(P_atm_0, T_init, b_ocean)

        _dT_probe = 5.0
        P_atm_p = self.P_background + P_CO2_est + august_roche_magnus_formula(T_init + _dT_probe) * self.relative_humidity
        P_CO2_p = get_P_CO2(P_atm_p, T_init + _dT_probe, b_ocean)
        _slope = (P_CO2_p - P_CO2_0) / _dT_probe   # Pa / K  (real PHREEQC sensitivity)

        def T_s_residual(T_guess):
            if not np.isfinite(T_guess):
                raise ValueError("T_guess is not finite")
            T_guess = np.clip(T_guess, T_min, 500)
            P_CO2 = max(0.0, P_CO2_0 + _slope * (T_guess - T_init))
            T_calc = get_T_surface(self.instellation, P_CO2, self.albedo, tidally_locked=self.tidally_locked)
            return T_guess - T_calc

        # Search for a root near T_init in progressively wider windows.
        # This keeps the solution on the local climate branch (warm or cold)
        # rather than letting Newton jump to a distant spurious root.
        # The linear P_CO2 model has a cold root at T_min (where P_CO2→0),
        # and Newton can land there if the warm root has an unstable derivative.
        T_s = None
        for half_width in [8.0, 18.0, 35.0, 70.0]:
            lo = max(T_min, T_init - half_width)
            hi = min(T_max, T_init + half_width)
            try:
                r_lo = T_s_residual(lo)
                r_hi = T_s_residual(hi)
                if r_lo * r_hi <= 0:
                    T_s = float(brentq(T_s_residual, lo, hi, xtol=0.1, maxiter=50))
                    break
            except Exception:
                pass

        if T_s is None:
            # No local root: check the global range for a genuine phase transition
            try:
                r_min = T_s_residual(T_min)
                r_max = T_s_residual(T_max)
                if r_min * r_max <= 0:
                    T_s = float(brentq(T_s_residual, T_min, T_max, xtol=0.1))
                elif r_max < 0:
                    T_s = float(T_max)  # Runaway Greenhouse
                else:
                    T_s = float(T_min)  # Snowball
            except Exception:
                T_s = float(T_max)

        P_H2O = august_roche_magnus_formula(T_s) * self.relative_humidity  # type: ignore
        P_atm = self.P_background + P_CO2_est + P_H2O
        P_CO2 = get_P_CO2(P_atm, T_s, b_ocean) # type: ignore

        assert ~np.isnan(T_s)

        return float(T_s), P_CO2 # type: ignore
    
    def get_seafloor_properties(self, T_s: float, P_CO2: float) -> tuple[float, float, float]:

        T_seafloor = 1.02 * T_s - 16.7

        P_H2O = august_roche_magnus_formula(T_s) * self.relative_humidity
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth

        T_seafloor = smooth_max(T_seafloor, 274)
        T_pore = T_seafloor + 9

        return T_seafloor, T_pore, P_pore
    
    def dY_dt(self, t, Y):

        self._n_eval = getattr(self, '_n_eval', 0) + 1
        print(f't = {t/YR:.4e} yr  evals={self._n_eval}', end='\r')

        try:

            # 1e-20 floor prevents PHREEQC from receiving exactly-zero C when
            # Alk > 0, which is an inconsistent charge state it cannot converge on.
            b_ocean = np.maximum(Y, 1e-20)

            # --- Atmospheric and Seafloor Properties ---
            # T and P_CO2 are algebraic outputs (quasi-static equilibrium).
            # Cache is used only as a warm-start; stale cache from Jacobian
            # probes is fine since solve_climate_from_chemistry is robust.

            T, P_CO2 = self.solve_climate_from_chemistry(
                b_ocean, T_init=self._T_cache, P_CO2_est=self._P_CO2_cache
            )
            self._T_cache = T
            self._P_CO2_cache = P_CO2

            T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T, P_CO2)

            # --- Outgassing ---

            F_vol = self.outgassing_flux / self.ocean_water_mass

            # --- Precipitation ---

            F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=carbonate_minerals + secondary_sink_minerals, precipitation_timescale=tau_prec)

            # --- Weathering ---

            _ocean_water_per_area = self.ocean_depth * 1000.0
            _F_carb = max(0.0, -F_prec[c_idx])
            _F_sil  = max(0.0, -F_prec[si_idx])

            S_sed = (_F_carb * 0.100 / 2710.0 + _F_sil * 0.060 / 2650.0) * _ocean_water_per_area
            weathering_flux, w_diag = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, self.alpha, self.crust_production_rate, clog=False, sedimentation_rate=S_sed)
            F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

            if self.land_area > 0:
                F_cont_per_area = get_continental_weathering_flux(T, P_CO2)
                F_cont = (F_cont_per_area * self.land_area) / self.ocean_water_mass
            else:
                F_cont = np.zeros(len(elements))

            # --- Total Flux ---

            F_net = F_vol + F_diss + F_cont + F_prec

            # Prevent over-depletion: cap the removal flux so no species can be
            # depleted faster than its characteristic timescale (tau_prec).  This
            # keeps concentrations positive during collocation probes and prevents
            # PHREEQC from receiving Alk > 0 with C = 0 (charge-balance failure).
            _max_removal = b_ocean / tau_prec
            F_net = np.where(F_net < -_max_removal, -_max_removal, F_net)

        except ChemistryError:
            # PHREEQC failed on a transient collocation probe; return outgassing
            # only (always valid) so Radau gets a correct-order Jacobian and can
            # try a smaller step rather than spiralling on wrong cached data.
            F_net = self.outgassing_flux / self.ocean_water_mass
            print(f'\nChemistryError at t={t/YR:.4e} yr; Y={Y}')

        return F_net
    
    def time_evolve(self):

        # State vector is b_ocean only — T and P_CO2 are algebraic.
        Y0 = np.zeros(elements.shape[0])

        self._T_cache = 300.0
        self._P_CO2_cache = 280e-6 * self.P_background
        self._last_F_net = np.zeros(elements.shape[0])
        self._n_eval = 0

        P_CO2_max = 1e5
        T_runaway = 350
        T_snowball = 260

        def _get_climate(Y):
            b = np.maximum(Y, 0.0)
            return self.solve_climate_from_chemistry(b, T_init=self._T_cache, P_CO2_est=self._P_CO2_cache)

        def runaway_event(t, Y):
            T, _ = _get_climate(Y)
            return T - T_runaway

        runaway_event.direction = 1
        runaway_event.terminal = True

        def snowball_event(t, Y):
            T, P_CO2 = _get_climate(Y)
            return min(T_snowball - T, P_CO2 - P_CO2_max)

        snowball_event.direction = 1
        snowball_event.terminal = True

        sol = solve_ivp(
            self.dY_dt,
            (0, 1e9 * YR),
            Y0,
            method='Radau',
            events=[runaway_event, snowball_event],
            max_step=1e7 * YR,
            rtol=1e-3,
            atol=1e-10,
        )

        print()
        print(f'T = {self._T_cache:.2f} K  P_CO2 = {self._P_CO2_cache:.3e} Pa')
        print(sol.y[:, -1])
    
    # def time_evolve_to_steady_state(self, dt=1e5 * YR, tol=5e-3, check_every=100, max_steps=30000, dt_max=1e7 * YR, dt_min=100*YR, output_dir=output_path, b_ocean_init: dict | None = None, verbose: bool=False, write_csv: bool=False):

    #     output_csv = os.path.join(output_dir, f'{self.name}_timeseries.csv')
    #     summary_file = os.path.join(output_dir, f'{self.name}_summary.json')

    #     Y = np.zeros(2 + elements.shape[0])
    #     Y[0] = 300
    #     Y[1] = 280e-6 * self.P_background

    #     _trace_defaults = {
    #         'Alkalinity': 1e-9,   # mol eq/kgw  — sets a valid pH baseline
    #         'C':          1e-9,   # mol/kgw
    #         'Ca':         1e-9,
    #         'Mg':         1e-9,
    #         'Si':         1e-9,
    #         'Al':         1e-9,
    #         'Fe':         1e-9,
    #         # 'Cl':         1e-9
    #     }

    #     for elem, val in _trace_defaults.items():
    #         idx = int(np.where(elements == elem)[0][0])
    #         Y[2 + idx] = val

    #     if b_ocean_init is not None:
    #         for elem, val in b_ocean_init.items():
    #             idx = int(np.where(elements == elem)[0][0])
    #             Y[2 + idx] = val

    #     t_current = 0
    #     convergence_status = 'max_steps'
    #     convergence_reason = f'Max steps ({max_steps}) reached without convergence'
    #     last_diagnostics = None

    #     # Setup lists before the loop
    #     t = []
    #     P_CO2 = []
    #     T = []
    #     concentrations = []
    #     fluxes = []
    #     pH = []

    #     if write_csv:
    #         _csv_fields = (
    #             ['t_Myr', 'dt_kyr', 'T', 'P_CO2', 'pH', 'rel_change', 'slowest', 'Da', 'A_reactive', 'supply_efficiency']
    #             + [f'Y_{e}' for e in elements]
    #             + [f'F_vol_{e}' for e in elements]
    #             + [f'F_diss_{e}' for e in elements]
    #             + [f'F_prec_{e}' for e in elements]
    #             + [f'F_cont_{e}' for e in elements]
    #             + [f'F_net_{e}' for e in elements]
    #         )
    #         _csv_file = open(output_csv, 'w', newline='')
    #         _csv_writer = csv.DictWriter(_csv_file, fieldnames=_csv_fields)
    #         _csv_writer.writeheader()

    #     T_current, P_CO2_current = Y[0], Y[1]

    #     _cfl_frac  = 0.05   # max fractional change per step (CFL target)
    #     _dt_grow   = 2      # growth factor when step was safe

    #     for i in tqdm(range(max_steps)):

    #         # --- Forward step with adaptive dt ---

    #         dYdt, diagnostics = self._compute_fluxes_and_derivatives(Y, dt)
    #         last_diagnostics = diagnostics

    #         # CFL: reduce dt so no concentration changes by more than _cfl_frac
    #         safe_b = np.maximum(np.abs(Y[2:]), 1e-15)
    #         frac_chem = np.max(np.abs(dYdt[2:]) * dt / safe_b)
            
    #         # NEW: Don't let Temperature change by more than 0.5 K per step
    #         # dYdt[0] is your dT_dt
    #         frac_T = np.abs(dYdt[0] * dt) / 0.5 
            
    #         # Use whichever restriction is stricter
    #         frac = max(frac_chem, frac_T)

    #         if frac > _cfl_frac and dt > dt_min:
    #             dt = max(dt * (_cfl_frac / frac), dt_min)

    #         dt_applied = dt
    #         Y += dYdt * dt_applied
    #         Y[2:] = np.maximum(Y[2:], 0.0)
    #         t_current += dt_applied

    #         # if Y[2 + c_idx] > 2.0:  # 2.0 mol/kgw threshold
    #         #     print(f'\nTerminated at t = {t_current / YR / 1e6:.2f} Myr.')
    #         #     print('Reason: Carbon concentration exceeded physical solubility limit (Clathrate/Liquid CO2 formation).')
    #         #     convergence_status = 'snowball'
    #         #     convergence_reason = 'Carbon saturation limit reached (Snowball state)'
    #         #     last_diagnostics = diagnostics # Preserve state for JSON
    #         #     break

    #         P_CO2_max_threshold = 1e5
    #         T_snowball_threshold = 260
    #         T_moist_greenhouse_threshold = 350
            
    #         if Y[1] >= P_CO2_max_threshold and Y[0] < T_snowball_threshold:
    #             print(f'\nTerminated at t = {t_current / YR / 1e6:.2f} Myr.')
    #             print(f'Reason: Max CO2 greenhouse limit ({P_CO2_max_threshold/1e5:.1f} bar) reached without escaping Snowball state.')
    #             convergence_status = 'snowball'
    #             convergence_reason = 'Maximum CO2 greenhouse limit reached'
    #             last_diagnostics = diagnostics 
    #             break

    #         if  Y[0] > T_moist_greenhouse_threshold:
    #             print(f'\nTerminated at t = {t_current / YR / 1e6:.2f} Myr.')
    #             print(f'Temperature reached moist greenhouse state.')
    #             convergence_status = 'runaway'
    #             convergence_reason = 'Temperature reached moist greenhouse state'
    #             last_diagnostics = diagnostics 
    #             break

    #         # Grow dt toward dt_max for next step
    #         dt = min(dt_applied * _dt_grow, dt_max)

    #         T_s, P_CO2_s = self.solve_climate_from_chemistry(Y[2:], T_init=T_current, P_CO2_est=P_CO2_current)
    #         # Half-step blend to damp the 2-step C_chem/T oscillation.
    #         # The Newton map f(T) has f'(T*)≈-1 at the fixed point, so plain
    #         # iteration oscillates.  The blended map g(T)=0.5*f(T)+0.5*T has
    #         # g'(T*)≈0, giving rapid convergence to the self-consistent T*.
    #         T_current = 0.05 * T_s + 0.95 * T_current
    #         P_CO2_current = 0.05 * P_CO2_s + 0.95 * P_CO2_current
    #         Y[0] = T_current
    #         Y[1] = P_CO2_current
    #         b_ocean = Y[2:]

    #         # --- Record values ---

    #         t.append(t_current)
    #         T.append(Y[0])
    #         P_CO2.append(Y[1])
    #         concentrations.append(Y[2:].copy())
    #         fluxes.append(dYdt[2:].copy())
    #         pH.append(diagnostics['pH'])

    #         # --- Calculate fluxes for convergence ---

    #         F_net = diagnostics['F_net']
    #         F_in = np.abs(diagnostics['F_diss']) + np.abs(diagnostics['F_vol']) + np.abs(diagnostics['F_cont'])
    #         F_out = np.abs(diagnostics['F_prec'])

    #         # Flux-balance convergence: |F_net[i]| / (F_in[i] + F_out[i]).
    #         # This is zero only when inputs exactly match outputs for every element.
    #         # It is independent of how large Y has grown, so it cannot give false
    #         # convergence the way a relative-drift-rate criterion can.
    #         F_throughput = F_in + F_out
    #         active = F_throughput > 1e-40
    #         if active.any():
    #             rel_imbalance = np.abs(F_net[active]) / F_throughput[active]
    #             max_flux_imbalance = float(np.max(rel_imbalance))
    #         else:
    #             max_flux_imbalance = 0.0

    #         rel_change = max_flux_imbalance   # written to CSV

    #         # --- Writes to CSV file ---

    #         row = {
    #             't_Myr': f'{t_current/YR/1e6:.4f}',
    #             'dt_kyr': f'{dt_applied/YR/1e3:.2f}',
    #             'T': f'{Y[0]:.4f}',
    #             'P_CO2': f'{Y[1]:.6e}',
    #             'pH': f"{diagnostics['pH']:.4f}",
    #             'rel_change': f'{rel_change:.6e}',
    #             # 'slowest': slowest_label,
    #             'Da': f'{diagnostics["Da"]:.6e}',
    #             'A_reactive': f'{diagnostics["A_reactive"]:.6e}',
    #             'supply_efficiency': f'{diagnostics["supply_efficiency"]:.4%}' # Formatted as a percentage!
    #         }
            
    #         if write_csv:
    #             for j, e in enumerate(elements):
    #                 row[f'Y_{e}']      = f'{Y[2+j]:.6e}'
    #                 row[f'F_vol_{e}']  = f'{diagnostics["F_vol"][j]:.6e}'
    #                 row[f'F_diss_{e}'] = f'{diagnostics["F_diss"][j]:.6e}'
    #                 row[f'F_prec_{e}'] = f'{diagnostics["F_prec"][j]:.6e}'
    #                 row[f'F_cont_{e}'] = f'{diagnostics["F_cont"][j]:.6e}'
    #                 row[f'F_net_{e}']  = f'{diagnostics["F_net"][j]:.6e}'
    #             _csv_writer.writerow(row)
    #             _csv_file.flush()

    #         # --- Termination checks ---

    #         if max_flux_imbalance < tol:
    #             print(f'\nNormal convergence at t = {t_current / YR / 1e6:.2f} Myr.')
    #             print(f'Max flux imbalance: {max_flux_imbalance:.2%}')
    #             break

    #         # if t_current > t_max:
    #         #     print('Ran out of time. No steady state likely.')

    #     else:
    #         print(f'Max flux imbalance: {max_flux_imbalance:.2%}')
    #         print('Not converged in maximum steps. Need more iterations.')
        
    #     if write_csv:
    #         _csv_file.close()
    #         print(f"\nTimeseries written to: {os.path.abspath(output_csv)}")

    #     # --- Write summary JSON ---
    #     final_ocean_chemistry = {str(e): float(Y[2 + j]) for j, e in enumerate(elements)}
    #     summary = {
    #         'input_parameters': {
    #             'name': self.name,
    #             'mass_kg': float(self.mass),
    #             'radius_m': float(self.radius),
    #             'surface_pressure_Pa': float(self.P_background),
    #             'instellation_solar': float(self._input_instellation),
    #             'tectonics': float(self._input_tectonics),
    #             'land_fraction': float(self.land_fraction),
    #             'ocean_depth_m': float(self.ocean_depth),
    #             'alpha': float(self.alpha),
    #             'tidally_locked': bool(self.tidally_locked),
    #             #'cl_chemistry': bool(self.cl_chemistry),
    #             #'hcl_co2_ratio': float(self.hcl_co2_ratio),
    #         },
    #         'convergence': {
    #             'status': convergence_status,
    #             'reason': convergence_reason,
    #             'time_to_converge_Myr': float(t_current / YR / 1e6),
    #             'snowball': bool(convergence_status == 'snowball'),
    #             'runaway': bool(convergence_status == 'stagnated' and abs(Y[0] - T_max) < 0.5),
    #             'stagnated': bool(convergence_status == 'stagnated'),
    #         },
    #         'final_state': {
    #             'T_K': float(Y[0]),
    #             'P_CO2_Pa': float(Y[1]),
    #             'pH': float(last_diagnostics['pH']) if last_diagnostics is not None else None,
    #             'ocean_chemistry_mol_kgw': final_ocean_chemistry,
    #         },
    #         'final_diagnostics': {
    #             'Da': float(last_diagnostics['Da']) if last_diagnostics is not None else None,
    #             'A_reactive_m2': float(last_diagnostics['A_reactive']) if last_diagnostics is not None else None,
    #             'supply_efficiency': float(last_diagnostics['supply_efficiency']) if last_diagnostics is not None else None,
    #         },
    #     }
    #     with open(summary_file, 'w') as _sf:
    #         json.dump(summary, _sf, indent=2)
    #     print(f"Summary written to:    {os.path.abspath(summary_file)}")

    #     # --- Plot results ---
    #     t = np.array(t)
    #     T = np.array(T)
    #     P_CO2 = np.array(P_CO2)
    #     concentrations = np.array(concentrations)
    #     fluxes = np.array(fluxes)
    #     pH = np.array(pH)

    #     fig, (ax1, ax2, ax4, ax5) = plt.subplots(4, 1, figsize=(15, 24), sharex=True)

    #     for j, elem in enumerate(elements):
    #         ax1.plot(t / YR, concentrations[:, j], label=elem)

    #     for j, elem in enumerate(elements):
    #         ax4.plot(t / YR, fluxes[:, j] * YR, label=elem)
        
    #     ax1.set_yscale('log')
    #     # ax1.set_xscale('log')
    #     ax1.set_ylabel('Concentration (mol/kgw)')
    #     ax1.set_title('Time Evolution of Ocean Chemistry')
    #     ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    #     ax1.grid(True, which="both", ls="--", alpha=0.5)
    #     ax1.set_ylim([1e-6, 2.1])

    #     color_co2 = 'tab:blue'
    #     ax5.set_xlabel('Time (years)')
    #     ax2.set_ylabel('$P_{CO2}$ (bar)', color=color_co2)
    #     ax2.plot(t / YR, P_CO2 / 1e5, color=color_co2, label='P_CO2')
    #     ax2.tick_params(axis='y', labelcolor=color_co2)
    #     ax2.set_yscale('log')
    #     ax2.grid(True, alpha=0.5)
    #     ax2.set_ylim([1e-7, 1])

    #     ax3 = ax2.twinx()  
    #     color_t = 'tab:red'
    #     ax3.set_ylabel('Temperature (K)', color=color_t)
    #     ax3.plot(t / YR, T, color=color_t, label='Temperature')
    #     ax3.tick_params(axis='y', labelcolor=color_t)
    #     ax3.set_ylim([T_min, T_max])

    #     ax4.set_yscale('symlog', linthresh=1e-12)
    #     ax4.set_ylabel('Fluxes (mol/kgw/yr)')
    #     ax4.grid(True, which="both", ls="--", alpha=0.5)

    #     ax5.plot(t / YR, pH)
    #     ax5.set_ylabel('pH')
    #     ax5.grid(True, which="both", ls="--", alpha=0.5)

    #     fig.tight_layout()  
    #     plt.savefig(os.path.join(output_dir, f'{self.name}_plot.pdf'))
    #     plt.savefig(os.path.join(output_dir, f'{self.name}_plot.png'))
    #     plt.close()