import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import json
import os

from kamino.constants import *
from kamino.chemistry import elements, get_P_CO2, c_idx, si_idx, alk_idx, ca_idx, mg_idx, ChemistryError
from kamino.weathering import get_weathering_flux, get_continental_weathering_flux, J_ref_normalised, rate_ref
from kamino.precipitation import get_precipitation
from kamino.mineral_info import *


from kamino.climate.clima_interpolator import get_T_surface
from kamino.utils import august_roche_magnus_formula

output_path = os.path.join(os.path.dirname(__file__), '../../output/')

# Mg→Ca exchange partition coefficient for HT hydrothermal circulation.
# Calibrated so that at Earth conditions ([Mg]=53 mM, crust_rate=1) the exchange gives
# ~1.4 Tmol/yr: KD_MG_HT = 1.4e12 × 0.7 / (0.053 × J_ref × YR), where 0.7 = A_seafloor/A_total.
KD_MG_HT = 0.7 * 1.4e12 / (0.053 * 1.4e15)  # ≈ 0.013

tau_prec = 1e5 * YR
tau_atm = 1e4 * YR
tau_r_avg = 3e7 * YR   # EMA timescale for convergence rate smoothing (~30 Myr)

class Planet:

    def __init__(
            self,
            mass: float,
            radius: float,
            background_pressure: float,
            instellation : float,
            crust_production_rate: float,
            outgassing: float,
            ocean_depth: float,
            land_fraction: float=0.0,
            crust_composition: dict[str, float]=basalt_composition,
            crust_carbonate_content: float=0.0,
            reverse_weathering: bool=True,
            f_HT: float=0.0,
            f_carb: float=0.0,
            verbose: bool=False,
            name: str='planet'
            ):
        
        planet_config = {
            "name": name,
            "mass": mass,
            "radius": radius,
            "background_pressure": background_pressure,
            "ocean_depth": ocean_depth,
            "land_fraction": land_fraction,
            "crust_composition": crust_composition.copy(),
            "instellation": instellation,
            "crust_production_rate": crust_production_rate,
            "outgassing": outgassing,
            "crust_carbonate_content": crust_carbonate_content,
            "reverse_weathering": reverse_weathering,
            "f_HT": f_HT,
            "f_carb": f_carb
        }

        self.name = name
        self.mass = mass
        self.radius = radius
        self.gravity = (G * self.mass) / (self.radius ** 2)
        self.surface_area = 4 * np.pi * self.radius ** 2
        self.P_background = background_pressure

        self.ocean_depth = ocean_depth
        self.ocean_water_mass = self.ocean_depth * self.surface_area * 1000

        self.crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * crust_production_rate
        self.hydrothermal_flux = EARTH_HYDROTHERMAL_FLUX_PER_AREA * crust_production_rate

        self.crust_composition = crust_composition.copy()

        if crust_carbonate_content > 0:
            for mineral in self.crust_composition:
                self.crust_composition[mineral] *= (1 - crust_carbonate_content)
            self.crust_composition['Calcite'] = crust_carbonate_content

        # Calcite excluded from pore space: carbonate veins in basalt form on Myr timescales,
        # much slower than tau_prec = 1e5 yr. Including it captures all basalt-derived Ca in
        # situ and starves the ocean of Ca entirely.
        self.pore_precipitating_minerals = [m for m in carbonate_minerals if m != 'Calcite'] + clay_minerals
        self.ocean_precipitating_minerals = carbonate_minerals + clay_minerals + silica_minerals
        self.shelf_depth = 1000.0  # m — representative continental shelf depth for carbonate burial
        self.shelf_precipitating_minerals = carbonate_minerals

        if reverse_weathering:
            # Reverse weathering (saponite authigenesis) occurs in marine sediment porewaters,
            # not in the seafloor basalt alteration zone — Isson & Planavsky (2018).
            self.ocean_precipitating_minerals += reverse_weathering_minerals

        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[1] = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing

        self.tidally_locked = False

        self.f_HT = f_HT
        self.f_carb = f_carb

        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area

        ocean_albedo = 0.3
        land_albedo = 0.3

        self.relative_humidity = 0.5

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

        self._output_filename = f"{output_path}{self.name}.json"
        with open(self._output_filename, 'w') as f:
            json.dump(planet_config, f, indent=0)

        self.verbose = verbose

    def dY_dt(self, t, Y):

        P_CO2 = Y[0]
        P_H2O = Y[1]
        b_ocean = Y[2:-1]   # Y[-1] is r_avg, excluded from chemistry

        # input safety

        P_CO2 = np.clip(P_CO2, 0, 1e6)
        P_H2O = np.maximum(0, P_H2O)
        b_ocean = np.maximum(b_ocean, 0.0)

        # atmosphere properties

        P_surface = self.P_background + P_CO2 + P_H2O
        T_surface = get_T_surface(self.instellation, max(P_CO2, 1.0), self.albedo, self.tidally_locked)
        P_H2O_new = august_roche_magnus_formula(T_surface)

        # seafloor properties

        T_seafloor = 1.02 * T_surface - 16.7
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth
        T_seafloor = np.maximum(T_seafloor, 274)
        T_pore = T_seafloor + 9

        # ocean fluxes

        F_vol = self.outgassing_flux / self.ocean_water_mass

        self._T = T_surface

        try:
            P_CO2_new = get_P_CO2(P_surface, T_surface, b_ocean)
            assert P_CO2_new > 0
            F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=self.ocean_precipitating_minerals, precipitation_timescale=tau_prec)
            self._pH = pH
            self._SI = SI

            ocean_water_per_area = self.ocean_depth * 1000.0
            F_carb, F_sil = max(0.0, -F_prec[c_idx]), max(0.0, -F_prec[si_idx])

            S_sed = (F_carb * 0.100 / 2710.0 + F_sil * 0.060 / 2650.0) * ocean_water_per_area
            J_total = J_ref_normalised * (self.crust_production_rate / rate_ref)
            J_LT = J_total * (1 - self.f_HT)
            J_HT = J_total * self.f_HT

            flux_lt, _ = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, 
                                             rate=self.crust_production_rate, J=J_LT, crust_composition=self.crust_composition,
                                             sedimentation_rate=S_sed, precipitating_minerals=self.pore_precipitating_minerals)
            
            if self.f_HT > 0:
                flux_ht, _ = get_weathering_flux(P_pore, 600, P_CO2, b_ocean, high_temperature=True,
                                                 rate=self.crust_production_rate, J=J_HT, crust_composition=self.crust_composition,
                                                 sedimentation_rate=S_sed, precipitating_minerals=self.pore_precipitating_minerals)
            else:
                flux_ht = np.zeros_like(flux_lt)

            weathering_flux = flux_lt + flux_ht

            F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

            F_cont = np.zeros_like(F_diss)
            F_shelf_prec = np.zeros_like(F_diss)
            if self.land_fraction > 0:
                F_sil_cont = get_continental_weathering_flux(T_surface, P_CO2).copy()
                F_sil_cont[c_idx] = 0.0  # C came from atmosphere; adding it here raises P_CO2 (double-counting)
                # Continental carbonate weathering: CaCO3 + CO2 + H2O -> Ca2+ + 2HCO3-
                # C here is geological (from land CaCO3), so it is NOT zeroed out.
                F_carb_w = np.zeros(elements.shape)
                F_carb_w[ca_idx] = self.f_carb * F_sil_cont[ca_idx]
                F_carb_w[alk_idx] = 2.0 * self.f_carb * F_sil_cont[ca_idx]
                F_carb_w[c_idx] = self.f_carb * F_sil_cont[ca_idx]
                F_cont = (F_sil_cont + F_carb_w) * self.land_area / self.ocean_water_mass

                # Shallow carbonate precipitation on continental shelves (above the CCD).
                # Uses T_seafloor (below thermocline) and 1000m pressure — higher Calcite SI
                # than the deep seafloor drives stronger carbonate burial.
                P_shelf = P_surface + 1000 * self.gravity * self.shelf_depth
                F_shelf_prec, _, _ = get_precipitation(P_shelf, T_seafloor, b_ocean,
                                                       precipitating_minerals=self.shelf_precipitating_minerals,
                                                       precipitation_timescale=tau_prec)

            # Simple HT Mg→Ca exchange: hot fluids strip Mg from seawater, Ca released from basalt.
            # Rate proportional to [Mg] and crust production; replaces the complex f_HT PHREEQC path.
            _ht_rate = KD_MG_HT * b_ocean[mg_idx] * J_total * self.surface_area / self.ocean_water_mass
            F_ht_exchange = np.zeros(elements.shape)
            F_ht_exchange[mg_idx] = -_ht_rate
            F_ht_exchange[ca_idx] = +_ht_rate

            F_net = F_vol + F_prec + F_shelf_prec + F_diss + F_cont + F_ht_exchange

        except (ChemistryError, AssertionError):
            # Chemistry has left the valid domain (typically high P_CO2 where PHREEQC cannot converge).
            # Return pure outgassing so LSODA gets a finite derivative and the acid_ocean event can terminate cleanly.
            dYdt = np.zeros_like(Y)
            dYdt[2:-1] = F_vol
            self._F_net = dYdt[2:-1]
            return dYdt

        dYdt = np.zeros_like(Y)

        dYdt[0] = (P_CO2_new - P_CO2) / tau_atm
        dYdt[1] = (P_H2O_new - P_H2O) / tau_atm

        F_net[b_ocean <= 0.0] = np.maximum(F_net[b_ocean <= 0.0], 0.0)

        dYdt[2:-1] = F_net
        self._F_net = F_net

        # Relaxation equation for smoothed convergence rate (r_avg = Y[-1])
        significant = b_ocean > 1e-7
        if np.any(significant):
            max_frac_rate = np.max(np.abs(F_net[significant]) / np.maximum(b_ocean[significant], 1e-6))
        else:
            max_frac_rate = 0.0
        dYdt[-1] = (max_frac_rate - Y[-1]) / tau_r_avg

        carbon_flux = dYdt[3]
        carbon = Y[3]
        calcite_SI = SI['Calcite']

        if self.verbose and not getattr(self, '_jac_active', False):
            self._step_count = getattr(self, '_step_count', 0) + 1
            dt_str = f'  dt={((t - self._last_t)/YR):.1e}yr' if hasattr(self, '_last_t') else ''
            print(f't = {t/YR:.4e} yr  T = {T_surface:.1f}  P_CO2 = {P_CO2 / 1e5:.1e} bar  pH = {pH:.1f}  Calcite SI = {calcite_SI:.1f}  C flux = {(carbon_flux / carbon) * 1e9 * YR:.1e} / Gyr  step={self._step_count}{dt_str}  ', end='\r')
            self._last_t = t

        self._pH = pH

        return dYdt

    def time_evolve(self, t_end=2e9 * YR, jac_epsilon=0.01):

        Y0 = np.zeros(elements.shape[0] + 3)  # +2 for P_CO2/P_H2O, +1 for r_avg

        Y0[0] = 1000
        Y0[1] = 1000
        Y0[-1] = 1.0 / (1e6 * YR)  # r_avg starts high (1/Myr) so convergence can't fire immediately

        T_runaway = 350
        T_snowball = 260
        P_CO2_acid_threshold = 5e5   # Pa (5 bar)
        P_CO2_floor = 1.0            # Pa — below this the planet is unambiguously going snowball
        min_time = 2e6 * YR
        convergence_rate = 0.1 / (1e9 * YR)

        self._T = np.nan
        self._pH = np.nan

        def event_snowball(t, Y):
            if t < min_time:
                return -1.0  # negative guard: direction=-1 won't trigger on the guard→actual transition
            P_CO2 = np.clip(Y[0], 0, 1e6)
            T_surface = get_T_surface(self.instellation, max(P_CO2, 1.0), self.albedo, self.tidally_locked)
            return T_surface - T_snowball
        event_snowball.terminal, event_snowball.direction = True, -1 # type: ignore

        def event_hothouse(t, Y):
            if t < min_time:
                return -1.0  # negative guard: direction=+1 won't trigger on the guard→actual transition
            T = get_T_surface(self.instellation, max(Y[0], 1.0), self.albedo, self.tidally_locked)
            return T - T_runaway
        event_hothouse.terminal, event_hothouse.direction = True, +1 # type: ignore

        def event_acid_ocean(t, Y):
            if t < min_time:
                return 1.0
            return P_CO2_acid_threshold - np.clip(Y[0], 0, None)
        event_acid_ocean.terminal, event_acid_ocean.direction = True, -1 # type: ignore

        def event_co2_floor(t, Y):
            if t < min_time:
                return -1.0  # negative guard: direction=-1 won't trigger on the guard→actual transition
            P_CO2 = np.clip(Y[0], 0, 1e6)
            T_surface = get_T_surface(self.instellation, max(P_CO2, 1.0), self.albedo, self.tidally_locked)
            if T_surface < T_snowball:
                return 1.0  # snowball handles this case
            return np.clip(Y[0], 0, None) - P_CO2_floor
        event_co2_floor.terminal, event_co2_floor.direction = True, -1 # type: ignore

        atol = np.ones_like(Y0) * 1e-6
        atol[0] = 1.0   # P_CO2 in Pa
        atol[1] = 1.0   # P_H2O in Pa
        atol[-1] = convergence_rate * 0.1  # r_avg: resolve to 10% of the convergence threshold

        chem_significant = 1e-7  # mol/kgw
        min_time_frozen = 1e8 * YR  # 100 Myr: if no carbonate cycle by here, ocean is frozen

        def event_converged(t, Y):
            if t < min_time:
                return 1.0
            P_CO2 = np.clip(Y[0], 0, 1e6)
            T_surface = get_T_surface(self.instellation, max(P_CO2, 1e-2), self.albedo, self.tidally_locked)
            if T_surface < T_snowball:
                return 1.0  # snowball conditions — defer to event_snowball
            b = Y[2:-1]
            if not np.any(b > chem_significant):
                return min_time_frozen - t
            return Y[-1] - convergence_rate  # r_avg vs threshold; pure function of (t, Y)
        event_converged.terminal, event_converged.direction = True, -1 # type: ignore

        N = len(Y0)

        def macro_jacobian(t, y):

            jac = np.zeros((N, N))

            eps_abs = np.empty(N)
            eps_abs[0] = 1e-3     # P_CO2  [Pa] — small so jac_epsilon*|P_CO2| dominates at all relevant values
            eps_abs[1] = 0.1      # P_H2O  [Pa]
            eps_abs[2:] = 1e-9   # b_ocean [mol/kgw]  (trace_approx threshold)

            delta = np.maximum(jac_epsilon * np.abs(y), eps_abs)

            self._jac_active = True # type: ignore
            try:
                for j in range(N):

                    if j == 1:
                        jac[1, 1] = -1.0 / tau_atm
                        continue

                    y_plus = np.copy(y)
                    y_minus = np.copy(y)

                    y_plus[j] += delta[j]
                    y_minus[j] -= delta[j]

                    try:
                        f_plus = self.dY_dt(t, y_plus)
                        f_minus = self.dY_dt(t, y_minus)
                    except ChemistryError:
                        continue

                    jac[:, j] = (f_plus - f_minus) / (2 * delta[j])
            finally:
                self._jac_active = False # type: ignore

            return jac

        start = time.time()

        sol = solve_ivp(
            self.dY_dt,
            (0, t_end),
            Y0,
            method='LSODA',
            max_step=1e7 * YR,
            rtol=1e-3,
            atol=atol,
            jac=macro_jacobian,
            events=[event_snowball, event_hothouse, event_acid_ocean, event_converged, event_co2_floor],
        )

        end = time.time()

        event_names = ['snowball', 'hothouse', 'acid_ocean', 'converged', 'co2_floor']

        sol.y = np.maximum(sol.y, 0.0)
        time_steps = sol.t.tolist()
        state_variables = sol.y.tolist()

        if sol.t[-1] >= t_end:
            termination = "timeout"
        elif sol.status == 1:
            termination = next(name for name, t_ev in zip(event_names, sol.t_events) if len(t_ev) > 0)
        else:
            termination = "solver_failure"

        if self.verbose:
            print()
            print(f'Terminated: {termination} at t = {sol.t[-1]/YR:.3e} yr')
            print(f'Simulation time: {end - start:.0f} s')
            print(f'Y: {sol.y[2:, -1]} mol/kgw')

        with open(self._output_filename, 'r') as f:
            output_data = json.load(f)

        output_data.update({
            "simulation_time_seconds": end - start,
            "termination": termination,
            "end_time_yr": sol.t[-1] / YR,
            "T": self._T,
            "P_CO2": sol.y[0, -1] / 1e5,
            "pH": self._pH,
            "data": {
                "time": time_steps,
                "y": state_variables
            }
        })

        with open(self._output_filename, 'w') as f:
            json.dump(output_data, f, indent=0)

        if self.verbose:
            print(f"Results successfully saved to {self._output_filename}")



if __name__ == '__main__':

    BACKGROUND_PRESSURE = 1e5   # Pa (1 bar)
    OCEAN_DEPTH = 3000          # m
    TECTONICS = 1.0

    instellation = [0.6, 0.8, 1.0, 1.2]

    # for rw in [False, True]:
    #     tag = 'rw' if rw else 'norw'
    for s in instellation:
        print(f's = {s}')
        p = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, s, 1.0, 3.0, 3000,
                    name=f'test_s_{s}', f_HT=0.005, verbose=True)
        p.time_evolve()
        print('')
