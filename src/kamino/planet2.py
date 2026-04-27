import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import json
import os

from kamino.constants import *
from kamino.chemistry.ocean_chemistry import *
from kamino.chemistry.mineral_info import *
from kamino.climate.clima_interpolator import get_T_surface
from kamino.utils import august_roche_magnus_formula

output_path = os.path.join(os.path.dirname(__file__), '../../output/')

tau_prec = 1e4 * YR
tau_atm = 1e4 * YR

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
            crust_carbonate_content: float=0.0,
            reverse_weathering: bool=False,
            name: str='planet'):        

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

        self.crust_composition = basalt_composition

        if crust_carbonate_content > 0:
            for mineral, fraction in self.crust_composition.items():
                self.crust_composition[mineral] *= (1 - crust_carbonate_content) * fraction
            self.crust_composition['Calcite'] = crust_carbonate_content

        self.precipitating_minerals = carbonate_minerals + secondary_sink_minerals
        if reverse_weathering:
            self.precipitating_minerals += reverse_weathering_minerals

        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[1] = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing

        self.tidally_locked = False

        self._input_instellation = instellation
        self._input_crust_production_rate = crust_production_rate
        self._input_outgassing = outgassing
        self._input_crust_carbonate_content = crust_carbonate_content
        self._input_reverse_Weathering = reverse_weathering

        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area

        ocean_albedo = 0.3
        land_albedo = 0.3

        self.relative_humidity = 0.5

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

        planet_config = {
            "name": self.name,
            "mass": self.mass,
            "radius": self.radius,
            "background_pressure": self.P_background,
            "ocean_depth": self.ocean_depth,
            "land_fraction": self.land_fraction,
            "instellation": self._input_instellation,
            "crust_production_rate": self._input_crust_production_rate,
            "outgassing": self._input_outgassing,
            "crust_carbonate_content": self._input_crust_carbonate_content,
            "reverse_weathering": self._input_reverse_Weathering
        }

        filename = f"{output_path}{self.name}_config.json"
        with open(filename, 'w') as f:
            json.dump(planet_config, f, indent=4)

    def dY_dt(self, t, Y):

        P_CO2 = Y[0]
        P_H2O = Y[1]
        b_ocean = Y[2:]

        # input safety

        P_CO2 = np.clip(P_CO2, 0, 1e6)
        P_H2O = np.maximum(0, P_H2O)
        b_ocean = np.maximum(b_ocean, 0.0)

        # atmosphere properties

        P_surface = self.P_background + P_CO2 + P_H2O
        T_surface = get_T_surface(self.instellation, P_CO2, self.albedo, self.tidally_locked)
        P_H2O_new = august_roche_magnus_formula(T_surface)
        P_CO2_new = get_P_CO2(P_surface, T_surface, b_ocean)

        assert P_CO2_new > 0

        # seafloor properties

        T_seafloor = 1.02 * T_surface - 16.7
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth
        T_seafloor = np.maximum(T_seafloor, 274)
        T_pore = T_seafloor + 9

        # ocean fluxes

        F_vol = self.outgassing_flux / self.ocean_water_mass

        try:
            F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=carbonate_minerals + secondary_sink_minerals, precipitation_timescale=tau_prec)
            self._pH = pH
            self._SI = SI

            ocean_water_per_area = self.ocean_depth * 1000.0
            F_carb, F_sil = max(0.0, -F_prec[c_idx]), max(0.0, -F_prec[si_idx])

            S_sed = (F_carb * 0.100 / 2710.0 + F_sil * 0.060 / 2650.0) * ocean_water_per_area
            weathering_flux, _ = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, None, self.crust_production_rate, crust_composition=self.crust_composition, clog=False, sedimentation_rate=S_sed)

            F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass
            F_net = F_vol + F_prec + F_diss
        except ChemistryError:
            # Chemistry has left the valid domain (typically high P_CO2 where PHREEQC cannot converge). 
            # Return pure outgassing so LSODA gets a finite derivative and the acid_ocean event can terminate cleanly.
            dYdt = np.zeros_like(Y)
            dYdt[2:] = F_vol
            self._F_net = dYdt[2:]
            return dYdt

        dYdt = np.zeros_like(Y)

        dYdt[0] = (P_CO2_new - P_CO2) / tau_atm
        dYdt[1] = (P_H2O_new - P_H2O) / tau_atm

        dYdt[2:] = F_net
        self._F_net = F_net

        carbon_flux = dYdt[3]
        carbon = Y[3]
        calcite_SI = SI['Calcite']

        print(f't = {t/YR:.4e} yr  T = {T_surface:.1f}  P_CO2 = {P_CO2 / 1e5:.1e} bar  pH = {pH:.1f}  Calcite SI = {calcite_SI:.1f}  C flux = {(carbon_flux / carbon) * 1e9 * YR:.1e} / Gyr  ', end='\r')

        self._T = T_surface
        self._pH = pH

        return dYdt

    def time_evolve(self, t_end=2e9 * YR, jac_epsilon=0.01, convergence_rate = 0.01 / (1e9 * YR)):

        Y0 = np.zeros(elements.shape[0] + 2)

        Y0[0] = 1000
        Y0[1] = 1000

        P_CO2_max = 1e6
        T_runaway = 350
        T_snowball = 260
        P_CO2_acid_threshold = 5e5   # Pa (5 bar)

        def event_snowball(t, Y):
            if Y[0] < P_CO2_max:
                return 1.0
            T = get_T_surface(self.instellation, max(Y[0], 1e-2), self.albedo, self.tidally_locked)
            return T - T_snowball
        event_snowball.terminal = True

        def event_hothouse(t, Y):
            T = get_T_surface(self.instellation, max(Y[0], 1e-2), self.albedo, self.tidally_locked)
            return T - T_runaway
        event_hothouse.terminal = True
        event_hothouse.direction = 1  # only fire when T crosses upward through T_runaway

        def event_acid_ocean(t, Y):
            return P_CO2_acid_threshold - Y[0]
        event_acid_ocean.terminal = True
        event_acid_ocean.direction = -1  # only fire when P_CO2 crosses upward through threshold

        atol = np.ones_like(Y0) * 1e-6
        atol[0] = 1.0   # P_CO2 in Pa
        atol[1] = 1.0   # P_H2O in Pa

        chem_significant = 1e-7  # mol/kgw

        def event_converged(t, Y):
            b = Y[2:]
            mask = b > chem_significant
            if not np.any(mask):
                return 1.0
            try:
                dYdt = self.dY_dt(t, Y)
            except ChemistryError:
                return 1.0
            F_net = dYdt[2:]
            max_fractional_rate = np.max(np.abs(F_net[mask]) / b[mask])
            return max_fractional_rate - convergence_rate
        event_converged.terminal = True
        event_converged.direction = -1

        N = len(Y0)

        def macro_jacobian(t, y):

            jac = np.zeros((N, N))

            eps_abs = np.empty(N)
            eps_abs[0] = 0.1      # P_CO2  [Pa]   (atol=1 Pa)
            eps_abs[1] = 0.1      # P_H2O  [Pa]
            eps_abs[2:] = 1e-9   # b_ocean [mol/kgw]  (trace_approx threshold)

            delta = np.maximum(jac_epsilon * np.abs(y), eps_abs)

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
            events=[event_snowball, event_hothouse, event_acid_ocean, event_converged],
        )

        end = time.time()

        print()

        event_names = ['snowball', 'hothouse', 'acid_ocean', 'converged']
        for name, t_ev in zip(event_names, sol.t_events):
            if len(t_ev) > 0:
                print(f'Terminated: {name} at t = {t_ev[0]/YR:.3e} yr')

        if sol.t[-1] == t_end:
            print(f'Simulation time out at t = {t_end/YR:.3e} yr')

        print(f'Simulation time: {end - start:.0f} s')
        print(f'Y: {sol.y[2:, -1]} mol/kgw')

        time_steps = sol.t.tolist()
        state_variables = sol.y.tolist()
        
        events_dict = {}
        event_names = ['snowball', 'hothouse', 'acid_ocean', 'converged']
        for name, t_ev in zip(event_names, sol.t_events):
            events_dict[name] = t_ev.tolist()

        results_data = {
            "simulation_time_seconds": end - start,
            "status": sol.status,
            "message": sol.message,
            "events_triggered": events_dict,
            "T": self._T,
            "P_CO2": sol.y[0, -1] / 1e5,
            "pH": self._pH,
            "data": {
                "time": time_steps,
                "y": state_variables
            }
        }

        results_filename = f"{output_path}{self.name}_results.json"
        with open(results_filename, 'w') as f:
            json.dump(results_data, f, indent=4)
            
        print(f"Results successfully saved to {results_filename}")



if __name__ == '__main__':

    BACKGROUND_PRESSURE = 1e5   # Pa (1 bar)
    OCEAN_DEPTH = 3000          # m
    TECTONICS = 1.0

    instellation = [0.4, 0.6, 0.8, 1.0, 1.2]

    for s in instellation:
        print(f's = {s}')
        p1 = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, s, 1.0, 1.0, 3000, name=f'test_s_{s}')
        p1.time_evolve()
        print('')