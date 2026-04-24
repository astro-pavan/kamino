import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from kamino.constants import *
from kamino.chemistry.ocean_chemistry import *
from kamino.climate.clima_interpolator import get_T_surface
from kamino.utils import august_roche_magnus_formula

tau_prec = 1e4 * YR
tau_atm = 1e4 * YR

class Planet:

    def __init__(self, mass: float, radius: float, background_pressure: float, instellation : float, tectonics: float, ocean_depth: float, land_fraction: float=0.0, name: str='planet'):        

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

        F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=carbonate_minerals + secondary_sink_minerals, precipitation_timescale=tau_prec)

        ocean_water_per_area = self.ocean_depth * 1000.0
        F_carb, F_sil = max(0.0, -F_prec[c_idx]), max(0.0, -F_prec[si_idx])

        S_sed = (F_carb * 0.100 / 2710.0 + F_sil * 0.060 / 2650.0) * ocean_water_per_area
        weathering_flux, w_diag = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, self.alpha, self.crust_production_rate, clog=False, sedimentation_rate=S_sed)

        F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

        F_net = F_vol + F_prec + F_diss

        # return derivatives

        dYdt = np.zeros_like(Y)

        dYdt[0] = (P_CO2_new - P_CO2) / tau_atm
        dYdt[1] = (P_H2O_new - P_H2O) / tau_atm

        dYdt[2:] = F_net

        carbon_flux = dYdt[3]
        carbon = Y[3]

        # print(f't = {t/YR:.4e} yr', end='\r')
        print(f't = {t/YR:.4e} yr  T = {T_surface:.1f}  P_CO2 = {P_CO2 / 1e5:.1e} bar  C flux = {(carbon_flux / carbon) * 1e9 * YR:.1e}/Gyr ', end='\r')
        # print(Y[2:])
        # print(dYdt[2:])

        return dYdt

    def time_evolve(self, t_end=2e9 * YR):

        Y0 = np.zeros(elements.shape[0] + 2)

        Y0[0] = 1000
        Y0[1] = 1000

        P_CO2_max = 1e5
        T_runaway = 350
        T_snowball = 260

        def event_runaway(t, Y):
            return Y[0] - P_CO2_max
        event_runaway.terminal = True

        def event_snowball(t, Y):
            T = get_T_surface(self.instellation, max(Y[0], 1e-2), self.albedo, self.tidally_locked)
            return T - T_snowball
        event_snowball.terminal = True

        def event_hothouse(t, Y):
            T = get_T_surface(self.instellation, max(Y[0], 1e-2), self.albedo, self.tidally_locked)
            return T - T_runaway
        event_hothouse.terminal = True

        atol = np.ones_like(Y0) * 1e-6
        atol[0] = 1.0   # P_CO2 in Pa
        atol[1] = 1.0   # P_H2O in Pa

        def macro_jacobian(t, y):
            # print('Using custom jacobian                                                                                ')
            N = len(y)
            jac = np.zeros((N, N))

            f0 = self.dY_dt(t, y)

            rel_rates = np.abs(f0 / (np.abs(y) + 1e-7))
            v_max = np.max(rel_rates[2:])
            
            # 3. Apply the inverse function
            eps_max = 0.5
            eps_min = 1e-3
            
            # Tuning knob: Because your dy/dt is per second over geologic time, 
            # v_max will be tiny. We need a massive k to scale it up. 
            # You will need to tune this! Start with 1e10 or 1e12.
            k = 1e9 * YR 
            
            dynamic_eps = eps_max / (1.0 + k * v_max)
            
            # Clamp to safety limits
            dynamic_eps = np.clip(dynamic_eps, eps_min, eps_max)
            
            # Choose a "macro" step size. 
            eps = dynamic_eps
            min_delta = 1e-5
            delta = min_delta + eps * np.abs(y)

            
            for j in range(N):
                # Create perturbed states
                y_plus = np.copy(y)
                y_minus = np.copy(y)
                
                y_plus[j] += delta[j]
                y_minus[j] -= delta[j]
                
                # Evaluate ODE at perturbed states
                f_plus = self.dY_dt(t, y_plus)
                f_minus = self.dY_dt(t, y_minus)
                
                # Calculate the macro-derivative for column j
                jac[:, j] = (f_plus - f_minus) / (2 * delta[j])
                
            return jac

        sol = solve_ivp(
            self.dY_dt,
            (0, t_end),
            Y0,
            method='LSODA',
            max_step=1e7 * YR,
            rtol=1e-3,
            atol=atol,
            jac=macro_jacobian
            # events=[event_runaway, event_snowball, event_hothouse],
        )

        print()
        print(sol.y[2:, -1])

        Alk = sol.y[2, :]
        C = sol.y[3, :]
        t = sol.t / YR

        plt.plot(t, C)
        plt.yscale('log')
        plt.show()

        dC = np.gradient(C)
        dAlk = np.gradient(Alk)
        dt = np.gradient(t)

        plt.plot(t, ((dC / dt) / C) * 1e9)
        plt.plot(t, ((dAlk / dt) / Alk) * 1e9)
        plt.axhline(0, color='black', linestyle='--')
        # plt.yscale('log')
        plt.yscale('symlog', linthresh=1e-3)
        plt.show()




if __name__ == '__main__':

    BACKGROUND_PRESSURE = 1e5   # Pa (1 bar)
    OCEAN_DEPTH = 3000          # m
    TECTONICS = 1.0

    s = 0.7

    print(f's = {s}')
    p1 = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, s, 1.0, 3000, name=f'test_s_{s}')
    p1.time_evolve()
    print('')