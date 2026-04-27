import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time

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
        self._pH = pH
        self._SI = SI

        ocean_water_per_area = self.ocean_depth * 1000.0
        F_carb, F_sil = max(0.0, -F_prec[c_idx]), max(0.0, -F_prec[si_idx])

        S_sed = (F_carb * 0.100 / 2710.0 + F_sil * 0.060 / 2650.0) * ocean_water_per_area
        weathering_flux, w_diag = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, self.alpha, self.crust_production_rate, clog=False, sedimentation_rate=S_sed)

        F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass

        F_net = F_vol + F_prec + F_diss

        # F_max_negative = - Y[2:] / (1e7 * YR)
        # F_net = np.maximum(F_net, F_max_negative)

        # return derivatives

        dYdt = np.zeros_like(Y)

        dYdt[0] = (P_CO2_new - P_CO2) / tau_atm
        dYdt[1] = (P_H2O_new - P_H2O) / tau_atm

        dYdt[2:] = F_net
        self._F_net = F_net

        carbon_flux = dYdt[3]
        carbon = Y[3]
        calcite_SI = SI['Calcite']

        # print(f't = {t/YR:.4e} yr', end='\r')
        print(f't = {t/YR:.4e} yr  T = {T_surface:.1f}  P_CO2 = {P_CO2 / 1e5:.1e} bar  pH = {pH:.1f}  Calcite SI = {calcite_SI:.1f}  C flux = {(carbon_flux / carbon) * 1e9 * YR:.1e} / Gyr  ', end='\r')
        # print(Y[2:])
        # print(dYdt[2:])

        return dYdt

    def time_evolve(self, t_end=2e9 * YR, jac_epsilon=0.01, convergence_rate = 0.01 / (1e9 * YR)):

        Y0 = np.zeros(elements.shape[0] + 2)

        Y0[0] = 1000
        Y0[1] = 1000

        P_CO2_max = 1e5
        T_runaway = 350
        T_snowball = 260
        SI_acid_threshold = -2.5   # calcite SI below which the ocean can't buffer CO2

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

        def event_acid_ocean(t, Y):
            if Y[0] < P_CO2_max:
                return 1.0
            return getattr(self, '_SI', {}).get('Calcite', 0.0) - SI_acid_threshold
        event_acid_ocean.terminal = True

        atol = np.ones_like(Y0) * 1e-6
        atol[0] = 1.0   # P_CO2 in Pa
        atol[1] = 1.0   # P_H2O in Pa

        def event_converged(t, Y):
            F_net = getattr(self, '_F_net', None)
            if F_net is None:
                return 1.0
            b = np.maximum(Y[2:], atol[2:])
            max_fractional_rate = np.max(np.abs(F_net) / b)
            return max_fractional_rate - convergence_rate
        event_converged.terminal = True

        N = len(Y0)

        def macro_jacobian(t, y):

            jac = np.zeros((N, N))

            # Per-variable absolute floors, set to roughly the solver atol / 10 for
            # pressures and just above the trace_approximation threshold (1e-9 mol/kgw)
            # for ocean concentrations.  A uniform floor of 1e-5 was far too large
            # for trace elements (Al, Fe ~ 1e-8): it caused y_minus to go negative and
            # get clipped to 1e-9, making the finite difference wildly inaccurate.
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

                f_plus = self.dY_dt(t, y_plus)
                f_minus = self.dY_dt(t, y_minus)

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
            events=[event_snowball, event_hothouse, event_acid_ocean],
        )

        end = time.time()

        print()

        event_names = ['snowball', 'hothouse', 'acid_ocean']
        for name, t_ev in zip(event_names, sol.t_events):
            if len(t_ev) > 0:
                print(f'Terminated: {name} at t = {t_ev[0]/YR:.3e} yr')

        print(f'Simulation time: {end - start:.0f} s')
        print(f'Y: {sol.y[2:, -1]} mol/kgw')

        t = sol.t / YR

        dt = np.gradient(t)

        print('FLUXES (fractonal change per Gyr)')

        for i in range(len(sol.y[2:, 0])):
            b = sol.y[i+2, :]
            db = np.gradient(b)
            flux = ((db / dt) / b) * 1e9

            print(flux[-20:-1])
            
            plt.plot(t, flux)

        plt.axhline(0, color='black', linestyle='--')
        plt.axhspan(-1e-2, 1e-2, color='blue', alpha=0.2)
        # plt.yscale('log')
        plt.yscale('symlog', linthresh=1e-3)
        plt.show()


if __name__ == '__main__':

    BACKGROUND_PRESSURE = 1e5   # Pa (1 bar)
    OCEAN_DEPTH = 3000          # m
    TECTONICS = 1.0

    instellation = [0.8, 1.0, 1.2]

    for s in instellation:
        print(f's = {s}')
        p1 = Planet(M_EARTH, R_EARTH, BACKGROUND_PRESSURE, s, 1.0, 3000, name=f'test_s_{s}')
        p1.time_evolve()
        print('')