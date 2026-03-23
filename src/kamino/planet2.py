import numpy as np
import numpy.typing as npt
from scipy.optimize import newton, bisect
from tqdm import tqdm
import matplotlib.pyplot as plt

from kamino.constants import *
from kamino.kamino_chem.ocean_chemistry import *
from kamino.speedy_climate.clima_interpolator import get_T_surface
from kamino.utils import *

T_min = 275
T_max = 350

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

        self.alpha = 0.5

        self.tidally_locked = False

        land_fraction = 0
        ocean_albedo = 0.3
        land_albedo = 0.3

        self.instellation = instellation * SOLAR_CONSTANT
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

    def solve_climate_from_chemistry(self, b_ocean: npt.NDArray[np.float64], T_init: float=288) -> tuple[float, float]:

        def T_s_residual(T_guess):
            T_guess = np.clip(T_guess, 150, 500)
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

    def dYdt(self, t, Y):

        T = Y[0]
        P_CO2 = Y[1]
        b_ocean = Y[2:]

        b_ocean = smooth_max(b_ocean, np.zeros_like(b_ocean))

        T_new, P_CO2_new = self.solve_climate_from_chemistry(b_ocean)

        tau_atm = 10 * YR

        dT_dt = (T_new - T) / tau_atm
        dP_CO2_dt = (P_CO2_new - P_CO2) / tau_atm

        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T_new, P_CO2_new)

        tau_prec = 1e4 * YR

        F_out = self.outgassing_flux / self.ocean_water_mass
        weathering_flux, _ = get_weathering_flux(P_pore, T_pore, P_CO2_new, b_ocean, self.alpha, self.crust_production_rate)
        F_diss = (weathering_flux * self.surface_area) / self.ocean_water_mass
        F_prec = get_precipitation_flux(P_pore, T_seafloor, b_ocean, precipitating_minerals=carbonate_minerals) / tau_prec

        F = F_out + F_diss + F_prec

        dYdt = np.zeros_like(Y)

        dYdt[0] = dT_dt
        dYdt[1] = dP_CO2_dt
        dYdt[2:] = F

        return dYdt

    def solve_steady_state(self):

        def residual(Y):
            return self.dYdt(0, Y)

        Y_initial = np.ones(2 + elements.shape[0])
        Y_initial[0] = 288
        Y_initial[1] = 280e-6 * self.P_surface

        Y_steady = newton(residual, Y_initial, maxiter=100, tol=1e-6)
        return Y_steady
    
    def time_evolve(self, dt=5*YR):

        Y = np.zeros(2 + elements.shape[0])
        Y[0] = 288
        Y[1] = 280e-6 * self.P_surface

        t_current = 0
        n = 50000

        t = np.empty(n)
        P_CO2 = np.empty(n)
        T = np.empty(n)
        Alk = np.empty(n)

        for i in tqdm(range(n)):

            Y += self.dYdt(0, Y) * dt
            t_current += dt

            t[i] = t_current
            T[i] = Y[0]
            P_CO2[i] = Y[1]
            Alk[i] = Y[2]

        plt.plot(t / YR, T)
        plt.show()

        plt.plot(t / YR, P_CO2 / self.P_surface)
        plt.yscale('log')
        plt.show()

        plt.plot(t / YR, Alk)
        plt.yscale('log')
        plt.show()
        