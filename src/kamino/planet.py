import numpy as np
np.set_printoptions(precision=3)
from scipy.optimize import brentq, root, least_squares
from scipy.integrate import solve_ivp
import pandas as pd

from kamino.constants import *
from kamino.speedy_climate.emulator import climate_emulator
from kamino.ocean_chemistry.co2 import get_P_CO2
from kamino.ocean_circulation.analytic import get_T_ocean
from kamino.seafloor_weathering.weathering import get_weathering_rate
from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate

class planet:

    def __init__(self, P_surface, ocean_depth, instellation):
        
        # # MAIN VARIABLES

        # self.C_ocean: np.ndarray
        # self.C_pore: np.ndarray
        # self.Alk_ocean: np.ndarray
        # self.Alk_pore: np.ndarray
        # self.Ca_ocean: np.ndarray
        # self.Ca_pore: np.ndarray

        # # SECONDARY VARIABLES

        # self.T: np.ndarray
        # self.P_CO2: np.ndarray

        # CONSTANTS

        self.mass: float = M_EARTH
        self.radius: float = R_EARTH
        self.gravity: float = (G * self.mass) / (self.radius ** 2)

        self.P_surface: float = P_surface
        
        self.ocean_depth: float = ocean_depth
        self.ocean_mass: float = self.ocean_depth * 4 * np.pi * self.radius ** 2

        self.outgassing: float = (0.0147 * 4 * np.pi * self.radius ** 2)

        self.pore_space_flux: float = 1e14
        self.pore_space_mass = PORE_DEPTH * 4 * np.pi * self.radius ** 2 * 1000 * POROSITY

        self.instellation: float = instellation * SOLAR_CONSTANT
        self.albedo: float = 0.3

        # CLIMATE EMULATOR

        self.climate_emulator = climate_emulator("earth_rapid_rotator", "helios_1000_runs_earth_rapid_rotator.csv")
        self.climate_emulator.make_temperature_pco2_interpolator(self.instellation, self.P_surface, self.albedo)

    def solve_climate(self, Alk, C, Ca):

        def T_s_residual(T_val):
            pco2 = np.maximum(get_P_CO2(self.P_surface, T_val, Alk, C, Ca), 1e-6)
            T_calc = self.climate_emulator.get_temperature_from_pco2(pco2)
            return T_val - T_calc
        
        T_s = float(brentq(T_s_residual, 273.2, 373.15)) # type: ignore
        pco2 = get_P_CO2(self.P_surface, T_s, Alk, C, Ca)

        return T_s, pco2

    def dY_dt(self, t, Y):

        Y = np.maximum(Y, 1e-9)
        Co, Cp, Ao, Ap, Cao, Cap = Y

        F_out = self.outgassing
        J = self.pore_space_flux
        Mo = self.ocean_mass
        Mp = self.pore_space_mass

        T_s, pco2 = self.solve_climate(Ao, Co, Cao)

        T_pore = get_T_ocean(T_s, self.ocean_depth) + 9

        x_CO2 = np.maximum(pco2 / self.P_surface, 1e-8)
        P_pore = self.P_surface + 1000 * self.gravity * self.ocean_depth

        F_diss = get_weathering_rate(P_pore, T_pore, x_CO2, 0.05, PORE_DEPTH, 50e6) * 4 * np.pi * self.radius ** 2 # mol / s
        F_prec_o = get_calcite_precipitation_rate(P_pore, T_pore, Ao, Co, Cao)[0] * Mo * YR
        F_prec_p = get_calcite_precipitation_rate(P_pore, T_pore, Ap, Cp, Cap)[0] * Mp * YR

        print(f'F_prec_o = {F_prec_o / Mo:.3e} mol/kgw/yr  F_prec_p = {F_prec_p / Mp:.3e} mol/kgw/yr')

        delta_C = Co - Cp
        delta_A = Ao - Ap
        delta_Ca = Cao - Cap

        dCo_dt = (- J * delta_C + F_out - F_prec_o) / Mo
        dCp_dt = (+ J * delta_C - F_prec_p) / Mp
        dAo_dt = (- J * delta_A - 2 * F_prec_o) / Mo
        dAp_dt = (+ J * delta_A - 2 * F_prec_p + 2 * F_diss) / Mp
        dCao_dt = (- J * delta_Ca - F_prec_o) / Mo
        dCap_dt = (+ J * delta_Ca - F_prec_p + 0.5 * F_diss) / Mp

        print(f't = {t:.3e} yr  Y = {Y} mol / kgw  T_s = {T_s:.0f} K  P_CO2 = {pco2:.3e} Pa (x_CO2 = {x_CO2:.2e})')

        return [float(dCo_dt), float(dCp_dt), float(dAo_dt), float(dAp_dt), float(dCao_dt), float(dCap_dt)]
    
    def run_simulation(self, t_end):
        
        Y0 = [0, 0, 0, 0, 0, 0]

        sol = solve_ivp(
            self.dY_dt,
            (0, t_end),
            Y0,
            method='Radau',
            atol=1e-7,
            rtol=1e-4,
            max_step=10.0
        )

        results = pd.DataFrame(sol.y.T, columns=['C_ocean', 'C_pore', 'Alk_ocean', 'Alk_pore', 'Ca_ocean', 'Ca_pore'])
        results['time'] = sol.t

        T_s = []
        pco2 = []

        for Co, Ao, Cao in zip(results['C_ocean'], results['Alk_ocean'], results['Ca_ocean']):
            T_s_val, pco2_val = self.solve_climate(Ao, Co, Cao)
            T_s.append(T_s_val)
            pco2.append(pco2_val)

        results['T_surface'] = T_s
        results['P_CO2'] = pco2

        return results
    
    def find_steady_state(self):

        initial_results = self.run_simulation(1000 * YR)
        initial_Y_guess = initial_results.iloc[-1][['C_ocean', 'C_pore', 'Alk_ocean', 'Alk_pore', 'Ca_ocean', 'Ca_pore']].values

        target_function = lambda Y: self.dY_dt(0, Y)

        # solution = root(target_function, initial_Y_guess, method='hybr')
        solution = least_squares(target_function, initial_Y_guess, bounds=(0, np.inf))
        Co, Cp, Ao, Ap, Cao, Cap = solution.x
        
        print(solution)

        T_s, pco2 = self.solve_climate(Ao, Co, Cao)

        print(T_s)
        print(pco2)