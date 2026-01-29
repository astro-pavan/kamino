import numpy as np
np.set_printoptions(precision=1)
from scipy.optimize import brentq, root, least_squares, newton, basinhopping, direct, differential_evolution
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt

import time

from kamino.constants import *
from kamino.speedy_climate.emulator import climate_emulator, august_roche_magnus_formula
# from kamino.speedy_climate.simple import get_T_surface
from kamino.speedy_climate.simple_v2 import get_T_surface
from kamino.ocean_chemistry.co2 import get_P_CO2, get_P_CO2_v2
from kamino.ocean_circulation.analytic import get_T_ocean, get_T_ocean_KT18
from kamino.seafloor_weathering.weathering import get_weathering_rate, get_weathering_rate_old, get_land_weathering_rate_old
from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate
from kamino.ocean_chemistry.aqueous_geochemistry import PHREEQCError
from kamino.utils import *

iter = 0

class planet:

    def __init__(self, P_surface, ocean_depth, instellation):

        # CONSTANTS

        self.mass: float = M_EARTH
        self.radius: float = R_EARTH
        self.gravity: float = (G * self.mass) / (self.radius ** 2)

        self.P_surface: float = P_surface
        
        self.ocean_depth: float = ocean_depth
        self.ocean_mass: float = self.ocean_depth * 4 * np.pi * self.radius ** 2 * 1000

        self.outgassing: float = (0.0147 * 4 * np.pi * self.radius ** 2) / 2

        self.pore_space_flux: float = 1e14
        self.pore_space_mass = PORE_DEPTH * 4 * np.pi * self.radius ** 2 * 1000 * POROSITY

        self.instellation: float = instellation * SOLAR_CONSTANT
        self.albedo: float = 0.3

        # CLIMATE EMULATOR

        # self.climate_emulator = climate_emulator("earth_rapid_rotator", "helios_1000_runs_earth_rapid_rotator.csv")
        # self.climate_emulator.make_temperature_pco2_interpolator(self.instellation, self.P_surface)

        self.use_KT18_weathering = False

        self.last_t = 0
        self.last_T_s = 288
        self.debug_counter = 0

    def temperature_profile(self, latitude, T_avg, delta_T):
        return T_avg - delta_T * 0.5 * (3 * np.sin(latitude) ** 2 - 1)

    def solve_climate_from_chemistry(self, t, Alk, C, Ca, T_init=288):

        try:

            def T_s_residual(T_val):
                ph2o = august_roche_magnus_formula(T_val) * 0.77
                pco2 = get_P_CO2(self.P_surface, T_val, Alk, C, Ca)
                T_calc = get_T_surface(self.instellation, pco2, ph2o, 0.5)
                return T_val - T_calc
            
            T_s, _ = newton(T_s_residual, T_init, full_output=True, disp=False)
            
            pco2 = get_P_CO2(self.P_surface, T_s, Alk, C)
            ph2o = august_roche_magnus_formula(T_s) * 0.77

            assert ~np.isnan(T_s)

            return T_s, pco2, ph2o

        except RuntimeError:

            print(f'{Alk}, {C}, {Ca}')

            raise RuntimeError
        
    def solve_climate_from_CO2(self, P_CO2, T_init=288):

        def T_s_residual(T_val):
            P_H2O = august_roche_magnus_formula(T_val) * 0.77
            T_calc = get_T_surface(self.instellation, P_CO2, P_H2O, 0.5)
            return T_val - T_calc
        
        T_s, _ = newton(T_s_residual, T_init, full_output=True, disp=False)
        P_H2O = august_roche_magnus_formula(T_s) * 0.77

        return T_s, P_H2O

    def dY_dt(self, t, Y, verbose=True, T_surface=None, P_CO2=None):

        dt = t - self.last_t
        self.last_t = t

        self.debug_counter += 1
        should_print = self.debug_counter % 100

        calc_time = time.time()

        # keeps Y above 1e-9 smoothly
        Y_calc = smooth_max(Y, 1e-9)

        Co, Cp, Ao, Ap, Cao, Cap = Y_calc

        F_out = self.outgassing
        J = self.pore_space_flux
        Mo = self.ocean_mass
        Mp = self.pore_space_mass

        climate_calc_time = time.time()
        if P_CO2 is None or T_surface is None:
            T_s, pco2, ph2o = self.solve_climate_from_chemistry(t, Ao, Co, Cao)
        else:
            T_s = T_surface
            pco2 = P_CO2
            ph2o = august_roche_magnus_formula(T_s)
        climate_calc_time = time.time() - climate_calc_time

        T_seafloor = get_T_ocean(T_s, self.ocean_depth)
        T_seafloor = get_T_ocean_KT18(T_s)
        T_seafloor = smooth_max(T_seafloor, 273.5)
        T_pore = T_seafloor + 9

        # T_seafloor = self.temperature_profile(90, T_s, 30)
        # T_seafloor = smooth_max(T_seafloor, 273.5)
        # T_pore = T_seafloor + 9

        # T_seafloor = self.temperature_profile(0, T_s, 30)
        # T_seafloor = smooth_max(T_seafloor, 273.5)
        # T_pore = T_seafloor + 9

        x_CO2 = pco2 / (self.P_surface + pco2 + ph2o)

        P_pore = (self.P_surface + pco2 + ph2o) + 1000 * self.gravity * self.ocean_depth

        weathering_calc_time = time.time() 
        if self.use_KT18_weathering:
            F_diss = get_weathering_rate_old(P_pore, T_pore, x_CO2) * (self.radius / R_EARTH) ** 2 # mol / s
        else:
            F_diss = get_weathering_rate(P_pore, T_pore, x_CO2, 0.05, PORE_DEPTH, 50e6) * (4 * np.pi * self.radius ** 2) # mol / s
        weathering_calc_time = time.time() - weathering_calc_time

        precipitation_calc_time = time.time()

        # rate_o, SI_o = get_calcite_precipitation_rate(P_pore, T_seafloor, Ao, Co, Cao)
        # rate_p, SI_p = get_calcite_precipitation_rate(P_pore, T_pore, Ap, Cp, Cap)

        # F_prec_o = rate_o * Mo * YR
        # F_prec_p = rate_p * Mp * YR

        # tau = 1 # in yrs
        # F_prec_o_max = (Co * Mo) / tau
        # F_prec_p_max = (Cp * Mp) / tau

        F_prec_o = 0.5 * F_diss
        F_prec_p = 0.5 * F_diss
        SI_o = 0
        SI_p = 0

        precipitation_calc_time = time.time() - precipitation_calc_time

        delta_C = Co - Cp
        delta_A = Ao - Ap
        delta_Ca = Cao - Cap

        dCo_dt = (- J * delta_C + F_out - F_prec_o) / Mo
        dCp_dt = (+ J * delta_C - F_prec_p) / Mp
        dAo_dt = (- J * delta_A - 2 * F_prec_o) / Mo
        dAp_dt = (+ J * delta_A - 2 * F_prec_p + 2 * F_diss) / Mp
        dCao_dt = (- J * delta_Ca - F_prec_o) / Mo
        dCap_dt = (+ J * delta_Ca - F_prec_p + 0.5 * F_diss) / Mp

        dYdt = np.array([float(dCo_dt), float(dCp_dt), float(dAo_dt), float(dAp_dt), float(dCao_dt), float(dCap_dt)])

        calc_time = time.time() - calc_time

        flux_diff = J * delta_C / Mp
        flux_prec = -F_prec_p / Mp
        flux_diss = F_diss / Mp
        
        dCp_dt = flux_diff + flux_prec + flux_diss

        if should_print and verbose:
            print(f't = {t:.1e} yr  Y = {Y_calc[0]:.1e}, {Y_calc[1]:.1e}, {Y_calc[2]:.1e}, {Y_calc[1]:.1e}, {Y_calc[4]:.1e}, {Y_calc[5]:.1e} mol/kgw  T_s = {T_s:.0f} K  T_f = {T_pore:.0f} K  P_CO2 = {pco2:.1e} Pa  Calcite SI = {SI_o:.3f}, {SI_p:.3f}')
            # print(f'dY/dt = {dYdt} mol/kgw/s')
            # print(f"\n--- DEBUG t={t:.2e} ---")
            # print(f"Cp: {Cp:.6e}")
            # print(f"dCp_dt: {dCp_dt:.2e} = (Diff: {flux_diff:.2e}) + (Prec: {flux_prec:.2e})")
            # print(f"Prec Rate (raw): {F_prec_p:.2e} | Max Allowed: {F_prec_p_max:.2e}")
            # print(f"T_pore: {T_pore:.4f} K")
            # print(f'Step time: {calc_time:.1e} s (Climate: {climate_calc_time:.1e} s  Weathering: {weathering_calc_time:.1e} s  Precipitation: {precipitation_calc_time:.1e} s)')

        return dYdt
    
    def run_simulation(self, t_end, make_plots, Y0):
        
        # Y0 = [0, 0, 0, 0, 0, 0]

        sol = solve_ivp(
            self.dY_dt,
            (0, t_end),
            Y0,
            method='LSODA',
            atol=1e-10,
            rtol=1e-2,
            first_step=1
        )

        res = np.maximum(sol.y.T, 1e-9)
        results = pd.DataFrame(res, columns=['C_ocean', 'C_pore', 'Alk_ocean', 'Alk_pore', 'Ca_ocean', 'Ca_pore'])
        results['time'] = sol.t

        T_s = []
        T_f = []
        pco2 = []
        dCodt = []
        dAodt = []
        dCaodt = []
        F_dissolution = []
        F_precipiation = []
        F_max = []
        SI = []

        Mo = self.ocean_mass
        Mp = self.pore_space_mass

        for t, Co, Ao, Cao, Cp, Ap, Cap in zip(results['time'], results['C_ocean'], results['Alk_ocean'], results['Ca_ocean'], results['C_pore'], results['Alk_pore'], results['Ca_pore']):
            
            T_s_val, pco2_val, ph2o_val = self.solve_climate_from_chemistry(t, Ao, Co, Cao)
            T_s.append(T_s_val)
            pco2.append(pco2_val)

            dYdt = self.dY_dt(0, np.array([Co, Cp, Ao, Ap, Cao, Cap]), verbose=False)
            dCodt.append(dYdt[0])
            dAodt.append(dYdt[2])
            dCaodt.append(dYdt[4])

            T_seafloor = get_T_ocean_KT18(T_s_val)
            T_seafloor = smooth_max(T_seafloor, 273.5)
            T_pore = T_seafloor + 9

            # T_seafloor = self.temperature_profile(90, T_s, 30)
            # T_seafloor = smooth_max(T_seafloor, 273.5)
            # T_pore = T_seafloor + 9

            # T_seafloor = self.temperature_profile(0, T_s, 30)
            # T_seafloor = smooth_max(T_seafloor, 273.5)
            # T_pore = T_seafloor + 9

            T_f.append(T_seafloor)

            x_CO2 = pco2_val / (self.P_surface + pco2_val + ph2o_val)
            P_pore = (self.P_surface + pco2_val + ph2o_val) + 1000 * self.gravity * self.ocean_depth

            F_diss_val = get_weathering_rate(P_pore, T_pore, x_CO2, 0.05, PORE_DEPTH, 50e6) * (4 * np.pi * self.radius ** 2)

            rate_o, SI_o = get_calcite_precipitation_rate(P_pore, T_seafloor, Ao, Co, Cao)
            rate_p, SI_p = get_calcite_precipitation_rate(P_pore, T_pore, Ap, Cp, Cap)
            SI.append(SI_o)
            F_prec_o = rate_o * Mo * YR
            F_prec_p = rate_p * Mp * YR
            tau = 1 # in yrs
            F_prec_o_max = (Co * Mo) / tau
            # F_prec_o = smooth_min(F_prec_o, F_prec_o_max)

            F_max = F_prec_o_max
            F_precipiation.append(F_prec_o + F_prec_p)
            F_dissolution.append(F_diss_val)


        results['T_surface'] = T_s
        results['T_seafloor'] = T_f
        results['P_CO2'] = pco2
        results['F_diss'] = F_dissolution
        results['F_max'] = F_max
        results['F_prec'] = F_precipiation
        results['SI'] = SI

        if make_plots:

            fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20, 20))
            plt.subplots_adjust(hspace=0)

            ax1.plot(results['time'], results['C_ocean'], label='C_ocean', color='black')
            ax1.plot(results['time'], results['Alk_ocean'], label='Alk_ocean', color='orange')
            ax1.plot(results['time'], results['Ca_ocean'], label='Ca_ocean', color='green')

            ax12 = ax1.twinx()
            ax12.plot(results['time'], np.array(dCodt), label='dC_ocean/dt', color='black', linestyle='--')
            ax12.plot(results['time'], np.array(dAodt), label='dA_ocean/dt', color='orange', linestyle='--')
            ax12.plot(results['time'], np.array(dCaodt), label='dCa_ocean/dt', color='green', linestyle='--')

            ax1.set_ylabel('Concentration (mol/kgw)')
            ax1.set_yscale('log')
            ax12.set_ylabel('Rate of change of concentration (mol/kgw/s)')
            ax12.set_yscale('log')

            lines1, labels1 = ax1.get_legend_handles_labels()
            lines12, labels12 = ax12.get_legend_handles_labels()
            ax1.legend(lines1 + lines12, labels1 + labels12, loc='lower right')

            # ax1.legend(loc='lower right')

            ax22 = ax2.twinx()
            ax22.plot(results['time'], results['T_surface'], 'k--', label='T_surface')
            ax22.plot(results['time'], results['T_seafloor'], 'b--', label='T_seafloor')
            ax2.plot(results['time'], results['P_CO2'], 'r-', label='P_CO2')

            ax2.set_ylabel('$P_{CO2}$ (Pa)')
            ax2.set_yscale('log')

            ax22.set_ylabel('T (K)')
            ax22.set_ylim([250, 340])

            # Combine legends from ax2 and ax3
            lines2, labels2 = ax2.get_legend_handles_labels()
            lines22, labels22 = ax22.get_legend_handles_labels()
            ax2.legend(lines2 + lines22, labels2 + labels22, loc='lower right')

            # 3. Set the x-axis label ONLY on the bottom plot
            ax3.set_xlabel('Time (yr)')

            ax3.plot(results['time'], results['F_diss'], 'b-', label='Dissolution rate')
            ax3.plot(results['time'], results['F_prec'], 'r-', label='Precipitation rate')
            ax3.axhline(self.outgassing, color='green', label='Outgassing rate')
            # ax3.plot(results['time'], results['F_max'], 'm', label='Max precipitaion rate')

            ax32 = ax3.twinx()

            ax32.plot(results['time'], results['SI'], 'k--', label='Calcite SI')
            ax32.set_ylabel('Saturation Index')

            ax3.set_ylabel('Rate (mol/yr)')
            ax3.set_yscale('log')

            lines3, labels3 = ax3.get_legend_handles_labels()
            lines32, labels32 = ax32.get_legend_handles_labels()
            ax3.legend(lines3 + lines32, labels3 + labels32, loc='upper right')

            # Optional: Ensure tick labels on the top plot are definitely off (sharex usually handles this)
            ax1.tick_params(labelbottom=False)

            # plt.tight_layout()
            plt.savefig('Evolution.png')
            plt.close()

        return results
    
    def find_steady_state(self, t_init: float=1000):

        print(f'Initializing planet by evolving for {t_init} yrs...')
        initial_results = self.run_simulation(t_init, True, [0.001, 0.001, 0.001, 0.001, 0.001, 0.001])
        initial_Y_guess = initial_results.iloc[-1][['C_ocean', 'C_pore', 'Alk_ocean', 'Alk_pore', 'Ca_ocean', 'Ca_pore']].values

        target_function = lambda Y: self.dY_dt(0, Y)

        print(f'Initial state: {initial_Y_guess} mol/kgw')

        # solution = root(target_function, initial_Y_guess, method='hybr')
        print('Finding steady state solution...')
        solution = least_squares(target_function, initial_Y_guess, bounds=(0, np.inf))
        Co, Cp, Ao, Ap, Cao, Cap = solution.x
        
        print('Solution found: ')
        print(solution)

        T_s, pco2, ph2o = self.solve_climate_from_chemistry(1e6, Ao, Co, Cao)

        jacobian = solution.jac

        eigval, _ = np.linalg.eig(jacobian)

        print(f'Surface Temperature  : {T_s:.0f} K')
        print(f'P_CO2                : {pco2:.1e} Pa')
        print(f'Jacobian Eigenvalues : {eigval}')

        stable = True

        if np.all(eigval <= 0):
            print(f'Climate is stable')
        else:
            instability_timescale = 1 / np.max(eigval)
            if instability_timescale > 1e15:
                print(f'Climate is stable')
            else:
                print(f'Climate will become unstable in {instability_timescale:.1e} yrs')
                stable = False

        return T_s, stable, Co
    
    def find_steady_state_simple(self):

        # pco2_range = np.logspace(-2, 4)
        # r = []
        # for pco2 in pco2_range:
        #     r.append(self.solve_climate_from_CO2(pco2)[0])

        # plt.plot(pco2_range, r)
        # plt.scatter(40, 288)
        # plt.xscale('log')
        # plt.show()
        # plt.close()

        def H(x, k=1e20):
            return 0.5 *((k*x) / (np.sqrt(1 + (k * x) ** 2)) + 1)
                                    
        def target_function_climate(pco2):

            T_s, ph2o = self.solve_climate_from_CO2(pco2)

            x_CO2 = pco2 / (self.P_surface + pco2 + ph2o)
            P_pore = (self.P_surface + pco2 + ph2o) + 1000 * self.gravity * self.ocean_depth

            T_seafloor = get_T_ocean_KT18(T_s)
            T_seafloor = smooth_max(T_seafloor, 273.5)

            # T_seafloor_cold = self.temperature_profile(90, T_s, 30)
            # T_seafloor_cold = smooth_max(T_seafloor_cold, 273.5)

            # T_seafloor_hot = self.temperature_profile(0, T_s, 30)
            # T_seafloor_hot = smooth_max(T_seafloor_hot, 273.5)

            # T_threshold = 300
            # T_seafloor = T_seafloor_cold * H(T_threshold - T_s) + T_seafloor_hot * H(T_s - T_threshold)
            
            T_pore = T_seafloor + 9

            # print(f'P_CO2 = {float(pco2):.2e} Pa')
            weathering = get_weathering_rate(P_pore, T_pore, x_CO2, 0.05, PORE_DEPTH, 50e6) * 4 * np.pi * self.radius ** 2
            # weathering = get_land_weathering_rate_old(self.P_surface, T_s, x_CO2) * 4 * np.pi * self.radius ** 2
            # weathering = get_weathering_rate_old(self.P_surface, T_s - 14, x_CO2) * 4 * np.pi * self.radius ** 2

            res = 0.5 * weathering - self.outgassing
            print(f'P_CO2 = {float(pco2):.2e} Pa  T_pore = {T_pore:.0f} K  Residual = {float(res):.2e}')

            return res / self.outgassing
        
        pco2_range = np.logspace(-2, 5)
        r = []
        for pco2 in pco2_range:
            r.append((target_function_climate(pco2) * self.outgassing + self.outgassing) * 2)

        plt.plot(pco2_range, r)
        plt.axhline(0)
        plt.xscale('log')
        plt.show()
        plt.close()
        
        print('Solving climate state...')
        sol_climate = least_squares(target_function_climate, 40, bounds=(1e-2, 1e5))
        # sol_climate = direct(target_function_climate, [(0, 1e4)])
        print('Solved.')

        print(sol_climate)

        P_CO2 = float(sol_climate.x)
        T_s, P_H2O = self.solve_climate_from_CO2(P_CO2)

        x_CO2 = P_CO2 / (self.P_surface + P_CO2 + P_H2O)
        P_pore = (self.P_surface + P_CO2 + P_CO2) + 1000 * self.gravity * self.ocean_depth

        T_seafloor = get_T_ocean_KT18(T_s)
        T_seafloor = smooth_max(T_seafloor, 273.5)
        T_pore = T_seafloor + 9

        weathering = get_weathering_rate(P_pore, T_pore, x_CO2, 0.05, PORE_DEPTH, 50e6) * 4 * np.pi * self.radius ** 2

        print(f'T_s = {T_s:.0f}')

        target_function_chemistry = lambda Y: self.dY_dt(0, Y, T_surface=T_s, P_CO2=P_CO2)

        def get_chemistry(Y):
            
            logC, logAlk, logCa = Y

            C, Alk, Ca = 10 ** float(logC), 10 ** float(logAlk), 10 ** float(logCa)

            rate_o, SI_o = get_calcite_precipitation_rate(P_pore, T_seafloor, Alk, C, Ca)

            precipitation = rate_o * self.ocean_mass * YR

            res = precipitation - 0.5 * weathering

            print(f'Y = {Y[0]:.3f} {Y[1]:.3f} {Y[2]:.3f} mol/kgw  SI = {SI_o:.2f}  Residual = {res / weathering}')

            return res / weathering


        init = 0.001
        x0 = [init, init, init, init, init, init]

        lower_bounds = [-7, -7, -7, -7, -7, -7]
        upper_bounds = [0, 0, 0, 0, 0, 0]

        bounds = [(-7, 0), (-7, 0), (-7, 0)]

        sol_chemistry = least_squares(
            target_function_chemistry, 
            x0, 
        )

        # sol_chemistry = differential_evolution(get_chemistry, bounds)

        Y_chem = sol_chemistry.x

        return T_s, P_CO2
        
