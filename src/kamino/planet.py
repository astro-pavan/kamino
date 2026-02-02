import numpy as np
np.set_printoptions(precision=1)
from scipy.optimize import least_squares, newton, bisect
from scipy.integrate import solve_ivp
import pandas as pd
import matplotlib.pyplot as plt

import time

from kamino.constants import *
# from kamino.speedy_climate.simple_v2 import get_T_surface
from kamino.speedy_climate.clima_interpolator import get_T_surface
from kamino.ocean_chemistry.co2 import get_P_CO2
from kamino.ocean_circulation.analytic import get_T_ocean, get_T_ocean_KT18
from kamino.seafloor_weathering.weathering import *
from kamino.ocean_chemistry.precipitation import get_calcite_precipitation_rate
from kamino.ocean_chemistry.aqueous_geochemistry import PHREEQCError
from kamino.ocean_chemistry.fast_chemistry import get_calcite_data_fast
from kamino.ocean_chemistry.precipitation_interpolator import *
from kamino.utils import *

iter = 0

class planet:

    def __init__(
            self,
            P_surface: float,
            ocean_depth: float,
            instellation: float,
            outgassing_rate: float,
            hydrothermal_flow_rate: float,
            seafloor_age: float,
            flow_path_length: float,
            tidally_locked: bool=False,
            albedo: float=0.05,
            use_KT18_weathering: bool=False,
            use_WHAK_weathering: bool=False,
            use_MAC_weathering: bool=False
            ):

        # CONSTANTS

        self.mass: float = M_EARTH
        self.radius: float = R_EARTH
        self.gravity: float = (G * self.mass) / (self.radius ** 2)
        self.surface_area: float = 4 * np.pi * self.radius ** 2

        self.tidally_locked = tidally_locked

        self.P_surface: float = P_surface
        
        self.ocean_depth: float = ocean_depth
        self.ocean_mass: float = self.ocean_depth * self.surface_area * 1000

        self.outgassing: float = 0.0147 * self.surface_area * outgassing_rate

        self.pore_space_flux: float = 1e14
        self.pore_space_mass = PORE_DEPTH * self.surface_area * 1000 * POROSITY

        self.instellation: float = instellation * SOLAR_CONSTANT
        self.albedo: float = albedo

        self.hydrothermal_flow_rate = hydrothermal_flow_rate
        self.seafloor_age = seafloor_age
        self.flow_path_length = flow_path_length

        self.use_KT18_weathering = use_KT18_weathering
        self.use_WHAK_weathering = use_WHAK_weathering
        self.use_MAC_weathering = use_MAC_weathering

        self.last_t = 0
        self.last_T_s = 288
        self.debug_counter = 0

    def temperature_profile(self, latitude: float, T_avg: float, delta_T: float):
        return T_avg - delta_T * 0.5 * (3 * np.sin(latitude) ** 2 - 1)

    def solve_climate_from_chemistry(self, t: float, Alk: float, C: float, Ca: float, T_init: float=288) -> tuple[float, float, float]:

        try:

            def T_s_residual(T_val):
                pco2 = get_P_CO2(self.P_surface, T_val, Alk, C, Ca)
                T_calc = get_T_surface(self.instellation, pco2, self.albedo, tidally_locked=self.tidally_locked)
                return T_val - T_calc
            
            T_s, r = newton(T_s_residual, T_init, full_output=True, disp=False)
            
            pco2 = get_P_CO2(self.P_surface, T_s, Alk, C, Ca)
            ph2o = august_roche_magnus_formula(T_s) * 0.5

            assert ~np.isnan(T_s)

            return float(T_s), pco2, ph2o

        except RuntimeError as e:

            print(e)

            raise RuntimeError
        
    def solve_climate_from_CO2(self, P_CO2: float, T_init: float=288) -> tuple[float, float]:

        def T_s_residual(T_val):
            T_calc = get_T_surface(self.instellation, P_CO2, self.albedo, tidally_locked=self.tidally_locked)
            return T_val - T_calc
        
        T_s = float(newton(T_s_residual, T_init))
        P_H2O = august_roche_magnus_formula(T_s) * 0.5

        return T_s, P_H2O

    def dY_dt(self, t, Y, full_precipitation_calculation=True, verbose=True):

        should_print = (self.debug_counter % 20) == 0
        self.debug_counter += 1

        # keeps Y above 1e-9 smoothly
        Y_calc = np.maximum(Y, 1e-9)

        # T, P_CO2, Co, Cp, Ao, Ap, Cao, Cap = Y_calc
        T, P_CO2, Co, Ao, Cao = Y_calc

        F_out = self.outgassing
        J = self.pore_space_flux
        Mo = self.ocean_mass
        Mp = self.pore_space_mass

        T_new, P_CO2_new, P_H2O = self.solve_climate_from_chemistry(t, Ao, Co, Cao)
        F_diss, F_prec_o, T_pore, SI_o = self.get_fluxes(t, Y, full_precipitation_calculation)

        # Calculate derivatives

        tau_atm = 10

        dT_dt = (T_new - T) / tau_atm
        dP_CO2_dt = (P_CO2_new - P_CO2) / tau_atm

        # delta_C = Co - Cp
        # delta_A = Ao - Ap
        # delta_Ca = Cao - Cap

        # dCo_dt = (- J * delta_C + F_out - F_prec_o) / Mo
        # dCp_dt = (+ J * delta_C - F_prec_p) / Mp
        # dAo_dt = (- J * delta_A - 2 * F_prec_o) / Mo
        # dAp_dt = (+ J * delta_A - 2 * F_prec_p + 2 * F_diss) / Mp
        # dCao_dt = (- J * delta_Ca - F_prec_o) / Mo
        # dCap_dt = (+ J * delta_Ca - F_prec_p + 0.5 * F_diss) / Mp

        # dYdt = np.array([float(dT_dt), float(dP_CO2_dt), float(dCo_dt), float(dCp_dt), float(dAo_dt), float(dAp_dt), float(dCao_dt), float(dCap_dt)])

        dCo_dt = (F_out - F_prec_o) / Mo
        dAo_dt = (- 2 * F_prec_o + 2 * F_diss) / Mo
        dCao_dt = (- F_prec_o + 1.0 * F_diss) / Mo

        dYdt = np.array([float(dT_dt), float(dP_CO2_dt), float(dCo_dt), float(dAo_dt), float(dCao_dt)])

        if should_print and verbose:
            print(f't = {t:.1e} yr  Y = {Y[2]:.2e} {Y[3]:.2e} {Y[4]:.2e} mol/kgw  T_s = {T_new:.0f} K  T_f = {T_pore:.0f} K  P_CO2 = {P_CO2_new:.1e} Pa  F_prec = {F_prec_o / F_out:.2e}  F_diss = {F_diss / F_out:.2e}  SI = {SI_o:.4f}')

        return dYdt
    
    def get_fluxes(self, t, Y, full_precipitation_calculation=True):

        Y_calc = smooth_max(np.array(Y), 1e-9)
        # T, P_CO2, Co, Cp, Ao, Ap, Cao, Cap = Y_calc
        T, P_CO2, Co, Ao, Cao = Y_calc

        P_H2O = august_roche_magnus_formula(T) * 0.5

        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T, P_CO2)

        x_CO2 = P_CO2 / (self.P_surface + P_CO2 + P_H2O)

        Mo = self.ocean_mass
        Mp = self.pore_space_mass

        F_diss = self.get_weathering(T, P_CO2)

        if full_precipitation_calculation:

            # rate_o, SI_o = get_calcite_precipitation_rate(P_pore, T_seafloor, Ao, Co, Cao)
            # rate_p, SI_p = get_calcite_precipitation_rate(P_pore, T_pore, Ap, Cp, Cap)

            rate_o, SI_o = get_calcite_data_interpolated(P_pore, T_seafloor, Ao, Co, Cao)
            # rate_p, SI_p = get_calcite_data_interpolated(P_pore, T_pore, Ap, Cp, Cap)

            kinetics_scaler = 1e-5
            smoothness = 0.01
            switch_o = smoothness * np.logaddexp(0, SI_o / smoothness)
            # switch_p = smoothness * np.logaddexp(0, SI_p / smoothness)
        
            rate_o *= kinetics_scaler
            # rate_p *= kinetics_scaler

            F_prec_o = rate_o * Mo * YR * switch_o
            # F_prec_p = rate_p * Mp * YR * switch_p

            tau_prec = 100

            F_Co_max = (Co * Mo) / tau_prec
            F_Ao_max = (Ao * Mo) / tau_prec
            F_Cao_max = (Cao * Mo) / tau_prec
            # F_prec_o_max = np.min([F_Co_max, F_Ao_max, F_Cao_max])

            # F_Cp_max = (Cp * Mp) / tau_prec
            # F_Ap_max = (Ap * Mp) / tau_prec
            # F_Cap_max = (Cap * Mp) / tau_prec
            # F_prec_p_max = np.min([F_Cp_max, F_Ap_max, F_Cap_max])

            # F_prec_o = smooth_min(F_prec_o_max, F_prec_o)
            F_prec_o = smooth_min(20 * self.outgassing, F_prec_o)
            # F_prec_p = smooth_min(F_prec_p_max, F_prec_p)

        else:

            F_prec_o = 0.5 * F_diss
            # F_prec_p = 0.5 * F_diss
        
        return F_diss, F_prec_o, T_pore, SI_o
    
    def get_weathering(self, T_s: float, P_CO2: float) -> float:

        P_H2O = august_roche_magnus_formula(T_s) * 0.5

        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T_s, P_CO2)

        x_CO2 = P_CO2 / (self.P_surface + P_CO2 + P_H2O)

        if self.use_MAC_weathering:
            weathering = get_weathering_rate_MAC(T_s, P_CO2) * self.surface_area
        elif self.use_WHAK_weathering:
            weathering = get_weathering_rate_WHAK(self.P_surface, T_s, x_CO2) * self.surface_area
        elif self.use_KT18_weathering:
            weathering = get_weathering_rate_KT18(self.P_surface, T_pore, x_CO2) * self.surface_area
        else:
            weathering = get_weathering_rate(P_pore, T_pore, x_CO2, self.hydrothermal_flow_rate, self.flow_path_length, self.seafloor_age) * self.surface_area

        return weathering

    def get_seafloor_properties(self, T_s: float, P_CO2: float) -> tuple[float, float, float]:

        P_H2O = august_roche_magnus_formula(T_s) * 0.5
        P_pore = (self.P_surface + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth

        T_seafloor = get_T_ocean_KT18(T_s)
        T_seafloor = smooth_max(T_seafloor, 273.5)
        T_pore = T_seafloor + 9

        return T_seafloor, T_pore, P_pore

    def run_simulation(self, t_end, make_plots, Y0):

        def check_equilibrium(t, Y):
            F_diss, F_prec_o, F_prec_p, T_pore = self.get_fluxes(t, Y)
            return F_prec_o - (0.5 * F_diss)
        
        check_equilibrium.terminal = True
        check_equilibrium.direction = 0

        sol1 = solve_ivp(
            self.dY_dt,
            (0, t_end),
            Y0,
            method='LSODA',
            atol=1e-10,
            rtol=1e-2,
            first_step=1,
            # events=check_equilibrium
        )

        # results = pd.DataFrame(res, columns=['T_surface', 'P_CO2',  'C_ocean', 'C_pore', 'Alk_ocean', 'Alk_pore', 'Ca_ocean', 'Ca_pore'])

        results = self.post_process_evolution(sol1.t, sol1.y)

        if make_plots:
            self.plot_evolution(results)

        return results
    
    def find_steady_state(self, t_init: float=1000, init_concentration: float=0.0001):

        print(f'Initializing planet by evolving for {t_init} yrs...')
        initial_results = self.run_simulation(t_init, True, [288, 1, init_concentration, init_concentration, init_concentration])
        initial_Y_guess = initial_results.iloc[-1][['T_surface', 'P_CO2', 'C_ocean', 'Alk_ocean', 'Ca_ocean']].values

        target_function = lambda Y: self.dY_dt(0, Y)

        print(f'Initial state: {initial_Y_guess} mol/kgw')

        print('Finding steady state solution...')
        solution = least_squares(target_function, initial_Y_guess, bounds=(0, np.inf))
        # T, P_CO2, Co, Cp, Ao, Ap, Cao, Cap = solution.x
        T, P_CO2, Co, Ao, Cao = solution.x
        
        print('Solution found: ')
        print(solution)

        jacobian = solution.jac

        eigval, _ = np.linalg.eig(jacobian)

        print(f'Surface Temperature  : {T:.0f} K')
        print(f'P_CO2                : {P_CO2:.1e} Pa')
        print(f'Chemistry            : {Co:.1e}, {Ao:.1e}, {Cao:.1e} mol/kgw')
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

        return T, stable, Co

    def find_steady_state_no_evolution(self, use_KT18_weathering=False, use_WHAK_weathering=False, use_MAC_weathering=False, diagnostic_plots: bool=False, solve_chemistry: bool=False) -> tuple[float, float, float, float, float, float]:
                                    
        def target_function_climate(pco2: float):
            T_s, _ = self.solve_climate_from_CO2(pco2)
            weathering = self.get_weathering(T_s, pco2)
            residual = weathering - self.outgassing
            return residual / self.outgassing
        
        if diagnostic_plots:
            pco2_range = np.logspace(-2, 5)
            r = []
            for pco2 in pco2_range:
                r.append(target_function_climate(pco2))

            plt.figure(figsize=(3,3))
            plt.title(f'S={self.instellation}')
            plt.plot(pco2_range, r)
            plt.axhline(0)
            plt.xscale('log')
            plt.show()
            plt.close()
        
        print('Solving climate state...')
        try:
            sol_climate = float(bisect(target_function_climate, 1e-2, 1e5)) # type: ignore
            print('Solved.')
        except ValueError:
            print('No solution.')
            return np.nan, np.nan, np.nan

        P_CO2 = float(sol_climate)
        T_s, P_H2O = self.solve_climate_from_CO2(P_CO2)
        T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T_s, P_CO2)

        if use_WHAK_weathering or use_MAC_weathering:
            T_weather = T_s
        else:
            T_weather = T_pore

        if solve_chemistry:

            k_conservative = -1.0e-4

            def target_function_chemistry(Y):
                Co, Ao = Y
                
                # We CALCULATE Calcium based on the conservative constant constraint
                # This reduces the degrees of freedom so the system is solvable.
                Cao = (Ao - k_conservative) / 2.0
                
                P_CO2_calc = get_P_CO2(self.P_surface, T_s, Ao, Co, Cao)
                res_pco2 = (P_CO2_calc - P_CO2) / P_CO2
                
                F_diss, F_prec_o, _, _  = self.get_fluxes(0, [T_s, P_CO2, Co, Ao, Cao])
                res_flux = (F_prec_o - self.outgassing) / self.outgassing

                return [res_pco2, res_flux]
            
            solution = least_squares(target_function_chemistry, [0.0001, 0.0001], bounds=(0, np.inf))
            Co, Ao = solution.x
            Cao = (Ao - k_conservative) / 2.0

            Y0 = [Co, Ao, Cao]

            def dY_dt_chemistry_only(t, Y):
                Co, Ao, Cao = Y
                dY = self.dY_dt(t, [T_s, P_CO2, Co, Ao, Cao])
                return dY[2], dY[3], dY[4]

            sol = solve_ivp(
                dY_dt_chemistry_only,
                (0, 1e5),
                Y0,
                method='LSODA',
                atol=1e-10,
                rtol=1e-2,
                first_step=1,
                # events=check_equilibrium
            )

            Y = sol.y[:, -1]
            Co, Ao, Cao = Y
        
        else:
            Co, Ao, Cao = np.nan, np.nan, np.nan

        print(f'Surface Temperature  : {T_s:.0f} K')
        print(f'P_CO2                : {P_CO2:.1e} Pa')
        print(f'Chemistry found      : {Co:.1e}, {Ao:.1e}, {Cao:.1e} mol/kgw')

        return T_s, P_CO2, Co, Ao, Cao, T_weather
    
    def stability_analysis(self, T_s, P_CO2, Co, Ao, Cao):

        # peturb instellation

        Y0 = T_s, P_CO2, Co, Ao, Cao

        t_switch_01 = 1e5
        t_switch_12 = 1e5
        t_end = 1e6

        sol0 = solve_ivp(
            self.dY_dt,
            (0, t_switch_01),
            Y0,
            method='LSODA',
            atol=1e-10,
            rtol=1e-2
        )

        Y1 = sol0.y[:, -1]

        self.instellation *= 1.1

        sol1 = solve_ivp(
            self.dY_dt,
            (t_switch_01, t_switch_01 + t_switch_12),
            Y1,
            method='LSODA',
            atol=1e-10,
            rtol=1e-2
        )

        # self.instellation /= 1.1

        Y2 = sol1.y[:, -1]

        sol2 = solve_ivp(
            self.dY_dt,
            (t_switch_01 + t_switch_12, t_end),
            Y2,
            method='LSODA',
            atol=1e-10,
            rtol=1e-2
        )

        sol_y = np.hstack([sol0.y, sol1.y, sol2.y])
        sol_t = np.hstack([sol0.t, sol1.t, sol2.t])

        results = self.post_process_evolution(sol_t, sol_y)

        self.plot_evolution(results)

    def post_process_evolution(self, sol_t, sol_y):
        
        res = np.maximum(sol_y.T, 1e-9)
        
        results = pd.DataFrame(res, columns=['T_surface', 'P_CO2',  'C_ocean', 'Alk_ocean', 'Ca_ocean'])
        results['time'] = sol_t

        T_f = []
        F_dissolution = []
        F_precipiation = []
        SI = []
        P_H2O_arr = []

        for t, T, P_CO2, Co, Ao, Cao in zip(results['time'], results['T_surface'], results['P_CO2'], results['C_ocean'], results['Alk_ocean'], results['Ca_ocean']):

            T_seafloor, T_pore, P_pore = self.get_seafloor_properties(T, P_CO2)
            F_diss, F_prec_o, T_pore, SI_o = self.get_fluxes(t, [T, P_CO2, Co, Ao, Cao])

            SI.append(SI_o)
            T_f.append(T_seafloor)
            F_precipiation.append(F_prec_o)
            F_dissolution.append(F_diss)
            P_H2O_arr.append(august_roche_magnus_formula(T) * 0.5)


        results['T_seafloor'] = T_f
        results['F_diss'] = F_dissolution
        results['F_prec'] = F_precipiation
        results['SI'] = SI
        results['P_H2O'] = P_H2O_arr

        return results

    def plot_evolution(self, results):
        
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(20, 20))
        plt.subplots_adjust(hspace=0)

        ax1.plot(results['time'], results['C_ocean'], label='C_ocean', color='black')
        ax1.plot(results['time'], results['Alk_ocean'], label='Alk_ocean', color='orange')
        ax1.plot(results['time'], results['Ca_ocean'], label='Ca_ocean', color='green')

        # ax12 = ax1.twinx()
        # ax12.plot(results['time'], np.array(dCodt), label='dC_ocean/dt', color='black', linestyle='--')
        # ax12.plot(results['time'], np.array(dAodt), label='dA_ocean/dt', color='orange', linestyle='--')
        # ax12.plot(results['time'], np.array(dCaodt), label='dCa_ocean/dt', color='green', linestyle='--')

        ax1.set_ylabel('Concentration (mol/kgw)')
        ax1.set_yscale('log')
        # ax12.set_ylabel('Rate of change of concentration (mol/kgw/s)')
        # ax12.set_yscale('log')

        # lines1, labels1 = ax1.get_legend_handles_labels()
        # lines12, labels12 = ax12.get_legend_handles_labels()
        # ax1.legend(lines1 + lines12, labels1 + labels12, loc='lower right')

        ax1.legend(loc='lower left')

        ax22 = ax2.twinx()
        ax22.plot(results['time'], results['T_surface'], 'k--', label='T_surface')
        ax22.plot(results['time'], results['T_seafloor'], 'g--', label='T_seafloor')
        ax2.plot(results['time'], results['P_CO2'], 'r-', label='P_CO2')
        ax2.plot(results['time'], results['P_H2O'], 'b-', label='P_H2O')

        ax2.set_ylabel('P (Pa)')
        ax2.set_yscale('log')

        ax22.set_ylabel('T (K)')
        ax22.set_ylim([250, 355])

        # Combine legends from ax2 and ax3
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines22, labels22 = ax22.get_legend_handles_labels()
        ax2.legend(lines2 + lines22, labels2 + labels22, loc='lower left')

        # 3. Set the x-axis label ONLY on the bottom plot
        ax3.set_xlabel('Time (yr)')

        ax3.plot(results['time'], results['F_diss'] / self.outgassing, 'b-', label='Dissolution rate')
        ax3.plot(results['time'], results['F_prec'] / self.outgassing, 'r-', label='Precipitation rate')
        # ax3.plot(results['time'], self.outgassing / self.outgassing, 'g-', label='Outgassing rate')
        # ax3.axhline(0.5, color='black')

        ax32 = ax3.twinx()

        ax32.plot(results['time'], results['SI'], 'k--', label='Calcite SI')
        ax32.set_ylabel('Saturation Index')

        ax3.set_ylabel('Rate (F_out)')
        # ax3.set_xscale('log')
        # ax3.set_ylim([1e12, 1e14])

        lines3, labels3 = ax3.get_legend_handles_labels()
        lines32, labels32 = ax32.get_legend_handles_labels()
        ax3.legend(lines3 + lines32, labels3 + labels32, loc='lower left')

        # Optional: Ensure tick labels on the top plot are definitely off (sharex usually handles this)
        ax1.tick_params(labelbottom=False)
        ax2.tick_params(labelbottom=False)

        # plt.tight_layout()
        plt.savefig('Evolution.png')
        plt.show()
        plt.close()
