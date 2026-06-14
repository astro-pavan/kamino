import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import json
import os

from kamino.constants import *
from kamino.chemistry import elements, get_P_CO2, get_pH, c_idx, si_idx, alk_idx, ca_idx, mg_idx, na_idx, cl_idx, so4_idx, ChemistryError
from kamino.weathering import get_weathering_flux, get_continental_weathering_flux, J_ref_normalised, rate_ref, EARTH_CONTINENTAL_WEATHERING_REF, ALPHA_REF
from kamino.precipitation import get_precipitation
from kamino.mineral_info import *


from kamino.climate.clima_interpolator import get_T_surface as _get_T_interp
from kamino.climate.analytic import get_T_surface_analytic as _get_T_analytic
from kamino.utils import august_roche_magnus_formula

output_path = os.path.join(os.path.dirname(__file__), '../../output/')

KD_MG_HT = 1.197244e-02
K_CL_SUBDUCTION = 1.373251e-04
K_NA_ALBITIZATION = 2.194806e-03
_M_OCEAN_REF = 1.4e21 # kg

# CaCO3: calibrated to net Ca burial rate matching continental + seafloor weathering inputs
# (~6 Tmol/yr at modern [Ca]=10.3 mmol/kg). NET geological burial
K_BIO_CA = 6e12 / (YR * 10.3e-3 * _M_OCEAN_REF)  # s⁻¹, residence time ~2.4 Myr

# Opal: ~6 Tmol/yr net burial (Tréguer et al. 2021) at [Si]_ref = 0.1 mmol/kg
K_BIO_SI = 6e12  / (YR * 1.0e-4  * _M_OCEAN_REF)  # s⁻¹, residence time ~233 kyr

_S_TERR_EARTH = 5 / (1e6 * YR)   # m/s at land_fraction = 0.3

_EARTH_CA_SOURCES  = 5.731e12 / YR   # continental + seafloor + HT at Earth [mol/s]
_ABIOTIC_CA_3MYR   = 0.321e12 / YR   # abiotic Ca sink at tau_prec=3 Myr, Earth [mol/s]
_TAU_PREC_REF_K    = 3e6 * YR

TAU_ATM = 1e4 * YR

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
            crust_composition: dict[str, float]=basalt_49,
            reverse_weathering: bool=True,
            alpha: float=1.43,
            f_HT: float=0.0,
            cl_outgassing_ratio: float=0.02,
            tau_prec: float=100e3 * YR,
            tau_rw: float=5e6 * YR,
            f_bio: float=0.0,
            verbose: bool=False,
            name: str='planet',
            climate_model: str='analytic',
            ):
        
        planet_config = {
            "name": name,
            "mass": mass,
            "radius": radius,
            "background_pressure": background_pressure,
            "ocean_depth": ocean_depth,
            "land_fraction": land_fraction,
            "crust_composition": crust_composition.copy(),
            "instellation": instellation(0.0) if callable(instellation) else instellation,
            "crust_production_rate": crust_production_rate,
            "outgassing": outgassing,
            "reverse_weathering": reverse_weathering,
            "alpha": alpha,
            "cl_outgassing_ratio": cl_outgassing_ratio,
            "tau_prec": tau_prec,
            "tau_rw": tau_rw,
            "f_bio": f_bio
        }

        self._get_T_surface = _get_T_analytic if climate_model == 'analytic' else _get_T_interp

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

        # if crust_carbonate_content > 0:
        #     for mineral in self.crust_composition:
        #         self.crust_composition[mineral] *= (1 - crust_carbonate_content)
        #     self.crust_composition['Calcite'] = crust_carbonate_content

        self.pore_precipitating_minerals = clay_minerals
        self.fast_ocean_precipitating_minerals = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals
        self.rw_ocean_precipitating_minerals = list(reverse_weathering_minerals) if reverse_weathering else []
        self.ocean_precipitating_minerals = self.fast_ocean_precipitating_minerals + self.rw_ocean_precipitating_minerals
        self.shelf_depth = 1000.0  # m — representative continental shelf depth for carbonate burial
        self.shelf_precipitating_minerals = carbonate_minerals

        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[c_idx]  = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing
        self.outgassing_flux[cl_idx] = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing * cl_outgassing_ratio
        # HCl outgassing acidifies the ocean: H+ + HCO3- → CO2 + H2O, consuming 1 eq Alk per mol Cl.
        self.outgassing_flux[alk_idx] -= self.outgassing_flux[cl_idx]

        self.cl_subduction_k = K_CL_SUBDUCTION
        self.na_albit_k = K_NA_ALBITIZATION
        self.tau_prec = tau_prec
        self.tau_rw = tau_rw
        self.f_bio = f_bio

        self.tidally_locked = False

        self.f_HT = f_HT

        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area
        self._s_terr = _S_TERR_EARTH * (land_fraction / 0.3)

        _bio_ca_target = _EARTH_CA_SOURCES - _ABIOTIC_CA_3MYR * (_TAU_PREC_REF_K / self.tau_prec)
        self._K_BIO_CA = max(_bio_ca_target, 0.0) / (10.3e-3 * self.ocean_water_mass)
        self._K_BIO_SI = 6e12 / (YR * 1.0e-3 * self.ocean_water_mass)

        ocean_albedo = 0.3
        land_albedo = 0.3

        self.relative_humidity = 0.5

        if callable(instellation):
            self._instellation_fn = lambda t: instellation(t) * SOLAR_CONSTANT
        else:
            self._instellation_fn = lambda t: instellation * SOLAR_CONSTANT
        self.instellation = self._instellation_fn(0.0)  # t=0 value for config/logging
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

        self.alpha = alpha

        self._output_filename = f"{output_path}{self.name}.json"
        with open(self._output_filename, 'w') as f:
            json.dump(planet_config, f, indent=0)

        self.verbose = verbose

    def dY_dt(self, t, Y):

        P_CO2 = Y[0]
        P_H2O = Y[1]
        b_ocean = Y[2:-1]   # Y[-1] is r_avg, excluded from chemistry

        # Input safety

        P_CO2 = np.clip(P_CO2, 0, 1e6)
        P_H2O = np.maximum(0, P_H2O)
        b_ocean = np.maximum(b_ocean, 0.0)

        # Atmosphere properties

        P_surface = self.P_background + P_CO2 + P_H2O
        T_surface = self._get_T_surface(self._instellation_fn(t), max(P_CO2, 1.0), self.albedo, self.tidally_locked)
        P_H2O_new = august_roche_magnus_formula(T_surface)

        # Seafloor physical properties

        T_seafloor = 1.02 * T_surface - 16.7
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth
        T_seafloor = np.maximum(T_seafloor, 274)
        T_pore = T_seafloor + 9
        self._T = T_surface

        # Volcanic outgassing fluxes

        F_vol = self.outgassing_flux / self.ocean_water_mass

        try:

            P_CO2_new = get_P_CO2(P_surface, T_surface, b_ocean)
            assert P_CO2_new > 0

            # Fast precipitation: carbonates, clays, silica, evaporites (tau_prec ~100 kyr)
            F_prec, pH, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=self.fast_ocean_precipitating_minerals, precipitation_timescale=self.tau_prec)
            # Slow precipitation: reverse weathering clays (tau_rw ~10-100 Myr)
            if self.rw_ocean_precipitating_minerals:
                F_prec_rw, _, SI_rw = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=self.rw_ocean_precipitating_minerals, precipitation_timescale=self.tau_rw)
                F_prec = F_prec + F_prec_rw
                SI.update(SI_rw)
            F_carb_abiotic = max(0.0, -F_prec[c_idx])
            F_sil_abiotic  = max(0.0, -F_prec[si_idx])

            # Ocean chemistry — pH at surface conditions (not pore conditions)
            self._pH = get_pH(P_surface, T_surface, b_ocean)
            self._SI = SI
            ocean_water_per_area = self.ocean_depth * 1000.0
            
            # Biogenic rates with sigmoid gate on Calcite SI
            _calcite_si = SI.get('Calcite', -1.0)
            _bio_calcite_gate = 1.0 / (1.0 + np.exp(-5.0 * _calcite_si))  # ~0 at SI=-1, ~1 at SI=+1
            _bio_ca = self.f_bio * self._K_BIO_CA * b_ocean[ca_idx] * _bio_calcite_gate if self.f_bio > 0 else 0.0
            _bio_si = self.f_bio * self._K_BIO_SI * b_ocean[si_idx] if self.f_bio > 0 else 0.0

            S_sed = ((F_carb_abiotic + _bio_ca) * 0.100 / 2710.0 + (F_sil_abiotic  + _bio_si) * 0.060 / 2650.0) * ocean_water_per_area + self._s_terr
            J_total = J_ref_normalised * (self.crust_production_rate / rate_ref)
            J_LT = J_total * (1 - self.f_HT)
            J_HT = J_total * self.f_HT

            # Off axis hydrothermal weathering
            flux_LT, _ = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, alpha=self.alpha, rate=self.crust_production_rate, J=J_LT, crust_composition=self.crust_composition, sedimentation_rate=S_sed, precipitating_minerals=self.pore_precipitating_minerals)
            
            # High temperature hydrothermal weathering

            if self.f_HT > 0:
                flux_HT, _ = get_weathering_flux(P_pore, 473, P_CO2, b_ocean, high_temperature=True, alpha=self.alpha, rate=self.crust_production_rate, J=J_HT, crust_composition=self.crust_composition, sedimentation_rate=S_sed, precipitating_minerals=self.pore_precipitating_minerals, water_rock_ratio=2.0)
            else:
                flux_HT = np.zeros_like(flux_LT)

            # Balanced-Na correction (Coogan 2022): at LT, 1 mol Na from Albite
            # dissolution comes with ~1 mol HCO3-, which is cancelled when Na
            # eventually leaves the ocean. Apply only to the LT flux.
            # At HT (200°C, pH~6.7), the Na/Alk ratio from Albite dissolution is
            # ~0.05, so PHREEQC already gives the correct, much smaller Alk; the
            # LT convention would over-subtract Alk by ~20x.
            F_diss_LT = (flux_LT * self.surface_area) / self.ocean_water_mass
            F_diss_LT[alk_idx] -= F_diss_LT[na_idx]
            F_diss_HT = (flux_HT * self.surface_area) / self.ocean_water_mass
            F_diss = F_diss_LT + F_diss_HT

            F_cont = np.zeros_like(F_diss)
            F_shelf_prec = np.zeros_like(F_diss)

            if self.land_fraction > 0:

                # Continental silicate weathering
                F_sil_cont = get_continental_weathering_flux(T_surface, P_CO2).copy()
                
                # Continetal carbonate weathering: CaCO3 + CO2 + H2O -> Ca2+ + 2HCO3-
                F_carb_w = np.zeros(elements.shape) # switched off
                # if self.f_carb > 0:
                #     F_carb_w[ca_idx]  = self.f_carb * F_sil_cont[ca_idx]
                #     F_carb_w[alk_idx] = 2.0 * self.f_carb * F_sil_cont[ca_idx]
                #     F_carb_w[c_idx]   = self.f_carb * F_sil_cont[ca_idx]

                F_cont = (F_sil_cont + F_carb_w) * self.land_area / self.ocean_water_mass

                # Shallow carbonate precipitation on continental shelves
                P_shelf = P_surface + 1000 * self.gravity * self.shelf_depth
                F_shelf_prec, _, _ = get_precipitation(P_shelf, T_seafloor, b_ocean, precipitating_minerals=self.shelf_precipitating_minerals, precipitation_timescale=self.tau_prec)

            # HT Mg->Ca exchange (scales with J_HT: Mg removal is a HT-only process)
            _ht_rate = KD_MG_HT * b_ocean[mg_idx] * J_LT * self.surface_area / self.ocean_water_mass
            F_ht_exchange = np.zeros(elements.shape)
            F_ht_exchange[mg_idx] = -_ht_rate
            F_ht_exchange[ca_idx] = +_ht_rate

            # Cl subduction
            F_cl_subduct = np.zeros(elements.shape)
            F_cl_subduct[cl_idx] = -self.cl_subduction_k * b_ocean[cl_idx] * J_total * self.surface_area / self.ocean_water_mass

            # Na sink: albitization in HT crust + LT reverse weathering (zeolites/clays).
            F_na_rw = np.zeros(elements.shape)
            F_na_rw[na_idx] = -self.na_albit_k * b_ocean[na_idx] * J_total * self.surface_area / self.ocean_water_mass

            # Biogenic burial
            F_bio = np.zeros(elements.shape)
            if self.f_bio > 0:
                F_bio[ca_idx]  = -_bio_ca
                F_bio[alk_idx] = -2.0 * _bio_ca
                # CaCO3 burial + organic carbon burial (f_bio fraction of C outgassing)
                F_bio[c_idx]   = -_bio_ca - self.f_bio * F_vol[c_idx]
                F_bio[si_idx]  = -_bio_si

            F_net = F_vol + F_prec + F_shelf_prec + F_diss + F_cont + F_ht_exchange + F_cl_subduct + F_na_rw + F_bio

        except (ChemistryError, AssertionError): # Chemistry has left the valid domain (typically high P_CO2 where PHREEQC cannot converge).

            dYdt = np.zeros_like(Y) # Return pure outgassing so LSODA gets a finite derivative and the acid_ocean event can terminate cleanly.
            dYdt[2:-1] = F_vol
            self._F_net = dYdt[2:-1]
            return dYdt

        dYdt = np.zeros_like(Y)

        dYdt[0] = (P_CO2_new - P_CO2) / TAU_ATM
        dYdt[1] = (P_H2O_new - P_H2O) / TAU_ATM

        F_net[b_ocean <= 0.0] = np.maximum(F_net[b_ocean <= 0.0], 0.0)
        F_net[so4_idx] = 0.0  # SO4 pinned to background; no ODE evolution

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

    def time_evolve(self, t_end=2e9 * YR, jac_epsilon=0.01, b0=None, initial_pco2=1000, convergence_threshold=0.1):

        Y0 = np.zeros(elements.shape[0] + 3)  # +2 for P_CO2/P_H2O, +1 for r_avg

        Y0[0] = initial_pco2
        Y0[1] = 1000
        Y0[-1] = 1.0 / (1e6 * YR)  # r_avg starts high (1/Myr) so convergence can't fire immediately
        if b0 is not None:
            Y0[2:-1] = np.asarray(b0)

        T_runaway = 350
        T_snowball = 260
        P_CO2_acid_threshold = 5e5   # Pa (5 bar)
        P_CO2_floor = 1.0            # Pa — below this the planet is unambiguously going snowball
        min_time = 2e6 * YR
        min_time_co2_floor = 1e6 * YR  # guard against initial CO2 transient: blank ocean recovers in ~300 kyr
        convergence_rate = convergence_threshold / (1e9 * YR)

        self._T = np.nan
        self._pH = np.nan

        def event_snowball(t, Y):
            if t < min_time:
                return -1.0  # negative guard: direction=-1 won't trigger on the guard→actual transition
            P_CO2 = np.clip(Y[0], 0, 1e6)
            T_surface = self._get_T_surface(self._instellation_fn(t), max(P_CO2, 1.0), self.albedo, self.tidally_locked)
            return T_surface - T_snowball
        event_snowball.terminal, event_snowball.direction = True, -1 # type: ignore

        def event_hothouse(t, Y):
            if t < min_time:
                return -1.0  # negative guard: direction=+1 won't trigger on the guard→actual transition
            T = self._get_T_surface(self._instellation_fn(t), max(Y[0], 1.0), self.albedo, self.tidally_locked)
            return T - T_runaway
        event_hothouse.terminal, event_hothouse.direction = True, +1 # type: ignore

        def event_acid_ocean(t, Y):
            if t < min_time:
                return 1.0
            return P_CO2_acid_threshold - np.clip(Y[0], 0, None)
        event_acid_ocean.terminal, event_acid_ocean.direction = True, -1 # type: ignore

        def event_co2_floor(t, Y):
            P_CO2 = np.clip(Y[0], 0, 1e6)
            T_surface = self._get_T_surface(self._instellation_fn(t), max(P_CO2, 1.0), self.albedo, self.tidally_locked)
            if T_surface < T_snowball:
                return 1.0  # snowball handles this case
            # Smooth guard: linearly decreasing offset keeps function positive during the initial transient (first min_time_co2_floor). This avoids the discontinuous guard that confused Brent's root finder and fired spuriously at the guard boundary.
            guard_offset = max(0.0, 1.0 - t / min_time_co2_floor) * 1000.0
            return np.clip(Y[0], 0, None) - P_CO2_floor + guard_offset
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
            T_surface = self._get_T_surface(self._instellation_fn(t), max(P_CO2, 1e-2), self.albedo, self.tidally_locked)
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
                        jac[1, 1] = -1.0 / TAU_ATM
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
