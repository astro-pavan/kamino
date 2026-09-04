import numpy as np
np.set_printoptions(precision=1)
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import json
import os

from kamino.constants import (
    G, YR, SOLAR_CONSTANT, M_EARTH, R_EARTH,
    EARTH_OUTGASSING, EARTH_CRUST_PRODUCTION_RATE_PER_AREA,
    EARTH_HYDROTHERMAL_FLUX_PER_AREA, EARTH_MANTLE_MG_SI, EARTH_DELTA_IW,
    A_SEAFLOOR_EARTH, EARTH_CL_OUTGASSING_RATIO
)
from kamino.chemistry import elements, get_ocean_state, c_idx, si_idx, alk_idx, ca_idx, mg_idx, na_idx, cl_idx, so4_idx, ChemistryError
from kamino.weathering import get_weathering_flux, get_continental_weathering_flux, ALPHA_REF
from kamino.precipitation import get_precipitation
from kamino.mineral_info import (
    carbonate_minerals, clay_minerals, silica_minerals,
    reverse_weathering_minerals, evaporite_minerals,
)

from kamino.crust_composition import mineral_composition

from kamino.climate.clima_interpolator import get_T_surface
from kamino.climate.analytic import get_T_surface_analytic, maximum_greenhouse
from kamino.utils import august_roche_magnus_formula

output_path = os.path.join(os.path.dirname(__file__), '../../output/')
os.makedirs(output_path, exist_ok=True)

KD_MG_HT = 1.394755e-02
K_CL_SUBDUCTION = 1.373251e-04
K_NA_CONT_REMOVAL = 4.234317e-03

_S_TERR_EARTH = 5 / (1e6 * YR)   # m/s at land_fraction = 0.3

_EARTH_CA_SOURCES  = 5.731e12 / YR   # continental + seafloor + HT at Earth [mol/s]
_ABIOTIC_CA_3MYR   = 0.321e12 / YR   # abiotic Ca sink at tau_prec=3 Myr, Earth [mol/s]
_TAU_PREC_REF_K    = 3e6 * YR

# Fast-precipitation relaxation timescale, referenced to an Earth-depth ocean. `tau_prec`
# defaults to None and is resolved in __init__ as TAU_PREC_REF * (ocean_depth / OCEAN_DEPTH_REF),
# so a 3 km ocean keeps the historical 100 kyr exactly and only deeper oceans change.
#
# Why scale it. The fast bucket (carbonates, clays, silica, evaporites) is not a kinetic rate --
# those phases reach saturation and stay there, so their flux is set by solute SUPPLY and the
# timescale barely enters. What `tau_prec` has to do is relax a reservoir, and the reservoir is
# volumetric while precipitation is areal, so the relaxation time grows with ocean depth. Holding
# it at 100 kyr in a 20 km ocean therefore buys nothing physical and costs a 3-6x stiffer system:
# it forces the integrator to resolve a mode far faster than anything the chemistry is doing.
#
# Measured at 20 km, S = 1.0 (100 kyr -> 667 kyr): Mg/Si 1.25 dT = +0.29 K, Ca +6.4%, steps
# 1232 -> 384;  Mg/Si 0.5 dT = +0.07 K, Ca +1.1%, steps 1477 -> 243. Calcite stays buffered at
# SI = +0.021 (baseline +0.004).
#
# The scaling must NOT be applied as a flat increase for all depths. Setting tau_prec = 5 Myr at
# 3 km breaks the carbonate buffer outright -- calcite runs to SI = +1.219 (16x supersaturated),
# ocean Ca rises 11x to 3.33 mM, pH falls 0.72 and T rises 7.2 K. Shallow oceans need the fast
# bucket to stay ON saturation, which is exactly what the reference value delivers.
#
# `tau_rw` is deliberately NOT scaled with it. Sepiolite(d) sits ~4 log units supersaturated
# (SI = +3.8) at every timescale tested, i.e. it never reaches equilibrium, so tau_rw IS the
# reverse-weathering flux rather than a relaxation constant. Scaling both together to preserve
# their ratio cools a 20 km / Mg/Si 1.25 world by 25 K (T 344.09 -> 319.26, pCO2 0.553 -> 0.034)
# because slowing reverse weathering frees alkalinity for calcite. See docs/development_history.md.
# Redox state (pe) of the ocean and seafloor pore fluid.
#
# This model is ABIOTIC by construction -- there is no oxygenic photosynthesis, so there is no
# source of free O2. Its oceans should therefore be reducing, like the Archean. Until 2026-08-27
# no pe was set anywhere, which is NOT the same as making no assumption: PHREEQC then uses its
# own default of pe = 4.0, an oxidising ocean in which ferric Goethite (FeOOH) is supersaturated
# and strips every mole of dissolved iron, taking 2 eq of alkalinity with it. Every result
# produced before this parameter existed silently assumed an oxygenated planet.
#
# Measured effect of removing that assumption (3 km, Mg/Si 1.25, dIW -2, alpha 2):
#
#     S     oxic (pe 4)              anoxic (no ferric Fe)
#     0.8   T 298.3, Fe 0            T 268.2, Fe 69 uM
#     1.0   T 309.8, Fe 0            T 290.5, Fe 1107 uM
#
# i.e. 19-30 K, and ocean Fe moves from exactly zero into the Archean range (~0.05-0.5 mM). It
# also roughly doubles the influence of the crust dIW axis, because iron only matters to the
# carbon cycle when it stays dissolved.
#
# Reference values for natural waters: modern oxic seawater ~ +12.5; anoxic marine pore water
# ~ -3 to -5; Archean ocean reconstructions ~ -3 to 0. Goethite goes undersaturated below
# pe ~ -5.8, and Siderite (ferrous FeCO3, the Archean iron sink) becomes stable below pe ~ 0.
#
# NOTE this is a distinct quantity from the two dIW values in the crust pipeline. Those set the
# oxygen fugacity of the MANTLE at core formation and of the melt; this is the ambient redox of
# the WATER-ROCK system. They must not be conflated.
#
# Set to None to restore the pre-2026-08-27 behaviour (PHREEQC's default pe = 4).
PE_DEFAULT = -3.0

TAU_PREC_REF = 100e3 * YR
TAU_RW_REF = 5e6 * YR 
OCEAN_DEPTH_REF = 3000.0   # m; the depth at which TAU_PREC_REF applies

TAU_ATM = 1e4 * YR

tau_r_avg = 3e7 * YR   # EMA timescale for convergence rate smoothing (~30 Myr)

WATER_ROCK_RATIO_LT = 3

# Surface albedo. Module-level so plot_results.py can recompute T from a stored final
# state using the exact same formula Planet.__init__ uses, instead of duplicating the
# numbers (which would silently drift if these ever change here).
OCEAN_ALBEDO = 0.3
LAND_ALBEDO = 0.3


class ChemistryFallbackLimitExceeded(RuntimeError):
    pass


class WallClockLimitExceeded(RuntimeError):
    """A run exceeded its wall-clock budget and was abandoned.

    The fallback cap cannot catch every runaway: a run whose chemistry keeps CONVERGING but
    whose integrator crawls (tiny timesteps on stiff chemistry) never spends fallback budget and
    never trips `chemistry_void` either, so it grinds indefinitely. In the fast_18 sweep the
    slowest single run took 3.2 hours and the top 5% of runs consumed 57% of total CPU time.
    A wall-clock cap is the only bound that holds regardless of WHY a run is slow.
    """
    pass


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
            reverse_weathering: bool=True,
            alpha: float=ALPHA_REF,
            water_rock_ratio: float=WATER_ROCK_RATIO_LT,
            f_HT: float=0.0,
            cl_outgassing_ratio: float=EARTH_CL_OUTGASSING_RATIO,
            mantle_mg_si: float=EARTH_MANTLE_MG_SI,
            delta_iw: float=EARTH_DELTA_IW,
            pe: float | None=PE_DEFAULT,   # ocean/pore redox; see PE_DEFAULT
            tau_prec: float | None=None,   # None -> scaled with ocean depth, see TAU_PREC_REF
            tau_rw: float=TAU_RW_REF,
            kd_mg_ht: float=KD_MG_HT,
            k_na_cont_removal: float=K_NA_CONT_REMOVAL,
            k_cl_subduction: float=K_CL_SUBDUCTION,
            verbose: bool=False,
            name: str='planet',
            climate_model: str='analytic',
            ):

        if tau_prec is None:
            tau_prec = TAU_PREC_REF * (ocean_depth / OCEAN_DEPTH_REF)

        planet_config = {
            "name": name,
            "mass": mass,
            "radius": radius,
            "background_pressure": background_pressure,
            "ocean_depth": ocean_depth,
            "land_fraction": land_fraction,
            "mantle_mg_si": mantle_mg_si,
            "delta_iw": delta_iw,
            "instellation": instellation(0.0) if callable(instellation) else instellation,
            "crust_production_rate": crust_production_rate,
            "outgassing": outgassing,
            "reverse_weathering": reverse_weathering,
            "alpha": alpha,
            "water_rock_ratio": water_rock_ratio,
            "f_HT": f_HT,
            "climate_model": climate_model,
            "cl_outgassing_ratio": cl_outgassing_ratio,
            "pe": pe,
            "tau_prec": tau_prec,
            "tau_rw": tau_rw,
            "kd_mg_ht": kd_mg_ht,
            "k_na_cont_removal": k_na_cont_removal,
            "k_cl_subduction": k_cl_subduction,
        }

        # Setting up climate model
        self._get_T_surface = get_T_surface_analytic if climate_model == 'analytic' else get_T_surface

        # Planet properties
        self.name = name
        self.mass = mass
        self.radius = radius
        self.gravity = (G * self.mass) / (self.radius ** 2)
        self.surface_area = 4 * np.pi * self.radius ** 2
        self.P_background = background_pressure

        # Ocean properties
        self.ocean_depth = ocean_depth
        self.ocean_water_mass = self.ocean_depth * self.surface_area * 1000

        # Tectonic properties
        # Crust mineralogy from the two composition axes. `mantle_mg_si` is the mantle's molar
        # Mg/Si (Earth 1.25) and `delta_iw` is the CORE-FORMATION oxygen fugacity relative to
        # iron-wustite (Earth -2), which sets mantle FeO and hence the crust's iron. Mantle
        # potential temperature is no longer an input: it is solved per composition by the
        # melting model and recorded in crust_compositions.csv.
        #
        # NOTE the rename from `mg_si_ratio`. That parameter was a MULTIPLIER on Earth's 1.23,
        # not a Mg/Si, so a config passing mg_si_ratio=1 meant Mg/Si=1.23. Reusing the name
        # would have silently changed the meaning of every existing sweep config.
        self.mantle_mg_si = mantle_mg_si
        self.delta_iw = delta_iw
        self.crust_composition = mineral_composition(mantle_mg_si, delta_iw)
        self.crust_production_rate = EARTH_CRUST_PRODUCTION_RATE_PER_AREA * crust_production_rate
        self.hydrothermal_flux = EARTH_HYDROTHERMAL_FLUX_PER_AREA * crust_production_rate

        # Mineral precipiation properties
        self.pore_precipitating_minerals = clay_minerals
        self.fast_ocean_precipitating_minerals = carbonate_minerals + clay_minerals + silica_minerals + evaporite_minerals
        self.rw_ocean_precipitating_minerals = list(reverse_weathering_minerals) if reverse_weathering else []
        self.ocean_precipitating_minerals = self.fast_ocean_precipitating_minerals + self.rw_ocean_precipitating_minerals

        # Outgassing properties
        self.outgassing_flux = np.zeros(elements.shape)
        self.outgassing_flux[c_idx]  = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing
        self.outgassing_flux[cl_idx] = (EARTH_OUTGASSING / YR) * self.surface_area * outgassing * cl_outgassing_ratio
        self.outgassing_flux[alk_idx] -= self.outgassing_flux[cl_idx] # HCl outgassing acidifies the ocean: H+ + HCO3- → CO2 + H2O, consuming 1 eq Alk per mol Cl.

        # Biogenic flux properties (doesn't do anything)
        self.k_biogenic = np.zeros(elements.shape)

        # Chemistry properties
        self.cl_subduction_k = k_cl_subduction
        self.pe = pe
        self.tau_prec = tau_prec
        self.tau_rw = tau_rw
        self.alpha = alpha
        self.water_rock_ratio = water_rock_ratio
        self.f_HT = f_HT
        self.kd_mg_ht = kd_mg_ht

        # Land properties
        self.land_fraction = land_fraction
        self.land_area = land_fraction * self.surface_area
        self._s_terr = _S_TERR_EARTH * (land_fraction / 0.3)
        self.shelf_depth = 1000.0  # m — representative continental shelf depth for carbonate burial
        self.shelf_precipitating_minerals = carbonate_minerals
        self.na_cont_k = k_na_cont_removal

        # Climate calculation parameters
        ocean_albedo = OCEAN_ALBEDO
        land_albedo = LAND_ALBEDO
        self.relative_humidity = 0.5

        if callable(instellation):
            self._instellation_fn = lambda t: instellation(t) * SOLAR_CONSTANT
        else:
            self._instellation_fn = lambda t: instellation * SOLAR_CONSTANT

        self.instellation = self._instellation_fn(0.0)  # t=0 value for config/logging
        self.albedo = land_albedo * land_fraction + ocean_albedo * (1 - land_fraction)

        # Output parameters
        self._output_filename = f"{output_path}{self.name}.json"
        with open(self._output_filename, 'w') as f:
            json.dump(planet_config, f, indent=0)
        self.verbose = verbose

    def dY_dt(self, t, Y):

        # Wall-clock budget. Checked here rather than after solve_ivp returns, for the same
        # reason as the fallback cap: a check on the finished run would save nothing. Jacobian
        # probes are included deliberately -- they cost wall time like any other evaluation.
        _deadline = getattr(self, '_wall_deadline', None)
        if _deadline is not None and time.time() > _deadline:
            self._abort_t = t
            self._abort_Y = np.asarray(Y).copy()
            raise WallClockLimitExceeded(
                f'wall-clock limit of {self._wall_limit:.0f} s exceeded at t = {t / YR:.3e} yr'
            )

        # Extract state vector
        P_CO2 = Y[0]
        P_H2O = Y[1]
        b_ocean = Y[2:-1]   # Y[-1] is r_avg, excluded from chemistry

        # Input safety
        P_CO2 = np.clip(P_CO2, 0, 1e6)
        P_H2O = np.maximum(0, P_H2O)
        b_ocean = np.maximum(b_ocean, 0.0)

        # Update atmosphere properties
        P_surface = self.P_background + P_CO2 + P_H2O
        T_surface = self._get_T_surface(self._instellation_fn(t), P_CO2, self.albedo, False)
        P_H2O_new = august_roche_magnus_formula(T_surface)

        # Update seafloor physical properties
        T_seafloor = 1.02 * T_surface - 16.7
        P_pore = (self.P_background + P_CO2 + P_H2O) + 1000 * self.gravity * self.ocean_depth
        T_seafloor = np.maximum(T_seafloor, 274)
        T_pore = T_seafloor + 9
        self._T = T_surface
        ocean_water_per_area = self.ocean_depth * 1000.0

        # Volcanic outgassing fluxes
        F_vol = self.outgassing_flux / self.ocean_water_mass

        try:

            # Ocean surface conditions
            P_CO2_new, pH_surface = get_ocean_state(P_surface, T_surface, b_ocean)
            assert P_CO2_new > 0
            self._pH_surface = pH_surface

            # Fast precipitation: carbonates, clays, silica, evaporites (tau_prec ~100 kyr)
            F_prec_fast, pH_seafloor, SI = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=self.fast_ocean_precipitating_minerals, precipitation_timescale=self.tau_prec, pe=self.pe)
            F_prec = F_prec_fast
            F_prec_rw = np.zeros(elements.shape)

            # Slow precipitation: reverse weathering clays (tau_rw ~10-100 Myr)
            if self.rw_ocean_precipitating_minerals:
                F_prec_rw, _, SI_rw = get_precipitation(P_pore, T_seafloor, b_ocean, precipitating_minerals=self.rw_ocean_precipitating_minerals, precipitation_timescale=self.tau_rw, pe=self.pe)
                F_prec = F_prec + F_prec_rw
                SI.update(SI_rw)

            # SI Diagnostics
            self._SI = SI
            
            # Sedimentation rate
            F_carb_abiotic = max(0.0, -F_prec[c_idx])
            F_sil_abiotic  = max(0.0, -F_prec[si_idx])
            S_sed = (F_carb_abiotic * 0.100 / 2710.0 + F_sil_abiotic * 0.060 / 2650.0) * ocean_water_per_area + self._s_terr

            # Hydrothermal flux
            J_total = EARTH_HYDROTHERMAL_FLUX_PER_AREA * (self.crust_production_rate / EARTH_CRUST_PRODUCTION_RATE_PER_AREA)

            # Off axis (low-temperature) hydrothermal weathering
            flux_LT, diag_LT = get_weathering_flux(P_pore, T_pore, P_CO2, b_ocean, alpha=self.alpha, rate=self.crust_production_rate, J=J_total, crust_composition=self.crust_composition, sedimentation_rate=S_sed, precipitating_minerals=self.pore_precipitating_minerals, water_rock_ratio=self.water_rock_ratio, pe=self.pe)
            F_diss = (flux_LT * self.surface_area) / self.ocean_water_mass

            # HT Mg->Ca exchange (parameterized)
            _ht_rate = self.kd_mg_ht * b_ocean[mg_idx] * J_total * self.surface_area / self.ocean_water_mass
            F_ht_exchange = np.zeros(elements.shape)
            F_ht_exchange[mg_idx] = -_ht_rate
            F_ht_exchange[ca_idx] = +_ht_rate

            # Na sink (parameterized)
            F_na_rw = np.zeros(elements.shape)
            F_na_rw[na_idx] = -self.na_cont_k * b_ocean[na_idx] * J_total * self.surface_area / self.ocean_water_mass
            F_na_rw[alk_idx] = F_na_rw[na_idx]

            # Land based fluxes
            F_cont = np.zeros_like(F_diss)
            F_shelf_prec = np.zeros_like(F_diss)

            if self.land_fraction > 0:

                # Continental silicate weathering
                F_sil_cont = get_continental_weathering_flux(T_surface, P_CO2).copy()

                F_cont = F_sil_cont * self.land_area / self.ocean_water_mass

                # Shallow carbonate precipitation on continental shelves
                P_shelf = P_surface + 1000 * self.gravity * self.shelf_depth
                F_shelf_prec, _, _ = get_precipitation(P_shelf, T_seafloor, b_ocean, precipitating_minerals=self.shelf_precipitating_minerals, precipitation_timescale=self.tau_prec, pe=self.pe)

            # Cl subduction
            F_cl_subduct = np.zeros(elements.shape)
            F_cl_subduct[cl_idx] = -self.cl_subduction_k * b_ocean[cl_idx] * J_total * self.surface_area / self.ocean_water_mass
            F_cl_subduct[alk_idx] = -F_cl_subduct[cl_idx]

            # Total final flux
            F_net = F_vol + F_prec + F_shelf_prec + F_diss + F_ht_exchange + F_cont + F_cl_subduct + F_na_rw

            # Diagnostic record
            self._flux_terms = {
                'outgassing':         F_vol,
                'seafloor LT':        F_diss,
                'HT Mg-Ca exchange':  F_ht_exchange,
                'continental':        F_cont,
                'precipitation':      F_prec_fast,
                'shelf carbonate':    F_shelf_prec,
                'reverse weathering': F_prec_rw,
                'Cl subduction':      F_cl_subduct,
                'Na albitization':    F_na_rw,
            }
            self._state = {
                'T_surface':   T_surface,
                'T_seafloor':  T_seafloor,
                'T_pore':      T_pore,
                'P_CO2':       P_CO2,
                'P_surface':   P_surface,
                'P_pore':      P_pore,
                'pH_surface':  pH_surface,
                'pH_seafloor': pH_seafloor,
                'Da':          diag_LT.get('Da', np.nan),
                'A_reactive':  diag_LT.get('A_reactive', np.nan),
                'supply_efficiency': diag_LT.get('supply_efficiency', np.nan),
                'pore_SI':     diag_LT.get('secondary_SI', {}),
                'ocean_SI':    SI,
                # Kept so time_evolve can evaluate the pore-fluid calcite saturation once on the
                # accepted final state. It is not in `pore_SI`: the pore space precipitates clays
                # only, so Calcite never enters that equilibrium and has no SI there.
                'b_pore':      diag_LT.get('b_pore'),
                'alk_flux_lt': float(flux_LT[alk_idx]),
            }

        except (ChemistryError, AssertionError): # Chemistry has left the valid domain (typically high P_CO2 where PHREEQC cannot converge).

            # LSODA needs a finite derivative. Count the event either way -- a run that spent
            # time in a degraded state must not look like a clean termination (see
            # 'chemistry_fallbacks' in the output JSON). Jacobian probes are excluded from the
            # count because they only perturb a Jacobian column, not the trajectory.
            if not getattr(self, '_jac_active', False):
                self._chem_fallbacks = getattr(self, '_chem_fallbacks', 0) + 1

                # Abandon the run once the budget is gone. Checked here rather than after
                # solve_ivp returns, because the whole point is to stop paying: a check on
                # the finished run would save nothing. Jacobian probes are outside this
                # branch, so they can neither spend the budget nor trip the limit.
                _limit = getattr(self, '_fallback_limit', None)
                if _limit is not None and self._chem_fallbacks > _limit:
                    self._abort_t = t
                    self._abort_Y = np.asarray(Y).copy()
                    raise ChemistryFallbackLimitExceeded(
                        f'{self._chem_fallbacks} chemistry fallbacks exceeded the limit of '
                        f'{_limit} at t = {t / YR:.3e} yr'
                    )

            # Reuse the last derivative that DID converge
            if getattr(self, '_dYdt_last_good', None) is not None:
                return self._dYdt_last_good.copy() # type: ignore

            # Nothing has converged yet (failure on the very first evaluation), so there is no
            # previous derivative to hold: fall back to pure outgassing.
            dYdt = np.zeros_like(Y)
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


        # Diagnostic output
        carbon_flux = dYdt[3]
        carbon = Y[3]
        calcite_SI = SI['Calcite']

        if self.verbose and not getattr(self, '_jac_active', False):
            self._step_count = getattr(self, '_step_count', 0) + 1
            dt_str = f'  dt={((t - self._last_t)/YR):.1e}yr' if hasattr(self, '_last_t') else ''
            print(f't = {t/YR:.4e} yr  T = {T_surface:.1f}  P_CO2 = {P_CO2 / 1e5:.1e} bar  pH = {pH_seafloor:.1f}  Calcite SI = {calcite_SI:.1f}  C flux = {(carbon_flux / carbon) * 1e9 * YR:.1e} / Gyr  step={self._step_count}{dt_str}  ', end='\r')
            self._last_t = t

        self._pH = pH_seafloor

        # Held for the ChemistryError path above, so a failed evaluation can return a
        # continuous derivative instead of a fabricated outgassing-only one. Recorded from
        # Jacobian probes too: those sit within jac_epsilon (~1%) of the trajectory, which is
        # far closer to the truth than switching every sink off.
        self._dYdt_last_good = dYdt.copy()

        # Successful trajectory evaluations, counted on the same basis as _chem_fallbacks
        # (Jacobian probes excluded) so the two are directly comparable. Their ratio is what
        # says how much of a trajectory is real -- see the void check in time_evolve.
        if not getattr(self, '_jac_active', False):
            self._chem_ok = getattr(self, '_chem_ok', 0) + 1

        return dYdt

    def _final_diagnostics(self, t_state, Y_state, ca_amount):
        """Weathering diagnostics for a state the trajectory has already reached.

        `dY_dt` computes Da, pH, the ocean SI dict and b_pore on every call, so this just
        re-evaluates it once more at (t_state, Y_state) and reads them off `self._state`. That
        re-evaluation is a diagnostic query on an already-reached state, not new integration
        work, so both abort mechanisms are suspended for the duration of this one call and
        restored immediately after -- exactly as a Jacobian probe must never be allowed to trip
        either of them.

        Used for both the accepted final state of a normal completion AND the abort state of a
        wall_timeout/fallback_limit run. The abort state is just as real a state as a converged
        one -- solve_ivp had already reached it, `dY_dt` had already evaluated it once, and this
        merely re-evaluates that same state a second time under a suspended deadline -- but
        before this method existed, the early `return` in the wall_timeout/fallback_limit except
        branch skipped this step entirely, so every such run recorded da=NaN regardless of how
        close to a real answer its abort state actually was. As a side effect this also fixes
        `self._T`/`self._pH` for that abort state, which previously could be left holding
        whichever intermediate call (possibly a Jacobian probe) happened to run last before the
        wall-clock limit fired.

        `ca_amount` is the ocean Ca inventory (mol/kgw) at (t_state, Y_state), used only to gate
        ocean_si (PHREEQC returns a spurious -inf once Ca is exhausted).
        """
        diagnostics = {"da": np.nan, "calcite_si": np.nan, "ocean_si": np.nan,
                       "alk_flux": np.nan, "pH_seafloor": np.nan}

        _saved = (self._fallback_limit, self._wall_deadline, self._chem_fallbacks,
                 self._chem_ok, self._dYdt_last_good)
        self._fallback_limit = None
        self._wall_deadline = None
        try:
            self.dY_dt(t_state, Y_state)
        except Exception:
            pass
        finally:
            (self._fallback_limit, self._wall_deadline, self._chem_fallbacks,
             self._chem_ok, self._dYdt_last_good) = _saved

        _final = getattr(self, '_state', None)
        if _final:
            diagnostics["da"] = float(_final.get('Da', np.nan))
            diagnostics["pH_seafloor"] = float(_final.get('pH_seafloor', np.nan))
            # Tmol eq/yr, normalised on A_SEAFLOOR_EARTH (0.7 of Earth's surface) -- the
            # same fixed normalisation plot_results uses, so the two agree exactly.
            diagnostics["alk_flux"] = (float(_final.get('alk_flux_lt', np.nan))
                                       * A_SEAFLOOR_EARTH * YR / 1e12)
            # Ocean calcite SI is only meaningful while there is calcium left to saturate with;
            # at the ODE floor PHREEQC returns a spurious -inf.
            _ocean_si = (_final.get('ocean_SI') or {}).get('Calcite', np.nan)
            if ca_amount > 1e-6:
                diagnostics["ocean_si"] = float(_ocean_si)
            _b_pore = _final.get('b_pore')
            if _b_pore is not None:
                try:
                    _, _, _si_p = get_precipitation(
                        _final['P_pore'], _final['T_pore'],
                        np.maximum(np.asarray(_b_pore, dtype=float), 0.0),
                        precipitating_minerals=['Calcite'],
                        precipitation_timescale=1e6 * YR, pe=self.pe)
                    diagnostics["calcite_si"] = float(_si_p.get('Calcite', np.nan))
                except Exception:
                    pass
        return {k: (None if v is None or not np.isfinite(v) else v)
                for k, v in diagnostics.items()}

    def time_evolve(self, t_end=2e9 * YR, jac_epsilon=0.01, b0=None, initial_pco2=1000,
                    convergence_threshold=0.05, max_chemistry_fallbacks=None,
                    void_fraction=0.5, max_wall_seconds=None):
        """Integrate to steady state, or until the state leaves the model's validity box.

        max_chemistry_fallbacks: abandon the run once this many trajectory derivative
            evaluations have fallen back to a held derivative (see FallbackLimitExceeded).
            None (the default) means no limit, so existing callers are unaffected. Sweeps
            should set it -- the fallback tail is what makes them slow.
        void_fraction: if more than this fraction of trajectory derivative evaluations fell
            back to a held derivative, the run is reported as 'chemistry_void' instead of its
            nominal termination. Catches the failure mode an absolute fallback cap cannot:
            chemistry that fails outright is cheap, so it never trips the cap (see the void
            check below). The nominal outcome is kept as 'termination_raw'.
        max_wall_seconds: abandon the run after this much wall-clock time, reporting
            'wall_timeout'. None (the default) means no limit. This is the only cap that bounds
            a run whose chemistry converges but whose integrator crawls -- see
            WallClockLimitExceeded. Sweeps should set it: in fast_18 the top 5% of runs took 57%
            of the CPU, and 2.9% of runs (fallback_limit + chemistry_void) took 32%.
        """

        self._wall_limit = max_wall_seconds
        self._wall_deadline = (time.time() + max_wall_seconds) if max_wall_seconds else None

        Y0 = np.zeros(elements.shape[0] + 3)  # +2 for P_CO2/P_H2O, +1 for r_avg

        Y0[0] = initial_pco2
        Y0[1] = 1000
        Y0[-1] = 1.0 / (1e6 * YR)  # r_avg starts high (1/Myr) so convergence can't fire immediately
        if b0 is not None:
            Y0[2:-1] = np.asarray(b0)

        # pCO2 ceiling: the maximum-greenhouse pCO2 for this planet's instellation, not a flat
        # 10 bar. Beyond this pCO2 more CO2 COOLS the planet (past maximum greenhouse) and the
        # OLR fit is extrapolating outside its 1-10 bar range, so this is the true edge of where
        # the climate model means anything -- and it IS the physical inner-habitable-zone limit,
        # rather than an artifact of the fit's domain.
        #
        # Computed ONCE here, not inside _domain_margins. maximum_greenhouse() scans ~60 pCO2
        # values (each a full get_T_surface_analytic root-find) plus an optimizer refinement --
        # call it every event evaluation (every solver step) and it dominates the run cost.
        # Uses self.instellation (the t=0 value), matching how the rest of _domain_margins
        # already treats a time-varying instellation_fn as fixed for the ceiling.
        P_CO2_HI, _ = maximum_greenhouse(self.instellation, self.albedo)
        T_LO, T_HI = 181.0, 389.0          # K, one degree inside the analytic scan bracket

        min_time = 2e6 * YR
        convergence_rate = convergence_threshold / (1e9 * YR)

        self._T = np.nan
        self._pH = np.nan
        self._chem_fallbacks = 0
        self._chem_ok = 0
        self._fallback_limit = max_chemistry_fallbacks   # read by dY_dt; None disables the cap
        self._dYdt_last_good = None   # no converged derivative to hold yet

        def _domain_margins(t, Y):
            P = max(float(np.clip(Y[0], 0.0, None)), 1e-12)
            T = self._get_T_surface(self._instellation_fn(t), P, self.albedo, False)
            return {
                'co2_high': np.log10(P_CO2_HI / P),
                'cold':     (T - T_LO) / 100.0,
                'hot':      (T_HI - T) / 100.0,
            }

        min_time_domain = 1e4 * YR

        def event_domain(t, Y):
            if t < min_time_domain:
                return 1.0  # same side as the real signal, so the crossing is always detectable
            return min(_domain_margins(t, Y).values())
        event_domain.terminal, event_domain.direction = True, -1 # type: ignore

        atol = np.ones_like(Y0) * 1e-6
        atol[0] = 1.0   # P_CO2 in Pa
        atol[1] = 1.0   # P_H2O in Pa
        atol[-1] = convergence_rate * 0.1  # r_avg: resolve to 10% of the convergence threshold

        def event_converged(t, Y):
            if t < min_time:
                return 1.0
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

        try:
            sol = solve_ivp(
                self.dY_dt,
                (0, t_end),
                Y0,
                method='LSODA',
                max_step=2e7 * YR,
                rtol=1e-3,
                atol=atol,
                jac=macro_jacobian,
                events=[event_domain, event_converged],
            )
        except (ChemistryFallbackLimitExceeded, WallClockLimitExceeded) as exc:

            _why = ('wall_timeout' if isinstance(exc, WallClockLimitExceeded)
                    else 'fallback_limit')
            end = time.time()
            if self.verbose:
                print()
                print(f'Aborted: {_why} -- {exc}')

            # The abort state is a real state the trajectory reached, not a discard -- recover
            # its weathering diagnostics (and re-settle self._T/self._pH onto it) the same way a
            # normal completion's accepted final state is. See _final_diagnostics.
            diagnostics = self._final_diagnostics(
                self._abort_t, self._abort_Y, self._abort_Y[2 + ca_idx])

            with open(self._output_filename, 'r') as f:
                output_data = json.load(f)
            output_data.update({
                "diagnostics": diagnostics,
                "simulation_time_seconds": end - start,
                "termination": _why,
                "domain_wall": None,
                "chemistry_fallbacks": self._chem_fallbacks,
                "fallback_limit": self._fallback_limit,
                "end_time_yr": getattr(self, '_abort_t', np.nan) / YR,
                "T": self._T,
                "P_CO2": float(np.clip(self._abort_Y[0], 0.0, None)) / 1e5,
                "pH": self._pH,
                "data": {"time": [], "y": []},
            })
            with open(self._output_filename, 'w') as f:
                json.dump(output_data, f, indent=0)
            return

        end = time.time()

        event_names = ['out_of_domain', 'converged']

        sol.y = np.maximum(sol.y, 0.0)
        time_steps = sol.t.tolist()
        state_variables = sol.y.tolist()

        if sol.t[-1] >= t_end:
            termination = "timeout"
        elif sol.status == 1:
            termination = next(name for name, t_ev in zip(event_names, sol.t_events) if len(t_ev) > 0)
        else:
            termination = "solver_failure"

        # Void check: what fraction of the trajectory is fabricated?
        _fab_total = self._chem_fallbacks + self._chem_ok
        fabricated_fraction = (self._chem_fallbacks / _fab_total) if _fab_total else 0.0
        termination_raw = termination
        if fabricated_fraction > void_fraction:
            termination = "chemistry_void"

        # Which wall was hit
        domain_wall = None
        if termination == "out_of_domain":
            _m = _domain_margins(sol.t[-1], sol.y[:, -1])
            domain_wall = min(_m, key=lambda k: _m[k])

        if self.verbose:
            print()
            print(f'Terminated: {termination} at t = {sol.t[-1]/YR:.3e} yr')
            print(f'Simulation time: {end - start:.0f} s')
            print(f'Y: {sol.y[2:, -1]} mol/kgw')

        # Recompute T, pH and the weathering diagnostics from the ACTUAL final accepted state.
        # self._T/self._pH are set as a side effect on EVERY call to dY_dt -- including Jacobian
        # finite-difference probes and solve_ivp's internal event-root-finding trials -- so
        # whichever call happened to run LAST is not guaranteed to be the accepted
        # (sol.t[-1], sol.y[:, -1]) state also used for "P_CO2" below. Near a domain wall this
        # matters: the climate response can be a genuine cliff (e.g. approaching the runaway
        # threshold at high instellation), so a routine 1% Jacobian perturbation is enough to
        # flip between a real ~389 K state and the analytic model's literal 400.0 "no
        # equilibrium found" sentinel -- corrupting the reported T/pH with a value that has
        # nothing to do with the reported P_CO2. Verified: at one S=1.15 'hot' termination this
        # reported T=400.0 while the true final P_CO2 (0.07644 bar) evaluates to T=388.99 K; a
        # +1% Jacobian probe of that P_CO2 lands exactly on the 400.0 rail, which is what got
        # left in self._T. `_final_diagnostics` re-evaluates dY_dt at the exact accepted state
        # (with both abort mechanisms suspended, so this extra call can never trip either) and
        # also records da/calcite_si/ocean_si/alk_flux/pH_seafloor from it -- one extra PHREEQC
        # solve per run (the pore-fluid calcite SI, the only part not already covered by an
        # equilibrium dY_dt performs), which is far cheaper than the plotting code reconstructing
        # this same state and re-running the weathering chemistry for every run it draws (about
        # 0.9 s each, ~30 minutes for a 2000-run sweep, every time a figure is regenerated).
        diagnostics = self._final_diagnostics(sol.t[-1], sol.y[:, -1], sol.y[2 + ca_idx, -1])

        with open(self._output_filename, 'r') as f:
            output_data = json.load(f)

        output_data.update({
            "diagnostics": diagnostics,
            "simulation_time_seconds": end - start,
            "termination": termination,
            "domain_wall": domain_wall,
            "chemistry_fallbacks": self._chem_fallbacks,
            "chemistry_ok": self._chem_ok,
            "fabricated_fraction": fabricated_fraction,
            "termination_raw": termination_raw,
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
