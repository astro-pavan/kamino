#!/usr/bin/env python3
"""
Iterative calibration of Earth-reference ocean chemistry constants.

Adjusts K_NA_CONT_REMOVAL, ALPHA_REF and KD_MG_HT so that Planet.time_evolve()
holds modern seawater at steady state, from a charge-balanced seawater seed.

Key design choices
------------------
K_cl is determined analytically from the Cl flux balance and never updated from
the simulation.  Cl has a residence time of ~10 Gyr so it cannot equilibrate
within a realistic integration time, and it is the dominant charge-balancing
anion, so it is seeded to its analytically-known steady-state value.

The initial condition is modern seawater with alkalinity set to the exact ion
charge sum -- NOT the "blank ocean" (Cl and SO4 only, cations zero) this script
used previously.  That seed violated Alk = ION_CHARGE.b by 592.9 mEq/kg at t=0
and the run could never repair it, poisoning every result.  See make_b0().

This makes the script a stability test ("does the model HOLD modern seawater?")
rather than a from-scratch assembly test.  That is the more useful question, and
the only tractable one: Na needs ~500 Myr to accumulate and Ca/Mg never reach
target from zero, so a blank start spends the whole integration far from Earth.

Runs typically terminate "converged"; Cl's residual drift is small enough at the
seeded value that it no longer blocks the convergence event.

Scope: abiotic only
-------------------
This script previously had a Phase 2 that scanned `f_bio` (biogenic CaCO3 +
organic C burial) to close the carbon budget.  That has been REMOVED, because
the capability it targeted no longer exists in the model:

  * `Planet.__init__` has no `f_bio` parameter (passing one raises TypeError).
  * `self.k_biogenic` is hard-zeroed in planet.py and never read by `dY_dt`.

So the model is abiotic by construction, and pCO2 should be expected to sit
ABOVE the 280 ppm pre-industrial target -- an abiotic ocean has no biogenic
carbonate pump.  That offset is a known scope limitation, not a calibration
failure.  Reinstating Phase 2 requires first reinstating a biogenic flux in
planet.py.

`f_HT` is likewise no longer calibrated here.  It is still accepted by the
Planet constructor and stored as `self.f_HT`, but nothing in the model reads it
(verified: the only occurrence in planet.py is the assignment itself), so
scanning it did nothing.  The Ca budget is now set by the LT seafloor source and
the carbonate sink alone.

Four constants are fitted: K_na (Na) and KD_MG_HT (the Ca:Mg split) by least squares against
Na, Ca and Mg at fixed alpha and tau_rw, alternating with alpha rescaled to the net seafloor
alkalinity flux and tau_rw rescaled to the reverse-weathering Mg sink -- see calibrate().
tau_prec is held at its reference value.

Usage:
    /data/pt426/big-venv/bin/python experiments/calibrate_earth.py
"""
import sys
import os
import json
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import kamino.planet as planet_module
from kamino.planet import Planet, WATER_ROCK_RATIO_LT
from kamino.chemistry import elements, ION_CHARGE, alk_idx, c_idx, si_idx, ca_idx, mg_idx, na_idx, cl_idx, so4_idx, k_idx
from kamino.chemistry import SEAWATER_SO4, SEAWATER_K
from kamino.weathering import get_weathering_flux, ALPHA_REF as ALPHA_REF_CODE
from kamino.mineral_info import clay_minerals
from kamino.constants import (
    EARTH_HYDROTHERMAL_FLUX_PER_AREA as J_ref_normalised,
    EARTH_CRUST_PRODUCTION_RATE_PER_AREA as rate_ref,
    A_SEAFLOOR_EARTH as A_seafloor,
)
from kamino.constants import M_EARTH, R_EARTH, YR, EARTH_OUTGASSING, EARTH_ATM, G, SEAFLOOR_T_FLOOR

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), '../output')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Targets — modern seawater (Millero 2013)
# ---------------------------------------------------------------------------
T_Cl  = 546e-3
T_Na  = 469e-3
T_Ca  = 10.3e-3
T_Mg  = 52.8e-3
T_Alk = 2.3e-3
T_C   = 2.1e-3
T_pCO2 = 280.0   # ppm

TARGETS = dict(Na=T_Na, Ca=T_Ca, Mg=T_Mg, Alk=T_Alk, C=T_C)  # Cl excluded: see module docstring

N_ELEM = len(elements)

# ---------------------------------------------------------------------------
# Fixed Earth input parameters
# ---------------------------------------------------------------------------
OCEAN_DEPTH      = 3700.0
LAND_FRAC        = 0.3
OUTGASSING       = 1.0
CRUST_RATE       = 1.0
INSTELLATION     = 1.0
CL_OUTGASSING_RATIO = 0.02   # matches Planet constructor default

# 4 Gyr: Na (tau~500 Myr) is the slowest calibrated species and needs many
# e-foldings from a blank start. The convergence audit of fast_17 showed
# unconverged runs needed only ~1.1-1.8x more time than 2 Gyr, so 4 Gyr is
# comfortably past that for a single run.
T_END = 4e9 * YR

# Abandon a run whose chemistry has collapsed rather than let it burn hours
# fabricating derivatives; also surfaces 'chemistry_void' in the output.
MAX_CHEM_FALLBACKS = 5000

GRAVITY = G * M_EARTH / R_EARTH**2

# ---------------------------------------------------------------------------
# K_CL: analytical steady-state balance.
#
# At SS:  Cl_outgassing = Cl_subduction
#   (EARTH_OUTGASSING/YR) * A_surf * cl_ratio / M  =  K_cl * [Cl] * J_total * A_seafloor / M
#
# M cancels but THE TWO AREAS DO NOT:
#   K_cl = (EARTH_OUTGASSING/YR * cl_ratio / ([Cl]_target * J_total)) * (A_surf / A_seafloor)
#
# J_total = J_ref_normalised for crust_rate=1.
#
# CORRECTED 2026-09-09 (development_history.md section 35.2). This previously cancelled the two
# areas against each other, which was right until section 33.3: outgassing is emitted over the
# WHOLE SPHERE while Cl subduction acts over the SEAFLOOR ONLY, and Planet.seafloor_area is
# (1 - land_fraction) * surface_area. At LAND_FRAC = 0.3 the missing factor is 1/0.7 = 1.43 --
# the same 1.43x KD_MG_CALIB and K_NA_CALIB silently absorbed, but recoverable in closed form
# here rather than through a fit.
#
# It matters because K_CL_ANALYTIC is an INPUT to the least-squares run (run_planet passes it
# straight into Planet), so a wrong Cl puts the wrong anion charge in the ocean that K_na, alpha
# and KD_mg are fitted against. Uncorrected it targets a 780 mM Cl steady state, not 546.
# ---------------------------------------------------------------------------
_AREA_RATIO_CL = 1.0 / (1.0 - LAND_FRAC)   # surface_area / seafloor_area, as Planet defines them
K_CL_ANALYTIC = ((EARTH_OUTGASSING / YR * CL_OUTGASSING_RATIO)
                 / (T_Cl * J_ref_normalised)) * _AREA_RATIO_CL

# ---------------------------------------------------------------------------
# ALPHA_REF diagnostic: the seafloor reactive-area scaling that would make the NET seafloor
# alkalinity flux (after pore kaolinite/goethite, at the model's pe) equal FLUX_TARGET_NET at Earth
# pore conditions with modern seawater as input. Before 2026-09-24 this used the primary flux at
# PHREEQC's default pe = 4, where ~90% of the "alkalinity" was Fe that the model never releases.
#
# This is a SEPARATE anchor from the alpha the loop below calibrates, and the two
# answer different questions. This one asks "what alpha reproduces a 1 Tmol/yr
# seafloor flux at fixed modern seawater?"; the loop asks "what alpha reproduces
# modern Ca+Mg at steady state?". They need not agree -- at Earth pore conditions
# the system is transport-limited, so the flux saturates and alpha has far more
# authority over the flux target than over the ocean. Reported for comparison;
# the loop's value is the one to paste into weathering.py.
#
# The residual now uses water_rock_ratio=WATER_ROCK_RATIO_LT to match what
# Planet actually passes to get_weathering_flux. It previously left it at None,
# which selects a different equilibrium branch in get_b_eq (the
# lt_equilibrium_buffer_minerals guard), so the old number was calibrated under
# conditions the model never runs in.
# T_surface=288K -> T_seafloor=277K -> T_pore=286K; depth=3000m.
# ---------------------------------------------------------------------------
_alpha_T_pore    = 286.0
_alpha_P_pore    = 1000.0 * 10.0 * 3000.0
_alpha_P_CO2     = EARTH_ATM * 280e-6
FLUX_TARGET_NET  = 0.9    # Teq/yr, net low-T seafloor alkalinity flux (Coogan & Dosso 2022, GCA 329, 22)
_alpha_flux_norm = (FLUX_TARGET_NET * 1e12 / YR) / A_seafloor   # per m² of seafloor

_alpha_b = np.zeros(len(elements))
_alpha_b[alk_idx]  = 2.3e-3
_alpha_b[ca_idx]   = 10.3e-3
_alpha_b[mg_idx]   = 52.8e-3
_alpha_b[na_idx]   = 480e-3
_alpha_b[cl_idx]   = 550e-3
_alpha_b[so4_idx]  = SEAWATER_SO4
_alpha_b[k_idx]    = SEAWATER_K
_alpha_b[si_idx]   = 0.1e-3
_alpha_b[c_idx]    = 2.0e-3

def _alpha_residual(a):
    flux, _ = get_weathering_flux(
        _alpha_P_pore, _alpha_T_pore, _alpha_P_CO2,
        _alpha_b, alpha=float(a[0]), rate=rate_ref, precipitating_minerals=clay_minerals,
        water_rock_ratio=WATER_ROCK_RATIO_LT, pe=planet_module.PE_DEFAULT,
    )
    return (flux[alk_idx] - _alpha_flux_norm) / _alpha_flux_norm

ALPHA_REF_FITTED = float(least_squares(_alpha_residual, [100.0]).x[0])

# Starting points for iterated constants
K_NA_INIT     = planet_module.K_NA_CONT_REMOVAL
# KD_MG_HT: first-order Mg-Ca exchange at HT (scales with J_HT × [Mg]).
# Iterated against the Ca:Mg split -- see calibrate().
KD_MG_INIT    = planet_module.KD_MG_HT
ALPHA_INIT    = ALPHA_REF_CODE   # start from the value the model ships with
# tau_prec is depth-scaled in the model (planet.TAU_PREC_REF, §26), so the calibration must use
# the value Planet would resolve at OCEAN_DEPTH -- otherwise the Earth anchor is fitted at a
# timescale the model never uses here. At 3700 m that is 123 kyr rather than the 100 kyr reference.
TAU_PREC_INIT = planet_module.TAU_PREC_REF * (OCEAN_DEPTH / planet_module.OCEAN_DEPTH_REF)
TAU_RW_INIT   = planet_module.TAU_RW_REF   # start from the value the model ships with

# Starting point OVERRIDES (2026-09-09). The module defaults put x0 on the calcite-COLLAPSED
# branch (Ca ~ 0.36 mM, cost ~11) now that Cl is correct at 546 mM, and the first 3-parameter
# attempt escaped it by running alpha away 1330x. These are the 2-parameter fit's converged
# values, which sit on the Ca-alive branch; alpha starts where the flux is ~1 Tmol/yr at that
# ocean (0.32 measured at 1.57, and the flux is ~linear in alpha, so 1.57/0.32 ~ 4.9).
# Set to None to fall back to whatever constants.py ships (the latest fit).
K_NA_START, KD_MG_START, ALPHA_START = None, None, None

print(f"K_CL (analytic)         = {K_CL_ANALYTIC:.4e}  "
      f"(current in constants.py: {planet_module.K_CL_SUBDUCTION:.4e})")
print(f"K_NA (starting)         = {K_NA_INIT:.4e}")
print(f"KD_MG_HT (starting)     = {KD_MG_INIT:.4e}  (Mg-Ca exchange)")
print(f"tau_prec                = {TAU_PREC_INIT/YR/1e6:.2f} Myr")
print(f"tau_rw (starting)        = {TAU_RW_INIT/YR/1e6:.1f} Myr  (reverse weathering)")
print(f"water/rock ratio        = {WATER_ROCK_RATIO_LT}")
print(f"t_end                   = {T_END/YR/1e9:.1f} Gyr")
print()
print(f"ALPHA_REF in code       = {ALPHA_REF_CODE:.6f}   (used by the runs below)")
print(f"ALPHA_REF refitted      = {ALPHA_REF_FITTED:.6f}   (diagnostic only, w/r={WATER_ROCK_RATIO_LT})")
if ALPHA_REF_CODE > 0 and abs(ALPHA_REF_FITTED / ALPHA_REF_CODE - 1) > 0.10:
    print(f"  ** these differ by {100*(ALPHA_REF_FITTED/ALPHA_REF_CODE - 1):+.0f}% -- the {FLUX_TARGET_NET:g} Teq/yr")
    print(f"     seafloor anchor no longer holds at the current w/r and mineral lists.")
print()


# ---------------------------------------------------------------------------
# Core simulation wrapper
# ---------------------------------------------------------------------------

# SO4 and K are pinned (F_net = 0 in planet.py), so they are set once here at modern seawater values.
print(f"Pinned backgrounds: SO4 = {SEAWATER_SO4*1e3:.1f} mM, K = {SEAWATER_K*1e3:.1f} mM (Millero et al. 2008)")
print()


def make_b0():
    """Initial ocean composition: modern seawater, exactly charge-balanced.

    This REPLACES the previous "blank ocean" seed (Cl and SO4 only, every cation
    and alkalinity zero), which was not a physical ocean and silently broke the
    run.  Seeding 546 mM Cl and 23.45 mM SO4 with no cations puts -592.9 mEq/kg
    of strong anions in the box while the tracked alkalinity starts at 0, so the
    invariant Alk = ION_CHARGE.b is violated by 592.9 mEq at t=0.

    Measured: that offset is 592.900 mEq at t=0 and 592.860 after 1.3 Gyr -- the
    flux terms are charge-perfect (the S20.3 fix works; 0.04 mEq drift per Gyr),
    so nothing in the run can ever repair it.  The carbonate system then sees
    ~+493 mM of alkalinity that no cation supports, giving pH 9.7, DIC 303 mM,
    and calcite supersaturation that strips Ca to 0.20 mM.  Every "Ca is 98%
    low" result from this script predates this fix and was an artifact of it.

    The old seed's rationale -- let species equilibrate from zero -- cannot work
    here: Na alone needs ~500 Myr, and Ca/Mg never reach target at all, so the
    ocean carries the full anion excess for the entire integration.

    Alkalinity is set to the exact ion charge sum rather than to 2.3 mM, so the
    invariant holds identically at t=0 (they agree to 3 decimal places at modern
    concentrations, but deriving it keeps the two definitions from drifting).
    """
    b = np.zeros(N_ELEM)
    b[cl_idx]  = T_Cl
    b[so4_idx] = SEAWATER_SO4
    b[k_idx]   = SEAWATER_K
    b[na_idx]  = T_Na
    b[ca_idx]  = T_Ca
    b[mg_idx]  = T_Mg
    b[c_idx]   = T_C
    b[si_idx]  = 0.1e-3
    b[alk_idx] = float(np.dot(ION_CHARGE, b))
    return b


def run_planet(K_na, KD_mg, alpha, tau_prec, tau_rw, name='calib'):
    """Run from the charge-balanced seawater seed with the given calibration."""
    p = Planet(
        mass=M_EARTH,
        radius=R_EARTH,
        background_pressure=1e5,
        instellation=INSTELLATION,
        crust_production_rate=CRUST_RATE,
        outgassing=OUTGASSING,
        ocean_depth=OCEAN_DEPTH,
        land_fraction=LAND_FRAC,
        reverse_weathering=True,
        alpha=alpha,
        tau_prec=tau_prec,
        tau_rw=tau_rw,
        k_cl_subduction=K_CL_ANALYTIC,
        k_na_cont_removal=K_na,
        kd_mg_ht=KD_mg,
        name=name,
    )

    p.time_evolve(
        t_end=T_END,
        b0=make_b0(),
        convergence_threshold=0.05,
        max_chemistry_fallbacks=MAX_CHEM_FALLBACKS,
    )

    out_path = os.path.join(OUTPUT_DIR, f'{name}.json')
    with open(out_path) as fh:
        data = json.load(fh)

    # data['data']['y'][i] = time series of state variable i
    # Layout: Y[0]=P_CO2, Y[1]=P_H2O, Y[2..N_ELEM+1]=b_ocean (older outputs add Y[-1]=r_avg)
    y = data['data']['y']

    def final(idx):
        return float(y[2 + idx][-1])

    # Drift of pCO2 over the last decade of the run -- the same settling metric
    # used to audit the sweeps. |slope| < 0.05 means the run is at steady state.
    t  = np.array(data['data']['time'], dtype=float)
    Pt = np.array(y[0], dtype=float)
    m  = (t > 0.9 * t[-1]) & (Pt > 0)
    slope = float(np.polyfit(np.log(t[m]), np.log(Pt[m]), 1)[0]) if m.sum() > 3 else float('nan')

    # data['P_CO2'] is in bar (planet.py stores sol.y[0,-1] / 1e5)
    return {
        'Cl':  final(cl_idx),
        'Na':  final(na_idx),
        'Ca':  final(ca_idx),
        'Mg':  final(mg_idx),
        'Alk': final(alk_idx),
        'C':   final(c_idx),
        'pCO2_ppm': float(data['P_CO2']) * 1e6,
        'T':   float(data.get('T', float('nan'))),
        'pH':  float(data.get('pH', float('nan'))),
        'term': data.get('termination', '?'),
        'fab': float(data.get('fabricated_fraction', 0.0)),
        'slope': slope,
        'alk_flux': (data.get('diagnostics') or {}).get('alk_flux'),   # net seafloor Alk, Teq/yr
        'rw_mg_flux': (data.get('diagnostics') or {}).get('rw_mg_flux'),   # RW Mg sink, Tmol/yr
    }


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def seafloor_alk_flux_tmol(result):
    """Net seafloor alkalinity flux (Teq/yr) that dY_dt applied at the run's final state.

    Read from the run's own diagnostics, so the pore clays, redox state and ocean composition are
    exactly the model's. It replaced a separate primary-only re-evaluation at PHREEQC's default pe
    (2026-09-24; development_history.md section 37.16).
    """
    f = result.get('alk_flux')
    return float('nan') if f is None else float(f)


def print_state(label, result, K_na, KD_mg, alpha, tau_prec, tau_rw):
    try:
        sf_alk = f"{seafloor_alk_flux_tmol(result):.3f}"
    except Exception as e:
        sf_alk = f"n/a ({type(e).__name__})"
    bar = '─' * 70
    print(f"\n{bar}")
    print(f"  {label}")
    print(f"  term={result['term']}  T={result['T']:.1f} K  pH={result['pH']:.2f}  "
          f"net seafloor Alk={sf_alk} Teq/yr  (target {FLUX_TARGET:g})")
    print(f"  reverse-weathering Mg sink={result.get('rw_mg_flux') or float('nan'):.4f} Tmol/yr  "
          f"(target {RW_MG_TARGET:g})")
    print(f"  settling: |dlnP/dlnt|={abs(result['slope']):.3f} "
          f"({'AT STEADY STATE' if abs(result['slope']) < 0.05 else 'STILL DRIFTING'})"
          f"   fabricated={result['fab']:.3f}")
    print(bar)
    print(f"  {'Species':6s}  {'Sim (mM)':>10s}  {'Target (mM)':>11s}  {'Error':>8s}  {'Note':s}")
    print(f"  {'-'*55}")
    all_species = [
        ('Cl',  result['Cl'],  T_Cl,  '(analytic K_cl; sim value not used)'),
        ('Na',  result['Na'],  T_Na,  ''),
        ('Ca',  result['Ca'],  T_Ca,  ''),
        ('Mg',  result['Mg'],  T_Mg,  ''),
        ('Alk', result['Alk'], T_Alk, ''),
        ('C',   result['C'],   T_C,   ''),
    ]
    for sp, s, t, note in all_species:
        err  = (s - t) / t * 100
        flag = '  <--' if abs(err) > 10 and sp != 'Cl' else ''
        print(f"  {sp:6s}  {s*1e3:>10.2f}  {t*1e3:>11.2f}  {err:>+7.1f}%{flag}  {note}")
    # Ca+Mg (what alpha controls) and the Ca:Mg split (what KD_MG_HT controls),
    # printed separately so it is visible which knob owns which error.
    s_sum = (result['Ca'] + result['Mg']) / (T_Ca + T_Mg)
    s_rat = ((result['Ca'] / T_Ca) / (result['Mg'] / T_Mg)) if result['Mg'] > 0 else float('inf')
    print()
    print(f"  Ca+Mg / target      = {s_sum:6.3f}   (alpha controls this)")
    print(f"  (Ca/Ca_t)/(Mg/Mg_t) = {s_rat:6.3f}   (KD_MG_HT controls this; 1.0 = correct split)")
    print()
    print(f"  pCO2 = {result['pCO2_ppm']:.1f} ppm  (target: {T_pCO2:.0f}; abiotic model, expect high)")
    print(f"  K_NA      = {K_na:.4e}  (init: {K_NA_INIT:.4e})")
    print(f"  KD_MG_HT  = {KD_mg:.4e}  (init: {KD_MG_INIT:.4e})")
    print(f"  ALPHA     = {alpha:.4e}  (init: {ALPHA_INIT:.4e})")
    print(f"  K_CL      = {K_CL_ANALYTIC:.4e}  (analytic, fixed)")
    print(f"  tau_prec  = {tau_prec/YR/1e6:.3f} Myr")
    print(f"  tau_rw    = {tau_rw/YR/1e6:.2f} Myr  (init: {TAU_RW_INIT/YR/1e6:.2f} Myr)")


def calibrated(result, tol=0.07):
    """Return True when Na, Mg, Ca are all within tol of targets.
    Cl is excluded (residence time >> integration time).
    Alk/C are checked loosely as they follow from the carbonate system."""
    return all(abs(result[sp] / TARGETS[sp] - 1) < tol for sp in ('Na', 'Mg', 'Ca'))


# ---------------------------------------------------------------------------
# Abiotic iteration
# ---------------------------------------------------------------------------

# KD_MG_HT is a rate constant, not a fraction, but a negative or runaway value is
# unphysical. The floor is well below any value that does anything (at 1e-6 the
# exchange is ~5 orders below the LT seafloor Mg flux), so pinning there is a
# meaningful result: it says the Ca:Mg split is NOT controlled by HT exchange.
_KD_LO, _KD_HI       = 1e-6, 10.0
_ALPHA_LO, _ALPHA_HI = 1e-3, 1e4


MAX_RUNS_PER_ROUND = 30   # cap on solver steps per ion fit (scipy excludes Jacobian probes from max_nfev)
ALPHA_ROUNDS = 6           # cap on alpha / tau_rw rescalings
FLUX_TOL = 0.03            # accept |ln(flux / FLUX_TARGET)| below this (3 %)
RW_TOL = 0.05              # accept |ln(rw / RW_MG_TARGET)| below this (5 %)
_TAU_RW_LO, _TAU_RW_HI = 1e4 * YR, 1e10 * YR
DIFF_STEP_LN = 0.05        # finite-difference step in ln(parameter), ~5 %; scipy scales diff_step by max(1, |x|)
X_SHIFT = 20.0             # x = ln(p) + X_SHIFT > 0, so scipy's probes step UP in K_na and KD_mg, away from Ca collapse

# ---------------------------------------------------------------------------
# Alternating fit (2026-09-24, development_history.md section 37.20).
#   inner: least_squares on (K_na, KD_mg) against (Na, Ca, Mg) at fixed alpha;
#   outer: alpha <- alpha * FLUX_TARGET / flux, since the flux is ~linear in alpha (exponent 0.99)
#          and alpha moves the ions by ~1 % over a 1.7x step;
#          tau_rw <- tau_rw * rw / RW_MG_TARGET, since Sepiolite stays far supersaturated, so the
#          reverse-weathering flux is ~ excess / tau_rw (section 26.3).
# The joint 3-parameter fit failed: diff_step=0.2 on ln(x) gave x0.36 / x0.47 trial steps in K_na
# and KD_mg that crossed onto the Ca-collapsed branch, and the flux residual dominated the cost, so
# it traded Mg (26 mM) for flux and returned a finite-difference probe (alpha 55.4, 0.40 Teq/yr).
# ALPHA_PINNED: skip the outer loop and fit the ions at this alpha. None = rescale alpha to the flux.
# ---------------------------------------------------------------------------
ALPHA_PINNED = None
FLUX_TARGET = FLUX_TARGET_NET   # Teq/yr, NET seafloor alkalinity flux
# Tmol Mg/yr into authigenic clays in typical deep-sea sediment (Dunlea et al. 2017, Nat. Commun. 8, 844);
# 0.4-0.8 if Si-rich sedimentation covered 50-100 % of the seafloor. None = hold tau_rw at TAU_RW_INIT.
RW_MG_TARGET = 0.02

_history = []   # (cost, K_na, KD_mg, alpha, result) for every successful evaluation
_alpha_now = None   # alpha held fixed during the current ion fit
_tau_rw_now = TAU_RW_INIT   # tau_rw held fixed during the current ion fit


def _residuals(x):
    """Log-space residuals in (Na, Ca, Mg) for x = [ln K_na, ln KD_mg] + X_SHIFT at alpha = _alpha_now."""
    K_na, KD_mg = (float(v) for v in np.exp(np.asarray(x) - X_SHIFT))
    alpha = _alpha_now
    name = f'calib_ls_{len(_history):03d}'
    try:
        r = run_planet(K_na, KD_mg, alpha, TAU_PREC_INIT, _tau_rw_now, name=name)
    except Exception as e:
        print(f"    [eval {len(_history):03d}] FAILED {type(e).__name__}: {str(e)[:60]}")
        return np.full(3, 5.0)   # finite penalty; keeps the solver moving

    res = np.array([np.log(max(r[s], 1e-12) / TARGETS[s]) for s in ('Na', 'Ca', 'Mg')])
    try:
        f_sf = seafloor_alk_flux_tmol(r)
    except Exception:
        f_sf = float('nan')

    cost = float(np.sum(res**2))
    _history.append((cost, K_na, KD_mg, alpha, r))
    print(f"    [eval {len(_history)-1:03d}] K_na={K_na:.3e} alpha={alpha:.3e} kd={KD_mg:.3e}"
          f" tau_rw={_tau_rw_now/YR/1e6:.3g}Myr"
          f"  ->  Na={r['Na']*1e3:7.1f} Ca={r['Ca']*1e3:7.2f} Mg={r['Mg']*1e3:7.2f}"
          f"  Alk={r['Alk']*1e3:6.2f} sfAlk={f_sf:6.3f} rwMg={r.get('rw_mg_flux') or float('nan'):.4f}"
          f"  cost={cost:.4f}")
    return res


def fit_ions(K_na, KD_mg, alpha):
    """least_squares on (K_na, KD_mg) against (Na, Ca, Mg) at fixed alpha; returns the best point.

    Ca is buffered by calcite saturation, so the ocean snaps between a Ca-alive and a Ca-collapsed
    branch (Ca ~0.7 mM, Na ~1500 mM). The ~5 % finite-difference steps keep the Jacobian probes on
    the branch the iterate is on; the old 0.2 x |ln x| steps did not. Even so, a 5 % DROP in K_na tips
    Earth onto the collapsed branch, so X_SHIFT makes every probe step upwards.
    """
    global _alpha_now
    _alpha_now = alpha
    n0 = len(_history)
    x0 = np.log([K_na, KD_mg]) + X_SHIFT
    lo, hi = np.log([1e-8, _KD_LO]) + X_SHIFT, np.log([1e2, _KD_HI]) + X_SHIFT
    try:
        least_squares(_residuals, x0, bounds=(lo, hi), diff_step=DIFF_STEP_LN / x0,
                      max_nfev=MAX_RUNS_PER_ROUND, xtol=1e-3, ftol=1e-3, gtol=1e-3)
    except Exception as e:
        print(f"\n  least_squares aborted ({type(e).__name__}: {e}); using best seen.")
    this_round = _history[n0:]
    if not this_round:
        raise RuntimeError(f"no successful evaluations at alpha = {alpha:g}")
    # least_squares can finish at a point worse than one it visited, so report the best.
    cost, K_na_b, KD_mg_b, _, res_b = min(this_round, key=lambda h: h[0])
    return K_na_b, KD_mg_b, res_b


def calibrate(K_na, KD_mg, alpha, tau_rw):
    """Alternate the ion fit with alpha and tau_rw rescaled to their flux targets until both hold."""
    global _tau_rw_now
    _tau_rw_now = tau_rw
    print(f"\n{'#'*70}")
    if ALPHA_PINNED is None:
        print(f"  Abiotic calibration: (K_na, KD_mg) vs (Na, Ca, Mg); alpha vs net seafloor Alk "
              f"-> {FLUX_TARGET:g} Teq/yr")
    else:
        print(f"  Abiotic calibration: (K_na, KD_mg) vs (Na, Ca, Mg), alpha PINNED at {ALPHA_PINNED:g}")
    if RW_MG_TARGET is not None:
        print(f"  tau_rw vs reverse-weathering Mg sink -> {RW_MG_TARGET:g} Tmol/yr")
    print(f"  targets: Na={T_Na*1e3:.0f} Ca={T_Ca*1e3:.1f} Mg={T_Mg*1e3:.1f} mM")
    print(f"{'#'*70}")

    if ALPHA_PINNED is not None:
        alpha = ALPHA_PINNED
    for k in range(ALPHA_ROUNDS):
        print(f"\n  -- round {k + 1}: alpha = {alpha:.4g}, tau_rw = {_tau_rw_now/YR/1e6:.4g} Myr")
        tau_fit = _tau_rw_now
        K_na, KD_mg, res = fit_ions(K_na, KD_mg, alpha)
        flux = seafloor_alk_flux_tmol(res)
        rw = res.get('rw_mg_flux')
        rw = float('nan') if rw is None else float(rw)
        print(f"  round {k + 1}: best K_na={K_na:.4e} KD_mg={KD_mg:.4e}  "
              f"net seafloor Alk={flux:.3f} Teq/yr  RW Mg sink={rw:.4f} Tmol/yr")
        flux_ok = ALPHA_PINNED is not None or abs(np.log(flux / FLUX_TARGET)) < FLUX_TOL
        rw_ok = RW_MG_TARGET is None or abs(np.log(rw / RW_MG_TARGET)) < RW_TOL
        if flux_ok and rw_ok:
            return K_na, KD_mg, alpha, tau_fit, res
        if ALPHA_PINNED is None:
            if not (np.isfinite(flux) and flux > 0):
                raise RuntimeError(f"non-positive seafloor flux at alpha = {alpha:g}")
            alpha = float(np.clip(alpha * FLUX_TARGET / flux, _ALPHA_LO, _ALPHA_HI))
        if RW_MG_TARGET is not None:
            if not (np.isfinite(rw) and rw > 0):
                raise RuntimeError(f"no reverse-weathering Mg sink at tau_rw = {_tau_rw_now/YR/1e6:g} Myr")
            _tau_rw_now = float(np.clip(_tau_rw_now * rw / RW_MG_TARGET, _TAU_RW_LO, _TAU_RW_HI))
    print(f"\n  *** alpha / tau_rw did not converge in {ALPHA_ROUNDS} rounds; the ions below are fitted at "
          f"the last alpha and tau_rw, whose fluxes are off target ***")
    return K_na, KD_mg, _alpha_now, tau_fit, res


K_na, KD_mg, alpha, tau_rw, result1 = calibrate(K_NA_START or K_NA_INIT,
                                                KD_MG_START or KD_MG_INIT,
                                                ALPHA_START or ALPHA_INIT,
                                                TAU_RW_INIT)
tau_prec = TAU_PREC_INIT

if not calibrated(result1):
    print("\n  *** Did not reach all three targets within tolerance ***")
if KD_mg <= _KD_LO * 1.001:
    print("\n  NOTE: KD_MG_HT pinned at its floor. The Ca:Mg split is not controlled by")
    print("        HT exchange -- the residual Mg deficit is set by another sink")
    print("        (reverse-weathering clays / LT seafloor). Calibrating it further is futile.")

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------

print("\n\n" + "="*70)
print("  CALIBRATION COMPLETE — FINAL CONSTANTS (ABIOTIC)")
print("="*70)
print_state("Best result", result1, K_na, KD_mg, alpha, tau_prec, tau_rw)

print("""
  ── Paste into src/kamino/constants.py ────────────────────────────────""")
print(f"  KD_MG_HT = {KD_mg:.6e}")
print(f"  K_NA_CONT_REMOVAL = {K_na:.6e}")
print(f"  K_CL_SUBDUCTION = {K_CL_ANALYTIC:.6e}")
print(f"  ALPHA_REF = {alpha:.6f}   (was {ALPHA_INIT:.6f})")
print(f"  TAU_RW_REF = {tau_rw/YR:.6e} * YR   (was {TAU_RW_INIT/YR:.6e} * YR)")
print("""
  ── Planet constructor defaults ───────────────────────────────────────""")
print(f"  tau_prec = {tau_prec/YR:.4e} * YR   # {tau_prec/YR/1e6:.3f} Myr")
print()
