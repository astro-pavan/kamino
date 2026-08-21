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

Three constants are iterated: K_na (Na balance), alpha (total Ca+Mg supply) and
KD_MG_HT (the Ca:Mg split).  They are close to orthogonal on those three targets
-- see calibrate().  tau_prec / tau_rw are held at their literature values.

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
from kamino.chemistry import elements, ION_CHARGE, alk_idx, c_idx, si_idx, ca_idx, mg_idx, na_idx, cl_idx, so4_idx
from kamino.weathering import get_weathering_flux, ALPHA_REF as ALPHA_REF_CODE
from kamino.constants import (
    EARTH_HYDROTHERMAL_FLUX_PER_AREA as J_ref_normalised,
    EARTH_CRUST_PRODUCTION_RATE_PER_AREA as rate_ref,
    A_SEAFLOOR_EARTH as A_seafloor,
)
from kamino.constants import M_EARTH, R_EARTH, YR, EARTH_OUTGASSING, EARTH_ATM, G

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
#   (EARTH_OUTGASSING/YR) * A * cl_ratio / M  =  K_cl * [Cl] * J_total * A / M
#
# A and M cancel, leaving:
#   K_cl = (EARTH_OUTGASSING/YR * cl_ratio) / ([Cl]_target * J_total)
#
# J_total = J_ref_normalised for crust_rate=1.
# ---------------------------------------------------------------------------
K_CL_ANALYTIC = (EARTH_OUTGASSING / YR * CL_OUTGASSING_RATIO) / (T_Cl * J_ref_normalised)

# ---------------------------------------------------------------------------
# ALPHA_REF diagnostic: the seafloor reactive-area scaling that would make the
# primary-dissolution alkalinity flux equal 1 Tmol/yr at Earth pore conditions
# with modern seawater as input.
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
_alpha_flux_norm = (1e12 / YR) / A_seafloor   # 1 Tmol/yr normalised per m²

_alpha_b = np.zeros(len(elements))
_alpha_b[alk_idx]  = 2.3e-3
_alpha_b[ca_idx]   = 10.3e-3
_alpha_b[mg_idx]   = 52.8e-3
_alpha_b[na_idx]   = 480e-3
_alpha_b[cl_idx]   = 550e-3
_alpha_b[so4_idx]  = 28e-3
_alpha_b[si_idx]   = 0.1e-3
_alpha_b[c_idx]    = 2.0e-3

def _alpha_residual(a):
    flux, _ = get_weathering_flux(
        _alpha_P_pore, _alpha_T_pore, _alpha_P_CO2,
        _alpha_b, alpha=float(a[0]), rate=rate_ref, precipitating_minerals=[],
        water_rock_ratio=WATER_ROCK_RATIO_LT,
    )
    return (flux[alk_idx] - _alpha_flux_norm) / _alpha_flux_norm

ALPHA_REF_FITTED = float(least_squares(_alpha_residual, [1.43]).x[0])

# Starting points for iterated constants
K_NA_INIT     = planet_module.K_NA_CONT_REMOVAL
# KD_MG_HT: first-order Mg-Ca exchange at HT (scales with J_HT × [Mg]).
# Iterated against the Ca:Mg split -- see calibrate().
KD_MG_INIT    = planet_module.KD_MG_HT
ALPHA_INIT    = ALPHA_REF_CODE   # start from the value the model ships with
TAU_PREC_INIT = 100e3 * YR   # literature alkalinity residence time ~100 kyr
TAU_RW_INIT   = 5e6 * YR     # reverse weathering timescale (secondary Mg control)

print(f"K_CL (analytic)         = {K_CL_ANALYTIC:.4e}  "
      f"(current in planet.py: {planet_module.K_CL_SUBDUCTION:.4e})")
print(f"K_NA (starting)         = {K_NA_INIT:.4e}")
print(f"KD_MG_HT (fixed)        = {KD_MG_INIT:.4e}  (Mg-Ca exchange; not iterated)")
print(f"tau_prec                = {TAU_PREC_INIT/YR/1e6:.2f} Myr")
print(f"tau_rw                  = {TAU_RW_INIT/YR/1e6:.1f} Myr  (reverse weathering)")
print(f"water/rock ratio        = {WATER_ROCK_RATIO_LT}")
print(f"t_end                   = {T_END/YR/1e9:.1f} Gyr")
print()
print(f"ALPHA_REF in code       = {ALPHA_REF_CODE:.6f}   (used by the runs below)")
print(f"ALPHA_REF refitted      = {ALPHA_REF_FITTED:.6f}   (diagnostic only, w/r={WATER_ROCK_RATIO_LT})")
if ALPHA_REF_CODE > 0 and abs(ALPHA_REF_FITTED / ALPHA_REF_CODE - 1) > 0.10:
    print(f"  ** these differ by {100*(ALPHA_REF_FITTED/ALPHA_REF_CODE - 1):+.0f}% -- the 1 Tmol/yr")
    print(f"     seafloor anchor no longer holds at the current w/r and mineral lists.")
print()


# ---------------------------------------------------------------------------
# Core simulation wrapper
# ---------------------------------------------------------------------------

# SO4 background: computed from charge balance at target concentrations.
# SO4 is pinned (F_net[so4_idx]=0) so it must be set in the initial condition.
# Real seawater SO4 is ~28 mM, but the model omits K+ (~10.2 mEq/kg cation),
# so the effective background is lower to give the correct alkalinity:
#   SO4_bg = (2[Ca]_t + 2[Mg]_t + [Na]_t - [Cl]_t - [Alk]_t) / 2
SO4_BG = (2*T_Ca + 2*T_Mg + T_Na - T_Cl - T_Alk) / 2   # ≈ 23.45 mM

print(f"SO4 background (computed) = {SO4_BG*1e3:.2f} mM  "
      f"(real seawater ~28 mM; lower because K+ is absent from model)")
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
    b[so4_idx] = SO4_BG
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
    # Layout: Y[0]=P_CO2, Y[1]=P_H2O, Y[2..N_ELEM+1]=b_ocean, Y[-1]=r_avg
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
    }


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def seafloor_alk_flux_tmol(result, alpha=None):
    """Primary seafloor alkalinity flux (Tmol/yr) at the run's final conditions.

    Uses the run's own alpha and the model's water/rock ratio so the number is
    comparable to what dY_dt actually saw (it previously used alpha=1.0 and
    w/r=None, which matched no configuration the model runs in).
    """
    if alpha is None:
        alpha = ALPHA_REF_CODE
    b = np.zeros(N_ELEM)
    b[alk_idx] = result['Alk']
    b[c_idx]   = result['C']
    b[si_idx]  = 0.1e-3
    b[ca_idx]  = result['Ca']
    b[mg_idx]  = result['Mg']
    b[na_idx]  = result['Na']
    b[cl_idx]  = result['Cl']
    b[so4_idx] = 28e-3

    T_sf   = max(1.02 * result['T'] - 16.7, 274.0)
    T_pore = T_sf + 9.0
    P_CO2_Pa = result['pCO2_ppm'] * 1e-6 * 1e5
    P_pore = 1e5 + P_CO2_Pa + 1000.0 * GRAVITY * OCEAN_DEPTH

    flux, _ = get_weathering_flux(
        P_pore, T_pore, P_CO2_Pa, b,
        alpha=alpha,
        rate=rate_ref,
        precipitating_minerals=[],  # primary dissolution only
        water_rock_ratio=WATER_ROCK_RATIO_LT,
    )
    return flux[alk_idx] * A_seafloor * YR / 1e12   # Tmol/yr


def print_state(label, result, K_na, KD_mg, alpha, tau_prec, tau_rw):
    try:
        sf_alk = f"{seafloor_alk_flux_tmol(result, alpha):.2f}"
    except Exception as e:
        sf_alk = f"n/a ({type(e).__name__})"
    bar = '─' * 70
    print(f"\n{bar}")
    print(f"  {label}")
    print(f"  term={result['term']}  T={result['T']:.1f} K  pH={result['pH']:.2f}  "
          f"seafloor Alk={sf_alk} Tmol/yr  (target ~1)")
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
    print(f"  tau_rw    = {tau_rw/YR/1e6:.1f} Myr")


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


MAX_RUNS = 60   # hard cap on planet integrations spent by the solver

_history = []   # (cost, K_na, KD_mg, alpha, result) for every successful evaluation


def _residuals(x):
    """Log-space residuals in (Na, Ca, Mg) for log-parameters x = [ln K_na, ln alpha, ln KD_mg].

    Logs on both sides. The parameters span decades, and the targets differ by a factor
    of 50 (Ca 10.3 mM vs Na 469 mM), so a linear residual would let Na dominate the fit
    entirely; log residuals weight all three by relative error, which is what we want.
    """
    K_na, alpha, KD_mg = (float(v) for v in np.exp(x))
    name = f'calib_ls_{len(_history):03d}'
    try:
        r = run_planet(K_na, KD_mg, alpha, TAU_PREC_INIT, TAU_RW_INIT, name=name)
    except Exception as e:
        print(f"    [eval {len(_history):03d}] FAILED {type(e).__name__}: {str(e)[:60]}")
        return np.array([5.0, 5.0, 5.0])   # finite penalty; keeps the solver moving

    res = np.array([np.log(max(r[s], 1e-12) / TARGETS[s]) for s in ('Na', 'Ca', 'Mg')])
    cost = float(np.sum(res**2))
    _history.append((cost, K_na, KD_mg, alpha, r))
    print(f"    [eval {len(_history)-1:03d}] K_na={K_na:.3e} alpha={alpha:.3e} kd={KD_mg:.3e}"
          f"  ->  Na={r['Na']*1e3:7.1f} Ca={r['Ca']*1e3:7.2f} Mg={r['Mg']*1e3:7.2f}"
          f"  Alk={r['Alk']*1e3:6.2f}  cost={cost:.4f}")
    return res


def calibrate(K_na, KD_mg, alpha):
    """Bounded least-squares solve for (K_na, alpha, KD_mg) against (Na, Ca, Mg).

    Replaces the previous hand-rolled ratio-scaling loop, which OSCILLATED rather than
    converged. Measured 2-cycle: kd 2.65e-2 -> 1.62e-2 -> 6.95e-2 -> 1.62e-2, with the
    Ca:Mg ratio swinging 2.7 -> 0.054 -> 18.4 -> 0.048. A 4x change in kd moves the
    ratio by 380x (log-log sensitivity ~4.1), so any damping exponent above ~0.24
    amplifies the error instead of shrinking it.

    That sensitivity is physical, not numerical: Ca is buffered by calcite saturation.
    Ca and alkalinity trade off along Ca.CO3 = Ksp (measured Ca=0.55 mM at Alk=13.5 mM
    versus Ca=55 mM at Alk=1.8 mM), so the system snaps between "calcite precipitates,
    Ca -> 0" and "calcite undersaturated, Ca accumulates". A gradient method with a
    proper trust region handles that; independent per-species ratio updates cannot,
    because each knob's correct step depends on where the others sit relative to the
    saturation boundary.

    Knob roles remain as before -- K_na sets Na (first-order sink), alpha sets the total
    Ca+Mg supply (seafloor reactive area), KD_mg sets the Ca:Mg split (HT exchange moves
    Mg to Ca mole-for-mole) -- but they are solved jointly rather than one-at-a-time.

    HISTORY: KD_mg used to be pinned here, on the measurement that HT exchange was 1.5%
    of Mg removal and 0.6% of the Ca source at the Earth steady state. That was measured
    when make_b0 violated Alk = ION_CHARGE.b by 592.9 mEq/kg, which pinned Ca near
    0.20 mM by calcite supersaturation and made every knob look inert. With the seed
    fixed, Ca responds and the pin no longer applies.
    """
    print(f"\n{'#'*70}")
    print(f"  Abiotic calibration — least_squares on (K_na, alpha, KD_mg)")
    print(f"  targets: Na={T_Na*1e3:.0f} Ca={T_Ca*1e3:.1f} Mg={T_Mg*1e3:.1f} mM   "
          f"max {MAX_RUNS} runs")
    print(f"{'#'*70}")

    x0 = np.log([K_na, alpha, KD_mg])
    lo = np.log([1e-8, _ALPHA_LO, _KD_LO])
    hi = np.log([1e2,  _ALPHA_HI, _KD_HI])

    try:
        least_squares(_residuals, x0, bounds=(lo, hi), diff_step=0.2,
                      max_nfev=MAX_RUNS, xtol=1e-3, ftol=1e-3, gtol=1e-3)
    except Exception as e:
        print(f"\n  least_squares aborted ({type(e).__name__}: {e}); using best seen.")

    if not _history:
        raise RuntimeError("no successful evaluations")

    # least_squares can finish at a point worse than one it visited, so report the best.
    cost, K_na_b, KD_mg_b, alpha_b, res_b = min(_history, key=lambda h: h[0])
    print(f"\n  Best of {len(_history)} evaluations: cost={cost:.4f}")
    return K_na_b, KD_mg_b, alpha_b, res_b


K_na, KD_mg, alpha, result1 = calibrate(K_NA_INIT, KD_MG_INIT, ALPHA_INIT)
tau_prec, tau_rw = TAU_PREC_INIT, TAU_RW_INIT

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
  ── Paste into planet.py ──────────────────────────────────────────────""")
print(f"  K_CL_SUBDUCTION         = {K_CL_ANALYTIC:.6e}")
print(f"  K_NA_CONT_REMOVAL       = {K_na:.6e}")
print(f"  KD_MG_HT                = {KD_mg:.6e}")
print("""
  ── Paste into weathering.py ──────────────────────────────────────────""")
print(f"  ALPHA_REF               = {alpha:.6f}")
print(f"    (was {ALPHA_INIT:.6f}; Planet's `alpha` default mirrors this and must match)")
print("""
  ── Planet constructor defaults ───────────────────────────────────────""")
print(f"  tau_prec = {tau_prec/YR:.4e} * YR   # {tau_prec/YR/1e6:.3f} Myr")
print(f"  tau_rw   = {tau_rw/YR:.4e} * YR   # {tau_rw/YR/1e6:.1f} Myr  (reverse weathering)")
print()
