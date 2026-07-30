#!/usr/bin/env python3
"""
Iterative calibration of Earth-reference ocean chemistry constants.

Adjusts K_CL_SUBDUCTION, K_NA_ALBITIZATION, f_HT, and tau_prec so
that Planet.time_evolve() reaches modern seawater concentrations at steady state,
starting from a blank ocean (b0=None) with no pre-seeding.

Key design choices
------------------
K_cl is determined analytically from the Cl flux balance and never updated from
the simulation.  Cl has a residence time of ~10 Gyr so it cannot equilibrate from
a blank ocean within a realistic integration time.  More importantly, Cl is the
dominant charge-balancing anion: if [Cl] is far below its target the carbonate
system compensates with ~400 mEq/kg of spurious alkalinity, making Alk/C/pCO2
calibration impossible.  The initial condition therefore seeds [Cl] to its
analytically-known steady-state value while leaving all other species at zero.
All other target species (Na τ~500 Myr, Mg τ~50 Myr, Ca τ~1 Myr, Alk τ~0.1 Myr)
do equilibrate within 2.5 Gyr from a blank start.

Because Cl's high fractional rate of change prevents the convergence event from
firing, simulations typically terminate with "timeout" rather than "converged".
This is expected and acceptable; the fast species are at steady state by t_end.

Phase 1: abiotic (f_carb=0, f_bio=0) — calibrate K_na, f_HT, tau_prec.
Phase 2: scan f_bio to close the C budget and any residual Ca error.
         f_bio controls both biogenic CaCO3 burial and organic C burial
         (f_bio × F_C_outgas removed as organic matter).  The literature
         value for the organic burial fraction is f_org ≈ 0.32–0.41
         (Kipp et al. 2021; Derry 2024), consistent with what is needed
         to balance the carbon budget at EARTH_OUTGASSING = 0.0147 mol/m²/yr.

Usage:
    /data/pt426/big-venv/bin/python tests/calibrate_earth.py
"""
import sys
import os
import json
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import kamino.planet as planet_module
from kamino.planet import Planet
from kamino.chemistry import elements, alk_idx, c_idx, si_idx, ca_idx, mg_idx, na_idx, cl_idx, so4_idx
from kamino.weathering import get_weathering_flux
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
# ALPHA_REF: calibrate seafloor weathering strength so that the alkalinity
# flux from primary basalt dissolution equals 1 Tmol/yr at Earth steady-state
# pore conditions with modern seawater as input.
# T_surface=288K → T_seafloor=277K → T_pore=286K; depth=3000m.
# Secondary precipitation excluded to isolate primary dissolution.
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
    )
    return (flux[alk_idx] - _alpha_flux_norm) / _alpha_flux_norm

ALPHA_REF = float(least_squares(_alpha_residual, [100.0]).x[0])

# Starting points for iterated constants
K_NA_INIT     = planet_module.K_NA_ALBITIZATION
# KD_MG_HT: first-order Mg-Ca exchange at HT (scales with J_HT × [Mg]).
# Balances Mg: larger KD → more Mg removal → lower steady-state [Mg].
KD_MG_INIT    = planet_module.KD_MG_HT
# f_HT: controls Ca from HT PHREEQC (Anorthite/Diopside dissolution).
F_HT_INIT     = 0
TAU_PREC_INIT = 100e3 * YR   # literature alkalinity residence time ~100 kyr
TAU_RW_INIT   = 5e6 * YR     # reverse weathering timescale (secondary Mg control)

print(f"K_CL (analytic)         = {K_CL_ANALYTIC:.4e}  "
      f"(current in planet.py: {planet_module.K_CL_SUBDUCTION:.4e})")
print(f"K_NA (starting)         = {K_NA_INIT:.4e}")
print(f"KD_MG_HT (starting)     = {KD_MG_INIT:.4e}  (Mg-Ca exchange, Mg balance)")
print(f"f_HT (starting)         = {F_HT_INIT:.4f}  (PHREEQC HT Ca source)")
print(f"tau_prec (starting)     = {TAU_PREC_INIT/YR/1e6:.2f} Myr")
print(f"tau_rw (starting)       = {TAU_RW_INIT/YR/1e6:.1f} Myr  (reverse weathering)")
print(f"ALPHA_REF (calibrated)  = {ALPHA_REF:.4f}  (target: 1 Tmol/yr seafloor Alk)")
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
    """Initial ocean composition: zeros everywhere except Cl and SO4.

    Cl is seeded to its analytically-derived steady-state value.
    SO4 is seeded to the background value computed from the charge balance
    at target concentrations (accounting for the absence of K+ from the
    model).  Both are physically motivated, not arbitrary: Cl was accumulated
    over ~4 Gyr and K_cl is analytically determined; SO4 is the constant
    background that makes the carbonate charge balance self-consistent.
    All other species start at zero (blank ocean).
    """
    b = np.zeros(N_ELEM)
    b[cl_idx]  = T_Cl
    b[so4_idx] = SO4_BG
    return b


def run_planet(K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio=0.0, name='calib'):
    """Run from a Cl-seeded blank ocean with the given calibration, return result dict."""
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
        tau_prec=tau_prec,
        tau_rw=tau_rw,
        f_bio=f_bio,
        f_HT=f_HT,
        k_cl_subduction=K_CL_ANALYTIC,
        k_na_albitization=K_na,
        kd_mg_ht=KD_mg,
        name=name,
    )

    # Run to 1 Gyr.  Ca (τ~1 Myr), Mg (τ~50 Myr), Na (τ~500 Myr) reach SS.
    # Cl stays near T_Cl throughout (K_cl is analytic; τ_Cl~10 Gyr so
    # drift is negligible).  Termination is usually "timeout" because
    # Cl's tiny residual fractional rate prevents the convergence event;
    # the fast species are at steady state regardless.
    p.time_evolve(
        t_end=1e9 * YR,
        b0=make_b0(),
        convergence_threshold=0.05,
    )

    out_path = os.path.join(OUTPUT_DIR, f'{name}.json')
    with open(out_path) as fh:
        data = json.load(fh)

    # data['data']['y'][i] = time series of state variable i
    # Layout: Y[0]=P_CO2, Y[1]=P_H2O, Y[2..N_ELEM+1]=b_ocean, Y[-1]=r_avg
    y = data['data']['y']

    def final(idx):
        return float(y[2 + idx][-1])

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
    }


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def seafloor_alk_flux_tmol(result):
    """Estimate primary seafloor alkalinity flux (Tmol/yr) at steady-state conditions."""
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
        alpha=1.0,
        rate=rate_ref,
        precipitating_minerals=[],  # primary dissolution only
    )
    return flux[alk_idx] * A_seafloor * YR / 1e12   # Tmol/yr


def print_state(label, result, K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio=0.0):
    sf_alk = seafloor_alk_flux_tmol(result)
    bar = '─' * 70
    print(f"\n{bar}")
    print(f"  {label}")
    print(f"  term={result['term']}  T={result['T']:.1f} K  pH={result['pH']:.2f}  "
          f"seafloor Alk={sf_alk:.2f} Tmol/yr  (target ~1)")
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
    print()
    print(f"  pCO2 = {result['pCO2_ppm']:.1f} ppm  (target: {T_pCO2:.0f})")
    print(f"  K_NA      = {K_na:.4e}  (init: {K_NA_INIT:.4e})")
    print(f"  KD_MG_HT  = {KD_mg:.4e}  (init: {KD_MG_INIT:.4e}; Mg balance)")
    print(f"  f_HT      = {f_HT:.4f}   (init: {F_HT_INIT:.4f}; Ca balance)")
    print(f"  K_CL      = {K_CL_ANALYTIC:.4e}  (analytic, fixed)")
    print(f"  tau_prec  = {tau_prec/YR/1e6:.3f} Myr")
    print(f"  tau_rw    = {tau_rw/YR/1e6:.1f} Myr")
    print(f"  f_bio     = {f_bio:.2f}")


def calibrated(result, tol=0.07):
    """Return True when Na, Mg, Ca are all within tol of targets.
    Cl is excluded (residence time >> integration time).
    Alk/C are checked loosely as they follow from the carbonate system."""
    return all(abs(result[sp] / TARGETS[sp] - 1) < tol for sp in ('Na', 'Mg', 'Ca'))


# ---------------------------------------------------------------------------
# Phase 1: abiotic iteration
# ---------------------------------------------------------------------------

def phase1(K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio=0.0, n_iter=6, tag='p1', fix_tau=True):
    """Iteratively update K_na, KD_mg, and f_HT until Na/Mg/Ca converge.

    Update rules (ratio-scaling):
      K_na   — Na  balance: albitization sink
      KD_mg  — Mg  balance: first-order HT Mg-Ca exchange (independent of Si)
      f_HT   — Ca  balance: HT PHREEQC Ca source from Anorthite/Diopside
      tau_rw — fixed; secondary Mg sink via reverse-weathering clays
    """
    print(f"\n{'#'*70}")
    print(f"  Phase 1  f_bio={f_bio:.2f}  tag={tag}")
    print(f"{'#'*70}")

    for i in range(n_iter):
        name   = f'calib_{tag}_i{i:02d}'
        result = run_planet(K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio=f_bio, name=name)
        print_state(f"Iter {i+1}/{n_iter}", result, K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio)

        if calibrated(result):
            print(f"\n  *** Converged after {i+1} iteration(s) ***")
            return K_na, KD_mg, f_HT, tau_prec, tau_rw, result

        K_na  = K_na  * (result['Na'] / T_Na)
        # KD_mg: too-high Mg → increase removal rate; too-low → decrease. Damped (0.6).
        KD_mg = max(KD_mg * (result['Mg'] / T_Mg) ** 0.6, 1e-6)
        # f_HT: too-low Ca → increase; too-high → decrease. Damped (0.5). Clamped.
        # f_HT  = min(max(f_HT * (T_Ca / max(result['Ca'], 1e-6)) ** 0.5, 1e-3), 0.5)

        if not fix_tau:
            tau_prec = tau_prec * (T_Ca / max(result['Ca'], 1e-9)) ** 0.3

    print(f"\n  *** Did not fully converge in {n_iter} iterations ***")
    return K_na, KD_mg, f_HT, tau_prec, tau_rw, result # type: ignore


K_na, KD_mg, f_HT, tau_prec, tau_rw, result1 = phase1(
    K_NA_INIT, KD_MG_INIT, F_HT_INIT, TAU_PREC_INIT, TAU_RW_INIT, f_bio=0.0, tag='abiotic', fix_tau=True
)

print("\n" + "="*70)
print("  PHASE 1 (ABIOTIC) COMPLETE")
print("="*70)
print_state("Abiotic final", result1, K_na, KD_mg, f_HT, tau_prec, tau_rw, f_bio=0.0)

# ---------------------------------------------------------------------------
# Phase 2: scan f_bio to close residual Ca/Alk error.
#
# Biogenic CaCO3 burial is the dominant Ca sink in the modern ocean; without
# it the abiotic model will likely over-predict [Ca].  f_bio adds a rate-
# limited Ca+Alk sink calibrated to the observed total Ca sources at Earth
# conditions.  We scan f_bio in steps, re-running a shortened Phase 1 (5
# iterations) from the abiotic-converged constants for each trial value.
# ---------------------------------------------------------------------------

def score(result):
    return (abs(result['Ca']  / T_Ca   - 1) +
            0.5 * abs(result['Alk'] / T_Alk  - 1) +
            0.5 * abs(result['pCO2_ppm'] / T_pCO2 - 1))


best = dict(score=score(result1), K_na=K_na, KD_mg=KD_mg, f_HT=f_HT,
            tau_prec=tau_prec, tau_rw=tau_rw, f_bio=0.0, result=result1)

# Enter Phase 2 if Ca or pCO2 are off — with f_bio=0 the C budget is open
# and pCO2 will typically be >> 280 ppm regardless of Ca.
if result1['Ca'] > T_Ca * 1.05 or result1['pCO2_ppm'] > T_pCO2 * 1.5:
    print(f"\n\n{'#'*70}")
    print("  PHASE 2 — scanning f_bio (organic C burial + biogenic CaCO3)")
    print(f"  Abiotic [Ca]  = {result1['Ca']*1e3:.2f} mM  (target {T_Ca*1e3:.2f} mM)")
    print(f"  Abiotic pCO2  = {result1['pCO2_ppm']:.0f} ppm  (target {T_pCO2:.0f} ppm)")
    print(f"{'#'*70}")

    for fbio in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]:
        print(f"\n  --- f_bio = {fbio:.1f} ---")
        K_na_t, KD_mg_t, f_HT_t, tau_t, tau_rw_t, res_t = phase1(
            K_na, KD_mg, f_HT, tau_prec, tau_rw,
            f_bio=fbio, n_iter=5,
            tag=f'bio{int(fbio * 10):02d}',
            fix_tau=True,
        )
        s = score(res_t)
        print(f"  f_bio={fbio:.1f}  Ca={res_t['Ca']*1e3:.2f} mM  "
              f"Alk={res_t['Alk']*1e3:.2f} mM  score={s:.4f}")

        if s < best['score']:
            best = dict(score=s, K_na=K_na_t, KD_mg=KD_mg_t, f_HT=f_HT_t,
                        tau_prec=tau_t, tau_rw=tau_rw_t, f_bio=fbio, result=res_t)

        if res_t['pCO2_ppm'] < T_pCO2 * 0.8:
            print("  pCO2 below 80% of target — stopping scan.")
            break

else:
    print(f"\n  Abiotic [Ca] within 5% of target — skipping Phase 2.")

# ---------------------------------------------------------------------------
# Final report
# ---------------------------------------------------------------------------

r  = best['result']
sf = seafloor_alk_flux_tmol(r)

print("\n\n" + "="*70)
print("  CALIBRATION COMPLETE — FINAL CONSTANTS")
print("="*70)
print_state("Best result", r, best['K_na'], best['KD_mg'], best['f_HT'], best['tau_prec'], best['tau_rw'], best['f_bio']) # type: ignore
print(f"\n  Seafloor Alk flux: {sf:.3f} Tmol/yr  (target ~1)")

print("""
  ── Paste into weathering.py ──────────────────────────────────────────""")
print(f"  ALPHA_REF               = {ALPHA_REF:.6f}")
print("""
  ── Paste into planet.py ──────────────────────────────────────────────""")
print(f"  K_CL_SUBDUCTION         = {K_CL_ANALYTIC:.6e}")
print(f"  K_NA_ALBITIZATION       = {best['K_na']:.6e}")
print(f"  KD_MG_HT                = {best['KD_mg']:.6e}")
print("""
  ── Planet constructor defaults ───────────────────────────────────────""")
# print(f"  f_HT     = {best['f_HT']:.4f}   # PHREEQC HT Ca source + Mg-Ca exchange")
print(f"  tau_prec = {best['tau_prec']/YR:.4e} * YR   # {best['tau_prec']/YR/1e6:.3f} Myr") # type: ignore
print(f"  tau_rw   = {best['tau_rw']/YR:.4e} * YR   # {best['tau_rw']/YR/1e6:.1f} Myr  (reverse weathering)") # type: ignore
print(f"  f_bio    = {best['f_bio']:.2f}")
print()
