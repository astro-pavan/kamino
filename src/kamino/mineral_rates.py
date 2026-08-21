"""
mineral_rates.py
================
Python implementations of mineral dissolution rate laws extracted from the
PHREEQC databases bundled with this project.

Sources
-------
Silicates / clays (acid–neutral–base TST framework):
    Kinec_v3.dat = basic_v2.dat
    Hermanska et al. (2022) Part I — primary silicate minerals
    Hermanska et al. (2023) Part II — secondary silicate minerals
    Oelkers & Addassi (2024*) Part III — non-silicate minerals

Pyrite:
    Williamson & Rimstidt (1994), GCA

General dissolution rate form (TST-based, Aagaard & Helgeson 1982)
-------------------------------------------------------------------
    r [mol m⁻² s⁻¹] = k_eff(T, pH) × (1 − Ω^(1/σ))

where
    k_eff = A_acid · exp(−E_acid / RT) · [H⁺]^n_acid
          + A_neut · exp(−E_neut / RT)
          + A_base · exp(−E_base / RT) · [H⁺]^n_base   (n_base < 0)

    Ω  = IAP / K_sp  (saturation ratio, dimensionless; Ω < 1 → undersaturated)
    σ  = Temkin stoichiometric coefficient (controls the shape of f(Ω))

Function naming convention
--------------------------
``mineral_k(T, pH)``       — effective rate constant k_eff [mol m⁻² s⁻¹];
                             temperature- and pH-dependent only; no Ω dependence.
``mineral_rate(T, pH, Ω)`` — full dissolution rate = k_eff × (1 − Ω^(1/σ))
                             [mol m⁻² s⁻¹].  Multiply by the total reactive
                             surface area S [m²] to get mol s⁻¹.

JAX acceleration
----------------
All scalar functions are JIT-compiled via @jax.jit.  After the first
call (which triggers tracing/compilation), subsequent calls skip all Python
overhead and run as compiled XLA kernels.

Vectorised variants (``vec_*``) accept JAX arrays instead of scalars — they
apply ``vmap`` over the leading batch axis so a sweep over N (T, pH, Ω)
conditions costs one compiled kernel launch:

    k_vals = vec_albite_k(T_arr, pH_arr)             # shape (N,)
    rates  = vec_albite_rate(T_arr, pH_arr, om_arr)  # shape (N,)

``vec_all_rates(T_arr, pH_arr, omega_arr)`` evaluates ALL TST minerals at once
and returns a (N, n_minerals) array.  Mineral order is given by
``TST_MINERAL_NAMES``.

Enable 64-bit floats (matching the original numpy behaviour) by calling::

    import jax; jax.config.update("jax_enable_x64", True)

This module sets that flag automatically at import time.

Gradient computation is free via ``jax.grad`` / ``jax.jacfwd`` — useful for
sensitivity analysis and inverse modelling.

Notes
-----
* K-Feldspar: the original Kinec_v3.dat RATES block contains a bug —
  the variable 'nC' is undefined in the base mechanism line.  The correct
  value is nb = −0.75.  This is fixed here.

* Forsterite / Fayalite: both terms in the database have positive n exponents
  (two "acid-type" terms); there is no true base mechanism for these minerals.

* Carbonate minerals (Calcite, Aragonite, Magnesite, Siderite, Dolomite) use a
  combined acid + carbonate-saturation rate law (Oelkers & Addassi 2024*).
  Their k functions accept optional carbonate activities:
      k = calcite_k(T, pH, act_HCO3=..., act_CO3=...)

* Pyrite uses the empirical log-rate form of Williamson & Rimstidt (1994);
  it does not follow the k × f(Ω) form and has no separate k function.

* Minerals with NO rate data in any searched database:
    Ferrosilite  — PHASES in basic_v2.dat but no kinetic data found
    Palygorskite — no PHASES and no independent data; use Sepiolite as proxy
    (Glauconite has kinetic data in phreeqc_rates.dat but no thermodynamic
    PHASES in any database; rate function provided as glauconite_rate() only)

* Additional rate sources (Kinec.v2.dat):
    Wollastonite, Illite, Sepiolite, SiO2(am)/Opal,
    Smectite-high-Fe-Mg, Smectite-low-Fe-Mg, Gypsum, Dolomite

* Goethite: neutral-only rate from Palandri & Kharaka (2004) Table 31.
  goethite_k(T) takes only T (no pH dependence).
"""

import numpy as np
import jax
import jax.numpy as jnp
from jax import jit, vmap

jax.config.update("jax_enable_x64", True)   # preserve float64 precision

R = 8.314  # J mol⁻¹ K⁻¹


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def _tst_factor(omega, sigma):
    """
    Thermodynamic driving force:  f(Ω) = 1 − Ω^(1/σ)

    Parameters
    ----------
    omega : saturation ratio (IAP / K_sp); dissolution when omega < 1
    sigma : Temkin coefficient
    """
    return 1.0 - omega ** (1.0 / sigma)


def _carbonate_k(T: float, pH: float,
                 Aa: float, Ea: float, na: float,
                 Ac: float, Eac: float, kc: float,
                 act_HCO3: float, act_CO3: float) -> float:
    """
    Generic k_eff for carbonate minerals (Oelkers & Addassi 2024*).

    k = A_acid · exp(−E_acid/RT) · aH^n  +  A_c · exp(−E_c/RT) · f_c

    where f_c = 1 − kc·(a_HCO3⁻ + a_CO3²⁻) / (1 + kc·(a_HCO3⁻ + a_CO3²⁻))
    """
    aH      = 10.0 ** (-pH)
    act_c   = act_HCO3 + act_CO3
    carb_term = 1.0 - kc * act_c / (1.0 + kc * act_c)
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ac * jnp.exp(-Eac / (R * T)) * carb_term)


# ══════════════════════════════════════════════════════════════════════════════
# PRIMARY SILICATES
# ══════════════════════════════════════════════════════════════════════════════

@jit
def nepheline_k(T: float, pH: float) -> float:
    """
    Nepheline (NaAlSiO4) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 1)
    Source: Hermanska et al. (2022), Kinec_v3.dat

    The feldspathoid counterpart of albite, produced by the CIPW desilication cascade in
    silica-undersaturated (high mantle Mg/Si) crusts. ~17x faster than albite at 343 K / pH 6.6,
    and — because its solubility product carries H4SiO4 to the first power rather than the
    third — it is not suppressed by the high dissolved silica that pins albite at saturation.

    Parameters
    ----------
    T   : temperature [K]
    pH  : solution pH  (−log₁₀ activity of H⁺)

    Returns
    -------
    k_eff : mol m⁻² s⁻¹
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 5.0e7,   63_000.0, 1.00
    An, En     = 1.0e-1,  58_500.0
    Ab, Eb, nb = 7.5e-5,  58_000.0, -0.40
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


def albite_k(T: float, pH: float) -> float:
    """
    Albite (NaAlSi3O8) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3)
    Source: Hermanska et al. (2022), Kinec_v3.dat

    Parameters
    ----------
    T   : temperature [K]
    pH  : solution pH  (−log₁₀ activity of H⁺)

    Returns
    -------
    k_eff : mol m⁻² s⁻¹
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 0.7,     58_000.0, 0.30
    An, En     = 2.05e-1, 60_000.0
    Ab, Eb, nb = 1.5e-5,  50_000.0, -0.30
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def albite_rate(T: float, pH: float, omega: float) -> float:
    """
    Albite (NaAlSi3O8) dissolution rate [mol m⁻² s⁻¹].

    rate = albite_k(T, pH) × (1 − Ω^(1/3))

    Parameters
    ----------
    T     : temperature [K]
    pH    : solution pH  (−log₁₀ activity of H⁺)
    omega : saturation ratio Ω = IAP / K_sp  (< 1 for dissolution)
    """
    return albite_k(T, pH) * _tst_factor(omega, 3.0)


@jit
def anorthite_k(T: float, pH: float) -> float:
    """
    Anorthite (CaAl₂Si₂O₈) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 2)
    Source: Hermanska et al. (2022), Kinec_v3.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 9.82e4,  58_000.0, 1.22
    An, En     = 1.5e-1,  60_000.0
    Ab, Eb, nb = 1.5e-5,  50_000.0, -0.35
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def anorthite_rate(T: float, pH: float, omega: float) -> float:
    """
    Anorthite (CaAl₂Si₂O₈) dissolution rate [mol m⁻² s⁻¹].

    rate = anorthite_k(T, pH) × (1 − Ω^(1/2))
    """
    return anorthite_k(T, pH) * _tst_factor(omega, 2.0)


@jit
def k_feldspar_k(T: float, pH: float) -> float:
    """
    K-Feldspar (KAlSi₃O₈) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3)
    Source: Hermanska et al. (2022), Kinec_v3.dat
    Bug fix: 'nC' in original RATES block → nb = −0.75
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 0.05,    51_700.0, 0.45
    An, En     = 1.08e-2, 60_000.0
    Ab, Eb, nb = 1.2e-10, 62_195.0, -0.75  # nb fixed (was 'nC' in source)
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def k_feldspar_rate(T: float, pH: float, omega: float) -> float:
    """
    K-Feldspar (KAlSi₃O₈) dissolution rate [mol m⁻² s⁻¹].

    rate = k_feldspar_k(T, pH) × (1 − Ω^(1/3))
    """
    return k_feldspar_k(T, pH) * _tst_factor(omega, 3.0)


@jit
def enstatite_k(T: float, pH: float) -> float:
    """
    Enstatite (MgSiO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral  (σ = 1)
    Source: Hermanska et al. (2022), Kinec_v3.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 0.574,  46_080.0, 0.50
    An, En     = 6252.0, 89_538.0
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T)))


@jit
def enstatite_rate(T: float, pH: float, omega: float) -> float:
    """
    Enstatite (MgSiO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = enstatite_k(T, pH) × (1 − Ω^(1/1))
    """
    return enstatite_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def diopside_k(T: float, pH: float) -> float:
    """
    Diopside (CaMgSi₂O₆) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral  (σ = 2)
    Source: Hermanska et al. (2022), Kinec_v3.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 8.55e-5, 32_654.0, 0.25
    An, En     = 4.30e-4, 43_866.0
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T)))


@jit
def diopside_rate(T: float, pH: float, omega: float) -> float:
    """
    Diopside (CaMgSi₂O₆) dissolution rate [mol m⁻² s⁻¹].

    rate = diopside_k(T, pH) × (1 − Ω^(1/2))
    """
    return diopside_k(T, pH) * _tst_factor(omega, 2.0)


@jit
def forsterite_k(T: float, pH: float) -> float:
    """
    Forsterite (Mg₂SiO₄) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: two acid-type terms (both n > 0; no neutral or base)  (σ = 1)
    Source: Hermanska et al. (2022), Kinec_v3.dat

    Note: The database uses labels 'acid' and 'base' but both exponents are
    positive, reflecting the H₂O-promoted mechanism at high pH rather than a
    true OH⁻ dependency.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 14.8e4, 70_400.0, 0.44   # acid mechanism
    Ab, Eb, nb = 220.0,  60_900.0, 0.22   # H₂O-promoted mechanism (nb > 0)
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def forsterite_rate(T: float, pH: float, omega: float) -> float:
    """
    Forsterite (Mg₂SiO₄) dissolution rate [mol m⁻² s⁻¹].

    rate = forsterite_k(T, pH) × (1 − Ω^(1/1))
    """
    return forsterite_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def fayalite_k(T: float, pH: float) -> float:
    """
    Fayalite (Fe₂SiO₄) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: two acid-type terms (both n > 0)  (σ = 1)
    Source: Hermanska et al. (2022), Kinec_v3.dat

    Note: same structure as Forsterite — see that docstring.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 1.20e6, 70_400.0, 0.44
    Ab, Eb, nb = 1.91e3, 60_900.0, 0.22
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def fayalite_rate(T: float, pH: float, omega: float) -> float:
    """
    Fayalite (Fe₂SiO₄) dissolution rate [mol m⁻² s⁻¹].

    rate = fayalite_k(T, pH) × (1 − Ω^(1/1))
    """
    return fayalite_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def quartz_k(T: float, pH: float) -> float:
    """
    Quartz (SiO₂) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + base (n_base < 0 → OH⁻ promoted)  (σ = 1)
    Source: Hermanska et al. (2022), Kinec_v3.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 4.03e-4, 45_600.0,  0.309
    Ab, Eb, nb = 0.105,   80_000.0, -0.41
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def quartz_rate(T: float, pH: float, omega: float) -> float:
    """
    Quartz (SiO₂) dissolution rate [mol m⁻² s⁻¹].

    rate = quartz_k(T, pH) × (1 − Ω^(1/1))
    """
    return quartz_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def wollastonite_k(T: float, pH: float) -> float:
    """
    Wollastonite (CaSiO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: two acid-type terms (both n > 0; no neutral or base)  (σ = 1)
    Source: Kinec.v2.dat

    Note: both terms have positive exponents, similar to Forsterite/Fayalite.
    The second term (Ab, nb=0.15) has a lower activation energy and reflects
    the water-promoted mechanism at near-neutral pH.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 700.0, 56_000.0, 0.40
    Ab, Eb, nb =  20.0, 52_000.0, 0.15   # nb > 0 → acid-type (not a base term)
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def wollastonite_rate(T: float, pH: float, omega: float) -> float:
    """
    Wollastonite (CaSiO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = wollastonite_k(T, pH) × (1 − Ω^(1/1))
    """
    return wollastonite_k(T, pH) * _tst_factor(omega, 1.0)


# ══════════════════════════════════════════════════════════════════════════════
# CLAY MINERALS
# ══════════════════════════════════════════════════════════════════════════════

@jit
def kaolinite_k(T: float, pH: float) -> float:
    """
    Kaolinite (Al₂Si₂O₅(OH)₄) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 2)
    Source: Hermanska et al. (2023), Kinec_v3.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 2.85,    73_000.0,  0.45
    An, En     = 4.15e-3, 67_000.0
    Ab, Eb, nb = 2.40e-11, 61_000.0, -0.76
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def kaolinite_rate(T: float, pH: float, omega: float) -> float:
    """
    Kaolinite (Al₂Si₂O₅(OH)₄) dissolution rate [mol m⁻² s⁻¹].

    rate = kaolinite_k(T, pH) × (1 − Ω^(1/2))
    """
    return kaolinite_k(T, pH) * _tst_factor(omega, 2.0)


@jit
def illite_k(T: float, pH: float) -> float:
    """
    Illite (K₀.₆Mg₀.₂₅Al₁.₈Al₀.₅Si₃.₅O₁₀(OH)₂) effective rate constant
    k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.5)
    Source: Kinec.v2.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 7.3e-4,  50_000.0,  0.55
    An, En     = 3.348e-3, 70_000.0
    Ab, Eb, nb = 6.0e-8,  74_000.0, -0.60
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def illite_rate(T: float, pH: float, omega: float) -> float:
    """
    Illite dissolution rate [mol m⁻² s⁻¹].

    rate = illite_k(T, pH) × (1 − Ω^(1/3.5))
    """
    return illite_k(T, pH) * _tst_factor(omega, 3.5)


@jit
def smectite_high_k(T: float, pH: float) -> float:
    """
    Smectite-high-Fe-Mg effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.5)
    Source: Kinec.v2.dat
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 1.66e-3, 50_798.0,  0.55
    An, En     = 9.0e-10, 30_000.0
    Ab, Eb, nb = 1.5e-9,  48_000.0, -0.30
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def smectite_high_rate(T: float, pH: float, omega: float) -> float:
    """
    Smectite-high-Fe-Mg dissolution rate [mol m⁻² s⁻¹].

    rate = smectite_high_k(T, pH) × (1 − Ω^(1/3.5))
    """
    return smectite_high_k(T, pH) * _tst_factor(omega, 3.5)


@jit
def smectite_low_k(T: float, pH: float) -> float:
    """
    Smectite-low-Fe-Mg effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.75)
    Source: Kinec.v2.dat

    Note: same kinetic parameters as Smectite-high-Fe-Mg; only σ differs.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 1.66e-3, 50_798.0,  0.55
    An, En     = 9.0e-10, 30_000.0
    Ab, Eb, nb = 1.5e-9,  48_000.0, -0.30
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def smectite_low_rate(T: float, pH: float, omega: float) -> float:
    """
    Smectite-low-Fe-Mg dissolution rate [mol m⁻² s⁻¹].

    rate = smectite_low_k(T, pH) × (1 − Ω^(1/3.75))
    """
    return smectite_low_k(T, pH) * _tst_factor(omega, 3.75)


@jit
def sepiolite_k(T: float, pH: float) -> float:
    """
    Sepiolite (Mg₄Si₆O₁₅(OH)₂·6H₂O) effective rate constant k_eff
    [mol m⁻² s⁻¹].

    Mechanism: acid + neutral  (σ = 6)
    Source: Kinec.v2.dat

    Note: phreeqc_rates.dat (2023 table) states that these parameters are
    also applicable to Palygorskite (attapulgite), which has no independent
    kinetic data.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 5.89e-3, 50_200.0, 0.248
    An, En     = 8.0e-7,  40_700.0
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T)))


@jit
def sepiolite_rate(T: float, pH: float, omega: float) -> float:
    """
    Sepiolite dissolution rate [mol m⁻² s⁻¹].

    rate = sepiolite_k(T, pH) × (1 − Ω^(1/6))
    """
    return sepiolite_k(T, pH) * _tst_factor(omega, 6.0)


@jit
def sio2am_k(T: float, pH: float) -> float:
    """
    Amorphous silica / SiO₂(am) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + base  (σ = 1)
    Source: Kinec.v2.dat

    Note: SiO₂(am) is used as a proxy for biogenic opal / amorphous silica.
    Parameters are very close to quartz (same na and nb), reflecting a shared
    silica dissolution mechanism but faster rates (lower activation energies).
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 4.563e-4, 41_610.0,  0.309
    Ab, Eb, nb = 0.0353,   73_000.0, -0.41
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def sio2am_rate(T: float, pH: float, omega: float) -> float:
    """
    Amorphous silica / SiO₂(am) dissolution rate [mol m⁻² s⁻¹].

    rate = sio2am_k(T, pH) × (1 − Ω^(1/1))
    """
    return sio2am_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def goethite_k(T: float) -> float:
    """
    Goethite (FeOOH) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: neutral only — no pH dependence  (σ = 1)
    Source: Palandri & Kharaka (2004) Table 31, via phreeqc_rates.dat.

    The pre-exponential A_neut is derived from the tabulated log₁₀(k) = −7.94
    at 25 °C and Ea = 86.5 kJ mol⁻¹:
        A_neut = 10^(−7.94) / exp(−86500 / (8.314 × 298.15)) ≈ 1.87 × 10⁷

    Parameters
    ----------
    T : temperature [K]
    """
    An, En = 1.87e7, 86_500.0
    return An * jnp.exp(-En / (R * T))


@jit
def goethite_rate(T: float, omega: float) -> float:
    """
    Goethite (FeOOH) dissolution rate [mol m⁻² s⁻¹].

    rate = goethite_k(T) × (1 − Ω^(1/1))

    Parameters
    ----------
    T     : temperature [K]
    omega : saturation ratio Ω = IAP / K_sp  (< 1 for dissolution)
    """
    return goethite_k(T) * _tst_factor(omega, 1.0)


@jit
def glauconite_k(T: float, pH: float) -> float:
    """
    Glauconite effective rate constant k_eff [mol m⁻² s⁻¹]
    (Python only — no PHREEQC PHASES entry).

    Mechanism: acid + neutral  (σ = 1, assumed)
    Source: phreeqc_rates.dat 2023 parameter table

    Note: Glauconite thermodynamic data is absent from all available PHREEQC
    databases (noted in Kinec_v3.dat README).  This function provides a rate
    estimate for use in pure-Python models; it cannot be used as a PHREEQC
    RATES block without a corresponding PHASES definition.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 9.55e-7, 32_300.0, 0.37
    An, En     = 1.10e-7, 37_500.0
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T)))


@jit
def glauconite_rate(T: float, pH: float, omega: float) -> float:
    """
    Glauconite dissolution rate [mol m⁻² s⁻¹]
    (Python only — no PHREEQC PHASES entry).

    rate = glauconite_k(T, pH) × (1 − Ω^(1/1))
    """
    return glauconite_k(T, pH) * _tst_factor(omega, 1.0)


# ══════════════════════════════════════════════════════════════════════════════
# PRECIPITATED / CARBONATE MINERALS
# ══════════════════════════════════════════════════════════════════════════════

@jit
def calcite_k(T: float, pH: float,
              act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Calcite (CaCO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + carbonate-saturation  (σ = 1)
    Source: Oelkers & Addassi (2024*), Kinec_v3.dat

    Parameters
    ----------
    T        : temperature [K]
    pH       : solution pH
    act_HCO3 : activity of HCO3⁻  (default 0 → carbonate term minimised)
    act_CO3  : activity of CO3²⁻  (default 0)
    """
    return _carbonate_k(T, pH,
                        Aa=5.625, Ea=16_000.0, na=1.0,
                        Ac=62.5,  Eac=48_000.0, kc=160.0,
                        act_HCO3=act_HCO3, act_CO3=act_CO3)


@jit
def calcite_rate(T: float, pH: float, omega: float,
                 act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Calcite (CaCO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = calcite_k(T, pH, act_HCO3, act_CO3) × (1 − Ω^(1/1))

    Parameters
    ----------
    T, pH, omega : see albite_rate docstring
    act_HCO3  : activity of HCO3⁻  (default 0)
    act_CO3   : activity of CO3²⁻  (default 0)
    """
    return calcite_k(T, pH, act_HCO3, act_CO3) * _tst_factor(omega, 1.0)


@jit
def aragonite_k(T: float, pH: float,
                act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Aragonite (CaCO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + carbonate-saturation  (σ = 1)
    Source: Oelkers & Addassi (2024*), Kinec_v3.dat
    """
    return _carbonate_k(T, pH,
                        Aa=11.025, Ea=16_000.0, na=1.0,
                        Ac=122.5,  Eac=48_000.0, kc=160.0,
                        act_HCO3=act_HCO3, act_CO3=act_CO3)


@jit
def aragonite_rate(T: float, pH: float, omega: float,
                   act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Aragonite (CaCO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = aragonite_k(T, pH, act_HCO3, act_CO3) × (1 − Ω^(1/1))
    """
    return aragonite_k(T, pH, act_HCO3, act_CO3) * _tst_factor(omega, 1.0)


@jit
def magnesite_k(T: float, pH: float,
                act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Magnesite (MgCO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + carbonate-saturation  (σ = 3.94)
    Source: Oelkers & Addassi (2024*), Kinec_v3.dat
    """
    return _carbonate_k(T, pH,
                        Aa=5.0e-4, Ea=16_000.0, na=0.66,
                        Ac=2.7e-2, Eac=45_000.0, kc=380.0,
                        act_HCO3=act_HCO3, act_CO3=act_CO3)


@jit
def magnesite_rate(T: float, pH: float, omega: float,
                   act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Magnesite (MgCO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = magnesite_k(T, pH, act_HCO3, act_CO3) × (1 − Ω^(1/3.94))
    """
    return magnesite_k(T, pH, act_HCO3, act_CO3) * _tst_factor(omega, 3.94)


@jit
def siderite_k(T: float, pH: float,
               act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Siderite (FeCO₃) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + carbonate-saturation  (σ = 4)
    Source: Oelkers & Addassi (2024*), Kinec_v3.dat

    Note: the PHASES block for Siderite in Kinec_v3.dat uses the reaction
    FeCO₃ = Fe²⁺ + CO₃²⁻ (from carbfix.dat), while the llnl.dat uses
    FeCO₃ + H⁺ = Fe²⁺ + HCO₃⁻.  The rate law is independent of this choice.
    """
    return _carbonate_k(T, pH,
                        Aa=2.0e-3, Ea=16_000.0, na=0.70,
                        Ac=0.2,    Eac=48_000.0, kc=160.0,
                        act_HCO3=act_HCO3, act_CO3=act_CO3)


@jit
def siderite_rate(T: float, pH: float, omega: float,
                  act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Siderite (FeCO₃) dissolution rate [mol m⁻² s⁻¹].

    rate = siderite_k(T, pH, act_HCO3, act_CO3) × (1 − Ω^(1/4))
    """
    return siderite_k(T, pH, act_HCO3, act_CO3) * _tst_factor(omega, 4.0)


@jit
def gypsum_k(T: float, pH: float) -> float:
    """
    Gypsum (CaSO₄·2H₂O) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid only  (σ = 1)
    Source: Kinec.v2.dat

    Note: the small n (0.11) means relatively weak pH dependence.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 1.8e4, 37_700.0, 0.11
    return Aa * jnp.exp(-Ea / (R * T)) * aH**na


@jit
def gypsum_rate(T: float, pH: float, omega: float) -> float:
    """
    Gypsum (CaSO₄·2H₂O) dissolution rate [mol m⁻² s⁻¹].

    rate = gypsum_k(T, pH) × (1 − Ω^(1/1))
    """
    return gypsum_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def dolomite_k(T: float, pH: float,
               act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Dolomite (CaMg(CO₃)₂) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + carbonate-saturation inhibition  (σ = 1.9)
    Source: Kinec.v2.dat
    """
    return _carbonate_k(T, pH,
                        Aa=1.2e-3, Ea=10_000.0, na=0.50,
                        Ac=650.0,  Eac=65_000.0, kc=160.0,
                        act_HCO3=act_HCO3, act_CO3=act_CO3)


@jit
def dolomite_rate(T: float, pH: float, omega: float,
                  act_HCO3: float = 0.0, act_CO3: float = 0.0) -> float:
    """
    Dolomite (CaMg(CO₃)₂) dissolution rate [mol m⁻² s⁻¹].

    rate = dolomite_k(T, pH, act_HCO3, act_CO3) × (1 − Ω^(1/1.9))
    """
    return dolomite_k(T, pH, act_HCO3, act_CO3) * _tst_factor(omega, 1.9)


# ══════════════════════════════════════════════════════════════════════════════
# CLAY MINERALS — CHLORITE GROUP
# ══════════════════════════════════════════════════════════════════════════════

@jit
def chlorite_k(T: float, pH: float) -> float:
    """
    Chlorite (Clinochlore / Daphnite) effective rate constant k_eff
    [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 1)
    Source: phreeqc_rates.dat, 2023 parameter table (Clinochlore row).
    The table note explicitly states these parameters are also valid for
    Chamosite and Daphnite (Fe-chlorite endmembers).

    Thermodynamic phases in ocean_chem.dat:
      Clinochlore-14A : Mg5Al2Si3O10(OH)8  (Mg-chlorite)
      Daphnite-14A    : Fe5AlAlSi3O10(OH)8 (Fe-chlorite)

    Notes
    -----
    sigma = 1 is used as a default (not tabulated in the source).
    The activation energies (15–30 kJ mol⁻¹) are lower than most silicates,
    reflecting relatively fast chlorite dissolution kinetics.
    """
    aH = 10.0 ** (-pH)
    Aa, Ea, na = 1.50e-4,  30_000.0,  0.74
    An, En     = 4.70e-11, 15_000.0
    Ab, Eb, nb = 2.00e-12, 15_000.0, -0.20
    return (Aa * jnp.exp(-Ea / (R * T)) * aH**na
          + An * jnp.exp(-En / (R * T))
          + Ab * jnp.exp(-Eb / (R * T)) * aH**nb)


@jit
def chlorite_rate(T: float, pH: float, omega: float) -> float:
    """
    Chlorite (Clinochlore / Daphnite) dissolution rate [mol m⁻² s⁻¹].

    rate = chlorite_k(T, pH) × (1 − Ω^(1/1))

    Parameters
    ----------
    T     : temperature [K]
    pH    : solution pH  (−log₁₀ activity of H⁺)
    omega : saturation ratio Ω = IAP / K_sp  (< 1 for dissolution)
    """
    return chlorite_k(T, pH) * _tst_factor(omega, 1.0)


# ── Reverse-weathering authigenic clays ──────────────────────────────────────

@jit
def greenalite_k(T: float, pH: float) -> float:
    """
    Greenalite (Fe₃Si₂O₅(OH)₄) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 1)
    Source: proxy — chlorite (Clinochlore/Daphnite) parameters from
    phreeqc_rates.dat 2023 table. No independent Greenalite kinetic data
    found in Kinec_v3.dat or phreeqc_rates.dat; chlorite is the closest
    structurally related Fe-phyllosilicate with tabulated rate data.
    Thermodynamic data: llnl.dat (LLNL database).
    """
    return chlorite_k(T, pH)


@jit
def greenalite_rate(T: float, pH: float, omega: float) -> float:
    """Greenalite dissolution/precipitation rate [mol m⁻² s⁻¹]."""
    return greenalite_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def minnesotaite_k(T: float, pH: float) -> float:
    """
    Minnesotaite (Fe₃Si₄O₁₀(OH)₂) effective rate constant k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.5)
    Source: proxy — Smectite-high-Fe-Mg parameters (Kinec_v3.dat / Kinec.v2.dat).
    No independent Minnesotaite kinetic data in any searched database.
    Smectite used as proxy for this Fe²⁺-talc; talc kinetics are broadly
    similar to smectite in the literature.
    Thermodynamic data: llnl.dat (LLNL database).
    """
    return smectite_high_k(T, pH)


@jit
def minnesotaite_rate(T: float, pH: float, omega: float) -> float:
    """Minnesotaite dissolution/precipitation rate [mol m⁻² s⁻¹]."""
    return minnesotaite_k(T, pH) * _tst_factor(omega, 3.5)


@jit
def chamosite_7a_k(T: float, pH: float) -> float:
    """
    Chamosite-7A (Fe₂Al₂SiO₅(OH)₄, Berthierine) effective rate constant
    k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 1)
    Source: proxy — chlorite parameters from phreeqc_rates.dat 2023 table.
    The table note states these parameters apply to Chamosite (Fe-chlorite);
    Chamosite-7A is the 7 Å polymorph (Berthierine), structurally related.
    Thermodynamic data: llnl.dat (LLNL database).
    """
    return chlorite_k(T, pH)


@jit
def chamosite_7a_rate(T: float, pH: float, omega: float) -> float:
    """Chamosite-7A (Berthierine) dissolution/precipitation rate [mol m⁻² s⁻¹]."""
    return chamosite_7a_k(T, pH) * _tst_factor(omega, 1.0)


@jit
def nontronite_k(T: float, pH: float) -> float:
    """
    Nontronite-Mg (Mg₀.₁₆₅Fe₂Al₀.₃₃Si₃.₆₇·nH₂O) effective rate constant
    k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.65)
    Source: Kinec_v3.dat (Catalano 2013; Hermanska et al. 2022).
    Identical kinetic parameters to Smectite-high-Fe-Mg; sigma differs (3.65
    vs 3.5) per Kinec_v3.dat Nontronite-Mg entry.
    Thermodynamic data: llnl.dat (LLNL database).
    """
    return smectite_high_k(T, pH)


@jit
def nontronite_rate(T: float, pH: float, omega: float) -> float:
    """Nontronite-Mg dissolution/precipitation rate [mol m⁻² s⁻¹]."""
    return nontronite_k(T, pH) * _tst_factor(omega, 3.65)


@jit
def saponite_k(T: float, pH: float) -> float:
    """
    Saponite-Mg (Mg₃.₁₆₅Al₀.₃₃Si₃.₆₇O₁₀(OH)₂) effective rate constant
    k_eff [mol m⁻² s⁻¹].

    Mechanism: acid + neutral + base  (σ = 3.65)
    Source: Kinec_v3.dat (Catalano 2013; Hermanska et al. 2022).
    Identical kinetic parameters to Smectite-high-Fe-Mg; sigma differs (3.65)
    per Kinec_v3.dat Saponite-Mg-Mg entry.
    Thermodynamic data: llnl.dat (LLNL database).
    """
    return smectite_high_k(T, pH)


@jit
def saponite_rate(T: float, pH: float, omega: float) -> float:
    """Saponite-Mg dissolution/precipitation rate [mol m⁻² s⁻¹]."""
    return saponite_k(T, pH) * _tst_factor(omega, 3.65)


# ══════════════════════════════════════════════════════════════════════════════
# SULFIDE MINERALS
# ══════════════════════════════════════════════════════════════════════════════

@jit
def pyrite_rate(T: float, m_O2: float, m_H: float,
                log10_area: float = 0.0,
                M: float = 1.0, M0: float = 1.0,
                parm2: float = 0.67,
                parm3: float = 0.5,
                parm4: float = -0.11) -> float:
    """
    Pyrite (FeS₂) oxidative dissolution rate.

    Empirical log-linear form of Williamson & Rimstidt (1994):
        log₁₀(r) = −8.19  +  parm3 · log₁₀(m_O2)  +  parm4 · log₁₀(m_H)

    where r is in mol m⁻² s⁻¹.

    This rate applies only under undersaturated conditions with dissolved
    oxygen present.  In PHREEQC, the surface area scaling is:
        log₁₀(A_total) = log10_area + log₁₀(M0) + parm2 · log₁₀(M/M0)

    Note: Pyrite uses an empirical log-linear form rather than the TST
    k × f(Ω) framework; there is no separate pyrite_k function.

    Parameters
    ----------
    T           : temperature [K]  (not used — rate is not thermally activated
                  in the W&R 1994 formulation, valid at ~25 °C)
    m_O2        : dissolved O₂ molality [mol kg⁻¹]
    m_H         : H⁺ molality / activity (dimensionless activity or mol kg⁻¹)
    log10_area  : log₁₀ of specific surface area [log₁₀(m² mol⁻¹)]
                  (PARM(1) in PHREEQC, default 0 → 1 m² mol⁻¹)
    M           : current moles of pyrite
    M0          : initial moles of pyrite
    parm2       : exponent for (M/M0) geometric correction (default 0.67)
    parm3       : exponent for O₂ (default 0.5, Williamson & Rimstidt 1994)
    parm4       : exponent for H⁺ (default −0.11, Williamson & Rimstidt 1994)

    Returns
    -------
    moles_per_second : mol s⁻¹  (total dissolution rate, not per m²)
                       Returns 0 if M ≤ 0 or m_O2 ≤ 0.
    """
    log_rate = -8.19 + parm3 * jnp.log10(m_O2) + parm4 * jnp.log10(m_H)
    log_area = log10_area + jnp.log10(M0) + parm2 * jnp.log10(M / M0)
    result = 10.0 ** (log_area + log_rate)
    valid = jnp.logical_and(M > 0, m_O2 > 0)
    return jnp.where(valid, result, 0.0)


# ══════════════════════════════════════════════════════════════════════════════
# VECTORISED (vmap) VARIANTS
# ══════════════════════════════════════════════════════════════════════════════
# vec_*_k  : accept JAX arrays of shape (N,); return k_eff array of shape (N,)
# vec_*_rate: accept JAX arrays of shape (N,); return rate array of shape (N,)
# All are JIT-compiled; the first call compiles once, subsequent calls run at
# full XLA speed.

# k_eff variants (no omega)
vec_albite_k        = jit(vmap(albite_k))
vec_anorthite_k     = jit(vmap(anorthite_k))
vec_k_feldspar_k    = jit(vmap(k_feldspar_k))
vec_enstatite_k     = jit(vmap(enstatite_k))
vec_diopside_k      = jit(vmap(diopside_k))
vec_forsterite_k    = jit(vmap(forsterite_k))
vec_fayalite_k      = jit(vmap(fayalite_k))
vec_quartz_k        = jit(vmap(quartz_k))
vec_wollastonite_k  = jit(vmap(wollastonite_k))
vec_kaolinite_k     = jit(vmap(kaolinite_k))
vec_illite_k        = jit(vmap(illite_k))
vec_smectite_high_k = jit(vmap(smectite_high_k))
vec_smectite_low_k  = jit(vmap(smectite_low_k))
vec_sepiolite_k     = jit(vmap(sepiolite_k))
vec_sio2am_k        = jit(vmap(sio2am_k))
vec_goethite_k      = jit(vmap(goethite_k))   # signature (T,) — no pH
vec_gypsum_k        = jit(vmap(gypsum_k))
vec_chlorite_k      = jit(vmap(chlorite_k))
vec_glauconite_k    = jit(vmap(glauconite_k))

# Carbonate k variants: vmap over (T, pH); carbonate activities kept scalar
vec_calcite_k   = jit(vmap(calcite_k,   in_axes=(0, 0, None, None)))
vec_aragonite_k = jit(vmap(aragonite_k, in_axes=(0, 0, None, None)))
vec_magnesite_k = jit(vmap(magnesite_k, in_axes=(0, 0, None, None)))
vec_siderite_k  = jit(vmap(siderite_k,  in_axes=(0, 0, None, None)))
vec_dolomite_k  = jit(vmap(dolomite_k,  in_axes=(0, 0, None, None)))

# Full rate variants (with omega)
vec_albite_rate        = jit(vmap(albite_rate))
vec_anorthite_rate     = jit(vmap(anorthite_rate))
vec_k_feldspar_rate    = jit(vmap(k_feldspar_rate))
vec_enstatite_rate     = jit(vmap(enstatite_rate))
vec_diopside_rate      = jit(vmap(diopside_rate))
vec_forsterite_rate    = jit(vmap(forsterite_rate))
vec_fayalite_rate      = jit(vmap(fayalite_rate))
vec_quartz_rate        = jit(vmap(quartz_rate))
vec_wollastonite_rate  = jit(vmap(wollastonite_rate))
vec_kaolinite_rate     = jit(vmap(kaolinite_rate))
vec_illite_rate        = jit(vmap(illite_rate))
vec_smectite_high_rate = jit(vmap(smectite_high_rate))
vec_smectite_low_rate  = jit(vmap(smectite_low_rate))
vec_sepiolite_rate     = jit(vmap(sepiolite_rate))
vec_sio2am_rate        = jit(vmap(sio2am_rate))
vec_goethite_rate      = jit(vmap(goethite_rate))   # signature (T, omega)
vec_gypsum_rate        = jit(vmap(gypsum_rate))
vec_chlorite_rate      = jit(vmap(chlorite_rate))
vec_glauconite_rate    = jit(vmap(glauconite_rate))

# Carbonate rate variants: vmap over (T, pH, omega); carbonate activities scalar
vec_calcite_rate   = jit(vmap(calcite_rate,   in_axes=(0, 0, 0, None, None)))
vec_aragonite_rate = jit(vmap(aragonite_rate, in_axes=(0, 0, 0, None, None)))
vec_magnesite_rate = jit(vmap(magnesite_rate, in_axes=(0, 0, 0, None, None)))
vec_siderite_rate  = jit(vmap(siderite_rate,  in_axes=(0, 0, 0, None, None)))
vec_dolomite_rate  = jit(vmap(dolomite_rate,  in_axes=(0, 0, 0, None, None)))


# ══════════════════════════════════════════════════════════════════════════════
# BATCH EVALUATION — ALL TST MINERALS AT ONCE
# ══════════════════════════════════════════════════════════════════════════════

#: Ordered list of TST minerals evaluated by vec_all_rates / _stack_all_tst.
#: Index i in the output array corresponds to TST_MINERAL_NAMES[i].
TST_MINERAL_NAMES: list[str] = [
    'Albite', 'Anorthite', 'K-Feldspar', 'Enstatite', 'Diopside',
    'Forsterite', 'Fayalite', 'Quartz', 'Wollastonite',
    'Calcite', 'Aragonite', 'Magnesite', 'Siderite', 'Dolomite', 'Gypsum',
    'Kaolinite', 'Illite', 'Smectite-high-Fe-Mg', 'Smectite-low-Fe-Mg',
    'Sepiolite', 'Palygorskite', 'SiO2(am)', 'Goethite',
    'Clinochlore-14A', 'Daphnite-14A',
]


@jit
def _stack_all_tst(T: float, pH: float, omega: float) -> jnp.ndarray:
    """
    Evaluate all TST minerals at a single (T, pH, Ω) point.

    Returns a 1-D JAX array of length ``len(TST_MINERAL_NAMES)``.
    Intended to be vmapped; call ``vec_all_rates`` for array inputs.
    Carbonate activities default to zero (acid mechanism dominates).
    """
    return jnp.array([
        albite_rate(T, pH, omega),
        anorthite_rate(T, pH, omega),
        k_feldspar_rate(T, pH, omega),
        enstatite_rate(T, pH, omega),
        diopside_rate(T, pH, omega),
        forsterite_rate(T, pH, omega),
        fayalite_rate(T, pH, omega),
        quartz_rate(T, pH, omega),
        wollastonite_rate(T, pH, omega),
        calcite_rate(T, pH, omega),
        aragonite_rate(T, pH, omega),
        magnesite_rate(T, pH, omega),
        siderite_rate(T, pH, omega),
        dolomite_rate(T, pH, omega),
        gypsum_rate(T, pH, omega),
        kaolinite_rate(T, pH, omega),
        illite_rate(T, pH, omega),
        smectite_high_rate(T, pH, omega),
        smectite_low_rate(T, pH, omega),
        sepiolite_rate(T, pH, omega),
        sepiolite_rate(T, pH, omega),      # Palygorskite proxy
        sio2am_rate(T, pH, omega),
        goethite_rate(T, omega),           # neutral only — no pH dependence
        chlorite_rate(T, pH, omega),       # Clinochlore-14A
        chlorite_rate(T, pH, omega),       # Daphnite-14A (same parameters)
    ])


vec_all_rates = jit(vmap(_stack_all_tst))
"""
Evaluate all TST minerals over arrays of conditions in one compiled kernel.

Parameters
----------
T_arr, pH_arr, omega_arr : jnp.ndarray of shape (N,)

Returns
-------
rates : jnp.ndarray of shape (N, len(TST_MINERAL_NAMES))
    rates[i, j] is the rate of mineral TST_MINERAL_NAMES[j] at condition i.

Example
-------
>>> import jax.numpy as jnp
>>> T   = jnp.linspace(273, 373, 100)
>>> pH  = jnp.full(100, 7.0)
>>> om  = jnp.full(100, 0.01)
>>> rates = vec_all_rates(T, pH, om)          # shape (100, 25)
>>> albite_rates = rates[:, TST_MINERAL_NAMES.index('Albite')]
"""


# ══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE DICTIONARIES
# ══════════════════════════════════════════════════════════════════════════════

from collections.abc import Callable

#: Map mineral name → JIT-compiled k_eff function (T, pH) → k [mol m⁻² s⁻¹].
#: Carbonate minerals accept optional (act_HCO3, act_CO3) keyword arguments.
#: Goethite k function has signature (T,) with no pH argument.
K_FUNCTIONS: dict[str, Callable] = {
    # Primary silicates
    "Albite":        albite_k,
    "Nepheline":     nepheline_k,
    "Anorthite":     anorthite_k,
    "K-Feldspar":    k_feldspar_k,
    "Enstatite":     enstatite_k,
    "Diopside":      diopside_k,
    "Forsterite":    forsterite_k,
    "Fayalite":      fayalite_k,
    "Quartz":        quartz_k,
    "Wollastonite":  wollastonite_k,
    # Carbonates / sulfates
    "Calcite":       calcite_k,
    "Aragonite":     aragonite_k,
    "Magnesite":     magnesite_k,
    "Siderite":      siderite_k,
    "Dolomite":      dolomite_k,
    "Gypsum":        gypsum_k,
    # Clay minerals / oxides
    "Kaolinite":             kaolinite_k,
    "Illite":                illite_k,
    "Smectite-high-Fe-Mg":  smectite_high_k,
    "Smectite-low-Fe-Mg":   smectite_low_k,
    "Sepiolite":             sepiolite_k,
    "Palygorskite":          sepiolite_k,    # proxy: same parameters as Sepiolite
    "SiO2(am)":              sio2am_k,
    "Goethite":              goethite_k,     # signature (T,) — no pH
    "Clinochlore-14A":       chlorite_k,      # Mg-chlorite
    "Daphnite-14A":          chlorite_k,      # Fe-chlorite (same parameters)
    # Reverse-weathering authigenic clays
    "Greenalite":            greenalite_k,    # Fe-serpentine; proxy: chlorite
    "Minnesotaite":          minnesotaite_k,  # Fe-talc; proxy: smectite-high
    "Chamosite-7A":          chamosite_7a_k,  # Berthierine/Fe-serpentine; proxy: chlorite
    "Nontronite-Mg":         nontronite_k,    # Fe³⁺-smectite; Kinec_v3 params
    "Saponite-Mg":           saponite_k,      # Mg-smectite; Kinec_v3 params
    # Pyrite: empirical log-linear form — no separate k function
    # Glauconite: no PHREEQC PHASES; call glauconite_k() directly
}

#: Map mineral name → JIT-compiled rate function (T, pH, omega) → rate [mol m⁻² s⁻¹].
#: Carbonate minerals accept optional (act_HCO3, act_CO3) keyword arguments.
#: Goethite rate function has signature (T, omega) with no pH argument.
RATE_FUNCTIONS: dict[str, Callable] = {
    # Primary silicates
    "Albite":        albite_rate,
    "Anorthite":     anorthite_rate,
    "K-Feldspar":    k_feldspar_rate,
    "Enstatite":     enstatite_rate,
    "Diopside":      diopside_rate,
    "Forsterite":    forsterite_rate,
    "Fayalite":      fayalite_rate,
    "Quartz":        quartz_rate,
    "Wollastonite":  wollastonite_rate,
    # Carbonates / sulfates
    "Calcite":       calcite_rate,
    "Aragonite":     aragonite_rate,
    "Magnesite":     magnesite_rate,
    "Siderite":      siderite_rate,
    "Dolomite":      dolomite_rate,
    "Gypsum":        gypsum_rate,
    # Clay minerals / oxides
    "Kaolinite":             kaolinite_rate,
    "Illite":                illite_rate,
    "Smectite-high-Fe-Mg":  smectite_high_rate,
    "Smectite-low-Fe-Mg":   smectite_low_rate,
    "Sepiolite":             sepiolite_rate,
    "Palygorskite":          sepiolite_rate,  # proxy: same parameters as Sepiolite
    "SiO2(am)":              sio2am_rate,
    "Goethite":              jit(lambda T, pH, omega: goethite_rate(T, omega)),
    "Clinochlore-14A":       chlorite_rate,     # Mg-chlorite
    "Daphnite-14A":          chlorite_rate,     # Fe-chlorite (same parameters)
    # Reverse-weathering authigenic clays
    "Greenalite":            greenalite_rate,   # Fe-serpentine; proxy: chlorite
    "Minnesotaite":          minnesotaite_rate, # Fe-talc; proxy: smectite-high
    "Chamosite-7A":          chamosite_7a_rate, # Berthierine/Fe-serpentine; proxy: chlorite
    "Nontronite-Mg":         nontronite_rate,   # Fe³⁺-smectite; Kinec_v3 params
    "Saponite-Mg":           saponite_rate,     # Mg-smectite; Kinec_v3 params
    # Pyrite has a different signature; call pyrite_rate() directly
    # Glauconite: no PHREEQC PHASES; call glauconite_rate() directly
}

#: Minerals in the model with NO rate law data in the searched databases
NO_RATE_DATA = [
    "Ferrosilite",   # PHASES in basic_v2.dat but no kinetic data in any db
]
#: Minerals with Python-only rate functions (no PHREEQC PHASES / RATES block)
PYTHON_ONLY_RATE = [
    "Glauconite",    # kinetic data available (phreeqc_rates.dat 2023 table)
                     # but thermodynamic PHASES absent from all databases
]


# ══════════════════════════════════════════════════════════════════════════════
# Kinetic parameters table (for reference / bulk calculations)
# ══════════════════════════════════════════════════════════════════════════════

#: Kinetic parameters extracted directly from Kinec_v3.dat.
#: Units: A in mol m⁻² s⁻¹, E in J mol⁻¹, n dimensionless, sigma dimensionless.
#: For carbonate minerals, 'kc' is the carbonate-inhibition constant [L mol⁻¹].
KINETIC_PARAMS = {
    "Albite": {
        "sigma": 3.0,
        "acid":  {"A": 0.700,   "E": 58_000, "n":  0.30},
        "neut":  {"A": 2.05e-1, "E": 60_000},
        "base":  {"A": 1.5e-5,  "E": 50_000, "n": -0.30},
    },
    "Anorthite": {
        "sigma": 2.0,
        "acid":  {"A": 9.82e4,  "E": 58_000, "n":  1.22},
        "neut":  {"A": 1.5e-1,  "E": 60_000},
        "base":  {"A": 1.5e-5,  "E": 50_000, "n": -0.35},
    },
    "K-Feldspar": {
        "sigma": 3.0,
        "note":  "nC bug in Kinec_v3.dat fixed to nb = -0.75",
        "acid":  {"A": 0.050,   "E": 51_700, "n":  0.45},
        "neut":  {"A": 1.08e-2, "E": 60_000},
        "base":  {"A": 1.2e-10, "E": 62_195, "n": -0.75},
    },
    "Enstatite": {
        "sigma": 1.0,
        "acid":  {"A": 0.574,  "E": 46_080, "n": 0.50},
        "neut":  {"A": 6252.0, "E": 89_538},
    },
    "Diopside": {
        "sigma": 2.0,
        "acid":  {"A": 8.55e-5, "E": 32_654, "n": 0.25},
        "neut":  {"A": 4.30e-4, "E": 43_866},
    },
    "Forsterite": {
        "sigma": 1.0,
        "note":  "two acid-type terms (no neutral or base mechanism)",
        "acid":  {"A": 14.8e4, "E": 70_400, "n": 0.44},
        "acid2": {"A": 220.0,  "E": 60_900, "n": 0.22},
    },
    "Fayalite": {
        "sigma": 1.0,
        "note":  "two acid-type terms (no neutral or base mechanism)",
        "acid":  {"A": 1.20e6, "E": 70_400, "n": 0.44},
        "acid2": {"A": 1.91e3, "E": 60_900, "n": 0.22},
    },
    "Quartz": {
        "sigma": 1.0,
        "acid":  {"A": 4.03e-4, "E": 45_600, "n":  0.309},
        "base":  {"A": 0.105,   "E": 80_000, "n": -0.41},
    },
    "Calcite": {
        "sigma": 1.0,
        "acid":  {"A": 5.625,  "E": 16_000, "n": 1.0},
        "carb":  {"A": 62.5,   "E": 48_000, "kc": 160.0},
    },
    "Aragonite": {
        "sigma": 1.0,
        "acid":  {"A": 11.025, "E": 16_000, "n": 1.0},
        "carb":  {"A": 122.5,  "E": 48_000, "kc": 160.0},
    },
    "Magnesite": {
        "sigma": 3.94,
        "acid":  {"A": 5.0e-4, "E": 16_000, "n": 0.66},
        "carb":  {"A": 2.7e-2, "E": 45_000, "kc": 380.0},
    },
    "Siderite": {
        "sigma": 4.0,
        "acid":  {"A": 2.0e-3, "E": 16_000, "n": 0.70},
        "carb":  {"A": 0.2,    "E": 48_000, "kc": 160.0},
    },
    "Kaolinite": {
        "sigma": 2.0,
        "acid":  {"A": 2.85,    "E": 73_000, "n":  0.45},
        "neut":  {"A": 4.15e-3, "E": 67_000},
        "base":  {"A": 2.40e-11,"E": 61_000, "n": -0.76},
    },
    "Wollastonite": {
        "sigma": 1.0,
        "note":  "two acid-type terms (both nb > 0; no neutral or base mechanism)",
        "acid":  {"A": 700.0, "E": 56_000, "n": 0.40},
        "acid2": {"A":  20.0, "E": 52_000, "n": 0.15},
    },
    "Illite": {
        "sigma": 3.5,
        "acid":  {"A": 7.3e-4,   "E": 50_000, "n":  0.55},
        "neut":  {"A": 3.348e-3, "E": 70_000},
        "base":  {"A": 6.0e-8,   "E": 74_000, "n": -0.60},
    },
    "Smectite-high-Fe-Mg": {
        "sigma": 3.5,
        "acid":  {"A": 1.66e-3, "E": 50_798, "n":  0.55},
        "neut":  {"A": 9.0e-10, "E": 30_000},
        "base":  {"A": 1.5e-9,  "E": 48_000, "n": -0.30},
    },
    "Smectite-low-Fe-Mg": {
        "sigma": 3.75,
        "note":  "same kinetic parameters as Smectite-high-Fe-Mg; only sigma differs",
        "acid":  {"A": 1.66e-3, "E": 50_798, "n":  0.55},
        "neut":  {"A": 9.0e-10, "E": 30_000},
        "base":  {"A": 1.5e-9,  "E": 48_000, "n": -0.30},
    },
    "Sepiolite": {
        "sigma": 6.0,
        "note":  "also used as proxy for Palygorskite (phreeqc_rates.dat 2023 table)",
        "acid":  {"A": 5.89e-3, "E": 50_200, "n": 0.248},
        "neut":  {"A": 8.0e-7,  "E": 40_700},
    },
    "SiO2(am)": {
        "sigma": 1.0,
        "note":  "Opal/biogenic silica proxy; same na/nb as Quartz but faster rates",
        "acid":  {"A": 4.563e-4, "E": 41_610, "n":  0.309},
        "base":  {"A": 0.0353,   "E": 73_000, "n": -0.41},
    },
    "Goethite": {
        "sigma": 1.0,
        "note":  "neutral only; Palandri & Kharaka (2004) Table 31",
        "neut":  {"A": 1.87e7, "E": 86_500},
    },
    "Gypsum": {
        "sigma": 1.0,
        "note":  "acid only; weak pH dependence (n=0.11)",
        "acid":  {"A": 1.8e4, "E": 37_700, "n": 0.11},
    },
    "Dolomite": {
        "sigma": 1.9,
        "acid":  {"A": 1.2e-3, "E": 10_000, "n": 0.50},
        "carb":  {"A": 650.0,  "E": 65_000, "kc": 160.0},
    },
    # Chlorite group — Clinochlore-14A and Daphnite-14A share these parameters
    # (phreeqc_rates.dat 2023 table; Clinochlore row; note: also valid for Chamosite, Daphnite)
    "Clinochlore-14A": {
        "sigma": 1.0,
        "note":  "Ea in J/mol (30 and 15 kJ/mol); sigma=1 assumed (not tabulated)",
        "acid":  {"A": 1.50e-4,  "E": 30_000, "n":  0.74},
        "neut":  {"A": 4.70e-11, "E": 15_000},
        "base":  {"A": 2.00e-12, "E": 15_000, "n": -0.20},
    },
    "Daphnite-14A": {
        "sigma": 1.0,
        "note":  "same parameters as Clinochlore-14A (phreeqc_rates.dat 2023 table)",
        "acid":  {"A": 1.50e-4,  "E": 30_000, "n":  0.74},
        "neut":  {"A": 4.70e-11, "E": 15_000},
        "base":  {"A": 2.00e-12, "E": 15_000, "n": -0.20},
    },
    # Python-only (no PHREEQC PHASES)
    "Glauconite": {
        "sigma": 1.0,
        "note":  "phreeqc_rates.dat 2023 table; no thermodynamic PHASES in any database",
        "acid":  {"A": 9.55e-7, "E": 32_300, "n": 0.37},
        "neut":  {"A": 1.10e-7, "E": 37_500},
    },
}


# ══════════════════════════════════════════════════════════════════════════════
# Example / demo
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import time

    T     = 298.15   # 25 °C
    pH    = 7.0
    omega = 0.01     # strongly undersaturated (dissolution)

    print("k_eff values at T=25 °C, pH=7  (no TST factor)")
    print("(mol m⁻² s⁻¹)\n")
    print(f"{'Mineral':<24}  {'k_eff (mol m⁻² s⁻¹)':>20}")
    print("-" * 48)

    tst_k = [
        ("Albite",            albite_k),
        ("Anorthite",         anorthite_k),
        ("K-Feldspar",        k_feldspar_k),
        ("Enstatite",         enstatite_k),
        ("Diopside",          diopside_k),
        ("Forsterite",        forsterite_k),
        ("Fayalite",          fayalite_k),
        ("Quartz",            quartz_k),
        ("Wollastonite",      wollastonite_k),
        ("Kaolinite",         kaolinite_k),
        ("Illite",            illite_k),
        ("Smectite-high-Fe-Mg", smectite_high_k),
        ("Smectite-low-Fe-Mg",  smectite_low_k),
        ("Sepiolite",         sepiolite_k),
        ("SiO2(am)",          sio2am_k),
        ("Gypsum",            gypsum_k),
        ("Clinochlore-14A",   chlorite_k),
        ("Daphnite-14A",      chlorite_k),
        ("Glauconite",        glauconite_k),
    ]
    for name, fn in tst_k:
        k = fn(T, pH)
        print(f"  {name:<22}  {k:>20.3e}")

    print(f"\n  {'Goethite':<22}  {goethite_k(T):>20.3e}  (T only, no pH)")

    act_HCO3, act_CO3 = 1e-3, 1e-6
    print()
    carb_k = [
        ("Calcite",   calcite_k),
        ("Aragonite", aragonite_k),
        ("Magnesite", magnesite_k),
        ("Siderite",  siderite_k),
        ("Dolomite",  dolomite_k),
    ]
    for name, fn in carb_k:
        k = fn(T, pH, act_HCO3=act_HCO3, act_CO3=act_CO3)
        print(f"  {name:<22}  {k:>20.3e}")

    print("\n\nDissolution rates at T=25 °C, pH=7, Ω=0.01")
    print("(mol m⁻² s⁻¹ — multiply by surface area S [m²] for mol s⁻¹)\n")
    print(f"{'Mineral':<24}  {'rate (mol m⁻² s⁻¹)':>18}")
    print("-" * 46)

    tst_minerals = [
        ("Albite",            albite_rate),
        ("Anorthite",         anorthite_rate),
        ("K-Feldspar",        k_feldspar_rate),
        ("Enstatite",         enstatite_rate),
        ("Diopside",          diopside_rate),
        ("Forsterite",        forsterite_rate),
        ("Fayalite",          fayalite_rate),
        ("Quartz",            quartz_rate),
        ("Wollastonite",      wollastonite_rate),
        ("Kaolinite",         kaolinite_rate),
        ("Illite",            illite_rate),
        ("Smectite-high-Fe-Mg", smectite_high_rate),
        ("Smectite-low-Fe-Mg",  smectite_low_rate),
        ("Sepiolite",         sepiolite_rate),
        ("SiO2(am)",          sio2am_rate),
        ("Goethite",          lambda T, pH, omega: goethite_rate(T, omega)),
        ("Gypsum",            gypsum_rate),
        ("Clinochlore-14A",   chlorite_rate),
        ("Daphnite-14A",      chlorite_rate),
        ("Glauconite",        glauconite_rate),
    ]
    for name, fn in tst_minerals:
        r = fn(T, pH, omega)
        print(f"  {name:<22}  {r:>18.3e}")

    print()
    carb_minerals = [
        ("Calcite",   calcite_rate),
        ("Aragonite", aragonite_rate),
        ("Magnesite", magnesite_rate),
        ("Siderite",  siderite_rate),
        ("Dolomite",  dolomite_rate),
    ]
    for name, fn in carb_minerals:
        r = fn(T, pH, omega, act_HCO3=act_HCO3, act_CO3=act_CO3)
        print(f"  {name:<22}  {r:>18.3e}")

    # Pyrite
    m_O2 = 2.5e-4  # ~8 mg/L at 25°C
    m_H  = 10**(-pH)
    r_py = pyrite_rate(T, m_O2=m_O2, m_H=m_H)
    print(f"\n  {'Pyrite':<22}  {r_py:>18.3e}  (mol s⁻¹, includes surface area)")

    print(f"\nMinerals with NO rate data: {', '.join(NO_RATE_DATA)}")
    print(f"Python-only rate functions:  {', '.join(PYTHON_ONLY_RATE)}")

    # ── Speed comparison: scalar (JIT warm-up) vs batch (vmap) ──────────────
    print("\n── Speed demo ──────────────────────────────────────────────────────")
    N = 10_000
    T_arr   = jnp.linspace(273.15, 373.15, N)
    pH_arr  = jnp.full(N, 7.0)
    om_arr  = jnp.full(N, 0.01)

    # First call: JIT compilation + execution
    t0 = time.perf_counter()
    rates = vec_all_rates(T_arr, pH_arr, om_arr).block_until_ready()
    t1 = time.perf_counter()
    print(f"vec_all_rates (compile+run):  {N} cond × {len(TST_MINERAL_NAMES)} minerals "
          f"in {(t1-t0)*1e3:.1f} ms")

    # Second call: execution only (compiled kernel)
    t0 = time.perf_counter()
    rates = vec_all_rates(T_arr, pH_arr, om_arr).block_until_ready()
    t1 = time.perf_counter()
    print(f"vec_all_rates (cached):       {N} cond × {len(TST_MINERAL_NAMES)} minerals "
          f"in {(t1-t0)*1e3:.2f} ms")
    print(f"  albite @ T=298 K: {rates[N//4, 0]:.3e} mol m⁻² s⁻¹")

    # k_eff demo
    t0 = time.perf_counter()
    k_vals = vec_albite_k(T_arr, pH_arr).block_until_ready()
    t1 = time.perf_counter()
    print(f"vec_albite_k (cached):        {N} conditions in {(t1-t0)*1e3:.2f} ms")
    print(f"  albite k_eff @ T=298 K: {k_vals[N//4]:.3e} mol m⁻² s⁻¹")
