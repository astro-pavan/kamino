# Crust composition — methods reference

How kamino turns two astronomically-motivated numbers into the mineral assemblage the weathering
model dissolves. Written to be lifted into a methods section: every constant below is traceable to
a line of code and a citation.

**Inputs:** mantle molar Mg/Si, and core-formation oxygen fugacity ΔIW.
**Output:** a `crust_composition` dict of mineral weight fractions.

---

## 0. What you actually need to read

The pipeline is ~1000 lines across two files, but the *method* is eight functions:

| function | file | lines | role |
|---|---|---|---|
| `feo_from_delta_iw` | both | 5 / 19 | ΔIW → mantle FeO |
| `mantle_composition` | `.jl` | 12 | build the bulk mantle on the two axes |
| `T_at_entropy` | `.jl` | 13 | track the isentrope (MAGEMin has no isentropic mode) |
| `isentropic_melt` | `.jl` | 12 | decompression melting path |
| `melt_oxides` | `.jl` | 11 | extract the liquid composition |
| `Tp_for_F` | `.jl` | 11 | solve T_p for the target melt fraction |
| `oxide_composition` | `.py` | 18 | interpolate the resulting table |
| `cipw_norm` | `.py` | ~40 of 179 | melt oxides → normative minerals |

The remaining ~850 lines are CLI/sharding/validation tooling (`check_crust_table.py`,
`merge_crust_slices.py`, `--probe`, `--calibrate`, `--points`) and rationale comments recording
approaches that were tried and rejected. Those exist so the dead ends are not re-derived; none of
them is method.

---

## 1. The two axes

### 1.1 Mantle Mg/Si

Molar Mg/Si of the silicate mantle; Earth = 1.25 (McDonough & Sun 1995 pyrolite as constructed
here). Grid spans **0.5–2.0**, the stellar range. This is the dominant control on the
olivine/orthopyroxene ratio (Guimond et al. 2024, §3.1.1), and mantles become olivine-free below
Mg/Si ≈ 0.8 and orthopyroxene-free above ≈ 1.6.

### 1.2 Core-formation oxygen fugacity ΔIW

Oxygen fugacity relative to the iron–wüstite buffer, in log₁₀ units, at the conditions of
metal–silicate equilibrium. It sets how much iron the mantle retained as FeO rather than losing to
the core as metal, through Fe + ½O₂ ⇌ FeO:

```
ΔIW = 2 log₁₀( a_FeO(silicate) / a_Fe(metal) )
```

with ideal mixing, i.e. activities replaced by mole fractions — the same convention as Young et al.
(2023), who write it identically. Inverting for the mantle FeO content, on a cation-mole basis
against the pyrolite non-Fe budget:

```
X_FeO = C · 10^(ΔIW/2)                                  [feo_from_delta_iw]
FeO(wt%) = 100 · X_FeO · k / [ (1 − X_FeO)/M_FeO + X_FeO · k ]
```

where `k` = 0.020237 cation mol per wt% of the renormalised non-Fe budget, and **C = 0.567957 is
calibrated so that ΔIW = −2 reproduces BSE FeO = 8.05 wt% exactly.**

The axis is **logarithmic in FeO** — this is the single most important thing to state in a methods
section, because it is not intuitive:

| ΔIW | −5 | −4 | −3 | **−2** | −1 |
|---|---|---|---|---|---|
| mantle FeO (wt%) | 0.26 | 0.82 | 2.59 | **8.05** | 24.1 |

The grid stops at −1 because mantle FeO then reaches ~24 wt%, at the ≈25 wt% ceiling above which
Guimond et al. (2024, §3.1.1) state the thermodynamic models are unreliable.

**Two caveats to declare.** (i) The activity constant carries ~0.4 log units of definitional slop
(γ_FeO = 1.0 gives Earth ΔIW = −2.35; γ_FeO = 1.5 gives −1.99), so only the Earth anchor is
meaningful — the absolute value is a label. (ii) Calibrating on Earth alone and extrapolating three
log units down mislabels the reduced end: Young et al. put FeO = 0.07 wt% at ΔIW ≈ −5, whereas this
mapping puts it at −6.1, a ~1 log unit offset. The *melt* is insensitive there (SiO₂ moves 0.08 wt%
between ΔIW −5 and −4.5), so this is a labelling issue, not a physics one — but do not quote the
reduced-end ΔIW labels against Young's without a two-point recalibration.

### 1.3 What is deliberately not an axis

- **T_p** — not observable, and not free (a mantle that cannot melt cannot transport heat
  magmatically). It is *solved* per composition, not prescribed.
- **C/O** — see §1.4; it is rejected on three independent grounds, and its real effect is
  atmospheric rather than crustal.
- **Ca/Al** — real, but Ca and Al are both refractory so stellar Ca/Al varies least of the candidate
  ratios. Held at pyrolite.
- **Na₂O** — the strongest un-swept lever, and should be reported as such: alkalis move the solidus
  more than MgO/FeO does (Guimond et al. 2024, §4.1), and albite-vs-nepheline is what sets ocean Na
  in this model. Excluded because Na is moderately volatile and devolatilisation is stochastic, so
  it does not belong on a grid indexed by observables. Held at pyrolite 0.36 wt%.

### 1.4 Why C/O is not an axis

The obvious third parameter, and the one most often suggested. It is excluded for three separate
reasons, any one of which would be sufficient.

**(i) The stellar route is genuinely rare.** Carbide and nitride condensation requires C/O ≳ 0.9
(Larimer & Bartholomay 1979). That is ~1% of FGK dwarfs (Guimond et al. 2024 §3.2, after Teske et
al. 2013; Brewer et al. 2016). For M dwarfs — the hosts that dominate temperate rocky planet
targets, so the population that actually matters here — the limits are **tighter**, not looser:

| constraint | value |
|---|---|
| Nakajima & Sorahana (2016), 46 M dwarfs, K-band CO + H₂O | **none with C/O > 0.8** |
| Gaidos (2015), CONCH-SHELL | C/O ≈ 1 in **< 1.2×10⁻³** (95% conf.) |
| Gaidos (2015), SDSS DR7 | **< 6×10⁻⁴** (99% conf.) |
| Gizis et al. (2016) | **< 1%** have 0.8 < C/O < 1 |

Two things are worth knowing about that table. M dwarf C/O is measured *more* reliably than FGK
C/O, not less: the CO + H₂O method is non-differential and needs no assumed solar C/O, whereas the
FGK optical C I / [O I] route is differential and the literature disagrees violently (some find
C/O > 0.8 in >20% of stars, others in none). And the historical "carbon-rich M dwarfs" were
artefacts — weak TiO reads as high C/O but is also what low metallicity looks like; Gaidos found
all his candidates were metal-poor stars or systematic errors. The metallicity trend also runs the
wrong way for high C/O, since C/O *rises* with [Fe/H].

**(ii) The non-stellar route exists, but produces a different population and an inert crust.**
Li et al. (2026) point out that the soot line (~500 K) lies inside the water-snow line (~160 K), so
anything forming beyond the snow line also formed beyond the soot line and should carry refractory
organic carbon — up to 40% of the solid mass, judging by comets. Their soot planet (74% rock,
26% soot) reaches a **bulk planetary C/O ≈ 1 at a solar stellar C/O of ~0.55**. Location, not
composition. But: that carbon is refractory organic CHON (C:H:O ≈ 100:77:14), not carbide; the
archetypes are dry (no ocean) or carry 25–50% H₂O by mass (orders beyond any `ocean_depth` here),
so neither is a kamino planet; and reduced carbon is kinetically inert in a cold abiotic ocean.
Graphite is the stable polymorph at seafloor pressures, but its abiotic oxidation is negligibly
slow — on Earth that flux is microbially mediated and O₂-dependent, and this model's ocean is
abiotic with `fO2 = 0`. So soot would act as a **diluent**, exactly as Quartz does at low Mg/Si
(§4): to first order a soot planet's crust weathers like a silicate crust diluted by ~26%.

*(Diamond does not arise at all. It needs ~1.76 GPa at 25 °C, against ~0.03 GPa at a 3 km ocean
floor and ~0.49 GPa at the deepest ocean swept — and liquid water gives way to ice VI near 1.0 GPa,
below diamond's field. There is no habitable pressure at which a liquid-water seafloor sits on
diamond-stable rock.)*

**(iii) Its real effect is on the source term, which is already a parameter — and one whose
speciation assumption fails first.** Crust composition sets the weathering **sink**, a feedback
bounded by crust production and land area. C/O sets the carbon **source**, a linear driver. Kamino
already carries that as `outgassing`. Two things follow:

- *Magnitude.* Li et al.'s soot is 79.9 wt% carbon, so a 26 wt%-soot planet holds ~208,000 ppm bulk
  carbon against Earth's mantle 10–100 ppm — a factor of 2,000–21,000. The sweep tops out at 10×,
  so such a planet sits **200–2,000× beyond the swept range**, and the model saturates long before
  that (§8 of the development history records 50 runs pinned at 5.00 bar with `acid_ocean`).
- *Speciation.* `outgassing_flux[c_idx]` is inorganic carbon, tied to alkalinity. A soot planet is
  reduced enough that volcanogenic carbon emerges as **CH₄ and CO** (Li et al.: soot/oxidised-iron
  ≈ 1.08, against their 0.3 threshold for neutralising ferric iron). CH₄ does not form carbonic
  acid, so the ocean chemistry barely registers it — but `get_T_surface(S, P_CO2, albedo)` has no
  methane axis. **A carbon-rich planet therefore breaks this model at the climate interface, not
  the crust interface.**

The one place C/O leaves a mark that *is* captured: a soot-rich planet is deeply reduced, which
puts it at the low-ΔIW end of the existing redox axis. That corner is thus motivated by two
independent formation routes — enstatite-chondrite-like accretion (Young et al. 2023) and soot
accretion — not one.

---

## 2. Building the mantle  [`mantle_composition`]

McDonough & Sun (1995) pyrolite, modified on the two axes **in this order** so they stay orthogonal:

1. **Iron first.** Set FeO from ΔIW; renormalise the non-Fe oxides to (100 − FeO) at their pyrolite
   proportions. Physically: iron that went to the core is no longer in the mantle, so the remaining
   budget rescales.
2. **Then Mg/Si.** Within that budget, re-split MgO/SiO₂ at fixed (MgO + SiO₂) mass to hit the
   target molar ratio, leaving Al, Ca, Na, Ti, Cr untouched.

Applied the other way round the axes would not be independent. Non-Fe pyrolite (wt%): SiO₂ 45.00,
Al₂O₃ 4.45, CaO 3.55, MgO 37.80, K₂O 0.029, Na₂O 0.36, TiO₂ 0.201, Cr₂O₃ 0.384.

**Ferric iron is off** — the MAGEMin `O` component is 0.0, i.e. all iron ferrous. Fe³⁺ is
second-order for melt major elements and enabling it would perturb the Earth anchor. Guimond et al.
run their pMELTS grids the same way.

---

## 3. Melting  [`isentropic_melt`, `T_at_entropy`, `Tp_for_F`]

**Code:** MAGEMin (Riel et al. 2022) with the Holland, Green & Powell (2018) igneous dataset.
Chosen over pMELTS because pMELTS fails above molar Mg/Si ≈ 1.6 (its solution models do not span the
ferropericlase-bearing assemblages stable there) and extrapolates badly below ≈ 0.8, returning
~69 wt% SiO₂ rhyolites. MAGEMin converges across 0.5–2.0 and stabilises nepheline and ferropericlase
on its own.

**Path:** isentropic decompression melting, 3.0 → 1.0 GPa in 0.2 GPa steps, batch (the melt stays
with the residue, so the final liquid is the pooled primary melt).

- The start temperature is set from a solid adiabat, T(P) = (T_p + 273.15)·exp(αP/ρc_p) − 273.15,
  with α = 3×10⁻⁵ K⁻¹, ρ = 3300 kg m⁻³, c_p = 1200 J kg⁻¹ K⁻¹.
- MAGEMin is a fixed-(P,T) Gibbs minimiser with **no isentropic mode**, so the isentrope is tracked
  by root-finding the temperature that holds entropy constant at each pressure step (secant
  iteration, clamped to ±60 K per step so a flat secant near a phase boundary cannot throw the
  iterate across the melting interval). This matters: prescribing the solid adiabat instead ignores
  latent heat and over-melts, by ~50 K at T_p = 1350.
- **1.0 GPa is the melt segregation pressure, not the base of the crust.** Under batch melting the
  liquid re-equilibrates at every step, so carrying it to 0.2 GPa yields an over-equilibrated
  andesite. 1 GPa is a representative mean segregation depth beneath a ridge.

**Closure — melt fraction fixed at F = 0.20.** This is the value Guimond et al. (2024) adopt, on the
grounds that for typical Earth mantle it is the degree at which clinopyroxene leaves the melting
assemblage, beyond which melting is much less productive (Katz et al. 2003). Using their number
makes the table directly comparable to their Figure 9.

T_p is then **solved** per composition, by bisection on F over [1150, 1900] °C to ±5 °C. Holding F
rather than T_p is deliberate: mantle temperature is self-regulated, not free, because a mantle that
cannot melt cannot transport heat by magmatism. Holding T_p fixed instead makes refractory
compositions look like planets that barely melt, which is an artifact of the closure, not a result.

**Validation of the closure.** Solving instead for clinopyroxene exhaustion at the Earth anchor
returns **F = 0.213** — a 6.5% independent confirmation of Guimond et al.'s choice of 0.20, from a
different thermodynamic dataset. It is not an assumption inherited on trust.

---

## 4. Melt → minerals  [`cipw_norm`]

The norm itself is **pyrolite's implementation of the CIPW norm** (Williams et al. 2020). This
pipeline restricts the input to the oxides the database can express, and then converts the
normative phases that have no counterpart in it. Cations are allocated in the standard order, then
the silica balance is closed by a desilication cascade:

1. Albite (NaAlSi₃O₈), then Anorthite (CaAl₂Si₂O₈) — feldspar
2. Diopside (CaMgSi₂O₆) — remaining Ca paired with Mg
3. All Fe → Fayalite (Fe₂SiO₄); Mg split between Enstatite and Forsterite by the silica balance
4. **5b — silica-deficient:** Albite → Nepheline (NaAlSiO₄), releasing 2 mol SiO₂ each
5. **5c — still deficient:** 2 Diopside → Akermanite (Ca₂MgSi₂O₇) + ½ Forsterite + 1.5 SiO₂,
   releasing 0.75 mol SiO₂ per mol diopside converted
6. **Silica-oversaturated:** excess retained as Quartz
7. **Still deficient after both cascades:** raise. A residual deficit means the norm would assign
   more SiO₂ to minerals than the rock contains. No grid cell triggers this.

Deliberate simplifications, all inherited from the database's phase list: TiO₂ and P₂O₅ are dropped;
all iron is routed to fayalitic olivine (the only Fe silicate available); Diopside and Enstatite are
pure Mg endmembers; K is not a tracked ion so K-feldspar is off (MORB is K-poor, ~0.1% K₂O).

### 4.1 Restriction and correction

**Oxides are restricted BEFORE the norm runs**, not after. TiO₂, K₂O, MnO, P₂O₅ and Cr₂O₃ are
removed and the remainder renormalised. Removing their products afterwards does not work: standard
CIPW allocates Fe to magnetite and ilmenite and K to orthoclase/leucite *ahead of* the
ferromagnesian phases, so deleting those products strands their cations and leaves the silica
balance — which sets the olivine/pyroxene split — wrong. Measured: post-hoc deletion puts normative
olivine 60% out and orthopyroxene 4× out at the Earth anchor; pre-filtering brings feldspar
agreement to <0.01 wt%.

**Four normative phases have no database counterpart** and are converted by balanced reactions
using only available phases:

| | reaction | silica |
|---|---|---|
| Fe-clinopyroxene | hedenbergite + ½ forsterite → diopside + ½ fayalite | neutral |
| *(olivine-free fallback)* | hedenbergite + enstatite → diopside + ferrosilite | neutral |
| Fe-orthopyroxene | ferrosilite + ½ forsterite → enstatite + ½ fayalite | neutral |
| *(olivine-free fallback)* | 2 ferrosilite → fayalite + SiO₂ | releases |
| larnite | larnite + 2 diopside → 2 åkermanite + SiO₂ | releases |
| resilication | nepheline + 2 SiO₂ → albite | reabsorbs |

The first two exist because the database has no Fe-pyroxene, so all iron must sit in fayalitic
olivine; normative hedenbergite is present in **every** composition. The fallbacks are needed
because silica-oversaturated melts (Mg/Si ≲ 0.7) carry no normative olivine to exchange iron with.
Larnite is routed through diopside rather than enstatite because **none of the 46 silica-deficient
compositions carries free enstatite**, and the only other balanced route produces wollastonite,
which floods Ca.

> ⚠️ **pyrolite invents ferric iron unless stopped, and the documented way does not stop it.**
> `Fe_correction=None` means *use the default* (`normative.py`: `if Fe_correction is None:
> Fe_correction = "LeMaitre"`), and that default assigns Fe₂O₃/FeO from a **TAS rock-type
> classification** — producing ~3.6 wt% normative magnetite from an all-ferrous melt and silently
> contradicting the ΔIW axis. The correction is skipped only where **both** FeO and Fe₂O₃ exceed
> zero, so `Fe2O3` must be a tiny **positive** sentinel (1e-9), never 0.0. This is undocumented and
> version-fragile: `cipw_norm` raises if any unexpected phase appears, naming magnetite explicitly,
> so a pyrolite upgrade that breaks the trick fails loudly rather than quietly.

**Cross-check against a hand-rolled norm.** The original bespoke implementation is retained as
`_cipw_norm_native`. Across all 153 grid cells the two agree to ≤0.01 wt% on quartz, anorthite,
enstatite and fayalite, and **107 of 153 cells agree to <0.5 wt% on every phase** — precisely the
107 cells with no silica deficit. The remaining 46 differ only in how that deficit is absorbed
(native converts diopside→åkermanite directly; pyrolite absorbs it into larnite first), reaching
6.2 wt% on diopside in the Mg/Si ≥ 1.7 corner already flagged as ultracalcic.

**Cost.** pyrolite's CIPW runs at ~1.6 s per call against 0.03 ms for the hand-rolled norm.
`mineral_composition` is therefore `lru_cache`d on (Mg/Si, ΔIW): a sweep holds composition fixed
while varying instellation, outgassing and crust production, so without the cache every run would
pay the full cost again for an identical answer.

**Why akermanite and not larnite.** Textbook CIPW desilicates to larnite (Ca₂SiO₄), and Kinec_v3.dat
has larnite complete — both thermodynamics and kinetics — so it was the zero-proxy option. It was
rejected because (i) its k_eff at 300 K / pH 6 is 606× wollastonite and ~10⁵× diopside, so at the
2–8 wt% the cascade produces it would supply essentially the entire dissolution flux, and (ii) it is
a cement clinker phase, rare in nature, whereas these melts are melilititic (Médard et al. 2004) and
melilitites crystallise melilite. Akermanite is the phase that is actually there.

**Akermanite's kinetics are a proxy — this must be stated.** Its thermodynamics are real (llnl
lineage, in `hybrid_ocean.dat`), but no measured rate exists in any available database (checked
Kinec_v3, Kinec.v2, llnl, core10, Thermoddem — PHASES only in all five). Melilite is a sorosilicate
(Si₂O₇ dimers), structurally between the orthosilicates and the chain silicates, which brackets the
rate: wollastonite 2.7×10⁻⁹ (fast) / forsterite 4.5×10⁻¹⁰ (default) / diopside 1.5×10⁻¹¹ (slow),
mol m⁻² s⁻¹ at 300 K, pH 6. Across that bracket, **Ca delivery spans 48×** while Mg and alkalinity
span only 1.8× and 1.6×; at the fast bound akermanite supplies 98% of all Ca. Results from
akermanite-bearing cells that depend on Ca or carbonate burial must sweep the bracket
(`set_akermanite_proxy`); salinity, alkalinity and Mg results are insensitive to it.

---

## 5. The table and its interpolation  [`load_crust_table`, `oxide_composition`]

The melting calculation runs offline and writes `crust_compositions.csv`: a full **17 × 9** grid
(Mg/Si 0.5–2.0 non-uniform, ΔIW −5 to −1 at 0.5 spacing — uniform in log FeO). At runtime the model
bilinearly interpolates it. FeO is interpolated in **log** space because it spans nearly two orders
of magnitude across the redox axis; a linear interpolant is up to 17% high mid-interval. Off-grid
requests raise rather than extrapolate.

---

## 6. Validation

### 6.1 Internal consistency

| check | result |
|---|---|
| Grid completeness | 153/153 cells, 0 failures |
| Mass balance, every oxide | **153/153 within 0.05 wt%** |
| Earth anchor FeO | 8.050 wt%, exact by calibration |
| Earth melt | T_p 1383 °C, F 0.201, SiO₂ 47.88, CaO/Al₂O₃ 0.85 — basaltic |
| Misfit to PRIMELT primary melt | 4.85 wt% (vs 2.22 at F = 0.117; the cost of adopting F = 0.20) |
| PHREEQC end-to-end | equilibrates at Mg/Si 0.5, 1.25, 2.0 |

Run `python src/kamino/data/check_crust_table.py` to reproduce.

### 6.2 Against Guimond et al. (2024)

| claim | theirs | ours | verdict |
|---|---|---|---|
| Mantle olivine-free below Mg/Si | ≲ 0.8 | **0.70** | ✓ (see note) |
| Mantle orthopyroxene-free above Mg/Si | ≳ 1.6 | **1.50** | ✓ (see note) |
| Excess SiO₂ / (Mg,Fe)O form their own phases | predicted | quartz at ≤ 0.6, ferropericlase at ≥ 1.7 | ✓ |
| Melt-vs-mantle relations (§4.2) | 5 relations | **20/20** across 4 compositions | ✓ |
| Melt variance vs mantle variance (Fig. 10) | SiO₂/MgO less, others greater | **6/7** | ✓ (see note) |
| Mantle FeO controls melting temperature (§4.1) | up to ~100 °C | **+53 °C mean** | ~ |
| F = 0.20 is where cpx is lost | assumed | **F = 0.213** solved | ✓ |

*Boundary note.* Their boundaries are subsolidus mantle mineralogy; ours is the residue after 20%
melting, which is depleted and therefore more olivine-rich. Olivine should therefore persist to
*lower* Mg/Si than in their mantle, and opx should vanish at *lower* Mg/Si — both offsets are in the
predicted direction and within 0.1.

*Variance note.* SiO₂ fails (1.11) on the full grid but passes (0.76) once restricted to the
ordinary olivine + opx field they sample; the failure is our uniform grid oversampling the exotic
corner. FeO remains marginally below 1 (0.88–0.98) in every subset — a genuine small disagreement,
attributable to our ΔIW axis spanning a 100× range in mantle FeO where their Hypatia population
varies it only modestly. Not like-for-like on that element.

*Boundary reproduced without being encoded.* Guimond name mantle (Mg+Fe)/Si as a better predictor of
olivine/orthopyroxene than Mg/Si, with the boundary near 0.8. Along Mg/Si = 0.5 our melt switches
from 71 wt% SiO₂ (quartz-normative) to 50 wt% (basaltic) exactly between the cells where (Mg+Fe)/Si
crosses 0.702 → 0.894. Nothing in the pipeline knows about that threshold.

### 6.3 Against experiment — Brugman et al. (2021)

Piston-cylinder melting of HEX1, a hypothetical exoplanet mantle at molar Mg/Si = 1.42.

**Bulk composition.** Constructing our mantle at Mg/Si 1.42 and FeO 8.23 wt% (ΔIW = −1.98)
reproduces their independently-designed starting material:

| oxide | HEX1 | ours | rel. |
|---|---|---|---|
| SiO₂ | 42.00 | 42.40 | +1.0% |
| MgO | 40.04 | 40.39 | +0.9% |
| FeO | 8.23 | 8.23 | 0.0% |
| molar Mg/Si | 1.421 | 1.420 | — |

Al₂O₃ (−8%), CaO (−6%) and Na₂O (+71%) differ because we hold them at pyrolite by construction.
The Mg–Si–Fe backbone — the part the two axes control — was never fitted to HEX1.

**Melt, at matched conditions** (their exact bulk, isobaric 1.5 GPa, F = 0.05): we get SiO₂ 45.37,
MgO 11.50, FeO 7.26, Al₂O₃ 20.30, CaO 12.65 against their Table 5 average of 48.09 / 9.17 / 11.16 /
15.65 / 13.92. Agreement is moderate: Al₂O₃ is 4.7 wt% high and FeO 3.9 wt% low. Their Table 5 is
averaged over 1.0–2.0 GPa and F = 0.004–0.056, so it is not a single-condition comparison; the FeO
deficit is the direction expected from our all-ferrous assumption, since Fe³⁺ is strongly
incompatible and would raise melt FeOtot. **This is the weakest link in the validation chain and
should be reported as such.**

Their solidus statement — HEX1's is *"the same as Earth's nominally anhydrous peridotite solidus"* —
is consistent with ours: T_p at Mg/Si 1.4 and 1.25 differ by only 12 °C.

### 6.4 Against the standard melting parameterisation — Katz et al. (2003)

Independent of both MAGEMin and Guimond. Pyrolite at ΔIW = −2, isobaric:

| P (GPa) | F | T MAGEMin | T Katz | diff |
|---|---|---|---|---|
| 1.0 | 0.10 | 1287 | 1286 | +1 |
| 1.5 | 0.05 | 1338 | 1316 | +22 |
| 1.5 | 0.10 | 1359 | 1341 | +18 |
| 1.5 | 0.20 | 1391 | 1381 | +10 |
| 2.0 | 0.10 | 1423 | 1394 | +29 |

**Mean offset +16 °C, sd 10 °C, max 29 °C** over 1.0–2.0 GPa and F = 0.05–0.20. MAGEMin is
consistently marginally hotter for a given melt fraction, i.e. slightly less melt-productive than
Katz's fit to the experimental compilation.

Katz's own cpx-out melt fraction at these pressures is F = 0.23–0.26. So three independent estimates
of where clinopyroxene leaves the assemblage — Guimond's assumed 0.20, our MAGEMin solve 0.213, and
Katz 0.23–0.26 — agree within about 20%. The closure is not an arbitrary inheritance.

---

## 7. Known limitations, for the caveats paragraph

1. **Ultracalcic melts at high Mg/Si.** 41/153 cells have CaO/Al₂O₃ > 1.2, reaching 1.87, against
   MORB's ~0.78. Médard et al. (2004) exclude such melts from a volatile-free fertile lherzolite
   source, which is what this is. This is **not** a closure artifact above Mg/Si ≈ 1.5: those mantles
   are already orthopyroxene-free, so the melts are ultracalcic however far melting proceeds
   (verified — solving for cpx-out instead still gives CaO/Al₂O₃ = 1.28 at Mg/Si 1.6).
2. **Silicic melts at low Mg/Si.** 27/153 cells give SiO₂ > 60 wt%, up to 76 wt% at Mg/Si 0.5 under
   reducing conditions — quartz-normative crusts. Physically coherent (those mantles are
   olivine-free and quartz-bearing, and pMELTS gives ~69 wt% at the same corner), but at the edge of
   what a peridotite melting model should be trusted for.
3. **The axes are coupled in nature.** Si is siderophile, so core formation raises mantle Mg/Si
   relative to bulk (Guimond et al. 2024, §3.2), and in Young et al.'s H₂-ingassing mechanism the
   same reaction that oxidises the mantle drives Si into the core. The grid treats Mg/Si and ΔIW as
   orthogonal; real planets do not populate it uniformly.
4. **Akermanite kinetics are proxied** — see §4.
5. Guimond et al. (2024) §24.6 note that ~89% of Hypatia stars fall in the ordinary olivine + opx
   field, so the exotic ends of both axes are rare.

---

## 8. Citations

- **MAGEMin:** Riel, Kaus, Green & Berlie (2022), *G³* **23**, e2022GC010427
- **Thermodynamic dataset:** Holland, Green & Powell (2018), *J. Petrol.* **59**, 881
- **Pyrolite:** McDonough & Sun (1995), *Chem. Geol.* **120**, 223
- **Melt fraction, mantle/melt mapping, mineralogical boundaries:** Guimond, Wang, Seidler, Sossi,
  Mahajan & Shorttle (2024), *Rev. Mineral. Geochem.* **90**, 259
- **cpx-out as the productive-melting limit:** Katz, Spiegelman & Langmuir (2003), *G³* **4**, 1073
- **ΔIW convention, Earth's core-formation fO2:** Young, Shahar & Schlichting (2023), *Nature*
  **616**, 306
- **Ultracalcic melt origins:** Médard, Schmidt & Schiano (2004), *Contrib. Mineral. Petrol.*
- **Fe partitioning as a mineralogy control:** Putirka & Rarick (2019), *Am. Mineral.*
- **Dissolution kinetics:** Hermanska et al. (2022, 2023), Kinec_v3.dat
- **Experimental exoplanet-mantle melts:** Brugman, Phillips & Till (2021), *JGR Planets*
  **126**, e2020JE006731
- **Peridotite melting parameterisation:** Katz, Spiegelman & Langmuir (2003), *G³* **4**, 1073
- **CIPW norm implementation:** Williams et al. (2020), *pyrolite*, JOSS 5, 2314
- **IW buffer:** Frost & McCammon (2008) / O'Neill (1987)
- **M dwarf C/O:** Nakajima & Sorahana (2016), *ApJ* **830**, 159; Gaidos (2015); Gizis et al. (2016)
- **Soot planets:** Li, Bergin, Hirschmann, Blake, Ciesla & Kempton (2026), *ApJ Lett.* **997**, L29
