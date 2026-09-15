# Earth calibration and the habitable zone

**2026-09-09.** Companion to `development_history.md` §35. This file records the recalibrated
Earth anchor and the three habitable-zone definitions the model now supports, with the
measurements behind each. Numbers are from the `sweep_output` continental baseline run of
2026-09-09 (48 runs, 24 instellations × land 0.3 and land 0, reducing ocean, 3 km, 1× outgassing
and crust production).

---

## 1. The calibration

### 1.1 What was wrong

Three defects, all found by pulling on one thread: DIC in the Earth-like baseline was 180× modern
seawater.

The cause was not the carbon cycle. Since §7 made alkalinity charge-derived, tracked `Alk` equals
the conservative-ion charge residual **exactly** — measured `Na 678.6 + 2(76.8) + 2(0.36) − 235.1
= 597.8 mEq` against a tracked 597.7. DIC is slaved to that, so a carbon excess is a charge-balance
error wearing a disguise. Decomposed against Earth:

| ion | model | Earth | Δcharge | share of excess |
|---|---|---|---|---|
| Cl | 235.1 | 550 | +315 | 54% |
| Na | 678.6 | 480 | +199 | 34% |
| SO₄ | 0.0 | 28 | +56 | 9% |
| Mg | 76.8 | 52.8 | +48 | 8% |
| Ca | 0.36 | 10.3 | −20 | −3% |
| K (untracked) | — | 10.2 | −10 | −2% |

1. **The ocean was never seeded.** `time_evolve` zero-fills `Y0` and no sweep passed `b0`, while
   `calibrate_earth.py` seeded modern seawater — so the calibration and the production sweeps were
   running different oceans. Two species make that fatal rather than merely slow. SO₄ has **no
   source term** (`planet.py` pins `F_net[so4_idx] = 0`), so a zero seed is permanent. Cl relaxes
   on **τ = 5571 Myr** against a 2 Gyr integration: from blank that reaches 30.2% of steady state,
   and 0.302 × 780 = 235 mM, the observed value to three digits. It is also why every baseline run
   terminated `timeout` rather than `converged`.

2. **`K_CL_ANALYTIC` had been wrong since §33.3.** It derived Cl's rate constant from a
   steady-state balance in which the source and sink areas cancel — true until the seafloor-area
   fix, false after, because outgassing is emitted over the whole sphere and Cl subduction acts
   over the seafloor only. Missing factor `1/(1 − land_fraction) = 1.43`. This is an **input** to
   the fit, so it biased every other parameter.

3. **The seafloor-area fix's 1.43× was absorbed** by `KD_MG_CALIB` and `K_NA_CALIB` (§33.3, known).

Seeding is also the physically correct choice, not a convenience: seawater Cl is an inherited
inventory from early degassing, with modern volcanic Cl being recycled subducted seawater Cl
rather than primordial (Kendrick et al. 2021, PNAS; Sharp & Draper 2013, EPSL). A blank ocean is
the unphysical case.

### 1.2 How the fit was restructured

Two failed attempts preceded the working one, and both failures are informative.

**Three parameters against Na/Ca/Mg — ill-posed.** Measured over a 100× scan at Earth, `alpha`
moves the seafloor alkalinity flux **94×** while moving every ocean concentration **<25%** and T
by 0.7 K. So the ocean residuals carry almost no gradient in `alpha`, and `least_squares` —
starting on the calcite-collapsed branch — used it as a free escape, driving `alpha` 1.100 → 1488
in a single trust-region step. It "converged" at 1462.75, giving a seafloor flux of 30.68 Tmol/yr
against the ~1 Tmol/yr anchor, at a cost (0.593) **worse than an unvisited point at `alpha` = 3**
(0.390).

**Pin `alpha`, fit the rest — better but not separable.** `alpha` barely moves the ocean, but the
ocean strongly moves the flux: the fit pulled Ca 4.47 → 10.04 mM, bringing the pore fluid closer
to saturation with the Ca-bearing primary phases, and the flux fell 1.00 → 0.32 Tmol/yr at
unchanged `alpha`. The anchor it had been pinned for no longer held once the fit had moved.

**Joint fit, four residuals — this works.** `(K_na, alpha, KD_mg)` against
`(Na, Ca, Mg, seafloor Alk → 1 Tmol/yr)`, all log-space, equally weighted. `alpha` gets a residual
that responds strongly to it (flux ∝ `alpha`^0.9866 measured over 3–300) while the others keep the
ones that respond to them. No outer loop and no separability assumption. Controlled by
`FLUX_TARGET` in `calibrate_earth.py`; `ALPHA_PINNED` retains the two-parameter mode.

> The equal weighting is a choice. The concentration targets are known to ~1%; the 1 Tmol/yr
> figure is a literature estimate with real spread (Coogan & Dosso). **If the flux ever fights the
> ocean, down-weight the flux** rather than assume the ocean is wrong.

### 1.3 The result

```
K_CL_SUBDUCTION   = 1.961786e-04      (planet.py)
K_NA_CONT_REMOVAL = 6.099720e-03      (planet.py)
KD_MG_HT          = 1.969604e-02      (planet.py)
ALPHA_REF         = 4.900000          (weathering.py; parameter_sweep mirrors all four)
```

Earth, S = 1, land 0.3, 3000 m — `converged`, T = 294.3 K, pCO₂ 686 ppm, pH 7.71:

| ion | model | Earth | error |
|---|---|---|---|
| Cl | 546.00 | 546 | **+0.0%** |
| Ca | 10.10 | 10.3 | **−1.9%** |
| Na | 427.46 | 469 | −8.9% |
| Mg | 74.40 | 52.8 | +40.9% |
| Alk | 3.56 | 2.3 | +54.8% |
| C | 3.19 | 2.1 | +51.7% |
| Si | 4.68 | 0.1 | +4575% |

Against the pre-calibration baseline: DIC 363 → 3.19 mM (**112×**), Alk 598 → 3.56, Ca 0.36 →
10.10, Cl 235 → 546. Runs now converge at ~1 Gyr instead of hitting the 2 Gyr timeout.

### 1.4 The two large residuals are both expected, and neither is a fit failure

**Mg +40.9% is structural.** §28.1 independently measured +37.0% and traced it to §27: removing
the hedenbergite correction reaction shifted the Earth assemblage (Diopside −25.7%, Forsterite
+14.9%, Ca/Mg supply ratio −19.4%). Two independent routes to the same number. `KD_mg` trades Mg
for Ca mole-for-mole, and the split diagnostic reads `(Ca/Ca_t)/(Mg/Mg_t) = 0.698` — with Ca
already at −2.4%, spending more Mg overshoots Ca. It needs a crust-composition change or a real Mg
sink (§12), not a constant. The alkalinity and DIC excess follows arithmetically: Mg's surplus
carries ~42 mEq and Na's deficit offsets only ~40.

**Si +4575% is the absence of a biosphere.** Modern seawater Si is drawn down to 0.1 mM by
diatoms; an abiotic ocean sits near amorphous-silica saturation, which is roughly where this is.
It falls on the "biotically controlled" side of the ion figure's own divider. Caption it, don't
fix it.

### 1.5 Why pCO₂ is 686 ppm, not 280

Measured carbon budget at the calibrated steady state (Tmol/yr):

```
outgassing        C +7.498
continental     Alk +18.830   Ca +5.462   Mg +2.124   Na +3.659
HT Mg-Ca exch.               Ca +2.037   Mg −2.037
seafloor LT     Alk +0.034
precipitation     C −3.320    Ca −3.320
shelf carbonate   C −4.178    Ca −4.178
Na albitization Alk −3.661                Na −3.661
```

Every mole of carbon is buried as calcite, so **carbon burial = Ca supply**: 5.462 + 2.037 = 7.499
in, 3.320 + 4.178 = 7.498 out, equal to outgassing. pCO₂ is not free — it rises until continental
weathering delivers 5.46 Tmol/yr of Ca, and 686 ppm / 294.3 K is what that costs (a 1.91× WHAK
enhancement over the reference).

The root cause is that two Earth constants do not agree. `EARTH_CONTINENTAL_WEATHERING_REF = 8
Tmol/yr` is an *alkalinity* flux; converted, that is Ca+Mg = 4 Tmol/yr and, at the model's 28% Mg
fraction, **2.88 Tmol/yr of Ca**. With HT's 2.04 the Earth reference supplies 4.9 Tmol/yr against
`EARTH_OUTGASSING = 7.5`. The model closes the 1.53× gap the only way it can — by warming.

**Options, in order of preference:**

1. **Accept it, and quantify it.** The model is abiotic. Schwartzman & Volk (1989) found that if
   modern weathering is 10× the abiotic rate an abiotic Earth would be **~15 °C warmer**; 100×
   gives ~30 °C. This model is **+6.3 K**, implying a biotic enhancement factor well under 10 —
   so 686 ppm is *conservative* against that literature, not excessive. State the implied factor.
2. **Raise `EARTH_CONTINENTAL_WEATHERING_REF` to ~15 Tmol/yr.** That balances at 288 K / 280 ppm
   by construction. Gaillardet et al. (1999) and later compilations put global silicate CO₂
   consumption at 11.7–17.9 Tmol/yr, which for silicate weathering equals the alkalinity flux in
   equivalents — so 15.2 is inside the range and 8 is at its bottom edge.
3. **Let Mg bury carbon** (dolomite / Mg-calcite). Relieves the Ca limitation *and* the Mg excess,
   same root cause — but §12 warns that picking Mg sinks off an SI ranking is the trap.
4. Lower `EARTH_OUTGASSING`. Within range (6–10 Tmol/yr) but it rescales every sweep's outgassing
   axis.

Note also that **~2.3 K of the +6.3 K is not the carbon cycle**: §33.4 measured this climate model
giving 290.3 K at S = 1, 280 ppm against Earth's 288.

---

## 2. The habitable zone

### 2.1 Three nested constraints

Ceiling 340 K throughout (the post-processing habitability cut). Continental edges are from
converged runs only.

| constraint | floor | outer S | inner S | width |
|---|---|---|---|---|
| Radiative only (climate model) | 260 K | 0.383 | 1.142 | 0.759 |
| Radiative only (climate model) | 273 K | 0.410 | 1.142 | 0.732 |
| Continental (land 0.3) | 260 K | 0.488 | 1.125 | 0.637 |
| Continental (land 0.3) | 273 K | 0.745 | 1.125 | 0.380 |
| Ocean world (land free) | either | — | — | **no converged runs** |

The radiative zone is the envelope of `T(S, pCO₂)` scanned over pCO₂ from 1 Pa to 10 bar with no
chemistry: the outer edge is the coldest instellation whose *best* pCO₂ still clears the floor
(maximum greenhouse), the inner edge the hottest whose *bare* atmosphere stays under the ceiling.
It brackets what any carbon cycle could achieve, so a chemistry-constrained zone must sit inside
it — and does.

The climate envelope, for reference:

| S | T at bare minimum pCO₂ | T at best pCO₂ | pCO₂ at max |
|---|---|---|---|
| 0.35 | 180.0 | 240.5 | 1.90 bar |
| 0.40 | 180.9 | 268.6 | 2.78 bar |
| 0.45 | 190.9 | 288.2 | 3.54 bar |
| 0.50 | 205.7 | 302.9 | 4.36 bar |
| 0.80 | 254.9 | 347.6 | 9.33 bar |
| 1.00 | 279.7 | 361.5 | 10.0 bar |

### 2.2 The temperature floor is not a detail

`T_SNOWBALL = 260 K` is a **global mean** in a model with no latitude, so it is not a liquid-water
criterion. Earth's equator runs ~12 K above the global mean (288 vs ~300) and that gradient
steepens as a planet cools, so a 260 K global mean puts the tropics near 272 K — right at
freezing. It is best read as a **waterbelt / Jormungand limit**: the coldest global mean still
permitting an ice-free tropical band. That state is real and studied (Abbot, Voigt & Koll 2011;
Rose 2015), though recent work questions its stability under cloud uncertainty (Hörner et al.
2022). It is **not** comparable to Kopparapu, whose outer edge is explicitly "the furthest point
where CO₂ can keep the surface temperature above 273 K."

**Recommendation: report both, labelled** — 273 K as the conservative Kopparapu-comparable edge,
260 K as the optimistic waterbelt edge.

> ⚠️ **The choice barely matters radiatively and matters enormously with chemistry.** The
> radiative outer edge moves 0.383 → 0.410 (7%); the continental outer edge moves **0.488 → 0.745
> (53%)**, cutting the width from 0.637 to 0.380. The reason is that a carbon-cycle-limited planet
> sits on a very shallow `T(S)` slope in the cold region — T rises only 255.6 → 273.3 K across S =
> 0.40 → 0.75 — so a 13 K change in the criterion sweeps a third of the instellation range. Any
> statement about zone width must name its floor.

### 2.3 Why the outer edge is not Kopparapu's 0.343

**Kopparapu et al. (2013) has no weathering, and therefore no runoff term.** It is a 1-D cloud-free
radiative–convective model; both edges are purely radiative (inner from moist/runaway greenhouse,
outer from maximum greenhouse). Its 0.343 is a **ceiling, not a prediction** — it answers "if a
planet had the optimal CO₂, how far out could it stay above 273 K?", and is silent on whether any
carbon cycle would deliver that CO₂. A model with a carbon cycle *should* land inside it. That is
the standard result (Abbot et al. 2012; Menou 2015; Haqq-Misra et al. 2016; Graham &
Pierrehumbert 2020).

The gap splits into two independent parts:

**(a) Radiative, ~0.07 in S — a genuine like-for-like discrepancy.** At Kopparapu's own 273 K
maximum-greenhouse criterion this climate model gives **S ≈ 0.410** against their 0.343. Both are
purely radiative and the OLR fit (Haqq-Misra et al. 2016) is itself a fit *to* Kopparapu's model,
so they should agree better. Two numbers locate the problem: at S = 0.343 this model peaks at
**1.78 bar reaching only 236 K**, where Kasting/Kopparapu reach 273 K at several bar. Both a lower
optimal pCO₂ and a much colder peak points at the **albedo**, not the OLR — `albedo_funtion`
(`climate/analytic.py:35`) uses `tau_ray = 0.19513 × pCO2_bar`, so by 4 bar the Rayleigh optical
depth is ~0.78 and the planet is already substantially brighter. **Untested; the cheap check is to
suppress the Rayleigh term and re-scan.**

**(b) Carbon supply, ~0.10 in S — real, but partly a missing mechanism.** At S = 0.45 the carbon
cycle delivers **0.71 bar** where 3.54 bar would be optimal, giving 258.3 K instead of 288.2 K —
about 30 K. CO₂ stalls near 1 bar because **weathering never shuts down**: evaluated at S = 0.40
(255.6 K, 1.11 bar) the WHAK law gives 12.0× from the pCO₂ term and 0.149× from temperature, i.e.
**1.79× modern** — §22.8's "near 2× across the entire zone", confirmed. Reaching 3.5 bar would
demand ~9× modern, far more than 1× outgassing can feed.

Physically it should shut down. At a 256 K global mean the continents are glaciated and silicate
weathering should nearly stop, which is how max-greenhouse calculations reach multi-bar
atmospheres. The model cannot represent that:

- `LAND_ALBEDO = OCEAN_ALBEDO = 0.3` (`planet.py:111-112`) — **no ice-albedo feedback**
- **no runoff term** — the law has only pCO₂ and T, but weathering needs liquid water on land
- no supply or thermodynamic limit

**Suggested fix:** add runoff dependence to `get_continental_weathering_flux`. Graham &
Pierrehumbert (2020) show the outer edge is "quite sensitive to parameters representing hydrology,
weathering thermodynamics, surface properties, and albedo"; the Maher–Chamberlain formulation they
use makes weathering kinetically, thermodynamically *or* runoff limited. A minimal runoff factor
going to zero as the surface freezes should move the outer edge from 0.488 toward ~0.39. It will
not reach 0.343 — that last 0.07 is item (a).

> The paper's framing — that a carbon-cycle-limited zone is narrower than the radiative one — is
> **sound and is the standard result**. What needs fixing before *quantifying* the narrowing is
> that part of it is currently a missing weathering shutdown rather than a genuine supply limit.

---

## 3. Open items

- 🔴 **The ocean-world arm produced zero converged runs** (13 `out_of_domain`, 5 `timeout`), so its
  HZ edges are not measured. Any figure quoting them is reading unconverged states. Likely cause:
  `seawater_seed()` is fixed-concentration, but the Cl steady state goes as `1/(1 − land_fraction)`
  through `K_CL_SUBDUCTION`'s area ratio, so **a land-free world equilibrates to 0.7 × 546 = 382
  mM**, not 546. Those runs were seeded ~43% above their own steady state and are draining toward
  it on the 5.6 Gyr timescale — the seeding trap mirrored. Resolve before any land-free sweep.
- ⚠️ **The seed's depth-scaling convention is undecided.** Fixed-concentration vs fixed-inventory
  differ by the full 167× ocean-mass range of the depth sweep (300 m – 50 km). Only fixed
  concentration at land 0.3 is justified so far.
- ⚠️ **`alpha`'s flux anchor measures a quantity `dY_dt` does not use.** The budget delivers 0.034
  Tmol/yr of seafloor alkalinity where `seafloor_alk_flux_tmol` reports 0.96, because the latter
  computes primary dissolution only (no pore clays) *and* passes `rate=rate_ref`
  (`calibrate_earth.py:372`) where `dY_dt` passes `rate=crust_production_rate`. Negligible for
  pCO₂ (seafloor is 0.034 of a 7.5 budget) but it means `alpha` is not yet anchored on the
  quantity intended.
- ⚠️ `plot_results.CONTINENTAL_HZ_OUTER = 0.480`, this sweep measures **0.488**. Update, or the HZ
  lines on every other figure are stale.
- ⚠️ Continental runs at S ≥ 1.20 hit `wall_timeout`; the inner edge is bracketed, not crossed.
- ⚠️ The calibration anchors at 3700 m (0.98× Earth's ocean mass) while the baseline runs at 3000 m
  (0.79×), and `tau_prec` is depth-scaled, so the fit is made at 123 kyr and the figure runs at
  100 kyr (§33.13).

---

## 4. References

- Abbot, D. S., Voigt, A. & Koll, D. (2011), *JGR Atmospheres* **116**, D18103 — Jormungand state.
- Coogan, L. A. & Dosso, S. E. — low-temperature seafloor alteration fluxes (the ~1 Tmol/yr anchor).
- Gaillardet, J. et al. (1999), *Chem. Geol.* **159**, 3 — global silicate weathering and CO₂ consumption.
- Graham, R. J. & Pierrehumbert, R. (2020), *ApJ* **896**, 115 — thermodynamic and energetic limits on continental weathering.
- Haqq-Misra, J. et al. (2016), *ApJ* **827**, 120 — OLR fit; limit cycles narrowing the zone.
- Kendrick, M. A. et al. (2021), *PNAS* **118** — bulk silicate Earth halogen budget.
- Kopparapu, R. K. et al. (2013), *ApJ* **765**, 131 — habitable zone limits.
- Menou, K. (2015), *ApJ* **812**, 36 — plate tectonic–climate coupling and exposed land area.
- Rose, B. E. J. (2015), *JGR Atmospheres* **120**, 1404 — stable waterbelt climates.
- Schwartzman, D. W. & Volk, T. (1989), *Nature* **340**, 457 — biotic enhancement of weathering.
- Sharp, Z. D. & Draper, D. S. (2013), *EPSL* **369-370**, 71 — the chlorine abundance of Earth.
