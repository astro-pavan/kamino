# Kamino — Consolidated Development Context

**Compiled 2026-08-17, extended through 2026-08-19** from all Claude Code sessions on this project
(2026-07-14 → 2026-08-19). Later sessions are weighted more heavily; where an early decision was later
reversed, the reversal is what's recorded, with the original noted as history.

Verified against the working copy on branch `NaCl-chemistry`. **§11** lists what is actually in the
code versus what was only proposed. **§18–§24 are the most recent work** and supersede earlier framing
where they conflict; §20 records the iron charge-leak fix, §21 removes reverse weathering from the pore
space (which invalidates the transition evidence in §20.4, see below), **§24 is the current authority** — it replaces pMELTS with MAGEMin and, in §24.6, records that
the Mg/Si crust question was already published by this group. §22 remains the authority on
CALIBRATION — it records the Earth calibration, shows `alpha` is not identifiable from Earth, and
establishes that the LT seafloor flux carries the wrong ion relative to Coogan & Dosso. §22.0 also
corrects two stale entries in §11/§15.

**§13 is the reference point for results.** The 17 June 2026 seminar
(`Ocean_Chemistry_Seminar-2.pdf`) is the last time the model produced a complete, self-consistent set
of results — the target behaviour.

**Reading §20 after §21.** §20.4 reported the June transition signature reappearing in `fast_15`.
That sweep still had reverse-weathering minerals precipitating in the pore space, which §21.3 shows
inverts the sign of the seafloor alkalinity flux over a broad region of (T, pCO₂) and manufactures a
runaway that looks like a transition. **The `fast_15` reversals should not be cited.** §21.4
(`fast_17`) re-establishes a reversal on one line with the bug fixed, but it is not yet settled.

---

## 1. What the project is

Kamino is a PHREEQC-based geochemical box model of **seafloor (basalt) weathering and ocean
chemistry evolution on exoplanets**, being written up as a paper. The distinguishing goal is that
it works across a **wide range of ocean chemistries**, not just Earth-like ones:

- Mg/Ca/C-dominated oceans (the focus is largely **land-free ocean worlds**)
- salinity up to ~100 g/kg → ionic strength I ≈ 2
- pH up to ~10
- explicit Na, Cl and SO₄ tracking (this is what most comparable models assume away)

The model couples: a climate solver (analytic OLR fit or a grid interpolator) → PHREEQC ocean
chemistry → a Maher & Chamberlain-style kinetic/transport weathering law → mineral precipitation
sinks → tectonic fluxes (outgassing, subduction, hydrothermal circulation), integrated to 2 Gyr
with LSODA.

### Environment and working conventions

- **Always use `/data/pt426/big-venv/bin/python`** to run anything in this project.
- Sweeps live in `/home/pt426/data/kamino_experiments_fast_N/` (also reachable as
  `/data/pt426/...` — same directory via symlink). `fast_4` … `fast_10` are referenced below.
- Sweeps run 28-way parallel with `ProcessPoolExecutor` using `spawn`, so **every worker
  re-imports from `src/`** — editing `src/` mid-sweep silently gives one sweep mixed code
  versions. Always check whether a sweep is running before touching `src/`.
- Because of the above, the established workflow is: develop and test on a **patched copy** in a
  scratch dir or worktree, then **copy the finished files into the working copy uncommitted** for
  the user to review. Do not commit; do not leave a worktree behind.
- The user prefers **simple, pre-existing solutions** over complex installs (this killed the
  pMELTS route, §5).
- Stated constraint: *"I've been avoiding adding K as that's another ion to manage and have a sink
  for."* K-Feldspar therefore defaults **off** in the CIPW norm. Sulfur is deliberately left
  decoupled/pinned for the same reason.
- `notebooks/` and the old `tests/` are stale; the README is out of date. Old tests were deleted
  by the user (`bd5fe3a`). `src/kamino/H21/chili/` is third-party (Hakim et al. CHILI) — **not the
  user's code, ignore it**.

---

## 2. Session timeline

| Dates | Focus |
|---|---|
| **Jun 10–15** *(pre-sessions)* | The **seminar version** (§13): Na–Cl chemistry added, Wollastonite Ca source, PHREEQC-HT *and* parameterized `F_ht_exchange` both active. Produced the last complete set of results, presented 17 June. Not a Claude session — reconstructed from git and the slide deck. |
| **Jul 14** | Code cleanup pass over `src/kamino`: constants → constructor args, `chemistry.py` readability, first look at excluded minerals and the `dissolve_only` flag (later backed out). Diagnostic multi-panel planet plot built in a second session. |
| **Jul 15–17** | Thermodynamic **database** construction (§4). Verifying the equilibrium step of the weathering law. Enstatite/Wollastonite → **Diopside** swap. HT weathering investigation begins: water/rock ratio, Ca–Mg exchange, albitization. |
| **Jul 28** | **Crust composition pipeline** (§5): CIPW norm + PRIMELT1 spreadsheet, T_p and Mg/Si as the two input axes. Then HT albitization resumed with the new self-consistent basalt. |
| **Jul 30–31** | Na sink search (§6). The unified LT+HT scheme. **Alkalinity charge-consistency fix** (§7) — the single most important correctness fix; first temperate Earth. HT Ca–Mg exchange debugging. |
| **Aug 3** | Sweep plumbing (`parameter_sweep.py`, `plot_results.py` rewritten for the new model). Damköhler-number diagnostic corrected (§8). **Decision to abandon PHREEQC-HT and return to the parameterized Mg→Ca exchange** (§6). |
| **Aug 4–6** | Termination-event redesign (§9) and **performance work** (§10): PHREEQC KNOBS, retry-on-non-convergence, continuous fallback, domain-guard fix. Sweeps `fast_6` → `fast_10`. |
| **Aug 17** | `kd_mg_ht` recalibration 0.012 → 0.07 (§18); fast_12; the Na-shadowing diagnosis (§19). |
| **Aug 17–19** | Step-1 buffer / water-rock investigation, the **iron charge-leak bug and its one-line fix** (§20), sweeps `fast_13` → `fast_15`. Transition visible for the first time. |
| **Aug 19** | The `self._T` side-effect corruption, retrograde-solubility confirmation, and the **reverse-weathering pore-space bug** (§21); `fast_16` → `fast_17`, the best-behaved sweep so far. |
| **Aug 19–20** | **Earth calibration** (§22): `K_na` and `kd_mg_ht` fitted and alpha-independent; `alpha` shown to be unidentifiable from Earth; the 1 Tmol/yr anchor found to be 88% Fe measured without pore precipitation; **Coogan & Dosso LT fluxes adopted as the literature target**, exposing the Mg/Ca composition mismatch. |
| **Aug 21** | **MAGEMin + Mg/Si sweep** (§24): pMELTS replaced, Earth basalt calibrated at T_p = 1325, constant-F closure adopted and homologous-temperature rejected, ultracalcic melts diagnosed and fixed by stopping at cpx-out — then **§24.6: most of it was already published by this group** (Guimond et al. 2024). |
| **Aug 20** | **Crust-composition pipeline** (§23): the flat T_p/Mg-Si trends traced to a broken oxide mapping plus CIPW clipping; **pMELTS brought online** (superseding §5) and `make_crust_compositions.py` written; **Nepheline added** to the database, rates, norm and build path. |

---

## 3. Model architecture as it now stands

### State vector

`elements = [Alkalinity, C, Si, Al, Fe, Ca, Mg, Na, Cl, S]`, plus `P_CO2`, `P_H2O` and a
convergence-metric variable `r_avg` (last element). Alkalinity **is** tracked as its own ODE
variable (Coogan & Dosso-compatible), but every flux's alkalinity component is *derived* from ion
charge — see §7.

`Y[0] = P_CO2`, `Y[1] = P_H2O`, `Y[2:-1]` = ocean concentrations, `Y[-1] = r_avg`.

### Flux terms in `dY_dt` (`planet.py:~200–300`)

```
F_net = F_vol            # outgassing: C, plus Cl at cl_outgassing_ratio (HCl acidifies: Alk -= Cl)
      + F_prec           # fast ocean precipitation, tau_prec = 100 kyr
                         #   carbonates + clays + silica + evaporites
      + F_prec_rw        # slow reverse weathering, tau_rw = 5 Myr
                         #   Sepiolite(d), Saponite-Na, Greenalite
      + F_shelf_prec     # shelf carbonate (land worlds)
      + F_diss           # LT seafloor weathering, PHREEQC, full_equilibrium=True, J = J_total
      + F_ht_exchange    # parameterized 1:1 Mg -> Ca swap (see below)
      + F_cont           # continental weathering (only if land_fraction > 0)
      + F_cl_subduct     # Cl removal ∝ crust production
      + F_na_rw          # Na sink, always on
```

Atmosphere: `dYdt[0] = (P_CO2_new − P_CO2) / TAU_ATM` with `TAU_ATM = 10 kyr`.

**Structural fact worth remembering:** there is **no air–sea exchange term in `F_net`**. The
atmosphere relaxes toward the pCO₂ implied by ocean DIC, but carbon leaving the atmosphere is
never added to the ocean. The atmosphere is a *lagged readout* of the ocean, not a mass reservoir.
This is what makes the initial condition inconsistent (§9.3).

### The weathering law

Maher & Chamberlain (2014) / Hakim et al. (2021) form:

```
F = A_r · (b_eq − b_in) / (b_eq/k + A_r/J)
Da = k · A_r / (J · b_eq)
```

Kinetic-limited when Da ≪ 1 (`F ∝ A·k`); transport/thermodynamic-limited when Da ≫ 1
(`F → J·(b_eq − b_in)`). `b_eq` and `k` both come from PHREEQC; `k` uses the pH from the
equilibrium calculation.

The LT path runs with `full_equilibrium=True`: the **full mineral assemblage, no exclusions, no
kinetic Mg offset**. This was an explicit goal — the user wanted the law to work "with less of
these exclusions and exceptions." Pore-space secondary precipitation (clays + reverse-weathering
minerals) happens inside `get_weathering_flux`.

### The HT path

`f_HT` is now **vestigial** (default `0.0`). High-temperature chemistry is a parameterized
exchange, not a PHREEQC call — see §6.4:

```python
_ht_rate = kd_mg_ht * b_ocean[mg_idx] * J_total * surface_area / ocean_water_mass
F_ht_exchange[mg_idx] = -_ht_rate
F_ht_exchange[ca_idx] = +_ht_rate
```

Strict 1:1, alkalinity- and charge-neutral by construction, and **Ca supply ∝ ocean [Mg] ∝
weathering ∝ climate** — that coupling is what carries the thermostat.

### Current tunable defaults (`planet.py`)

| Parameter | Value | Role |
|---|---|---|
| `KD_MG_HT` | **`0.07`** (was `1.197244e-02`) | HT Mg→Ca exchange; recalibrated 2026-08-17, see §18 |
| `K_CL_SUBDUCTION` | `1.373251e-04` | Cl subduction |
| `K_NA_CONT_REMOVAL` | `2.194806e-03` | Na sink (always-on, `J_total`-scaled) |
| `alpha` | `1.43` | base reactive area per unit crust area |
| `f_HT` | `0.0` | vestigial |
| `tau_prec` | 100 kyr | fast precipitation |
| `tau_rw` | 5 Myr | reverse weathering |
| `TAU_ATM` | 10 kyr | atmosphere relaxation |
| `P_CO2_CLIMATE_FLOOR` | 1.0 Pa | climate-input clamp (§9.2) |
| `convergence_threshold` | 0.05 /Gyr | `event_converged` |

> ⚠️ **Memory conflict.** `memory/project_calibration_state.md` records `KD_MG_HT=0.233`,
> `K_NA=3.33e-3`, `f_HT=0.039`, and a convergence threshold of 0.5. **None of these match the
> current code** (values above). That memory predates the Aug 3 switch-back and should be treated
> as stale history, not current calibration.

---

## 4. Thermodynamic databases

The original databases were assembled ad hoc from the PHREEQC installation's bundled files
(`/home/pt426/Code/ocean_chem/databases/database`). The Jul 15–17 work was about finding a
**defensible, citable** basis for the paper.

**The chemistry regime that drives the choice:** a divalent-carbonate brine at I ≈ 2, pH up to 10,
Mg/Ca/C-dominated. That rules out plain Debye–Hückel and makes the activity model the central
question.

- **SIT** (Specific ion Interaction Theory) — valid to ~3–4 mol/kg, the right fit. Sourced from
  **ThermoChimie `sit.dat` v12a** (ANDRA / NWS / ONDRAF). A citation for this was specifically
  requested for the paper.
- **SUPCRTBL** (Zimmer et al. 2016) — extended Debye–Hückel / b-dot; ionic-strength headroom was
  the concern that prompted the SIT move.
- **Kinec_v3** (Hermanská et al. 2022/2023) — kinetic rate data.
- **llnl.dat** (SUPCRT92 lineage), **Thermoddem** (Blanc et al. 2012) — Thermoddem was judged more
  aimed at cement/waste systems than basalt weathering.

Runtime databases: `lt_weathering_sit.dat` (the new SIT low-temperature DB), `hybrid_ocean.dat`
(LT ocean), `hydrothermal.dat` (HT).

`src/kamino/data/make_database.py` builds these: `make_database()` imports the base model,
`retrieve_mineral_thermodynamic_data()` pulls log K, `retrieve_kinetic_data()` pulls rates.
Decisions taken: `parse_stoichiometry` **stays** in `chemistry.py` (not moved); SUPCRTBL is read
from **local files**, not fetched from the SupPHREEQC GitHub repo.

**Naming trap:** silica is `H4(SiO4)` in `sit.dat`, `SiO2` in llnl/Kinec, and `H4SiO4` in
`hybrid_ocean.dat`. Cross-database mixing has to reconcile this.

**Smectite-Na** does not exist in the PHREEQC databases; it was replaced by **Saponite-Na**
(trioctahedral Mg-smectite, a real basalt-alteration clay) and added to `hybrid_ocean.dat`.

**Policy on database size:** the user had deliberately kept the mineral list minimal. The
conclusion was that a *comprehensive* database is safer for thermodynamics, with **kinetic
filtering done explicitly** by the model's exclusion/precipitation lists rather than implicitly by
omitting phases. See §12 for why this matters.

---

## 5. Crust composition pipeline

Originally the basalt composition was hardcoded, eyeballed from
`/home/pt426/Pictures/Mineralogy_igneous_rocks_EN.svg.webp` (the `basalt_49` composition). Replaced
by a systematic bulk-chemistry → mineralogy generator in `src/kamino/crust_composition.py`:

```
T_p, Mg/Si  →  oxide_composition()  →  mineral_composition()  →  crust_composition dict
                    (PRIMELT1 interp)      (CIPW norm)
```

- **`import_primelt_spreadsheet()`** reads `src/kamino/data/ggge967-sup-0002-primelt1.xls`
  (Herzberg & O'Hara 2002 / Herzberg et al. 2007 PRIMELT supplement) and builds the interpolator.
- **`oxide_composition()`** maps mantle potential temperature → degree of melting → primary melt
  oxides, with an Mg/Si modification as a second, orthogonal axis (silica saturation:
  olivine-normative ↔ quartz-normative). Stellar/planetary Mg/Si ≈ 0.7–1.5; Earth ≈ 1.0.
- **`mineral_composition()`** applies the **CIPW norm** (`cipw_norm`, moved out of
  `mineral_info.py`), emitting only phases valid in the PHREEQC database, with K dropped by default.

**Two input parameters** — `mantle_potential_temperature` (note: currently passed in **°C**, not K)
and `mg_si_ratio` — replace the full oxide vector. Earth-like basalt is `T_p = 1350`:

```python
mineral_composition(1350) = {Anorthite 0.361, Albite 0.204, Diopside 0.183,
                             Fayalite 0.111, Forsterite 0.134, Enstatite 0.006}
```

**Routes explored and rejected:** ~~pMELTS/alphaMELTS via PetThermoTools / alphaMELTS-for-Python /
ENKI ThermoEngine (install complexity, GSL soname ABI issues — GSL 2.4→`libgsl.so.23`,
2.6→`.25`, 2.7→`.27`)~~ — **superseded, see §23.5: pMELTS now runs.** alphaMELTS 2.3.1 was already
installed; the only blocker was `libgsl.so.27`, and building GSL 2.7.1 from source takes ~1 minute.
`src/kamino/data/make_crust_compositions.py` uses it. Also rejected:
hand-implementing Katz (2003) + Niu (1997) or the Langmuir/Klein/Plank
(1992) parameterization. `pyMelt` is installed (Katz 2003) but the PRIMELT interpolation was
simpler and sufficient. `pyrolite` was installed for its CIPW implementation.

> **Side effect that cost time:** `pyrolite` writes `pyrolite.mplstyle` into
> `~/.config/matplotlib/stylelib/` on first import and applies it globally. Its line 27
> (`legend.bbox_to_anchor`) is not a valid rcParam, so **every** matplotlib import in **every**
> sweep worker printed a "Bad key" warning — even though nothing in the repo imports pyrolite.
> The file was deleted. A future `import pyrolite.plot` will re-create it.

### The Diopside change and its consequences

`mineral_info.py` compositions were swapped from **Enstatite/Wollastonite → Diopside**. This turned
out to be one of the most consequential changes in the project:

- Wollastonite (June model): fast-dissolving, **Ca-only**, high solubility → an easy, strong Ca source.
- Diopside: **Ca–Mg coupled**, ~40× slower kinetics.

Much of the Jul 30 – Aug 3 "why is Ca starved / why did the June version work and this one doesn't"
debugging traces back to exactly this. The user identified it themselves:
*"Wait, the June version had the Ca in Wollastonite, not Diopside. Wollastonite dissolves very fast."*

---

## 6. The Na, Cl and HT-exchange saga

This is the longest arc in the project and ended in a deliberate retreat from PHREEQC self-consistency.

### 6.1 The Na sink problem

On a land-free ocean world, **nothing in the literature clearly identifies the Na sinks**. Na must
still reach steady state or the ocean runs away. Candidates examined:

- **Halite** — sources were sought for it being Earth's dominant Na sink; evaporites are
  land-dependent, so this doesn't transfer to a water world.
- **Reverse-weathering clays** — PHREEQC tests showed **Saponite-Na barely takes up Na**. A search
  for Na-richer reverse-weathering clays found that the only Al-free Na sinks are **Na-carbonates**,
  which were then added to `make_database.py`.
- **Albitization** — Na into precipitating albite. This became the main line of attack.

### 6.2 What PHREEQC-HT actually did

Isolated HT tests (473 K, water/rock = 2.0, self-consistent crust) gave a clean result at
**seawater-like Na (~470 mM)**: Mg stripped into Clinochlore (b_eq[Mg] 53 → 0.02 mM), Ca released
from Diopside/Anorthite (b_eq[Ca] 10 → 74 mM), **dCa/−dMg ≈ 1.1–1.2**, Na near-neutral. Exactly the
intended exchange.

**But at low ocean Na it collapsed.** At the hothouse state (Na = 27 mM): F[Ca] ≈ 0 or negative,
and HT instead dumped **F[Na] = +10 Tmol/yr**. The controlled test (strip Albite from the crust)
proved the mechanism — with Albite removed, `b_eq[Ca]` pins at ~63 mM and `F[Ca] ≈ +5.7 Tmol/yr`
**regardless of ocean Na**.

**Albite is the switch.** Low ocean Na drives Albite dissolution, which floods Na *and*
simultaneously collapses Ca release (63 → 7 mM). Same failure shape as epidote earlier: an
Al-bearing phase suppressing Ca.

The literature check supported the mechanism: albitization is Na-in/Ca-out, and the Ca/Na
partitioning is set by T- and P-dependent plagioclase equilibrium (Berndt & Seyfried plagioclase
exchange experiments; Mottl & Holland; Seyfried & Bischoff seawater–basalt experiments, assemblage
albite–chlorite–epidote–actinolite). **Caveat:** the *total* Ca collapse at low Na is exaggerated,
because real calcic plagioclase and diopside still dissolve congruently and release some Ca. The
model's collapse comes from treating normative Albite as a freely-equilibrating two-way phase.

### 6.3 Fixes tried, in order, and why each failed

| Attempt | Result |
|---|---|
| **600 K** instead of 473 K | **Worse.** b_eq[Ca] collapses at *every* Na (74 → 8.8 mM at Na=470), F[Ca] goes negative, Na flooding intensifies. Temperature is the wrong lever. |
| **`dissolve_only` on all primaries** | Ca robust (b_eq[Ca] ≈ 228 mM, Na-independent) and albitization a clean one-way Na sink — but **it broke the thermostat**. F[Ca] became flat vs ocean Mg (dCa/−dMg swung 210 → 1.9), i.e. a climate-blind constant +31 Tmol/yr alkalinity pump. Sweep result: **bistable** — f_HT ≤ 0.005 hothouses, f_HT ≥ 0.02 collapses to the CO₂ floor. Also drove Na and Cl to ~0. |
| **Non-`dissolve_only` + Albite precipitate-only** | Chemically the best PHREEQC result: dCa/−dMg ≈ 0.8–1.1, F[Ca] ∝ ocean Mg, robust at low Na, moderate b_eq[Ca] 14–110 mM at pH 6.1–6.6. **But the full sweep still hothoused everything** — Ca stuck at 0.1–0.3 mM while Alk (145–203 mM) and C (171–227 mM) piled up as *Mg*-alkalinity. Na still drained to 0 (Albite-precip-only removed the only Na source on a landless world); Cl collapsed to ~4 mM. |

The user's diagnosis at the pivotal moment was the correct one:
*"Surely the Ca supply from the HT system is controlled by the Mg supply, which is controlled by the
weathering?"* — a true exchange caps Ca release stoichiometrically at the Mg removed, and since ocean
Mg is set by climate-sensitive weathering, **that coupling is what carries the thermostat**. The
`dissolve_only` variant violated it; that is precisely why it was bistable.

The honest read that emerged: the exchange is **alkalinity-neutral** (Ca in ≈ Mg out), so it adds no
drawdown capacity — it only unlocks the Mg fraction of weathering alkalinity for calcite burial, and
its throughput is capped by `f_HT·J`. A pure water world's seafloor exchange is a **weak thermostat**.

### 6.4 The Aug 3 decision: switch back to the parameterization

User: *"We could switch back to the Ca-Mg parameterization? The full PHREEQC machinery is causing
endless problems."* The June version (`12d082a`, 2026-06-14) was recovered from git and had two
components the newer code had lost:

1. **`F_ht_exchange`** — `KD_MG_HT · [Mg] · J`, strict 1:1 Mg→Ca, alkalinity-neutral and
   charge-consistent by construction, Ca ∝ [Mg] ∝ weathering.
2. **`F_na_rw`** — a parameterized Na sink scaling with `J_total` **not land**, which is why the
   June model gave water worlds a Na sink at all.

Implemented: dropped the `get_weathering_flux(high_temperature=True)` PHREEQC call and the J_LT/J_HT
split (all hydrothermal flux now goes through the LT path with `J = J_total`), added
`F_ht_exchange`, and made the Na sink **always-on** rather than land-gated. `chemistry.py`'s HT
`get_b_eq` edit was reverted to a plain equilibrium (kept clean, now unused).

**Immediate result:** still Ca-limited. `WW o=1.0` hothoused with Ca = 0.25 mM, Mg = 100 mM.
Quantified: outgassing at o=1.0 is 7.5 Tmol/yr C, the exchange supplies only **2.4 Tmol/yr** at
Mg = 100 mM; balancing 7.5 would need ocean Mg ≈ **313 mM**. Diagnosis: `KD_MG_HT` was calibrated in
June *with the PHREEQC HT dissolution also present*, so it is now ~3× too weak carrying the HT Ca
load alone. **Recalibrating `kd_mg_ht` upward is the outstanding action from this arc.**

### 6.5 Cl

Cl enters from outgassing at `cl_outgassing_ratio` (default 0.02) of the C flux, with HCl outgassing
correctly debited from alkalinity, and leaves by subduction ∝ crust production. **Cl only reaches
~165 mM against Earth's ~550 mM.** This matters more than it looks (§7).

---

## 7. Alkalinity charge-consistency — the key correctness fix

### The bug: phantom Al alkalinity

The weathering flux formula mixes a **kinetic** `k` with an **equilibrium** `b_eq` *per element*.
For alkalinity vs Al those two disagree badly:

- `k[Alk]` counts the H⁺ consumed to release Al³⁺ from Anorthite/Albite — which dissolve fast.
- `b_eq[Al]` is a trace (~0.0005 mM) because at equilibrium those phases are nearly inert.

So the flux delivered Al's **alkalinity** while pinning Al's **mass flux** at ~0: alkalinity with no
matching ocean cation. Measured inflation was **~3× the real cation charge**, growing with
temperature (2.37× at T=335 vs 1.08× at T=300).

The user asked the right question — *shouldn't the secondary Al clay precipitation remove that
alkalinity, since PHREEQC is self-consistent?* Tested directly: Kaolinite **does** fire (SI = +1.05)
and does remove some alkalinity, but only closes ~20% of the leak, because **the Al never reaches
the pore fluid** to precipitate. PHREEQC is self-consistent *within* `get_b_eq` and *within*
`get_precipitation`; the kinetic-flux layer between them is not.

### The fix (user's formulation)

Rather than dropping alkalinity from the state vector (my first proposal — closer to the textbook
"track conservative ions + DIC, derive Alk"), the user proposed keeping it tracked but **defining
each function's alkalinity output as the charge its ion fluxes actually deliver**. This stays
aligned with Coogan & Dosso (2022), who also track alkalinity.

An `ION_CHARGE` vector was added to `chemistry.py` (Ca +2, Mg +2, Na +1, Fe +2, Al +3, Cl −1, S −2)
and applied as `flux[alk_idx] = ION_CHARGE · flux`, computed **after** secondary precipitation so the
ion fluxes are final. Applied in:

- `get_weathering_flux` (covers LT and HT at once — both leaked)
- `get_continental_weathering_flux` (dropping the hard-coded balanced-Na assumption)
- `F_cl_subduct[alk] = −F_cl_subduct[cl]`
- `F_na_rw[alk] = F_na_rw[na]`

**Fe and Al are included**, on the user's call: under anoxic conditions Fe²⁺ is soluble and
genuinely conservative, so excluding it would create the opposite imbalance. Redox gating is handled
upstream (Fayalite is excluded when fO₂ > 0), and deriving from the *actual* `F[Al]` rather than
`k[Al]` is what kills the phantom, so including Al is exact and harmless.

### Outcome

- Per-process charge audit clean: seafloor LT, seafloor HT, Cl subduction, continental and Na terms
  all at leak ≈ 0 (were +7.5e-17, +7.3e-18, −2.8e-20, …). Only a ~1e-22 precipitation residual from
  the `tanh(SI)` smoothing.
- **Na fixed 1220 → 475 mM** at land = 0.3 — essentially Earth's 480.
- **First temperate, stable Earth from this scheme:** `timeout` at 2 Gyr, **T = 297 K,
  pCO₂ ≈ 1200 ppm, pH ≈ 9.15**, Ca = 0.26, Mg = 12.7, Na = 475, Cl = 165.

### Two known residual inconsistencies

1. **The clamps in `dY_dt` break tracked-Alk = ion charge.** After summing fluxes,
   `F_net[b_ocean<=0] = max(F_net, 0)` (positivity guard) and `F_net[so4_idx] = 0` (S pin) modify
   individual *ion* fluxes without adjusting alkalinity. On Earth, Ca hovers near 0 and gets
   clamped, so per-flux consistency doesn't survive into the ocean state.
2. **Cl is too low (165 vs ~550 mM)** — and this must be fixed *first*. Real seawater alkalinity
   (2.3 meq) is a tiny residual of near-cancelling large ions: Na (480) balanced almost entirely by
   Cl (550). With Cl at 165, Na is left unbalanced and the true charge balance is a large 336 meq
   (tracked Alk = 64). **Enforcing exact charge consistency before calibrating Cl would collapse
   Earth's CO₂.**

---

## 8. The Damköhler diagnostic

The user noticed the plotted Da was built from the **primary-dissolution alkalinity**
(`Da_primary[alk_idx]`) — i.e. the Damköhler number of the alkalinity that is then **thrown away**
and overwritten by the charge-derived value.

Consequence: `b_eq[alk]` carries the phantom Al alkalinity, so the plotted Da is **under-reported by
~2× right in the thermodynamic regime**, and the bias **grows with T** — it distorts the shape, not
just the level. The `Da = 1` contour was misplaced toward higher instellation.

The user proposed `Da = (Σk) · A / (J · Σb)` over the delivered cations. Verified: the unweighted sum
is within ~2% of the charge-weighted version at the kinetic end and ~9% at the thermo end, and
**both put the `Da = 1` crossing at T ≈ 315 K where the rigorous effective Da crosses it** — while
the old diagnostic said 0.92 at T = 315, still calling the system kinetic when it was already past
transition. The two things that actually matter are **sum the delivered cations** and **exclude Al**.

Two bugs found in the user's first implementation, both fixed:

1. **`Da` was `nan` everywhere.** `np.dot(ION_CHARGE, k_primary)` ran on the `inf`-substituted
   `k_primary` (zero rates → `inf`); C has `k = inf` and `ION_CHARGE = 0`, so `0*inf = nan` poisoned
   the dot product.
2. **Al was back in**, contributing 47% of `b_alk` — reintroducing exactly the contamination the
   exercise was meant to remove.

Current form (`weathering.py:146–152`):

```python
q_alk = np.where(ION_CHARGE > 0, ION_CHARGE, 0.0)   # drop anions
q_alk[al_idx] = 0.0                                 # drop Al (precipitates as clay)
k_finite = np.where(k_nonzero, k_primary, 0.0)      # avoid 0*inf -> nan
Da_alk = (np.dot(q_alk, k_finite) * A_reactive) / (J * np.dot(q_alk, b_eq_primary))
```

There isn't really one scalar Da: at T = 335 the per-ion values diverge widely (Mg 7.5, Ca 2.0,
Na 0.32) — Na stays kinetic while Mg is deep in the thermo limit. The carbon-relevant number is a
charge-weighted blend dominated by Ca/Mg.

### Why "P_CO2 doesn't rise when Da > 1" was a false alarm

The trend *is* in the model. Two things masked it:

1. **Da is temperature-dominated, not supply-dominated.** Correlation of `log Da` with inputs:
   instellation **+0.825**, crust production −0.410, outgassing +0.085. Along the dominant axis
   (instellation) you're raising T, the temperature feedback is healthy and stabilizing, so pCO₂
   *falls* as Da rises → net Da–pCO₂ correlation −0.20. Along the **supply** axis (crust production
   → hydrothermal J) the expected rise is exactly there: c = 10 → 0.01 gives Da 9.3e-4 → 2.14 and
   P = 0.52 → 5.0 bar.
2. **The runs that would populate the rising branch terminated instead.** When seafloor weathering
   is supply-limited its flux caps at `F_max ≈ J·b_eq(T)`; if outgassing exceeds that cap there is
   **no steady state**, pCO₂ rises without bound and hit the old 5-bar `acid_ocean` event. The 50
   runs piled at exactly P = 5.00 *were* the rising branch, censored.

Also measured: the direct CO₂-restoring feedback is **nearly zero** — `d log F_alk / d log P_CO2` ≈
0.00–0.13, against a classic Walker exponent of ~0.3. The equilibrium alkalinity basalt buffers to
is nearly independent of pCO₂ (1.33e-3 → 1.47e-3 while pCO₂ goes 0.01 → 3 bar). So once
supply-limited, the rise is near-vertical, not gentle. Root cause of the zero-feedback finding:
**all 160 runs in that sweep had `land_fraction = 0.0`**, so the only silicate sink was seafloor
weathering; continental weathering (the kinetic Walker-type sink) is gated on `land_fraction > 0`.

---

## 9. Termination events: six → two

### 9.1 The redesign

The user's instruction was decisive: *"I ultimately think there are too many events. Converged never
fires and a lot of the other ones are just triggering before the final state. All the events need to
do is end runs when the planet is outside the valid P_CO2 or T range for the climate solver."*

Confirmed: `converged` fired **1 time in 98** runs, and `snowball`/`hothouse`/`frozen` were outcome
*classifiers* dressed as stopping rules.

The key realization: **the climate solver already knows its own domain.** `get_T_surface_analytic`
returns exactly `180.0` or `400.0` when its scan finds no root in [180, 390] — a self-reported "no
equilibrium here". So `planet.py` had been re-deriving 260 K / 350 K / 5 bar thresholds that already
live in `climate/`. Both climate models top out at exactly **10 bar** CO₂ (the analytic OLR fit is
explicitly only fitted below 10 bar).

`event_snowball`, `event_hothouse`, `event_co2_ceiling`, `event_co2_floor` and a short-lived
`event_frozen` were replaced by a single **domain box** returning the smallest *normalised* margin,
with `direction = −1`:

```python
T_LO, T_HI = 181.0, 389.0    # K, one degree inside the analytic scan bracket
# margins: log10 decades for pCO2, (T - T_LO)/100 for temperature
```

Normalising matters: `min()` over raw Pa and K is dominated by whichever has bigger units, so log for
pCO₂ and /100 for T puts all walls at O(1) for the root finder. Stopping *just inside* the rails
(181/389) rather than at them keeps the event off the solver's flat region.

Which wall was hit is recovered afterwards from the final state as `domain_wall` in the JSON
(`cold` / `hot` / `co2_high`), so **adding a wall never means adding an event**. `event_converged`
survives purely as a compute optimisation, stripped of the T-gate and solute-free special case that
were why it never fired.

Prior to this, the 5-bar `acid_ocean` event was found to be "a numerical guard wearing a planetary
costume" — it fired on pCO₂ > 5 bar but was reported alongside `snowball`/`hothouse` and labelled
"Acid Ocean" even though nothing in the trigger looked at pH. Worse, **the 5-bar cut landed in the
middle of the maximum-greenhouse peak** (4 bar at S=0.5, 7 bar at S=0.7). With the ceiling moved to
10 bar, `ww_S0.5` stopped dying "acid" at a habitable 302 K and instead ran on to terminate as a
genuine **snowball at 141 Myr**.

Old `acid_ocean` JSONs still render via a legacy label in `plot_results.py`.

### 9.2 The domain-guard bug and the `co2_low` wall removal

This was a three-step sequence worth remembering as a unit:

1. **`min_time` was carried over onto the domain event (2 Myr) — wrong.** `min_time` exists so a
   *planetary verdict* isn't declared during startup; leaving the model's validity box is a fact
   about the state, not a verdict. The cost was extreme: `s_1.2_out_10_crust_1` exceeded 389 K at
   **t = 16.6 kyr** but couldn't terminate until 2 Myr, so it chattered across a **26.5 K
   discontinuity** at pCO₂ = 1.0829 bar **2,618 times** over 13,018 steps — then mislabelled itself
   `co2_high` at 10 bar. Fixing the guard to `min_time_domain = 1e4 * YR`: **1884 s → 1.1 s
   (~1700×)**, terminating honestly as `hot` at 16.3 kyr. Both controls bit-identical.
   A **short** guard rather than zero is required: `solve_ivp` fires on a sign *change*, so a planet
   already outside the box at t = 0 (T > 389 K at the initial pCO₂, from S ≈ 1.3 up) would have a
   margin that never changes sign and would never terminate. That is the same trap that had
   disabled `event_snowball`.
2. **The short guard then killed 13 healthy runs** as spurious `co2_low`. Cause: 10 kyr is far
   shorter than the blank-ocean CO₂ drawdown transient. Measured on the runs where it was allowed
   to continue: the pCO₂ dip **bottoms out at 105–210 kyr** (median worst margin −1.02 decades, two
   runs hitting exactly 0) and **recovers above 0.1 Pa only at 133 kyr – 1.47 Myr**. So even a 1 Myr
   guard would still have killed several.
3. **The user reframed it correctly** — the low-pCO₂ state isn't itself a problem, and if the
   numerics there were fine the wall wouldn't be needed. Investigation found: the oscillation is
   real (`Y[0]` reaches −0.097 Pa, negative on 1.67% of trajectory calls, `dP/dt` flipping sign 19
   times) **but costs essentially nothing** — `atol[0] = 1.0 Pa` is ten times larger than the whole
   excursion, so the region sits below the solver's noise floor. It is actually *cheaper* per step
   there (10.9 derivative calls per accepted step inside the sub-1 Pa window vs 16.5 outside); total
   cost across two sweeps was **0.07 h**. And crucially, the clamp `max(P_CO2, 1.0)` was **already
   applied consistently** at both call sites, so T never sees the region where the OLR fit
   collapses. **The `co2_low` wall was protecting against nothing.** It was removed, the magic `1.0`
   was promoted to `P_CO2_CLIMATE_FLOOR` with the real justification recorded (both OLR fits are
   polynomials in `log10(pCO2)`; below ~1e−2 Pa they extrapolate outside their range and collapse to
   353–359 K for *every* instellation, so 1 Pa is two decades of margin), and a false comment
   claiming a low-pCO₂ planet is "unambiguously frozen" was corrected — T at 1 Pa actually spans
   193–365 K, and one run had been labelled "CO₂ depleted" at 365 K.

**Cost of the removal: +1.87 h on a matched set** (my estimate of +0.25 h, revised to +1.10 h, was
still ~40% low both times — recovered runs are well above the median timeout cost). It bought back
**22 wrongly-killed runs** and removed a whole class of artifact. Habitable fraction went 26% → 40%.

### 9.3 Initial conditions

`Y0[0] = initial_pco2 = 1000` Pa sits above a chemically **blank** ocean (`Y0[2:-1] = 0`). Because
there's no air–sea exchange term (§3), the 1000 Pa is **fictitious carbon that simply evaporates**
over a few `TAU_ATM` while ocean carbon independently builds from outgassing. The pCO₂ U-shape has
**no physical content** — it is an inconsistent initial condition unwinding, not a young planet
degassing.

Two options were examined:

- **Equilibrated `b0`** — choose `b0` so the ocean already implies `initial_pco2`, making
  `dYdt[0] ≈ 0` at t = 0. Verified well-posed: pCO₂ is smooth and monotonic in `b[C]`, and `brentq`
  on `log10(b[C])` over [1e−7, 1e−1] converges in a handful of PHREEQC calls
  (`b[C] = 4.6552e−04 mol/kgw → pCO₂ = 999.96 Pa, pH 4.88`).
- **Start at pCO₂ = 1 Pa** (user's suggestion) — safe and slightly better, but **does not save
  time**, because spin-up has two parts and this only removes the cheap one. The atmosphere crash is
  a smooth exponential LSODA walks in a few steps; the expensive part is **ocean carbon filling from
  outgassing over ~1 Myr**, which is real physics and unaffected by where the atmosphere starts.
  Steps in the first 1.5 Myr barely moved (three configs went slightly *up*).

  **Valuable side result:** all five configs gave the **identical outcome** — same termination, same
  wall, ΔT ≤ 0.47 K, same final pCO₂ to 4 s.f. So the attractor does not depend on initial pCO₂; the
  model is well-posed in that respect. (Largest shift was Na at 1.83e−07 → 2.35e−07 mol/kgw — a
  trace that hasn't converged in 2 Gyr, Na having the longest residence time in the system.) Bonus:
  at 1 Pa, S = 1.3–1.5 start *in* domain where at 1000 Pa they did not.

Neither is implemented; `initial_pco2 = 1000` is still the default.

---

## 10. Performance work

91% of wall time is inside PHREEQC, and **81% of that is the Jacobian** (12 columns × 2 central
differences = 24 `dY_dt` calls, each costing ~5 PHREEQC solves; Jacobian:trajectory call ratio
4.3:1).

### What worked

| Change | Effect | Where |
|---|---|---|
| **PHREEQC `KNOBS` block** | `chemistry.py` had **no KNOBS at all** — stock solver defaults. Every failure was `Maximum iterations exceeded`, never a thermodynamic error. Replaying 200 captured failing inputs: baseline 144/200 → `-iterations 1000 -diagonal_scale -step_size 5 -pe_step_size 5` **200/200**. End-to-end S=0.7: **2785 → 393 errors**. **Answer-preserving** — worst relative difference 1.8e-9 across pCO₂, pH, all ten flux components and Calcite SI. Also makes calls *faster* (10.1 → 4.9 ms) by avoiding PHREEQC's retry cascade. | `chemistry._knobs_block()` |
| **Retry once on non-convergence** | The single biggest win. Escalating KNOBS further is **exhausted** (400 already-failing inputs: more iterations, smaller steps, no diagonal-scale all recover the same 235/400; only relaxing tolerance to 1e-6 helped, at 250/400, which trades accuracy). But replaying with **identical settings recovers 58.8%** — the failures are largely **warm-start path artifacts**. A second identical attempt: **78.5%**. End-to-end: **fallbacks −91–96%, runs 3.5–4× faster, final state bit-identical.** | `chemistry.py:406–418` |
| **Continuous fallback** | On `ChemistryError`, return `self._dYdt_last_good` instead of zeroing every sink (which was a discontinuity the solver had to crawl across). Combined with the retry: **5.1–5.9×** on affected runs, fallbacks 115 → 3. Resets per `time_evolve` so a reused `Planet` can't carry stale state. | `planet.py:341`, `:388` |
| **Short domain guard** | See §9.2 — up to ~2400× on the worst runs. | `planet.py:477` |

State-reset-then-retry recovers marginally more (80.5% vs 78.5%) for added complexity; **resetting
before every call is actively harmful** (62.2%).

### What did *not* work (measured, not assumed)

- **Forward instead of central differences for the Jacobian: 7.3× SLOWER** (36.7 s → 268 s). Halving
  the Jacobian cost degrades it enough that LSODA needs 6× more steps (341 → 1984). The answer is
  unchanged, so it's purely stiffness. The central-difference Jacobian is earning its price.
- **Raising `max_step`:** not binding. Only 10.6% of steps sit at the 1e7 yr cap; median dt is
  5e5 yr.
- **Widening the OLR blend (0.1 → 0.5 bar):** the 78% correlation between small steps and
  pCO₂ ∈ [0.9, 1.1] bar was **real but not causal**. Step counts moved ±5%, wall time was equal or
  worse. The real cause was the pre-`min_time` rail chattering (§9.2), and fixing *that* removed the
  hotspot entirely (78.2% → 4.2%, i.e. zero enrichment).

  For the record the blend kink is real: `analytic.OLR` blends two polynomial fits with
  `tanh((pCO2 − 1.0)/0.1)`; they agree in **value** at 1 bar to 0.000 W/m² but their **slopes**
  differ by 12–32 W/m² per decade, making T(pCO₂) non-monotonic across 0.90–1.20 bar (a spurious
  local max *and* min, ~0.03 K). The climate feedback gain passes through zero and reverses sign
  twice, so the thermostat is marginally stable there.

### Costed but NOT implemented

- **Trim the saturation-index list: ~8–11%.** PHREEQC computes an SI for **84 phases** on every
  call; only ~17 are ever read. Verified end-to-end at **1.08× with bit-identical results**.
  One-line change at `chemistry.py:363`/`262-263`. **The trap:** `get_ocean_state` computes pCO₂ as
  `EARTH_ATM * 10 ** output['si_CO2(g)']`, so the kept list must retain **all `(g)` phases** —
  dropping them breaks every run immediately. Everything else is safe because the SI dict is already
  filtered downstream by `precipitating_minerals`. *(Still not applied: `available_mineral_string`
  is `' '.join(minerals)` — all 84.)*
- **Skip two dead Jacobian columns: ~13%, exact.** `S`/SO₄ is identically zero across all sweep runs
  checked *and* pinned (`F_net[so4_idx] = 0`), so its Newton update is always zero. `r_avg`'s column
  is analytic.

---

## 11. Current code state — verified

Present in the working copy (uncommitted, on top of `f46de52`):

- ✅ `ION_CHARGE` + charge-derived alkalinity in `get_weathering_flux`,
  `get_continental_weathering_flux`, `F_cl_subduct`, `F_na_rw`
- ✅ Corrected Da diagnostic (`q_alk`, Al excluded, finite `k`)
- ✅ `F_ht_exchange` parameterization; HT PHREEQC call dropped; `f_HT` vestigial (default 0.0)
- ✅ Na sink always-on (`F_na_rw`, `J_total`-scaled)
- ✅ `_knobs_block()` and retry-once in `solve_solution`
- ✅ Continuous fallback via `_dYdt_last_good`
- ✅ `event_domain` (3 walls: cold/hot/co2_high) + `event_converged` only; `domain_wall` in JSON
- ✅ `min_time_domain = 1e4 * YR`, separate from `min_time`
- ✅ `P_CO2_CLIMATE_FLOOR = 1.0` Pa at both call sites; `co2_low` wall and `P_CO2_LO` removed
- ❌ **REVERSED (§21.3).** Reverse-weathering minerals are **no longer** in the pore precipitation
  list: `pore_precipitating_minerals = clay_minerals`. They were double-counted (also applied to the
  ocean via `rw_ocean_precipitating_minerals`) and, acting on a pore fluid loaded to ~100 mM Mg by
  dissolution, Sepiolite(d) cancelled the entire primary alkalinity flux and flipped `F_net` negative
  over a broad region of (T, pCO₂). Added in fast_10, removed 2026-08-19.
- ✅ **`b_pore = b_input + F_primary/J`** (§20.3) — pore fluid derived from the flux, not a second
  interpolation. The correctness fix that closed the iron charge leak and killed the Ca runaway.
- ✅ **`kd_mg_ht = 0.07`** default (§18); **`WATER_ROCK_RATIO_LT = 3`** — low, on the feedback
  plateau (§20.2). Briefly 600 during fast_13; now 3.
- ✅ **`lt_equilibrium_buffer_minerals = []`** — Kaolinite step-1 buffer wired but OFF (§20.2).
- ✅ **`chemistry_void` termination** + `fabricated_fraction`/`chemistry_ok` fields; buffer
  self-disables without a w/r (§20.4 guards).
- ✅ `crust_composition.py` pipeline; `plot_results.py` rewritten for the new model
- ❌ SI-list trim — **not applied**
- ❌ Dead Jacobian columns — **not applied**
- ❌ Equilibrated `b0` / low `initial_pco2` — **not applied**
- ⚠️ `convergence_threshold` default is **0.05**; `parameter_sweep.py` passes nothing so sweeps use
  it. Note (§20.4) it barely matters — the metric is set by slow Ca/Cl, and most timeout runs are
  climate-settled regardless of the threshold.
- ~~⚠️ `experiments/calibrate_earth.py` still passes `f_bio=`~~ — **stale, see §22.0.** It has been
  rewritten (charge-balanced seawater seed, `least_squares` on `log(K_na, alpha, KD_mg)`, Phase 2
  removed) and runs clean.
- ✅ **CO₂ ceiling is the instellation-dependent maximum greenhouse**, not a flat 10 bar
  (`maximum_greenhouse` in `climate/analytic.py`, called once in `time_evolve`). This closes §15
  item 1; no sweep has been run with it yet.

---

## 12. The Mg problem — the live chemistry failure

The remaining PHREEQC failures are **not numerical**; the ocean drifts somewhere PHREEQC cannot
solve. **100% of failures are one call**: the pore-space clay step inside `get_weathering_flux`.
`get_ocean_state` and the weathering `get_b_eq` never fail.

A representative failing state (58 °C, 296 bar, pH 10.04) is **10³⁴⁴× supersaturated in serpentine**:

| phase | SI | allowed to precipitate? |
|---|---|---|
| Antigorite | **+343.7** | no |
| Anthophyllite | +45.2 | no |
| Saponite-Mg | +26.5 | no |
| Talc | +25.0 | no |
| Chrysotile | +19.4 | no |
| Kaolinite | −0.31 | yes (but *under*saturated) |
| Goethite | +11.5 | yes (but Fe-limited at 3e−5 mol/kgw) |

Every Mg-silicate that would relieve the supersaturation is absent from the precipitating list, so
PHREEQC is asked to equilibrate two trace-limited phases inside a solution held wildly out of
equilibrium in a dozen others. (A charge-imbalance hypothesis was tested and **disproved** —
correlation with fallbacks r = 0.029, and the worst-imbalance runs have zero fallbacks.)

**Adding the reverse-weathering minerals to the pore list helped but could not fix it**, for two
independent reasons:

1. **The pore pathway is throttled by exactly the thing causing the problem.** `weathering.py:133`
   is `flux += J * d_b_secondary`, and `J_total ∝ crust_production_rate`. At `crust = 0.01` the
   whole pore contribution is scaled by **1%** — in precisely the runs that need it. The same
   applies to the HT Mg→Ca exchange, which is also ∝ `J_total`.
2. **Sepiolite(d) is starved of silicon.** Amorphous silica precipitates on `tau_prec` = 100 kyr and
   holds `SiO2(am)` at SI = 0.00 (exactly saturation); the reverse-weathering list runs on
   `tau_rw` = 5 Myr, **50× slower**, so the silica is gone before it gets a look. Sepiolite(d)
   doesn't precipitate at all.

> **Reconciling this with §21.3, which measured Sepiolite removing 34 Tmol/yr.** No contradiction —
> they are different calls. Point 2 above concerns the **ocean** RW step (`tau_rw` = 5 Myr, kinetically
> throttled, silica already consumed by `SiO2(am)` on the 50× faster `tau_prec`). §21.3 concerns the
> **pore** step inside `get_weathering_flux`, which runs at *instantaneous* equilibrium on a fluid that
> primary dissolution has just loaded with both Mg and Si. Same mineral, inert in one reservoir and
> dominant in the other. That asymmetry is exactly why the RW minerals belonged in the ocean list and
> not the pore list, and it is what made the defect hard to see: §12 correctly established Sepiolite
> was doing nothing, in the reservoir being looked at.

Measured Mg budget at the final state of `s_0.9_out_0.1_crust_0.01`: seafloor LT source **+0.342
Tmol/yr**, HT exchange −0.049, ocean precipitation **0**, reverse weathering −2.7e−06 → **net
+0.293**, residence time ~1056 Myr. Mg simply accumulates.

**The user's kinetic objection is correct and stands:** dolomite and magnesite are thermodynamically
favoured and kinetically absent (modern seawater is supersaturated in dolomite and it still doesn't
form abiotically; magnesite is inhibited by Mg²⁺'s dehydration energy). Picking sinks off an SI
ranking is the trap. The model's exclusion list is *already* doing implicit kinetic filtering —
correctly refusing Antigorite at 10⁶⁷ — and the honest problem is that **after filtering there is no
accessible Mg sink left at all**. The kinetically plausible low-temperature hydrated Mg carbonates
(nesquehonite −0.40, artinite −0.91) and brucite (−3.72) are **undersaturated**, so they wouldn't
help even if added. One consideration left open: these states are at **T_seafloor = 72.7 °C,
T_pore = 81.7 °C**, not 25 °C, which does weaken the magnesite kinetic argument.

---

## 13. The June 2026 seminar baseline — the last consolidated results

**`Ocean_Chemistry_Seminar-2.pdf`** (repo root) — *"The Ocean Chemistry of Rocky Waterworlds: Do
ocean worlds have stable climates?"*, Tanna & Shorttle, 17 June 2026, 42 slides.

**This is the last time the model produced a complete, self-consistent set of results.** Everything
in §6–§12 happened *after* it, and the current `NaCl-chemistry` branch has not yet reproduced this
quality of output. Treat these figures as the target behaviour, not as history.

The code behind it is the 2026-06-10 → 06-15 commit cluster — `4059acc` "Improved plotting"
(06-15, two days before the talk) and `12d082a` (06-14), the same "June version" the Aug 3 session
kept referring back to. Na–Cl chemistry was already in (`925c514`, `82c435a`, `489dc7f`, 06-10), so
these plots include salinity. Critically, this version still had **Wollastonite** as the Ca source
and both `F_ht_exchange` *and* the PHREEQC HT dissolution call (§6.4).

### Plot grammar (verified against `4059acc:experiments/plot_results.py`)

Every results figure shares one structure. Two encodings are easy to misread:

- **x-axis is always instellation** S/S₀ ∈ [0.3, 1.4]. Rows are T and P_CO₂, or pH and salinity.
  Columns are crust production (0.1×, 1×, 10×); colour is outgassing (0.1×–10×, log scale).
- **Line style is the weathering regime:** solid = Da < 1 (kinetic), dashed = Da ≥ 1
  (thermodynamic), open circles mark the **Da = 1 crossing**.
- **Dotted = "T_sf at floor (274 K)".** More meaningful than it looks:
  `T_seafloor = max(1.02·T_surface − 16.7, 274)`, so the floor engages whenever
  **T_surface ≤ 285 K**. Dotted means the seafloor has stopped tracking the surface.
- **Grey dash-dot = "Equilibrium temperature"** is
  `equilbrium_temperature(S, albedo=0.3, greenhouse=0.5) = ((1−A)·S·S₀/(4σg))^0.25` — a **fixed-
  greenhouse, no-carbon-cycle reference** scaling as S^0.25. Values: 224 K (S=0.3), 254.5 (0.5),
  276.9 (0.7), **302.7 (1.0)**, 316.8 (1.2), 329.3 (1.4).
- Shaded bands are the outcome thresholds: blue 235–260 K (`T_SNOWBALL = 260`), red 340–360 K
  (`T_RUNAWAY = 360`).

### The central result (slide 26)

The thermostat, and the grey reference curve is what makes it legible. Between S = 0.5 and 1.0 the
model sits on a near-**isothermal plateau at ~300–305 K**, rising ~8 K. The fixed-greenhouse
reference rises **48.2 K** across the same interval. **The carbon cycle suppresses the temperature
response by a factor of ~6**, paid for with roughly **three decades of CO₂ drawdown** (~3 bar at
S=0.5 → ~2×10⁻³ bar at S=1.05). Textbook Walker feedback.

At S ≈ 1.05 the line goes dashed at the open circle and the behaviour **inverts**: P_CO₂ starts
*rising* with instellation while T shoots vertically into the red band. This is the positive
feedback of slide 11 — in the thermodynamic limit the flux is `J·b_eq(T)` and `b_eq` *falls* with
warming, so hotter → less drawdown → hotter. The thermostat does not weaken gradually; **it reverses
sign.**

**The inner edge of the habitable zone here is set by the Da = 1 transition, not by the classical
runaway greenhouse.** That is the paper's most distinctive claim and it appears as a single crossing
point.

### What each parameter does (slides 27–28)

| Parameter | Effect |
|---|---|
| **Outgassing** ↑ | Sets plateau *height*: ~295 K (0.3×), ~305 K (1×), ~340 K (3×); 0.1× sits in/near the snowball band. T and CO₂ shift up together. Also moves the **Da = 1 crossing to *lower* instellation** (~1.10 at 0.1× → ~1.03 at 3×) — more carbon to dispose of reaches the thermodynamic limit sooner, so high-outgassing worlds have a *closer-in* inner edge. |
| **Crust production** ↑ | The supply term, and the mirror image: 0.1× → 10× drops the plateau ~320 K → ~300 K and the CO₂ minimum by ~1 decade. More new seafloor + more hydrothermal flux = more weathering capacity = lower steady-state CO₂ for the same carbon input. |

### The parameters collapse to one (slide 30)

Plotting outcome against **outgassing / crust production ratio** vs instellation, the habitable band
is a clean vertical stripe at **ratio ≈ 0.3–10**, spanning S = 0.4–1.2. Left of it (tectonically
active, volatile-poor) → snowball, weathering capacity swamps carbon supply. Right (inactive,
volatile-rich) → CO₂ piles up. Top (S ≳ 1.2–1.3) → hothouse regardless.

This is what licenses the conclusion that **stagnant-lid planets can still be habitable**: what
matters is the *ratio*, so 0.1× crust production is fine provided outgassing scales down with it.
The two-parameter tectonic problem is really one-parameter.

### Chemistry (slides 31–34, 36, 38)

pH and salinity both **track P_CO₂** rather than anything else, which makes every panel readable
once noticed:

- **pH** rises with instellation (~5.7 → 7.5), falls with outgassing, rises with crust production —
  in every case the mirror of CO₂. Carbonic acidification dominates the pH budget.
- **Salinity** falls with instellation (~25 g/kg at S=0.45 → ~6 g/kg at S=1.05), rises with
  outgassing, falls with crust production. The saltiest oceans are the cold, high-CO₂,
  low-tectonic-activity ones — an acidic ocean holds more dissolved solute.

That yields the coupling on slide 41, which is a genuine result: **the planets nearest the outer edge
are also the saltiest**, so their freezing point is most depressed and they are least vulnerable to
ice-albedo runaway. A second, independent negative feedback on the cold side, arguing the outer edge
extends further than a pure-water treatment gives.

### Second-order parameters

- **Ocean depth (slides 35–36):** 300 m → 30,000 m, two decades, and the T and P_CO₂ tracks
  **essentially collapse onto each other** in the kinetic regime — depth only separates them past
  Da = 1. Salinity varies strongly and inversely (shallower = saltier), close to pure dilution.
- **Crust composition (slides 37–38):** 44–51% SiO₂ barely moves climate (most mafic runs ~15 K
  warmer on the plateau) but moves chemistry a lot — **more mafic → higher pH and much higher
  salinity** (~100 g/kg at 44% SiO₂ vs ~20 at 51%). Mafic crust is richer in divalent cations and
  poorer in silica, so dissolution delivers more cation charge per unit rock.

Hence the conclusion that ocean depth and crust composition "generally do not affect the climate
state" while both strongly affect ocean chemistry.

### Earth calibration (slides 24–25)

Slide 24 puts an Earth marker on the T/P_CO₂/pH/salinity-vs-instellation curves at S = 1 (≈287 K,
pH ≈ 8.1, ≈33 g/kg) with the model curve passing close to it.

Slide 25 is the calibration argument, as **% model difference** split into biotically- and
abiotically-controlled ions:

| Biotically controlled | | Abiotically controlled | |
|---|---|---|---|
| Alk | ~+65% | Mg | ~+18% |
| C | ~+72% | Na | ~+55% |
| **Si** | **~+1235%** | **Cl** | **~+2%** |
| Ca | ~−55% | | |

The argument is that the abiotic ions are matched while the biotic ones are off for a known reason —
Si by 12× because there is no biogenic silica sink (diatoms), Ca low, Alk/C ~70% high.

### Three caveats to carry into the paper

1. **"Acid Ocean (>5 bar CO₂)" on slide 30 is not a planetary outcome.** It is the `acid_ocean`
   event later established to be a numerical guard mislabelled as physics (§9.1) — its trigger is
   atmospheric pCO₂ while its name refers to ocean pH, and the 5-bar cut lands *inside* the
   maximum-greenhouse peak (4 bar at S=0.5, 7 bar at S=0.7). With the ceiling at 10 bar at least one
   such run turned out to be a **snowball**, not a hothouse. The right-hand edge of the habitable
   stripe is a censoring artifact and the orange region means "runs we stopped", not a state.
2. **The dashed post-transition branches are drawn partly through non-steady-state points.** The
   June code has the same behaviour as the current version (`_plot_group_on_axes`: the line is drawn
   through the whole `group`, while `hab`/`non_hab` only control markers). Beyond the circles,
   plotted values are wherever each trajectory was when its event fired — one such run was still
   changing pCO₂ by **+611%** in its final 20% (§14, the `fast_4` discontinuity analysis). The
   *direction* of the inversion is sound physics; the magnitudes on that branch are not steady states.
3. **The outer edge is set by the 274 K seafloor floor, not by CO₂ exhaustion.** Below
   T_surface ≈ 285 K the seafloor is temperature-locked, so the weathering sink **stops weakening as
   the planet cools** — the cold-side negative feedback saturates. That is why the low-S branch
   stalls at ~8 bar CO₂ and then plunges into the snowball band instead of stabilising. The
   outer-edge location is therefore directly sensitive to that floor value, which is a modelling
   choice and should be stated as such.

### Two regressions this baseline exposes

1. **The Ca source.** These plots have a working thermostat and a clean habitable band because June
   got its Ca from **Wollastonite** — fast, Ca-only, high-solubility — *plus* the PHREEQC HT
   dissolution call. The Diopside swap (Ca–Mg coupled, ~40× slower, §5) is what put the current
   branch into Ca-starved hothouses, and `KD_MG_HT` is now ~3× too weak carrying the HT Ca load
   alone (§6.4). **Slide 26 is the specific behaviour to reproduce.**
2. **The Cl budget was calibrated in June and has since regressed.** Slide 25 shows **Cl at ~+2%**;
   the current branch sits at Cl = 165 mM against ~550 mM real, i.e. **~−70%** (§7). The open-items
   list treats Cl calibration as an unsolved problem — this slide says it was previously solved, so
   the June parameters are worth recovering rather than re-deriving.

---

## 14. Sweep history

Grid: instellation `s` × outgassing `out` × crust production `crust`, ~155–160 runs, 2 Gyr each.

| Sweep | What changed going in | Result |
|---|---|---|
| **fast_4** | New model (CIPW crust, salinity) | The Da/pCO₂ investigation (§8). Zero `converged`; everything pinned at a boundary. 50 runs piled at exactly 5.00 bar (`acid_ocean`). All runs `land_fraction = 0`. |
| **fast_5** | baseline for comparison | — |
| **fast_6** | 6 events → domain box | **1.20× slower overall**, but the *typical* run was not slower (median per-run ratio 1.01×). The entire slowdown was **10 runs (7%)** accounting for **104%** of the extra time — runs that used to exit early on `snowball`/`hothouse`/`co2_floor` and now sat inside the box for the full 2 Gyr. Confirmed **genuinely unconverged** (0.2–0.7 /Gyr drift vs a 0.1 threshold), i.e. the old events had been cutting off unsettled runs. `converged` fired **0/137**. |
| **fast_7** | short domain guard | **17.31 h → 8.00 h**; steps 279,438 → 86,910. 121/153 runs step-for-step identical, **0 runs with more steps** — exactly the intended signature. But introduced the **13 spurious `co2_low`** (§9.2). |
| **fast_8** | retry + continuous fallback | **15.18 h → 2.75 h (5.5×)** with **all 153 matched terminations byte-identical**. The `co2_low` regression persisted → habitable 25% vs fast_6's 34%. |
| **fast_9** | `co2_low` wall removed | **Zero `co2_low`**; all 22 recovered to `timeout`. `cold` (23) and `hot` (14) **identical** to fast_8 — the control that matters. Habitable 26% → 40%. Matched cost 3.94 h → 5.80 h. Top 12 runs = 45% of all wall time. |
| **fast_10** | RW minerals in pore list; convergence threshold raised | **`converged` fires for the first time — 23 runs (14.5%)**, up from zero in every previous sweep. Total wall **5.80 h → 4.24 h (−27%)**, max run 1509 s → 658 s. Mean fallbacks 64 → 23; Antigorite SI +344 → +81; pH 10.0 → 6.4. `s_0.8_out_0.1_crust_0.01`: 6117 fallbacks → **7**. |

> **fast_10's headline change was later reversed (§21.3).** Putting the RW minerals in the pore list
> did produce the convergence and fallback improvements tabulated above, but it also inverted the sign
> of the seafloor alkalinity flux over a broad region of (T, pCO₂), because the pore fluid carries
> ~100 mM Mg from dissolution and Sepiolite(d) removes essentially all of it. Every sweep from
> fast_10 through fast_16 carries this defect. The fast_10 caveats immediately below — the
> `crust_0.01` runs going the wrong way, and the 5 runs pushed to `co2_high` — read in hindsight as
> early symptoms of it.

**fast_10 caveats:** three `crust_0.01` runs still misbehave (2092, 856, 562 fallbacks — one with
more than one fallback per step, so its trajectory is substantially fabricated), and
`s_1.0_out_0.1_crust_0.01` went the **wrong way** (1 → 856 fallbacks) — the extra pore phases are
not uniformly beneficial. Also **5 runs went `timeout` → `co2_high`** from the pore mineral list,
not the threshold change.

**Two things that make plots look worse than the physics is:** `plot_results` classifies
`co2_high` (56–58 runs, ~36%) as *Unknown* rather than as the supply-limited runaway branch it
actually is; and until fast_10 there was no `converged` population at all, so the only positive
label was `timeout` ("survived 2 Gyr without leaving the box").

### The discontinuity artifact in the line plots

Real, and worth not re-diagnosing: the connecting line is drawn through the **entire** group
(`plot_results.py:525-528`, `group[col].values`), while the habitable/failed split (`hab` /
`non_hab`, defined at `:514-517`) is only used for the markers. So every
hothouse / acid_ocean / snowball terminal value — which are **transient snapshots, not steady
states** — gets wired into the curve. The s=0.90→0.95 jump (P: 0.54 → 3.16 bar) is a hothouse caught
mid-runaway, its pCO₂ still changing **+611%** in the final 20%. The converged branch itself is
genuinely converged (pCO₂ drifts <0.3% over the final 20% of 2 Gyr, `dP/dt → ~1e-14/yr`), which
confirms the user's position that **timeout runs have essentially converged**.

---

## 15. Open items, roughly in priority order

**The overall target is to reproduce §13 (the June seminar results) on the current chemistry.**
**Reordered 2026-08-19 by §20.** The iron charge-leak fix (§20.3) closed the correctness problem that
was blocking convergence; the transition is now visible (§20.4). The binding constraint has moved from
chemistry correctness to the **CO₂ ceiling** censoring the thermodynamic branch.

1. **Raise or taper the CO₂ ceiling** to the maximum-greenhouse peak (§9.2/§13), instellation-
   dependent rather than a flat 10 bar. The thermodynamic branch *is* the rising-pCO₂ branch, so it
   hits the wall almost immediately and only ~1 point past each transition survives. This is now the
   main limit on settled fraction (~26%) and on seeing each transition line in full (§20.5).
2. **Run a full-res sweep with the iron fix** (§20.3–20.4) to see the transition across many lines,
   not just the 4 that cleared ≥5 settled points in low-res fast_15.
3. **Seed Cl analytically in the sweep** (§20.4). Cl's residence time ∝ 1/crust reaches ~10 Gyr, so
   it never equilibrates from blank and drags many runs off "settled". `calibrate_earth.py` already
   does this; the formula generalises to any (out, crust). Note Cl is *not* what sets the convergence
   metric (Ca is) — this is about not discarding climatically-settled runs on a slow ion.
4. **Judge steadiness on pCO₂/T, not all species.** Many "drifting" timeouts are climate-settled
   (dP < 2%) while slow ions (Ca, Cl) equilibrate over >2 Gyr (§20.4). This is also why `converged`
   fires ~never (§14).
5. **The Da diagnostic** (§20.5, §8 note) — `Da_alk` reads 2–3 decades below the transition on
   crust=10 lines and `Da_eff` is no better; the reversal needs Da>1 *and* steep retrograde
   d b_eq/dT, so a two-condition regime map is the honest diagnostic. Not implemented.
6. **Na is still ~1 mM** (§20.5). The source is kinetically capped; no equilibrium tuning raises it
   (§20.2). Adopt Na as a swept O(0.1–1) mol/kg inventory (Kite & Ford, §19.6) rather than trying to
   generate it.
7. **Re-examine `kd_mg_ht`** (§18, currently 0.07). ~~With charge now consistent and Ca bounded, the
   HT exchange's true strength can finally be assessed on clean runs.~~ **Done, §22.1/§22.4** — the
   Earth fit gives **0.019** and Coogan's HT Ca flux implies **0.005–0.009**, against §18's 0.07.
   Three anchors, still disagreeing; §18's "Earth cannot constrain it" no longer holds.
8. **The Cl budget** — the problem is the Na:Cl *ratio*, not Cl's absolute value (§19.3); Cl is
   already 3× below Earth's. Do not "fix" it by lowering Cl.
9. **Cheap performance wins:** SI-list trim (~8–11%, verified bit-identical), dead Jacobian
   columns (~13%, exact).
10. ~~**`calibrate_earth.py` is broken** (`f_bio=`), deliberately deprioritized.~~ **Fixed and run,
   §22.** The ill-posed `kd_mg_ht` ratio update (§18) was replaced by a bounded `least_squares`,
   which is what the calcite-saturation bistability requires. **New top-priority item in its place:
   decide `alpha` (§22.9) — it is not identifiable from Earth but carries the thermostat on the
   land-free worlds the sweeps target.**
11. Consider the **equilibrated `b0`** to remove the unphysical spin-up U (§9.3).

**Closed / superseded:**
- **Iron charge-leak → Ca runaway** (§20.3) — fixed, the root cause of the convergence failure.
- **The two flooding prerequisites** (was item 1) — settled differently than expected (§20.1–20.2):
  the shadow is Al not Si, the buffer is Kaolinite-only and was *reverted*, and w/r=3 (low, not the
  600 first tried). Na turned out to be kinetically capped, not equilibrium-limited.
- **`dY_dt` clamps charge-neutral** (was item 3) — the clamp was *exposing* the b_pore bug, not
  itself the cause; with the flux self-consistent the clamp is a no-op (net Fe ≈ 0) and the leak is
  gone. A general clamp-charge correction is no longer needed for this failure mode.
- **The Mg flooding / "Mg sink problem" (§12)** — the w/r=3 + no-buffer configuration keeps b_eq[Mg]
  bounded; the runaway Ca that dominated it is fixed.
- **Clogging** (§19.8 — keep `clog=False`); the **"residual carbon-sink deficiency"** (§19.2 — does
  not exist).

---

## 16. Cross-cutting lessons

- **Never benchmark `dY_dt` on a fixed `Y`.** PHREEQC warm-starts, so results are path-dependent and
  a naive A/B inverts. Compare with `time_evolve()`. This is also why a "baseline" replay of failing
  inputs recovers 58.8% rather than 0% — the fair comparison is always knob-vs-knob under identical
  replay.
- **Microbenchmarks mislead here.** Forward differences looked like a 2× saving and were 7×
  slower; the SI list was the rare case where the microbenchmark held up end-to-end.
- **Correlation found the wrong culprit twice** — the OLR blend window (78% of small steps, not
  causal) and the charge imbalance (r = 0.029). Both times the fix came from reading an actual
  failing input rather than theorising.
- **`sol.t` only records accepted steps.** Solver thrashing (rejected steps, Jacobian re-forms) is
  invisible in it; measure derivative calls per accepted step. Likewise `sol.y = np.maximum(sol.y, 0)`
  masks negative excursions in the stored output — log the unclipped `Y` LSODA actually passes in.
- **Cost estimates for un-terminating runs were ~40% low, twice.** Recovered runs are well above the
  median timeout cost; don't estimate from the median.
- **PHREEQC is self-consistent within a call, not across the kinetic-flux layer between calls.**
  That gap is the origin of the phantom-Al bug and is worth suspecting whenever a `k`-derived and a
  `b_eq`-derived quantity are combined per-element.

---

## 17. Key references used

- **Maher & Chamberlain (2014)**, *Science* — fluid-transport weathering model, Damköhler
  coefficient `Dw`, `C_eq`. `/data/pt426/Downloads/science.1250770-1.pdf`
- **Hakim et al. (2021)** — CHILI model; primary-mineral equilibrium, γ = 1 ideal assumption, β_th.
  `/home/pt426/Documents/Hakim et al...`; code in `src/kamino/H21/chili/` (third-party)
- **Coogan & Dosso (2022)**, *GCA* **329**, 22–37 — Cenozoic seawater chemistry.
  `/home/pt426/Documents/Coogan and Dosso - 2022 - Controls on the evolution of Cenozoic seawater che.pdf`
  Tracks alkalinity as its own variable, but pins it at ~2.4 mmol/kg via a buffering carbonate sink,
  enforces charge balance every step with **Ca as the free balancing ion**, and does **not track Na,
  Cl or SO₄ at all**. That last assumption is exactly what Kamino cannot make.
  **§2.3–2.5 + Table 1 are the source of the LT/HT flux targets in §22.4–22.5.**
- **Coogan & Dosso (2026)**, *EPSL* **677**, 119811 — adds a carbon cycle to the above; Fig. 6 gives
  the Ca\* fluxes per process used to cross-check §22.4.
  `/home/pt426/Documents/Coogan and Dosso - 2026 - A model for Cenozoic seawater chemistry and carbon.pdf`
- **Coogan et al. (2019)**, *EPSL* **508**, 41–50 — source of `ALK_LT` (alkalinity per kg rock
  altered) and the Arrhenius temperature dependence of the altered-basalt fraction; Troodos ophiolite.
- **Coogan & Gillis (2013, 2018)** — low-temperature off-axis alteration as a silicate-carbonate
  weathering pathway, and the temperature dependence of chemical exchange during seafloor weathering.
  Also the source of the `t_clog_ref = 20 Myr` in `seafloor_reactive_area` (§19.8).
- **Dunlea et al. (2017)**, *Nat. Commun.* **8**, 844 — authigenic Mg sink in deep-sea sediments
  (0.02 Tmol/yr), the lower bound on Coogan's diagenetic ("reverse weathering") Mg flux.
- **ThermoChimie `sit.dat` v12a** (ANDRA/NWS/ONDRAF); **SUPCRTBL** (Zimmer et al. 2016);
  **Kinec_v3** (Hermanská et al. 2022/2023); **Thermoddem** (Blanc et al. 2012)
- **PRIMELT1/2/3** — Herzberg & O'Hara (2002), Herzberg et al. (2007), Herzberg & Asimow (2008)
- **pMELTS** — Ghiorso, Hirschmann, Reiners & Kress (2002), *G3* **3**(5), 1030. The MELTS
  calibration for peridotite melting at 1–3 GPa; used by `make_crust_compositions.py` (§23.5–23.6)
  via **alphaMELTS for Python** (Antoshechkina & Ghiorso 2018).
- **Brugman, Burney & Walter (2021)**, *JGR Planets* **126**, e2020JE006731 — experimental solidi and
  partial melts for two hypothetical exoplanet mantles. ⚠️ **HEX1 is a Mg/Si end-member (1.42);
  HEX2 is a Ca/Al end-member (1.07).** They are not two ends of one axis — see §23.4.
- **Spaargaren et al. (2023)**, *ApJ* **948**, 53 — plausible bulk compositions of terrestrial
  exoplanets in the solar neighbourhood; mantle Mg/Si 0.8–1.6 spans quartz-bearing to
  ferropericlase-bearing mantles (Earth 1.23), which bounds the useful sweep range.
- **McDonough & Sun (1995)**, *Chem. Geol.* **120**, 223 — pyrolite bulk mantle, the starting
  composition for the melting grid.
- **Guimond, Wang, Seidler, Sossi, Mahajan & Shorttle (2024)**, *Rev. Mineral. Geochem.* **90**, 259
  — "From stars to diverse mantles, melts, crusts and atmospheres of rocky exoplanets"
  (arXiv:2404.15427). **The prior work for the whole crust pipeline** — §4 is melts and crusts,
  Fig. 5 is the Mg/Si mineralogy sweep, Fig. 8 the pMELTS solidi. Shorttle is a co-author. See §24.6.
- **Putirka & Rarick (2019)**, *Am. Mineral.* — >4,000 Hypatia stars; core formation controls half
  or more of mantle mineralogy variation; classification by (FeO+MgO)/SiO2.
- **Riel, Kaus, Green & Berlie (2022)**, *G3* **23**, e2022GC010427 — MAGEMin; with
  **Holland, Green & Powell (2018)**, *J. Petrol.* **59**, 881 — the igneous thermodynamic dataset.
- **Médard et al. (2004)**, *Contrib. Mineral. Petrol.* — ultracalcic primitive melts; CaO/Al2O3 > 1
  requires volatiles or a cpx-rich source, NOT volatile-free fertile lherzolite (§24.4).
- **Nature Comms Earth Environ (2022) 3, 261** — mantle buffering at near-constant homologous
  temperature; **Korenaga (2016)**, *Sci. Adv.* — the counterargument that self-regulation is too
  slow (§24.3).
- **Robie & Hemingway (1995)** — nepheline molar volume (Vm = 54.16) for the database entry.
- **Katz (2003)**, **Niu (1997)**, **McKenzie & Bickle (1988)**, **Hirschmann (2000)** solidus,
  **Langmuir/Klein/Plank (1992)** — melting parameterizations (explored, not used)
- Seafloor alteration: **Berndt & Seyfried** plagioclase-exchange experiments; **Mottl & Holland**;
  **Seyfried & Bischoff** seawater–basalt experiments

---

## 18. `kd_mg_ht` recalibration (2026-08-17): 0.012 → 0.07

### The structural finding: Earth cannot constrain this parameter

Budget breakdown at the Earth steady state (`land_fraction=0.3`, depth 3700 m, S=1, out=1, crust=1),
in Tmol/yr:

| term | Mg | Ca |
|---|---|---|
| continental | **+2.78** | **+7.16** |
| seafloor LT | −2.30 | +0.04 |
| reverse weathering | −0.44 | 0 |
| **HT Mg–Ca exchange** | **−0.04** | **+0.04** |
| precipitation + shelf carbonate | 0 | −7.24 |

The exchange is **1.5% of Mg removal and 0.6% of the Ca source** on Earth — continental weathering
dominates both. Switch land off and its Ca share jumps to **35.5%**.

Two consequences:

1. **`calibrate_earth.py`'s update rule for it was ill-posed.** `KD_mg *= (Mg/T_Mg)**0.6` scales the
   constant by the Mg error, but the constant governs 1.5% of Mg removal, so the loop was fitting it
   to the residual of the *other* sinks. Measured directly: kd 0.012 → 0.073 moves Earth by **0.2 K**
   and leaves Ca at 0.20 mM either way. That is very likely how the June value survived the removal
   of the HT PHREEQC call without anyone noticing it had become meaningless.
   **The rule is now commented out**, with the reasoning recorded in place.
2. **It must be calibrated on land-free worlds**, or from first principles. It was set from the latter
   and *verified* on the former.

### The value

Because `surface_area` cancels in `_ht_rate`, the exchange reduces to
`tau_Mg = ocean_depth·1000 / (kd · J_total)`, so **kd is the fraction of circulating Mg removed per
pass** — a physical quantity with literature bounds. Three independent anchors:

| anchor | implied kd |
|---|---|
| Ca supply balances Earth outgassing (7.5 Tmol/yr) at Earth's Mg (52.8 mM) | 0.073 |
| axial hydrothermal flux incl. diffuse (~1e14 kg/yr) ÷ model's 1.4e15 kg/yr, Mg fully stripped at high T | 0.071 |
| literature high-T axial Mg removal (4e12 kg/yr × 53 mM = 0.21 Tmol/yr) vs the model's 0.041 | 0.061 |
| **chosen** | **0.07** |

**Honest caveat:** counting only *focused* high-T venting (4e12 kg/yr) gives **0.003**, four times
*smaller* than the old value. The range is therefore wide and the choice depends on crediting diffuse
axial flow with substantial Mg stripping; the top is preferred because the functional and
literature-flux anchors agree there independently.

### Verification

Water world, out=1, crust=1, depth 3000 m:

| S | kd = 0.012 | kd = 0.07 |
|---|---|---|
| 0.4 | `cold`, T=181.0, 9.97 bar | `cold`, T=180.9, 9.97 bar — **unchanged, no new snowball** |
| 0.5 | co2_high, Mg=14.8, Ca=41.8, Alk=68 | co2_high, Mg=2.0, Ca=51.2, Alk=62 |
| 0.6 | co2_high, Mg=157, Ca=1.37, Alk=300 | co2_high, Mg=26.9, Ca=13.8, Alk=67 |
| **0.7** | **co2_high**, 10 bar, Mg=200, Ca=0.54, Alk=388 | **timeout**, 8.34 bar, **Mg=53.9**, Ca=42.0, Alk=25 |
| 0.9 | co2_high, Mg=0, Alk=0 | **byte-identical** — exchange inert (see below) |
| 1.0 | co2_high, Mg=191, Ca=0.38, Alk=369 | co2_high, Mg=28.2, Ca=1.88, Alk=57 |
| 1.1 | co2_high, Mg=158, Ca=0.34, Alk=305 | co2_high, Mg=54.0, Ca=8.24, Alk=27 |

**What it fixed:** the Mg/Ca inversion across the whole grid (Mg falls 3–7×, Ca rises 4–30×, runaway
alkalinity falls from 300–390 mM to 25–67 mM), and the S=0.7 run stops running away. At S=0.7 the
ocean settles at **Mg = 53.9 mM against Earth's 52.8** — not fitted to that, which is the strongest
single piece of evidence for the value. Earth is unaffected, and the cold side gains no new snowballs.

**What it did not fix — three limits to carry forward:**

1. **Inert wherever Mg is clamped at zero.** `_ht_rate ∝ [Mg]`, so at S=0.9 all of kd = 0.012–0.12
   give byte-identical results. Fixing the clamp masking (§9.2 mechanism, fast_11 artifact D) is a
   prerequisite for this recalibration to act across the grid.
2. **pCO₂ still ~8 bar**, even at kd = 0.12 (T ≈ 337 K). June's S=0.7 point was ~300 K at ~1 bar.
   So the Ca supply was *not* the binding constraint on carbon burial. **§19.2/§19.3 identify what
   is — and it is NOT a carbon-sink deficiency; that framing is retracted.** With Na at 0, Cl
   consumes 86% of the cation charge and the ocean retains almost no alkalinity to hold CO₂. The HT
   exchange is alkalinity-neutral by construction, which is exactly why raising kd could not help.
3. **Costs some solver time at high instellation.** S=1.1 went from 0 fallbacks / 15 s to 142 / 350 s;
   S=1.0 from 3 to 20. Bounded by the new fallback cap, but worth watching in the next sweep.

---

## 19. fast_12, and the Na chemistry investigation (2026-08-17)

### 19.1 fast_12 — the kd recalibration sweep

Two deliberate changes from fast_11: `kd_mg_ht` 0.012 → 0.07 (§18) and the new
`MAX_CHEMISTRY_FALLBACKS = 5000` cap. Grid otherwise identical, so it is a clean paired A/B on
**921 matched configs** (analysed at 921 of 931 complete; full detail and `summary.csv` in the
sweep's own `CONTEXT.md`).

**94.8% of configs landed in the same state.** Net gain of 33 runs, dominated by **30 recovered
from the CO₂ wall**; habitable 25.4% → 29.0%.

| | fast_11 | fast_12 |
|---|---|---|
| median Mg | 5.67 mM | **1.87** (0.33×) |
| median Ca | 0.67 mM | **2.06** (3.07×) |
| Ca < 1 mM (starved) | 55% | 44% |
| Mg > 100 mM (runaway) | 11% | 4% |
| **Alk == 0 (clamp-pinned)** | **23%** | **23%** |
| mean fallbacks | 698 | **20** |
| runs > 10,000 fallbacks | **6** | **0** |
| slowest run | 33,924 s | **3,310 s** |
| elapsed | 10.30 h | **1.87 h** |
| trustworthy (steady) | 70 | **122** |

**The cap is not what fixed the chemistry — the recalibration is.** Only 2 runs hit the cap and the
count with any fallbacks barely moved (361 → 357), but the >10,000 bucket went 6 → 0. This
causally confirms §12: the Mg runaway to 200–320 mM was what drove Antigorite to 10⁸¹×
supersaturation with no permitted sink, producing unsolvable states. Fixing the Ca sink stopped the
Mg runaway, which fixed the solver. The cap is now cheap insurance.

**Versus June (§13):** the ratio degeneracy is substantially restored (at fixed ratio = 1 the
habitable fraction across (out, crust) went from `0.68, 0.53, 0.26, 0, 0, 0, 0` to
`0.68, 0.63, 0.53, 0.37, 0.11, 0.05, 0`) and the band's upper edge extended (ratio = 3: 0% → 21%).
But the **thermostat is untouched** — median dT/dT_noFB 1.36 → 1.33 against June's 0.17 — the band
median is still 1.23 decades low, and no line shows a plateau. Exactly as §18 predicted.

### 19.2 There is no carbon-sink deficiency — this supersedes §18's framing

**Correction to §18 and to both sweep CONTEXT.md files.** At the fast_12 S=0.7 steady state:

| | Tmol/yr |
|---|---|
| C outgassing in | 7.498 |
| C buried by precipitation | **7.511** |
| Alk needed to bury it as CaCO₃ (2 eq/mol) | 15.00 |
| Alk supplied by seafloor weathering | **15.16** |
| Calcite SI | **+0.06** (at saturation) |

Carbon in balances carbon out, the alkalinity supply meets the requirement, and calcite is
precipitating at saturation. **Nothing is failing to be buried.** The "residual carbon-sink
deficiency" was a misdiagnosis.

### 19.3 The real mechanism: Cl eats the cation charge because Na is dead

Alkalinity is the leftover positive charge after the negative ions are balanced, and it is what
lets an ocean hold CO₂. Charge accounting at that state:

| ion | conc | charge |
|---|---|---|
| Ca | 42.04 mM | +84.08 mEq |
| Mg | 53.92 mM | +107.84 mEq |
| **Na** | **0.00 mM** | **0.00** |
| Cl | 164.60 mM | −164.60 mEq |
| | cations **191.93** | **Alk = 27.3 mEq** |

**Cl consumes 86% of the available cation charge.** Alk/C is 0.169, so the DIC is overwhelmingly
dissolved CO₂ — pH 5.26, pCO₂ 8.34 bar.

**Decisive test** — hold C, Ca, Mg, Cl fixed and add only Na:

| Na (mM) | Alk (mEq) | pCO₂ (bar) | pH |
|---|---|---|---|
| 0 | 27.3 | **11.94** | 5.34 |
| 100 | 127.3 | 2.28 | 6.71 |
| **165** (Cl balanced) | 192.3 | **0.031** | 8.51 |

**165 mM of Na drops pCO₂ by a factor of 385**, with no change to carbon, weathering or burial. At
that Na, calcite is +3.94 supersaturated, so carbon would keep precipitating and drive pCO₂ lower.

This also explains why `kd_mg_ht` couldn't help: **the HT exchange is alkalinity-neutral by
construction** (Mg²⁺ out, Ca²⁺ in, 1:1). It cannot change the charge balance — it only relabels
which cation carries the alkalinity that already exists. Hence Ca rose 3× and pCO₂ didn't move.

And it explains the inverted thermostat: pCO₂ is set by the Na/Cl charge residual, which temperature
does not affect, so the weathering–climate feedback is largely *disconnected* rather than weak.

> **Note on framing.** Earth's 470 mM Na is *not* the target — ocean-world chemistry may be nothing
> like Earth's and the Na cycle there is essentially unknown. The defect is structural and holds on
> any planet: the model outgasses Cl as a strong acid with **no monovalent cation source of
> comparable magnitude**, so the charge balance is forced onto Ca and Mg.

### 19.4 Why Na is zero — and it is not what it looks like

- **Na is at a genuine steady state at ~1.3 mM**, not a numerical collapse: source ≈ 0.006 Tmol/yr
  against the parameterized sink (0.0044 Tmol/yr at 1 mM). Matches the observed 1.6 mM maximum
  across 921 runs. The balance is real; the source is ~1000× too small to matter.
- **The source is tiny because Albite is saturated, not absent.** The rock is 20% Albite, but
  b_eq[Na] = **7.6×10⁻⁶ mM** in the full assemblage versus **523 µM for Albite alone in pure
  water** — a 69,000× suppression. The other primaries flood the pore fluid (b_eq[Si] = **170 mM**,
  b_eq[Mg] = **382 mM**), and albite saturation depends on [Na][Al][H₄SiO₄]³, so the Si³ term alone
  drives the fluid to saturation.
- **`dissolve_only=True` is already set on the LT path and does not help.** It forbids
  precipitation; it does not force dissolution. A saturated phase contributes nothing either way.
  This is the subtle trap.
- **The pore precipitation is a weak Na *source* (+0.006 Tmol/yr), not a sink.** Saponite-Na is
  Na₀.₃₄Mg₃Al₀.₃₄Si₃.₆₆O₁₀(OH)₂ — 0.34 Na against 3 Mg, i.e. a Mg sink that takes a trace of Na, and
  Al-limited besides.
- **Al-bearing Na sinks are all dead ends.** Analcime carries 1 Na per Al — the same ratio as
  Saponite-Na — and Al is immobile (b_eq[Al] ≈ 2.5×10⁻⁴ mM). Zeolites also live only in
  `lt_weathering_sit.dat`, not the runtime `hybrid_ocean.dat`. Not worth pursuing.
- **Al-free sinks already exist** (Nahcolite in `carbonate_minerals`; Halite when land > 0), so the
  sink side is not the blocker. **The problem is entirely the source.**

### 19.5 Root cause: one incomplete change, never committed

| | June (`4059acc`, working) | current working copy |
|---|---|---|
| primary exclusions at LT | **Anorthite, Forsterite, Enstatite excluded unconditionally** | none (`full_equilibrium=True` → `exclude_primary=False`) |
| `full_equilibrium` | did not exist | `True` on the LT call |
| water/rock ratio at LT | n/a | **not passed** |
| silica in pore list | n/a | **absent** |
| Mg | added back kinetically, Mg index only | in the equilibrium, floods to 382 mM |
| Ca source | **Wollastonite** (fast, Ca-only) | Diopside (slower, Ca–Mg coupled) |

June never had Forsterite or Enstatite in the LT equilibrium, so there was no Si/Mg flooding, albite
was undersaturated, and Na flowed. The comment at `chemistry.py:462-467` states the prerequisites
for including them:

> *"Including them needs a realistic water_rock_ratio at LT plus secondary phases that can buffer
> Mg/Ca/Si; with those, it converges cleanly. Until that is settled, keep them excluded."*

**The exclusions were dropped and neither prerequisite was added.** That single half-completed change
explains the Na shadowing, the Ca starvation, the Mg flooding, the PHREEQC failures, and by extension
the pCO₂ ceiling and the inverted thermostat. `full_equilibrium=True` is **not at HEAD** — it exists
only in the uncommitted working copy. The comment is now actively misleading and should be corrected.

### 19.6 The adopted fix, and options rejected

**Rejected — uniform kinetic primary dissolution.** Proposed as the "consistent" alternative to
special-casing Albite: dose all primaries kinetically and let only secondary phases set the ceiling.
**Correctly rejected**, because the equilibrium formulation is what supplies a thermodynamic limit to
minerals with *no secondary counterpart* — a pyroxene is capped by its own saturation whether or not
a corresponding clay exists. Dosing it kinetically removes the affinity term, and with it the Da > 1
regime and the inner-edge feedback that are the paper's central novelty.

**Rejected — charge-balanced Cl outgassing as a parameterization.** The *mechanism* has strong
precedent (the Rubey/Holland acid-gas picture: degassed HCl leaches rock, yielding a chloride ocean
carrying cations in rock proportion). But read carefully that says the acid is neutralised *by rock
dissolution* — which is the weathering path itself, not a separate co-outgassing term. No precedent
was found for co-outgassing NaCl with a tunable neutral fraction.

**Adopted — cap the flooding, keep the formulation.** Two changes, no per-mineral special-casing, and
`b_eq` stays a genuine assemblage equilibrium so Da and the transport limit keep working:

1. **Add a silica phase (`SiO2(am)`) to `pore_precipitating_minerals`** — it exists in
   `silica_minerals` but only in the *ocean* list.
2. **Pass a realistic `water_rock_ratio` to the LT `get_weathering_flux` call** — the parameter
   exists and is already used on the HT path.

Measured at a fast_12 pore state:

| | production now | + silica cap + w/r = 10 |
|---|---|---|
| b_eq[Si] | 170 mM | **3.1 mM** (near silica saturation) |
| **b_eq[Na]** | **7.6×10⁻⁶ mM** | **15.7 mM** |
| b_eq[Ca] | 42.5 mM (≈ ocean, no net source) | **259 mM** (net source) |
| b_eq[Mg] | 382 mM (net source) | **1.4 mM** (net sink) |
| PHREEQC | **errors** at w/r = none | converges cleanly |

This resolves a design tension that looked unresolvable — June got no-flooding by hand-picking three
minerals; `full_equilibrium` got no-special-casing by accepting the flooding; uniform kinetic dosing
would have got both by sacrificing the thermodynamic limit. The two prerequisites get all three.

**Also worth considering (Kite & Ford 2018,** [arXiv:1801.00748](https://arxiv.org/abs/1801.00748)**):**
their waterworld model releases Na⁺ and Ca²⁺ from seafloor basalt, quotes **O(0.1–1) mol/kg** as
geologically reasonable ocean-world Na with no reference to Earth, and treats the cation inventory as
acquired early and thereafter fixed. That supports adding **ocean Na as a swept input parameter** over
that range, so results don't hinge on an unknown cycle. (Their fixed-inventory *justification* rests
on deep-waterworld seafloor pressures and is weak at the sweep's 3000 m.) Encouragingly, a rate-based
estimate for kinetic albite release also lands near 0.1 mol/kg, inside their range.

### 19.7 Consequences and constraints found along the way

- **Silica is 78% of the sediment blanket** (43.7 of 51.2 Tmol/yr; 1.94 of 2.49 m/Myr). Since
  `clog=False`, sediment cover alone sets the area, so moving silica removal into the pore space
  raises `A_reactive` **0.819 → 1.252 (+53%)** and increases weathering. Not contained — but the
  blanket it displaces is itself artifact-driven (b_eq[Si] is 55× amorphous-silica saturation), and
  the direction is helpful given how marginal the alkalinity supply is.
- **Pore-space silica is defensible** on nucleation surfaces (the inhibition is a nucleation
  barrier, and pore space is rock-lined), higher supersaturation, and long residence time; the pore
  list already contains Kaolinite and Sepiolite(d) at instantaneous equilibrium, which are *more*
  kinetically demanding. `SiO2(am)` is the right phase, not Quartz/Chalcedony. The honest way to
  honour the kinetics is a **positive `precipitation_SI`** — which exists (`chemistry.py:319`,
  `:350`), defaults to 0, and is **not exposed by `get_precipitation`** (small plumbing job).
- **Water/rock ratio becomes a sensitive parameter** (Na 4.1 → 15.7 mM, Mg 49 → 1.4 mM between
  w/r = 100 and 10). It needs a literature value, strongly water-dominated for off-axis flank
  circulation. ⚠️ **w/r = 2.0 returned results identical to `none`** — the parameter may not bite at
  low values; check before relying on it.

### 19.8 Clogging: investigated, and it should stay off

- **The exponential parameterisation reverses the feedback sign.** `t_clog = 20 Myr ·
  exp(−(T−280)/7)` gives d ln A/dT = **−0.143 /K** (Eₐ ≈ 93 kJ/mol) against the rate constants'
  **+0.075 to +0.090 /K** (Eₐ ≈ 63 kJ/mol). It doesn't cancel the thermostat, it overwhelms it:
  measured d ln F_alk/dT goes from **+0.09/+0.24** (off) to **−0.047/−0.052** (on). At T_pore = 320 K
  the alkalinity flux falls **8.49 → 0.010 Tmol/yr (850×)**, and Da is crushed to **5.4×10⁻⁵**,
  removing the transport limit the paper is about. It also dominates `t_cover` above T_pore =
  **275.1 K**, i.e. always, and is extrapolated **11 e-foldings** beyond its 280 K reference.
- **It is half of the Coogan & Gillis mechanism** — their T-dependent crustal precipitation is
  *stabilising* because the carbonate is itself the carbon sink. The clog term takes the area
  reduction without the carbon burial. Completing it (Calcite in the pore list) doesn't rescue it:
  the carbon sink appears only above ~305 K and the alkalinity flux goes **negative**.
- **Self-consistent clogging is the right formulation but is negligible.** Deriving
  `t_clog = POROSITY·PORE_DEPTH / Σ_m J·n_m·MW_m/ρ_m` (reusing the `S_sed` machinery) is
  parameter-free and self-limiting. Across **5 ocean states × 6 pore assemblages**, `t_clog` ≥ 10⁵ Myr
  in **26 of 30 cells**, and F_alk and d ln F/dT were identical to four decimals with and without it.
  **Cause: the pore phases are Al/Fe-limited.** Four of six need Al or Fe, and pore Al ≈ 3×10⁻⁷,
  Fe ≈ 1×10⁻¹² mol/kgw. Verified: Kaolinite forms 3.0×10⁻⁸ mol/kgw and remains at **SI +0.20** —
  element-limited, not equilibrated. Only Calcite gives finite timescales (61.5 Myr temperate), and
  it makes F_alk negative.
- **Verdict: keep `clog=False`**, now with a quantitative justification rather than an assumption —
  the model's own precipitation rates imply t_clog ~10⁷ Myr. The exponential should be deleted or
  clearly marked unusable so nobody enables it expecting physics.

### 19.9 Two incidental findings that matter more

**The intrinsic feedback collapses with temperature** — measured on the weathering law itself, with
no clogging involved:

| state | T_pore | d ln F_alk/dT |
|---|---|---|
| cold | 283 K | **+0.091 /K** |
| temperate | 292 K | +0.078 |
| warm | 338 K | +0.019 to +0.034 |
| hot | 389 K | **+0.001** |
| stripped | 349 K | **−0.017 to +0.007** |

The thermostat weakens monotonically and effectively vanishes above ~340 K. This localises the
anti-thermostatic behaviour of §14/§19.1 to the weathering law rather than inferring it from
trajectories.

**The pore assemblage is a 3.4× lever on the alkalinity flux.** Warm state:

| pore set | F_alk (Tmol/yr) |
|---|---|
| clays + silica | 64.1 |
| silica only | 62.6 |
| clays | 59.6 |
| **clays + RW (production default)** | **17.7** |

Adding the reverse-weathering minerals cuts the alkalinity supply 3.4×. Since that supply is the
binding constraint on carbon burial (15.16 vs 15.00 needed), this single choice matters far more for
the CO₂ ceiling than clogging does — and `reverse_weathering=True` is the sweep default. **Untested:
whether `reverse_weathering=False` changes the ceiling behaviour.**

---

## 20. The step-1 buffer, water/rock ratio, and the iron charge-leak bug (2026-08-17 → 08-19)

This section covers the arc after §19: attempts to un-shadow Na via the LT equilibrium, the design
tension that surfaced, and — the payoff — a genuine **correctness bug** whose fix (one line) closed
the charge leak, killed the Ca runaway, and made the kinetic→thermodynamic transition visible for the
first time.

### 20.1 The two-step weathering calculation, and where Na was being shadowed

`get_weathering_flux` runs in two steps: (1) `get_b_eq` — dissolve the primary rock to an equilibrium
ceiling `b_eq`; then the M&C flux formula; (2) `get_precipitation` — let secondary phases precipitate
out of the resulting pore fluid. **Step 1 had no secondary phases in it** (the HT path always did, via
`ht_secondary_minerals`; LT never did). With nothing buffering the fluid, the primaries flood it —
b_eq[Si] to 170 mM, b_eq[Mg] to 382 mM — and Albite (NaAlSi₃O₈) sits at saturation, releasing no Na:
**b_eq[Na] = 7.6×10⁻⁶ mM** versus 523 µM for Albite alone in pure water, a 69,000× suppression.

Key correction to the §19.6 diagnosis: the shadow is set by **aluminium**, not silica. Al is the
trace species in Albite's solubility product (~10⁻⁷ mol/kgw), so the Al sink gates feldspar
dissolution. Tested by adding phases to step 1: silica-only did nothing (F[Na] 0.34 → 0.34 at the
hot state); **Kaolinite alone** (the Al sink) recovered Na and was the minimal, best-converging
choice. The Mg sinks (reverse-weathering clays) must **not** go in step 1 — they strip Mg from the
through-flowing seawater, and each Mg²⁺ removed takes 2 eq of alkalinity with it (measured F[Alk]
+2.5 → −207 Tmol/yr).

`dissolve_only=True` is already set on the LT path and does **not** help — it forbids precipitation,
it does not force a saturated phase to dissolve. This was the subtle trap.

### 20.2 The design tension: Na source vs thermodynamic feedback (the water/rock ratio)

Kaolinite in step 1 needs a `water_rock_ratio` to converge (without one PHREEQC fails on every call —
see fast_13 below). But w/r turned out to control something more important than Na, and the two pull
opposite ways:

- **b_eq[Na] barely depends on w/r** (3.1–3.2 mM across w/r 3–600 at the near-Da=1 states; the
  feldspar-derived ions are kinetically limited, so `F → k·A`, independent of b_eq). Na is a
  *kinetic-source* problem, not an equilibrium one — no w/r or buffer choice raises the *delivered*
  Na flux meaningfully near the transition (measured −0.013 to +0.009 Tmol/yr across the whole range).
- **The temperature dependence of `b_eq` does depend on w/r, strongly.** `d ln charge(b_eq)/dT` is
  **−0.044/K at w/r ≤ 10** (steep, retrograde — the feedback that drives June's pCO₂ reversal) but
  **collapses to ~−0.002/K at w/r = 600**. At high w/r there is too little rock for solubility to
  bind: Mg is rock-*supply*-limited (b_eq[Mg] flat vs T) rather than solubility-limited.

So high w/r (fast_13) gave sweep-wide Na but killed the feedback; low w/r restores the feedback but
loses Na. **These are not jointly satisfiable from the equilibrium alone**, which is the same pattern
noted at the end of §19 and is worth stating in the paper: a self-consistent equilibrium chemistry and
the June transition may not both come from the same `get_b_eq` call. Resolution adopted: **w/r = 3**
(on the feedback plateau; also closest to the instantaneous pore geometry, w/r ≈ 0.04, and formally
correct because M&C already applies the supply limit separately via `A/J` — starving the rock in
`b_eq` double-counts it). The Kaolinite buffer was **reverted** (`lt_equilibrium_buffer_minerals = []`)
because it does not help Na where the transition is and it costs convergence; the wiring is kept, off.

`WATER_ROCK_RATIO_LT` was briefly 600 (derived as the *integrated* flow-through ratio
`J·t_exposure/rock_mass`); the **integrated** ratio is the wrong quantity here — b_eq should use the
instantaneous rock-dominated ratio, so the constant is now **3**.

### 20.3 The iron charge-leak bug — the actual root cause (THE finding)

Chasing why runs would not converge led to Ca accumulating **linearly at ~119 mM/Gyr, forever**, with
no steady state. Tracing it:

1. **Tracked alkalinity diverged from the ion charge by up to 253×** (charge sum 491 mEq while tracked
   Alk carried 1.9). This is the §7 "tracked Alk = ion charge" invariant failing — but not as a clamp
   residual. It grew monotonically, the signature of a **persistent per-step flux error integrating
   over 2 Gyr**.
2. **The flux audit was clean** — every term's alkalinity component equalled `ION_CHARGE·flux`
   exactly. But the *clamped* `F_net` leaked at a steady **+0.35 Tmol eq/yr**, and it was **iron**:
   net seafloor Fe was −0.18 Tmol/yr (more removed than supplied), the clamp blocked the removal
   (ocean Fe = 0), but alkalinity had already been debited 2×0.18 for it. 0.35 × 2 Gyr = 7.5×10⁸ Tmol
   eq — matching the observed gap to two significant figures.
3. **Why Fe removed more than existed** (the user's insight): `weathering.py` computed the same
   dissolution two different ways —
   - flux: `F = A(b_eq−b_in)/(b_eq/k + A/J)` = `J(b_eq−b_in)·Da/(1+Da)`
   - pore fluid: `b_pore = b_in + (b_eq−b_in)·(1−exp(−Da))`

   `Da/(1+Da)` and `1−exp(−Da)` **agree as Da→0 and Da→∞ but differ by up to 30% near Da≈1**. Ca, Mg,
   Si sit at Da~10⁻⁴ (kinetic); Al at Da~175 (transport); **Fe is the only species near Da≈1** (5.07),
   so b_pore overstated dissolved Fe by 1.190×. Goethite then stripped the phantom excess. Observed
   overstatement 1.1896× — matching Da/(1+Da) vs 1−exp(−Da) at Da=5.07 exactly.

**The fix (one line, `weathering.py`):** derive the pore fluid from the flux, so the two cannot
disagree:
```python
b_pore = b_input + F_primary / J        # was: b_in + (b_eq - b_in)*(1 - exp(-Da))
```
This makes `J·(b_pore − b_in) ≡ F_primary` for every species at any Da — the secondary step can never
remove what the primary step did not supply. The M&C flux law is the model's basis, so it is the
authority; `b_pore` follows from it. (It also makes the `Calcite`-not-in-crust special case for carbon
automatic, since F_primary is already zeroed there.)

**Verified:** the live charge audit leak went +0.35 → ~0 Tmol eq/yr at every trajectory point; on full
runs Ca went 248 → 6 mM (flagship config), median dCa 0.40 → 0.03, and Alk/charge 0.004 → 0.97+.

### 20.4 Sweep history: fast_13, fast_14, fast_15

- **fast_13** (kd=0.07, Kaolinite buffer, w/r absent). **First run was void** — the buffer needs a w/r
  and had none, so PHREEQC failed on 100% of calls, every run fell back to outgassing-only, and the
  sweep came back *looking cleaner than ever* (no CO₂ ceiling, 95% timeouts) while containing **no
  chemistry at all**. Caught only because the JSONs lacked new fields and pCO₂ was pinned at
  `initial_pco2`. This motivated two guards: the buffer self-disables without a w/r, and a
  **`chemistry_void` termination** fires when `fabricated_fraction > 0.5` (an absolute fallback cap
  cannot catch total failure — a run doing no chemistry is too cheap to reach it; fast_13 failed 100%
  on ~203 fallbacks, far below the 5000 cap). A *second* fast_13 (buffer + w/r=600, genuine) gave
  sweep-wide Na (median 6.6×10⁻⁷ → 0.911 mM) but no transition (w/r=600 flattened the feedback).
  **A resume-check trap also bit here:** `RERUN=False` skipped all 931 configs because broken JSONs
  from the first attempt already had a `termination`; the "rerun" returned the old results unchanged.
- **fast_14** (w/r=3, buffer reverted, but *before* the iron fix). Only **16% settled** — dominated by
  the charge/Ca corruption. Also exposed: the out=10 row has **no steady state by construction**
  (carbon input 75 Tmol/yr vs max alkalinity supply 45 Tmol/yr — a genuine supply-limited runaway,
  36 runs learning one bit); the drifting timeouts are **slow-converging on Ca, not stuck**; and Cl
  never equilibrates (residence time ∝ 1/crust, up to ~10 Gyr) — the sweep should seed it analytically
  as `calibrate_earth.py` already does, though Cl is *not* what sets the convergence metric (Ca is).
- **fast_15** (w/r=3, iron fix in). The payoff, paired against fast_14 on the same grid:

  | | fast_14 (before fix) | fast_15 (after fix) |
  |---|---|---|
  | Alk/charge within 10% of 1 | 49% | **97%** (median 1.000) |
  | max ocean Ca | 585 mM | **70 mM** |
  | runs with Ca > 100 mM | 17 | **0** |
  | ion-settled (dC < 5%) | 25/160 | **41/156** |

  **And the transition is visible for the first time:** of the 4 lines with ≥5 settled points, **3
  show the clean June signature** — pCO₂ falls, bottoms out, rises again, with a 22–45 K temperature
  jump (out=0.03/crust=0.1 reverses at S=1.00; out=0.1/crust=10 and out=1/crust=10 at S=0.90). No
  prior sweep showed this cleanly, because the charge/Ca corruption was destroying the settled states
  the transition needs.

### 20.5 What remains

- **The CO₂ ceiling is now the binding constraint on settled fraction** (still ~26%; out_of_domain
  91/156 essentially unchanged — the iron fix was charge/Ca correctness, orthogonal to the ceiling).
  The thermodynamic branch *is* the rising-pCO₂ branch, so it hits the 10-bar wall almost immediately;
  raising or tapering the ceiling (to the maximum-greenhouse peak, §9.2/§13) is the next lever to see
  more of each transition line.
- **Da diagnostic still unresolved** (§8 note): the plotted `Da_alk` (ratio of charge-weighted sums)
  reads two–three decades below where the transition actually is on the crust=10 lines; a rigorously
  flux-weighted `Da_eff` was no better (2/7 vs 1/7 match to observed reversals). The reversal needs
  **both** Da>1 *and* a steep retrograde `d b_eq/dT`, so a single scalar Da cannot locate it — a
  two-condition regime map is the honest diagnostic. Not yet implemented.
- **Na is still ~1 mM**, far short of Cl (~tens of mM). The Na *source* is kinetically capped and no
  equilibrium tuning raises it (20.2); the Kite & Ford route (Na as a swept O(0.1–1) mol/kg inventory)
  remains the pragmatic answer.

### 20.6 Cross-cutting lessons from this arc

- **A cleaner-looking sweep can mean less chemistry, not more.** fast_13's null run had the best
  headline stats of any sweep and contained nothing. Always check `fabricated_fraction` / that the
  ions actually moved before trusting terminations.
- **Two formulas for one physical quantity is a latent bug.** The b_pore/flux mismatch was invisible
  for the entire project because the two functions agree in both asymptotic limits; only a species
  parked at Da≈1 (iron) exposed it. Derive dependent quantities from the authority, don't recompute.
- **A monotonically growing conservation-law violation localises to a persistent per-step flux
  error** — audit the flux terms *and the clamps* separately; here the terms were clean and the clamp
  exposed the upstream b_pore bug.

---

## 21. The reverse-weathering pore-space bug, and the retrograde-solubility confirmation (2026-08-19)

This session started as an investigation of a "temperature falls at high instellation" trend in
`fast_16` and ended by finding two independent defects: a **reporting** bug that corrupted the
recorded `T`/`pH` of every run, and a **model** bug that inverted the sign of the seafloor alkalinity
flux over a broad region of parameter space. Fixing the second produced `fast_17`, the best-behaved
sweep the project has had.

### 21.1 The `self._T` / `self._pH` side-effect corruption

`Planet.dY_dt` sets `self._T` and `self._pH` as **side effects**, and `time_evolve` reads them after
`solve_ivp` returns. But `solve_ivp` calls `dY_dt` for far more than accepted trajectory steps:
Jacobian finite-difference probes (`jac_epsilon = 0.01`, i.e. ±1% state perturbations) and the
internal bisection trials of event root-finding. Whichever call happens to run **last** is what gets
reported — and near a domain wall the climate response is a near-discontinuity, so a routine 1%
perturbation can flip between a real state and the analytic climate model's literal `400.0`
"no equilibrium found" sentinel.

Verified case (S=1.15, out=0.01, crust=0.01): the run reported **T = 400.0 K** while its own final
`P_CO2 = 0.076442 bar` evaluates to **T = 388.99 K**. A +1% Jacobian probe of that `P_CO2` lands
exactly on the 400.0 rail. `P_CO2` itself was never affected — it is read directly from
`sol.y[0, -1]`, not through a side-effect channel.

**Fixes, both in place:**

- `planet.py` — `time_evolve` now re-evaluates `dY_dt(sol.t[-1], sol.y[:, -1])` after `solve_ivp`
  returns, so `_T`/`_pH` correspond to the same state as the reported `P_CO2`. The fallback counters
  (`_chem_fallbacks`, `_chem_ok`, `_dYdt_last_good`) and `_fallback_limit` are snapshotted and
  restored around the call, so this extra evaluation cannot perturb `fabricated_fraction` or trip the
  abort budget. The `fallback_limit` abort path does **not** need this: `self._T` is assigned before
  the `try` block that can raise, so it is already consistent with `self._abort_Y`.
- `plot_results.py` — `_recompute_T(instellation, P_CO2_bar, land_fraction)` recovers T from the
  stored final `P_CO2` via `get_T_surface_analytic`, and is wired into `load_data`,
  `_pore_conditions`, and `_get_mineral_si`; `_diag_from_json` likewise recomputes pH from the
  corrected `T_seafloor`/`P_pore`. This retroactively corrects sweeps already on disk. Albedo
  constants were promoted to module level in `planet.py` (`OCEAN_ALBEDO`, `LAND_ALBEDO`) so the
  plotting path uses the same numbers as `Planet.__init__` rather than duplicating them.

Verified across all 1078 `fast_16` rows: T range 193.8–389.0 K, **zero** rows at the 400.0
contamination value. The specific case from the screenshot (out=3, crust=0.1) went from a false
decline (S=1.15: 400.0 K → S=1.20: 389.9 K) to correctly flat (389.0 K at both).

**This fixed the reporting, but a real T-decline remained** — which is what led to §21.3.

### 21.2 Retrograde solubility — confirmed, after a methodological error

A first attempt to test whether `b_eq` is retrograde evaluated it **at each run's own final state**
and concluded it was *prograde*. That was wrong: along the trajectory `P_CO2` varies by five orders
of magnitude, so this measures `b_eq(P_CO2)` and mislabels it `b_eq(T)`.

Redone properly — `b_eq` at **fixed** `P_CO2` and **fixed** `b_input` (P_pore = 3.1e7 Pa, w/r = 3,
`b_input` = the `fast_16` S=1.00 ocean), scanning T alone:

| charge-weighted `b_alk` | 280 K | 320 K | 350 K | 380 K |
|---|---|---|---|---|
| P_CO2 = 0.01 bar | 237.5 mM | 24.3 | 10.1 | 9.4 |
| P_CO2 = 0.1 bar | 1071.9 mM | 91.3 | 24.2 | 10.7 |
| P_CO2 = 1.0 bar | 1311.5 mM | 417.1 | 94.7 | 28.9 |

**Strongly retrograde**, one to two orders of magnitude over the range. `b_eq[Ca]` is essentially
flat (~4 mM, buffered); the entire effect is `b_eq[Mg]`.

Two consequences:

- Since `Da = kA/(J·b_eq)`, the retrograde `b_eq` sits in the **denominator** and drives Da up
  alongside the Arrhenius rise in `k`. The two push the same way, which is why Da crosses 1 as
  sharply as it does. This supports §20.5's conclusion that a single scalar Da cannot locate the
  reversal — Da>1 and steep retrograde `d b_eq/dT` are not independent conditions here.
- The **primary** flux behaves exactly as the transport-limited law requires: it rises kinetically,
  peaks where Da≈1, then declines. At P_CO2 = 0.01 bar, `F_prim` (Tmol_eq/yr) over
  T = 280/300/320/340/350/360 K is 0.33 / 0.84 / 2.46 / **4.29** / 2.74 / 1.89, peaking at 340 K
  where Da = 3.35. Weathering does *not* keep increasing in the thermodynamic limit.

**Lesson:** when testing the T-dependence of an equilibrium quantity, hold every other state variable
fixed. Sampling along a trajectory conflates all of them, and the dominant covariate wins.

### 21.3 The real defect: reverse-weathering minerals in the pore space

`planet.py:133` read:

```python
self.pore_precipitating_minerals = clay_minerals + reverse_weathering_minerals if reverse_weathering else clay_minerals
```

The RW minerals were **also** applied to the ocean (`rw_ocean_precipitating_minerals`, line 135, used
at line 225) — so they were double-counted, and the pore-space copy was acting on the wrong
reservoir. The pore fluid has just been loaded to **~100+ mM Mg** by primary dissolution; ocean Mg in
these runs is ~0.5 mM. Sepiolite(d) then removes essentially all of it back out.

Per-mineral breakdown at T = 350 K, P_CO2 = 1 bar (Tmol_eq/yr):

```
F_prim alk       = +32.503
  Kaolinite        +0.000
  Goethite         -0.078
  Sepiolite(d)    -34.099    (Mg -17.05, Si -25.57)
  Saponite-Na      -0.022
  Greenalite       -0.078
F_net alk        =  -1.679
```

**`F_net` is a small residual between two large opposing fluxes**, so its sign is numerically fragile.
Mapped over (T, P_CO2) it is jagged and flips negative across a broad region:

```
F_net (Tmol_eq/yr)
   T\pCO2    1e-03   1e-02   1e-01   1e+00   5e+00
     290     0.135   0.178   0.231   0.334   0.485
     310     0.609   1.017   0.813  -1.032  -1.712
     330     0.120   3.116  -1.819  -2.053  -2.066
     350     0.318   1.650   5.652  -1.679  -2.037
     370     0.338   0.482   2.609  23.918  37.474
```

`F_prim` alone is positive and monotonic everywhere on the same grid.

A negative seafloor alkalinity flux is an alkalinity **sink**, hence a CO₂ **source** — a positive
feedback: alk removed → CO₂ up → hotter → still inside the negative patch. This is what generated the
"massive temperature jump" that had been read as the thermodynamic transition. The
`fast_16` S=0.90 / out=0.01 / crust=1 point sits at T_pore = 352 K, P_CO2 = 6 bar, squarely in the
negative region, and is **not a steady state**: `dlnP/dlnt = 0.890` at the 2 Gyr cutoff (still
climbing, 3.2 → 6.0 bar over the second half of the run) with **31% of its derivative evaluations
fabricated**.

`mineral_info.py:96–100` already carried a warning about exactly this pathology — measured
`F[Alk] +2.5 → -207 Tmol/yr` — but scoped it to `lt_equilibrium_buffer_minerals` (inside `get_b_eq`),
concluding the RW minerals were safe in "the post-equilibrium precipitation step". They were not:
that step is still the *pore* fluid, not the ocean.

**Fix:** `self.pore_precipitating_minerals = clay_minerals`. RW remains an ocean sink only, where it
was already being applied.

### 21.4 fast_17 — the sweep with the fix

Coarse grid (10 S × 4 outgassing × 4 crust, depth 3000, rw=True, mt=1350), an exact subset of
`fast_16`'s grid, so 150 runs match one-to-one:

| | fast_16 (RW in pore) | fast_17 (no RW in pore) |
|---|---|---|
| out_of_domain | 98 | **87** |
| timeout (survived to end) | 50 | **62** |
| fallback_limit aborts | 1 | **0** |
| converged (\|dlnP/dlnt\| < 0.05) | 33.7% | **36.4%** |
| median \|dlnP/dlnt\| | 0.761 | **0.529** |
| fabricated, mean | **0.023** | 0.050 |
| runs > 10% fabricated | **9** | 20 |

**T(S) is monotonic in every (out, crust) pair** — the jump-then-decline artifact is gone. With the
pore assemblage reduced to Kaolinite+Goethite the secondary term is ≈ −0.08 Tmol/yr instead of −34,
so `F_net` can no longer change sign. At out=1/crust=1 `fast_16` had every point pinned at the 10-bar
ceiling; `fast_17` converges S=0.7–1.1 as real timeouts.

The two regimes are now cleanly separable:

*Kinetic line* (out=0.03, crust=1) — Da never reaches 1, textbook negative feedback, well converged:

```
S:     0.6    0.7    0.8    0.9    1.0    1.1
Da:  0.001  0.001  0.001  0.004  0.038  0.662
pCO2: 0.563  0.254  0.086 0.0157 9.3e-4 1.9e-4     |slope| < 0.03
T:     286    288    289    292    296    315 K
```

*Transition line* (out=0.1, crust=0.1) — Da crosses 1 near S≈1.05, and pCO₂ reverses:

```
S:     0.7    0.8    0.9    1.0    1.1     1.2
Da:  0.008    nan  0.057  0.809  2.028   590
pCO2: 0.716  0.358  0.192  0.551  1.156  0.042
T:     308    311    318    344    355    389 K
                      ^min   ^^^^^^^^^^ reversal
```

The pCO₂ minimum sits immediately before Da→1 and pCO₂ rises with instellation through the crossing —
the June signature (§13), now arising on its own rather than being manufactured by the Sepiolite
cancellation.

### 21.5 What remains

- **Fabricated fraction went up**, not down (mean 0.023 → 0.050; runs >10% went 9 → 20). Unexplained
  and worth tracking down — removing a precipitation step should not make the chemistry harder.
- **The reversal line is not settled.** Slopes on S=0.9/1.0/1.1 are −0.07 / −0.23 / −0.08, and the
  S=0.9 point is **26% fabricated**. Needs a targeted rerun at finer S spacing and longer integration
  before it can carry a paper claim.
- **out=3 (both crusts) and out=1/crust=0.1 remain entirely pinned at the 10-bar ceiling**, unchanged
  from `fast_16`. The ceiling is still the binding constraint (§20.5).
- **The S=1.2 column is the hot domain wall (389 K) everywhere**, not a steady state; it should not be
  read as a trend endpoint.

### 21.6 Cross-cutting lessons

- **A flux that is the difference of two large opposing terms has an untrustworthy sign.** `F_net`
  was +32.5 − 34.1; nothing about the magnitudes was wrong, but the residual carried all the physics
  and none of the precision. Check dominant-term cancellation before trusting a small net flux.
- **A warning comment scoped to one call site can apply to another.** The RW-strips-Mg pathology was
  documented in `mineral_info.py` and still shipped, because the note reasoned about `get_b_eq` and
  the same minerals were live in the post-equilibrium step. Ask *which reservoir* a sink acts on, not
  which function it is called from.
- **Diagnostics written as side effects will be read at the wrong time.** `self._T` was correct on
  every accepted step and wrong in the output, because the solver calls the RHS for its own purposes.
  Recompute reported quantities from the accepted state.
- **The user's physical intuition caught two wrong conclusions this session** (that RW could remove
  carbon — it consumes alkalinity and is a CO₂ *source*; and that `b_eq` was prograde). Both were
  resolved by measurement, not argument. When a chemically-trained objection contradicts a model
  result, instrument the model.

---

## 22. Earth calibration, the `alpha` anchor, and the Coogan LT fluxes (2026-08-19 → 08-20)

Run-up to a large sweep campaign: the goal was simply to get calibrated values of `kd_mg_ht`, `k_na`
and `alpha` before launching. Two of the three came out cleanly and are **alpha-independent**. The
third, `alpha`, turned out not to be identifiable from Earth at all, and chasing its literature anchor
exposed a **composition mismatch between the model's low-temperature seafloor flux and the
literature's**: the model delivers the right *amount* of alkalinity as the wrong *ion*.

**Nothing was patched into the working copy.** The `alpha` decision is still open — see §22.9.

### 22.0 Code state going in

Two things had moved past §11/§15 before this session and are worth recording:

- **The CO₂ ceiling is implemented** (§15 item 1, which §21.5 still lists as the binding constraint).
  `time_evolve` now computes `P_CO2_HI, _ = maximum_greenhouse(self.instellation, self.albedo)` once
  and uses it for the `co2_high` margin (`climate/analytic.py:111`, `planet.py:421`). No sweep has
  been run with it — `parameter_sweep.py` still points at `fast_17`.
- **`calibrate_earth.py` is no longer broken.** §11's ⚠️ (`f_bio=`) and §15 item 10 are stale. It has
  been rewritten: Phase 2 removed (the model is abiotic by construction), a charge-balanced modern
  seawater seed, Cl seeded analytically, and a bounded `least_squares` over
  `log(K_na, alpha, KD_mg)` replacing the hand-rolled ratio loop that used to oscillate.

### 22.1 The calibration run

28 evaluations, terminated on its own `xtol`/`ftol` (cap was 60), best cost 0.0024. It reproduced an
earlier same-day run to the digit, so the result is deterministic.

```
planet.py:      K_NA_CONT_REMOVAL = 3.904380e-03   (was 2.194806e-03)
                KD_MG_HT          = 1.898657e-02   (was 7.0e-02)
                K_CL_SUBDUCTION   = 1.373251e-04   (unchanged — analytic, already in code)
weathering.py:  ALPHA_REF         = 0.908383       (was 1.43)
```

At Earth (land = 0.3, 3700 m, S = 1, out = 1, crust = 1): `converged`, **T = 294.3 K, pH 7.77,
pCO₂ 679 ppm**, `|dlnP/dlnt| = 0.001`, zero fabricated derivatives.

| | Na | Ca | Mg | Alk | C |
|---|---|---|---|---|---|
| sim (mM) | 464.8 | 10.60 | 54.90 | 2.86 | 2.63 |
| error | **−0.9%** | **+2.9%** | **+4.0%** | +24% | +25% |

Ca:Mg split 0.990 (1.0 = correct), so the three-target fit is genuinely solved rather than traded off.
Alk/C ~25% high with pCO₂ at 679 vs 280 ppm is the expected abiotic offset — there is no biogenic
carbonate pump in the model.

**What moved:** dropping `kd_mg_ht` 0.07 → 0.019 is what unlocked Ca; `K_na` rose 2.19e-3 → 3.9e-3.
**The calcite snap is visible in the trace:** wherever the solver probed `kd` too low (evals 007, 011,
017, 022, 027) Ca collapses to ~0.2 mM and Mg/Alk run away to 100+ mM. This is the bistability the
script's own docstring documents — Ca and alkalinity trade along Ca·CO₃ = Ksp — and it is why a
gradient method with a trust region is required and independent per-species ratio updates cannot work.

> **This supersedes §18's conclusion that Earth cannot constrain `kd_mg_ht`.** That measurement (HT
> exchange = 1.5% of Mg removal, 0.6% of the Ca source) was made when `make_b0` violated
> Alk = ION_CHARGE·b by 592.9 mEq/kg, which pinned Ca near 0.20 mM by calcite supersaturation and made
> every knob look inert. With the seed fixed, Ca responds sharply to `kd`. §18's *method* stands; its
> "Earth is blind to this parameter" conclusion does not.

### 22.2 `alpha` is not identified by the Earth fit

Paired evaluations differing only in `alpha` return **identical oceans to four significant figures**:

```
eval 008  K_na=3.824e-3  alpha=0.9230  kd=1.861e-2  ->  Na=475.1  Ca=5.06  Mg=55.66   cost 0.5068
eval 010  K_na=3.824e-3  alpha=0.9083  kd=1.861e-2  ->  Na=475.1  Ca=5.06  Mg=55.66   cost 0.5071
eval 024  K_na=3.904e-3  alpha=0.9230  kd=1.899e-2  ->  Na=464.8  Ca=10.60 Mg=54.91   cost 0.0025
eval 026  K_na=3.904e-3  alpha=0.9084  kd=1.899e-2  ->  Na=464.8  Ca=10.60 Mg=54.90   cost 0.0024
```

The reported 0.908383 is wherever the trust region stopped, not a measured value. Two causes: at Earth
pore conditions the seafloor flux is **transport-limited**, so it saturates in `alpha`; and at
`land_fraction = 0.3` **continental weathering supplies most of the Ca+Mg** anyway (§18's budget table).

Confirmed by two fixed-`alpha` refits (§22.8): across a **220× change in `alpha`**, `K_na` moves 3%
and `kd_mg_ht` 2.5%. **Those two constants can be adopted regardless of how `alpha` is settled.**

### 22.3 The 1 Tmol/yr anchor is measured in a configuration the model never runs

`weathering.py`'s `ALPHA_REF` docstring says it is "calibrated to give ~1 Tmol/yr seafloor Alk flux at
modern Earth pore conditions with modern seawater composition". The diagnostic that produces that
number calls `get_weathering_flux` with `precipitating_minerals=[]` — **primary dissolution only**. At
Earth pore conditions with modern seawater, at the `alpha = 19.92` that anchor implies:

```
Fe   0.438 Tmol/yr  -> charge +0.875   (88% of the "alkalinity")
Mg   0.051          -> charge +0.102
Al   0.008          -> charge +0.023
Ca   0.0002         -> charge +0.000
```

**88% of it is dissolved Fe²⁺ from Fayalite.** But `dY_dt` passes
`pore_precipitating_minerals = ['Kaolinite', 'Goethite']`, and Goethite exists precisely to strip that
Fe — which under the §20.3 charge convention takes 2 eq of alkalinity with it.

| alpha | F_alk primary-only | F_alk net (clay pore list) |
|---|---|---|
| 0.908 (fitted) | 0.093 | **0.005** |
| 19.92 (anchor) | **1.000** | **0.103** |
| 202.25 | 2.590 | **0.990** |

So `alpha = 19.92` delivers 1 Tmol/yr of primary dissolution and **0.10 Teq/yr net to the ocean**.
Reaching a *net* 1 Teq/yr needs `alpha ≈ 202`.

### 22.4 Coogan & Dosso: the literature low-temperature fluxes

**Coogan & Dosso (2022)**, *GCA* 329, 22–37, §2.4 + Table 1 parameterizes seafloor weathering as

```
J_Mg^LT  = M_L · f_basalt · X_Mg^SW · Kd_Mg          Kd_Mg  = 5 (±3),  prior bounds 0–100
J_Alk^LT = M_L · f_basalt · ALK_LT                   ALK_LT = 1.7 (±1.4) eq/kg,  bounds 0–10
J_Ca^LT  set by charge balance                       M_L = 4.42e12·A_oc kg/yr, f_basalt = 0.12 (0.01)
```

Evaluated at modern values:

| low-T seafloor weathering | flux | note |
|---|---|---|
| Mg | **−0.14 Tmol/yr** | sink |
| Ca | **+0.59 Tmol/yr** | source |
| Alk | **+0.90 Teq/yr** | source; 1σ range **0.16–1.64** |

Charge accounting is the same convention as our `ION_CHARGE·flux`: +1.18 eq from Ca, −0.28 eq from Mg
uptake, net +0.90. **Cross-checked** against **Coogan & Dosso (2026)**, *EPSL* 677, 119811, Fig. 6:
Ca\* = Ca − (Alk/2)·f_calcite = 0.59 − 0.41 = **0.18**, matching the ~+0.2 Tmol/yr plotted for seafloor
weathering.

**High-temperature** (2026 Fig. 6): Ca\* ≈ **+0.5 to +0.9 Tmol/yr**; §2.3 assumes *all* seawater Mg and
alkalinity are consumed by fluid-rock reaction with Ca released for charge balance — structurally
identical to our `F_ht_exchange`.

**Implication for `kd_mg_ht`.** Stripping all Mg from kamino's hydrothermal flux would give
**105.6 Tmol/yr** Ca at Earth (i.e. `J_total` ≈ 2e15 kg/yr, the *total* axial + flank circulation), so
Coogan's 0.5–0.9 implies **`kd_mg_ht` = 0.005–0.009**. The fitted 0.019 is 2–4× high; **§18's 0.07 is
~10× high**. Notably this lands almost exactly on §18's own "focused venting only" bound of 0.003 —
the bound §18 considered and chose against.

**No constraint on `k_na`.** Coogan explicitly excludes it: *"The alkalinity fluxes associated with
other major ions (largely Na, SO₄ and K) are balanced"* (2022 §2), and the 2026 paper tracks a
univalent-ion alkalinity proxy rather than Na and K. Keep the Earth-concentration fit.

### 22.5 The Mg-removal mechanism in Coogan — a partition coefficient, not a saturation

Mg is removed by **uptake into clay minerals, parameterized as an empirical partition coefficient**,
in three places:

1. **Low-T off-axis basalt alteration** — `Kd_Mg` = 5 is Table 1's *"partition coefficient for Mg into
   clays"*, the Mg concentration of altered lava relative to seawater, calibrated on altered ocean
   crust and the Troodos ophiolite (Coogan & Gillis 2013, 2018; Coogan et al. 2019). **The temperature
   dependence is not in `Kd_Mg`** but in `f_basalt`, the *fraction of lavas altered*, through an
   Arrhenius term with `E_a^SFW` ≈ 97 kJ/mol. ≈ 0.14 Tmol/yr.
2. **Marine sediment diagenesis — explicitly "reverse weathering"** — the *same* `Kd_Mg` applied to
   `M_dia` = 190e9 kg/yr of authigenic minerals (*"Both Li and Mg are largely taken up in clay minerals
   during diagenesis, as they are during low-temperature basalt alteration, and so the same partition
   coefficients are used"*), anchored on Dunlea et al. (2017)'s 0.02 Tmol/yr deep-sea authigenic Mg
   sink. Alkalinity consumption follows by charge balance. ≈ 0.05 Tmol/yr.
3. **High-T axial hydrothermal** — complete stripping per pass, magnitude set by fluid flux ∝ crust
   production.

**Why the functional form matters more than the number.** A partition sink is *first-order in seawater
Mg and capped by the mass of rock altered*; it cannot strip a fluid to zero. `get_precipitation` fires
whenever SI > 0 and can. That is exactly the difference between Coogan's LT term removing **24%** of
the Ca-derived alkalinity and fast_16's equilibrium Sepiolite(d) removing **105%** of it (§21.3). The
§21.3 fix (drop RW from the pore list) and the fast_16 bug are the two ends of a spectrum whose
literature-supported middle is a *rate-limited partial* uptake.

Second structural point: **in Coogan's LT term the Ca release is not a dissolution calculation at
all** — `ALK_LT` and `J_Mg` are independent parameterizations and Ca is the charge-balance residual.
There is no equilibrium step that can pin Ca at zero, which is what ours does (§22.6).

**kamino already contains a Coogan-style partition sink:** `F_ht_exchange = kd_mg_ht · [Mg] · J_total`
has precisely the form `Kd × [Mg] × throughput`. It is the **LT path that lacks one**.

### 22.6 The composition mismatch — measured

`b_eq[Ca]` tracks the input ocean Ca to within **0.01 mM at every ocean Ca tested (0.1 → 30 mM)**, and
stays pinned at ~10.3–10.5 mM across T = 286–350 K and pCO₂ = 280 ppm → 1 bar. **The LT equilibrium
releases essentially no Ca, ever** — the Ca-bearing primaries never dissolve. Meanwhile Mg is a strong
source (b_eq[Mg] up to 704 mM at pCO₂ = 1 bar).

At the `alpha = 180` that reproduces Coogan's 0.90 Teq/yr:

| | Mg | Ca | Alk |
|---|---|---|---|
| Coogan & Dosso 2022 | **−0.14** | **+0.59** | +0.90 |
| kamino | **+0.437** | **+0.001** | +0.900 |

Right total, wrong ion, wrong sign on Mg. **`alpha` scales magnitude, not composition, so no value of
it can fix this.** Root cause is the Diopside swap (§5, §13) plus the absent pore Mg sink (§21.3).

**Independent confirmation from the fits:** at `alpha = 202` the concentration fit degrades
*specifically in Mg* (+32.9%), because raising `alpha` scales a Mg **source**. Had the LT sign been
right, pushing `alpha` toward the literature alkalinity flux would have *improved* the Mg fit. `kd`
cannot absorb the excess either — reducing Mg through the HT exchange runs into the calcite snap that
kills Ca.

### 22.7 Incidental: `WATER_ROCK_RATIO_LT` is inert at Earth pore conditions

`b_eq` is **identical** for w/r = 0.1, 1, 3, 10, 100 and `None` (T = 286 K, pCO₂ = 280 ppm, modern
seawater, `dissolve_only=True`, empty buffer list). The equilibrium dissolves less rock than even the
most rock-starved setting provides, so the constraint never binds.

This promotes §19.7's warning ("w/r = 2.0 returned results identical to `none` — the parameter may not
bite at low values") from a suspicion to a measurement *at these conditions*. §20.2's choice of w/r = 3
was made at the hotter near-Da=1 states, where `d ln charge(b_eq)/dT` did respond strongly, so the
parameter is presumably still live there — but it does nothing at Earth.

### 22.8 Three calibrations

`alpha` fixed, `(K_na, kd_mg_ht)` refitted at each — same seed, targets, 4 Gyr integration:

| alpha | K_na | kd_mg_ht | Na | Ca | Mg | cost | net LT Alk | LT Mg |
|---|---|---|---|---|---|---|---|---|
| **0.908** (free fit) | 3.904e-3 | 1.899e-2 | −0.9% | +2.9% | +4.0% | 0.0024 | 0.005 | +0.002 |
| **19.92** (primary anchor) | 3.894e-3 | 1.894e-2 | −2.1% | +7.8% | +8.1% | 0.0121 | 0.184 | +0.090 |
| **202.25** (Coogan net Alk) | 3.765e-3 | 1.852e-2 | −7.0% | −5.8% | **+32.9%** | 0.0896 | 1.215 | +0.586 |
| Coogan & Dosso 2022 | — | — | — | — | — | — | **0.90** | **−0.14** |

**Adopt regardless of the `alpha` decision:**

```
K_NA_CONT_REMOVAL = 3.9e-3   (range 3.77e-3 – 3.90e-3)
KD_MG_HT          = 1.9e-2   (range 1.85e-2 – 1.90e-2)
```

### 22.9 What remains — the open decision

**`alpha` is undecided.** Four options, none free:

1. **19.92** — concentrations still good (Na −2.1%, Ca +7.8%, Mg +8.1%) and it matches the ~1 Tmol/yr
   figure under the measurement convention `weathering.py` already documents. Middle of the range.
2. **202** — closest to Coogan's net 0.90 Teq/yr, but Mg overshoots +33% and the ion is still wrong.
3. **0.908** — best Earth fit, but the net seafloor alkalinity flux is 0.005 Teq/yr, ~180× below
   Coogan. Effectively switches seafloor weathering off, which is fatal on land-free worlds.
4. **Fix the composition first** — give the LT path a first-order Mg uptake of Coogan's form
   (`Kd_Mg × [Mg] ×` altered rock mass) plus an alkalinity/Ca source term, then recalibrate. Makes
   `alpha` a meaningful parameter instead of one fitted to a flux carried by the wrong ion.

**The Earth calibration structurally cannot settle this.** `alpha` is nearly free on Earth *because*
continental weathering dominates at land_fraction = 0.3; on the land-free worlds the sweeps are about,
seafloor weathering is the only silicate sink and `alpha` carries the entire thermostat. This is the
same shape of finding as §18 reached for `kd_mg_ht`: a parameter that must be set on land-free worlds
or from first principles.

Also outstanding: `kd_mg_ht` now has **three** disagreeing anchors — Earth fit 0.019, Coogan HT
0.005–0.009, §18 first-principles 0.07.

### 22.10 Cross-cutting lessons

- **A diagnostic anchor must be evaluated in the configuration the model actually runs.** The
  `ALPHA_REF` anchor omitted the pore precipitation step, so it was matching a flux that is 88% Fe²⁺
  which Goethite removes on the very next line. Ten-fold error, invisible for the life of the project.
- **Check which ion carries a matched flux, not just its magnitude.** Kamino and Coogan agree on
  0.90 Teq/yr of seafloor alkalinity and disagree on essentially everything about how it gets there.
- **A parameter can be unidentifiable in the calibration and load-bearing in the application.**
  `alpha` is invisible to Earth and carries the thermostat on water worlds. Fitting it on Earth and
  reporting a 6-digit value would have been false precision.
- **Partition-coefficient sinks and saturation-driven sinks are not interchangeable.** One is capped by
  rock supply and first-order in the ocean; the other goes to completion. Substituting the second for
  the first is what produced the fast_16 alkalinity sign inversion.
- **A calibration that reproduces to the digit is worth the second run.** The re-run matched the
  previous session exactly, which is what made it safe to treat the constants as deterministic rather
  than as a solver artifact.

---

## 23. The crust-composition pipeline: pMELTS, and Nepheline (2026-08-20)

> ⚠️ **Partly superseded by §24.** pMELTS is replaced by MAGEMin (§24.1) and
> `make_crust_compositions.py` by the `.jl` version. The Nepheline work (§23.7) STANDS and is
> independently confirmed by MAGEMin. §24.6 records that the Mg/Si → crust question was already
> published by this group, which reframes §23.1–23.6.

Started from a sweep result — mantle temperature barely moves the chemistry, unlike June's crust
composition series — and ended by replacing the reason. The T_p and Mg/Si axes were weak because
**the mapping from mantle composition to melt was wrong**, and because the CIPW norm clipped away
whatever signal survived. Both are now fixed: pMELTS is running (superseding §5), and Nepheline is
in the database.

### 23.1 Why the mantle-temperature sweep looked flat (§22 follow-up)

The `fast_18` T_p sweep (out=0.1, crust=1, alpha=2, T_p 1350→1600) gives **T −5.8 K, salinity −3.4%,
Mg 0.1%, Alk 1.4%** — against June's slides 37–38, where 44–51% SiO₂ moved salinity from ~20 to
~100 g/kg. Three reasons, in order of importance:

1. **The compositional lever is ~10× smaller.** Cation charge delivered per kg of rock:
   June's named series spans **1.33×** (25.79 → 19.43 eq/kg); T_p 1350–1600 spans **1.03×**; the
   full 2-D (T_p, Mg/Si) pipeline grid spans only **1.14×**.
2. **The pipeline cannot make June's end-members.** `komatiite_42` is 100% olivine with *zero*
   plagioclase and zero Al. Every PRIMELT+CIPW composition carries 25–36% normative Anorthite,
   because a primary mantle melt always has substantial Al₂O₃ and CaO. June varied the assemblage;
   T_p only re-weights a fixed one.
3. **An internal cancellation.** Hotter mantle raises Ca+Mg (4.96 → 7.03 mol/kg, +42% — a lever
   comparable to June's 1.38×) but simultaneously drops Na (0.78 → 0.32) as normative Albite falls
   0.20 → 0.08. In *total* cation charge the two nearly cancel. June's komatiites had Na ≈ 0, so
   their Ca+Mg change passed straight through.

Downstream this is reinforced: ocean Mg sits at **18.9 mM for S ≥ 0.5 at every T_p, identical to
three significant figures**, because steady-state Mg ≈ (LT source)/(`kd`·J) and the source moves 3%.

### 23.2 The CIPW clipping mechanism

`mineral_composition` returned *identical* crusts for Mg/Si 0.5–0.8, and again for 1.1–1.5. The
cause is in `cipw_norm` step 5, the only silica-sensitive step: after Albite/Anorthite/Diopside/
Fayalite are allocated, the Mg remainder is split between Enstatite (1 Si per Mg) and Forsterite
(0.5 Si per Mg). **The assemblage can therefore only absorb silica over a factor of 2**, so the
responsive window is exactly `Si_rem/Mg_rem ∈ (0.5, 1.0)` — measured as **Mg/Si ≈ 0.87–1.04** at
T_p = 1350, about ±10% around Earth.

Outside it, two asymmetric failure modes:

- **Oversaturated:** excess silica hits step 6 and is discarded (`emit_quartz=False`). At Mg/Si 0.5
  that is **0.48 mol of SiO₂ thrown away**. Wasteful, mass-safe.
- **Undersaturated:** `Si` goes **negative** and step 6 only fires on `Si > 1e-9`, so the deficit was
  **silently dropped** — the norm assigned more SiO₂ to minerals than the rock contained
  (−0.32 mol at Mg/Si 1.5). Not mass-conservative. Those were never "high-Mg/Si crusts"; they were
  "the most olivine-rich crust the norm can build", relabelled.

Root cause is a missing phase: standard CIPW handles undersaturation by desilicating feldspar
(albite → nepheline), and there was no feldspathoid in the database.

> **`emit_quartz=True` is not the fix, and its stated justification was obsolete.** The docstring
> kept Quartz out of the crust because it "doubles as an HT precipitation buffer" — but **nothing
> sets `high_temperature=True` anywhere**; `f_HT` is stored and never read, and `ht_secondary_minerals`
> is unreachable. Testing it anyway: Quartz is *not* in `primary_minerals`, so it precipitates,
> clamping pore b_eq[Si] 12.4 → 0.39 mM (32×) with Mg 65 → 120 mM and Al 10 → 317 µM in tow — a
> low-temperature quartz buffer the model deliberately avoids elsewhere (`silica_minerals` is
> `SiO2(am)`, §19.7). Made dissolve-only it is **bit-identical to `emit_quartz=False`**, because the
> pore fluid is already 32× supersaturated in quartz (SI +1.50) so it cannot dissolve, and it carries
> no cation charge. Its only effect is diluting the reactive phases by its weight fraction (−12.7%
> on `k_charge` at Mg/Si 0.7). Not a lever in either configuration.
>
> ⚠️ **Do not delete `hydrothermal.dat` when cleaning up the dead HT path.** It is still load-bearing
> for *stoichiometry*: `STOICHIOMETRY_SOURCE` takes `Anorthite` from it deliberately (Al³⁺ basis,
> the one Kaolinite is written on, so Anorthite + Kaolinite nets to the textbook Alk = 2). The solver
> path is dead; the parsed database is not.

### 23.3 The Mg/Si → melt mapping was wrong (measured)

`oxide_composition` implemented Mg/Si as `oxides['SiO2'] *= MG_SI_REF / mg_si_ratio` followed by
renormalisation to 100 wt%. Consequences, all measured:

| defect | reality | old mapping |
|---|---|---|
| CaO vs Mg/Si | ~flat | **9.85 → 14.22 wt% (+44%)** — pure renormalisation artifact |
| Al₂O₃ | modest | 14.11 → 20.36 wt% (20 wt% is a gabbroic cumulate, not a melt) |
| FeO | falls with Mg/Si | rises (renormalisation) |
| endpoints | 46–53 wt% SiO₂ | **38.8–57.6 wt%** — sub-komatiite to andesite; not primary melts |

And the input is mislabelled: the docstring says "the planet's mantle molar Mg/Si", but it is a
*relative multiplier* on Earth's 1.23, so `mg_si_ratio` 0.7–1.5 is really mantle Mg/Si ≈ 0.86–1.85.

**The norm was clipping because it was being handed compositions that shouldn't exist.**

### 23.4 A retracted result — HEX1 and HEX2 are different axes

**Brugman et al. (2021)** melted two hypothetical exoplanet mantles. It is tempting, and wrong, to
treat them as two ends of a Mg/Si axis:

- **HEX1** — Earth's undepleted mantle adjusted to molar **Mg/Si = 1.42** (Earth 1.23). A Mg/Si
  end-member.
- **HEX2** — adjusted to molar **Ca/Al = 1.07** (Earth 0.72). A **Ca/Al** end-member; its Mg/Si of
  0.95 is a side effect.

Their 68% melt-Al₂O₃ difference (9.31 vs 15.65 wt%) is the *Ca/Al* axis. An interpolation between
them in Mg/Si would have hard-coded a spurious Al trend into the crust pipeline. **A conclusion
derived that way — "k_charge falls 0.83× with Mg/Si" — was retracted.** (It later reappeared, 0.82×,
from a clean orthogonal axis in §23.6, but that is a coincidence of magnitude, not a rescue.)

Second trap in the same comparison: these are experimental partial melts at fixed P–T and melt
fraction, while PRIMELT gives primary melts. Comparing directly is not like-for-like — the melt-MgO
difference is largely melt-fraction driven. **Match melt fraction, or the comparison is meaningless.**

### 23.5 pMELTS works — §5's rejection is superseded

§5 recorded pMELTS/alphaMELTS as rejected for install complexity and GSL soname ABI issues. That is
no longer true and the blocker was a *single* missing shared library:

- **alphaMELTS for Python 2.3.1 is already installed** at `/data/pt426/alphamelts/`.
- Its `libalphamelts.so` needs `libgsl.so.27`; the system has GSL 2.6 (`libgsl.so.25`). The
  pre-existing `/data/pt426/melts-deps/lib/libgsl.so.27` is unusable — built against GLIBC 2.35,
  this box has 2.34.
- **Building GSL 2.7.1 from source takes ~1 minute** (`./configure && make -j8 && make install`)
  and works cleanly. That is the whole fix.

> ⚠️ **Never symlink the system GSL 2.6 as `libgsl.so.27`.** All 38 `gsl_` symbols alphaMELTS needs
> are present, so it loads *and runs* — then corrupts the heap (`realloc(): invalid pointer`,
> SIGABRT). Symbol presence is not ABI compatibility, and the failure mode is silent numerical
> corruption before it is a crash.

Three API facts that cost time and are not documented anywhere obvious:

1. **`MELTSdynamic(2)` is pMELTS** (not 4), and the calculation mode can only be set **once per
   process**.
2. **`calcEquilibriumState` does not raise on failure.** It prints, leaves the library unrecoverable
   ("Could not re-initialize MELTS library after failure"), and returns — so every later calculation
   in that process silently degrades. Read **`engine.status.failed`** / `status.message`, and run one
   process per grid point.
3. **fO₂ must be buffered** (`setSystemProperties("Log fO2 Path", "FMQ")`). Without it pMELTS drives
   the Fe₂O₃ liquid component negative and returns a **36 wt% SiO₂** "melt" — plausible-looking and
   wrong.

Conditions must be marched with `addNodeAfter()`, not by mutating one engine.

### 23.6 `make_crust_compositions.py`

New: `src/kamino/data/make_crust_compositions.py`. pMELTS adiabatic decompression melting over a
(T_p, Mg/Si) grid → `crust_compositions.csv`, to replace the PRIMELT spreadsheet + SiO₂-rescale proxy
that `oxide_composition` interpolates.

- Mantle Mg/Si is set by **re-splitting MgO/SiO₂ at fixed total mass**, leaving Al/Ca/Fe/Na at their
  pyrolite values — so the axis is genuinely orthogonal, which it was not before.
- Path: 3.0 → **1.0 GPa** isentropic (runMode 3), batch melting, so the final liquid is the pooled
  primary melt. **`P_END` is the melt segregation pressure, not the base of the crust** — this is
  batch melting, so the liquid re-equilibrates at every step, and carrying it to 0.2 GPa gives
  **57 wt% SiO₂ at F = 0.29**, an over-equilibrated andesite. This is the most consequential choice
  in the script.

**Validation.** At T_p = 1350, Mg/Si = 1.25 it gives F = 0.075, SiO₂ 49.92, MgO 11.80, Al₂O₃ 17.01,
CaO 9.43 — against PRIMELT's 48.76 / 11.27 / **17.05** / 11.91 at the same T_p. Independent
thermodynamics landing that close to the existing interpolation is good evidence for both. And the
artifact is gone: across Mg/Si 1.0→1.5 at fixed T_p, **CaO moves +9% and Al₂O₃ +10%**, versus the
proxy's spurious +44% each.

**Not yet generated for production** — the grid extent is a decision, and Nepheline had to land first.

### 23.7 Nepheline added

The pMELTS melts are **genuinely nepheline-normative** — every one at T_p = 1350 clipped
silica-deficient, because low-F batch melts are alkali-rich (Na₂O 3.6–3.8 wt% at F = 0.075). So this
was a real missing phase, not an artifact of the bad mapping. Checked first: **Nepheline alone clears
the deficit on every melt** (0.003–0.048 mol converted of 0.07–0.20 mol albite available), so no
leucite (K is dropped) or larnite (absent from the database) is needed.

Changes:

- **`hybrid_ocean.dat`** — Nepheline PHASES block from `llnl.dat` (SUPCRT lineage, as for the other
  primary silicates), rewritten onto this database's H4SiO4 basis via SiO₂ + 2 H₂O = H4SiO4
  (log_k 0), exactly as the Albite entry was: `NaAlSiO4 + 4 H+ = Al+3 + Na+ + H4SiO4`,
  log_k 13.8006, plus delta_H, the analytic expression and Vm = 54.16.
- **`mineral_rates.py`** — `nepheline_k` from the Kinec_v3 Hermanska parameters
  (Aa 5e7 / Ea 63 kJ, An 0.1 / En 58.5 kJ, Ab 7.5e-5 / Eb 58 kJ, na 1.0, nb −0.4), + `K_FUNCTIONS`.
- **`mineral_info.py`** — `MINERAL_MOLAR_MASS['Nepheline'] = 142.05`; added to `primary_minerals`
  so `dissolve_only` applies (a phase this soluble must not be allowed to precipitate).
- **`crust_composition.py`** — step **5b**, the desilication cascade: albite → nepheline releases
  2 mol Si each; a residual deficit now *warns* instead of being silently discarded.
- **`make_database.py`** — added to `ADDED_PHASES` and `RATES_MINERALS` so the documented build path
  stays complete.

> ⚠️ **The PHASES block must go INSIDE the existing PHASES section** (before `PITZER`), not appended
> at the end of the file. PHREEQC stops reading at `END`, but `chemistry.parse_stoichiometry` does
> not — so a block appended after `END` parses fine in Python, appears in `stoichiometry` and
> `minerals`, and PHREEQC still reports "Phase not found in database". Caught only by an end-to-end
> `get_b_eq` call.

**Why this is the right phase, and why it matters beyond the norm.** Albite's solubility product
carries [H₄SiO₄]³; nepheline's carries it to the **first** power, and its log K is 11 orders higher
(13.80 vs 2.76). So `SI(Neph) − SI(Alb) = −2·log[SiO₂] − 11.04`, which at the pore fluid's 12.4 mM
silica is **−7.2 log units** — they would only be equally saturated below **3 µM**. Nepheline is
therefore immune to the exact mechanism that kills Na in this model (§19.4: b_eq[Na] = 7.6×10⁻⁶ mM,
a 69,000× suppression caused by the Si³ term). It is also ~17× kinetically faster than albite at
343 K / pH 6.6.

**Verified end-to-end:** the cascade clears every deficit with no warnings, PHREEQC equilibrates
nepheline-bearing crusts, and across Mg/Si 1.0→1.5 at T_p = 1350 **`k[Na]` rises 2.8×**
(7.4e-11 → 2.0e-10) with Nepheline reaching 0.047 weight fraction. `k_charge` falls 0.82×.
Regression-checked: `basalt_49` and `mineral_composition(1350)` are unchanged.

### 23.8 What remains

- **Generate the production table** and wire `oxide_composition` to interpolate the CSV instead of
  the PRIMELT spreadsheet. Grid extent is a decision; **Spaargaren et al. (2023)** put the
  mineralogically meaningful mantle range at **Mg/Si 0.8–1.6** (quartz-bearing to
  ferropericlase-bearing), Earth 1.23.
- **Re-check the responsive window** with the cascade in place — it should widen well beyond the old
  0.87–1.04, but the new limit has not been measured.
- **Watch for Na flooding.** A kinetic-limit estimate from a 10 wt% nepheline crust gives ~7 Tmol/yr
  Na against the current model's ~0.006, implying steady-state ocean Na of 250–1260 mM at
  `k_na` = 0.004–0.02. That brackets Earth's 469 mM and sits in Kite & Ford's O(0.1–1) mol/kg
  waterworld range — but §19.3 measured that just 165 mM of Na drops pCO₂ by **385×**, so this can
  over-correct into a frozen soda ocean exactly as alpha = 20 did through a different ion. `k_na`
  was calibrated on Earth, where Na comes from continents, and will need re-deriving.
- **Batch vs fractional melting.** PRIMELT models accumulated fractional melts; the script does batch.
  Worth revisiting if the melt compositions carry a headline result.
- `pyrolite.mplstyle` has **re-created itself** in `~/.config/matplotlib/stylelib/` (§5) and is again
  printing "Bad key" from every sweep worker.

### 23.9 Cross-cutting lessons

- **Two parsers reading one file will disagree.** PHREEQC stops at `END`; `parse_stoichiometry` does
  not. The Python side reported the phase as present while PHREEQC could not see it — a green unit
  check and a broken model. Validate database edits with an end-to-end solve, never with the parser.
- **Check what an end-member actually varies before using it as an axis.** HEX1/HEX2 differ in Mg/Si
  *and* Ca/Al; treating them as a Mg/Si pair produced a confident, wrong trend.
- **A clipped input looks exactly like an insensitive model.** Three T_p values returning identical
  mineralogy read as "crust composition doesn't matter", when the norm was discarding the difference.
  When a sweep comes back flat, check the input range is actually reaching the model.
- **Symbol presence is not ABI compatibility.** All 38 GSL symbols resolved and the library ran before
  corrupting its heap.
- **A justification can outlive the thing it justified.** `emit_quartz=False` was defended by an HT
  buffer that no longer exists anywhere in the executed code.

---

## 24. MAGEMin, the Mg/Si sweep, and a literature refocus (2026-08-21)

Replaced pMELTS with MAGEMin, calibrated Earth's basalt, swept Mg/Si across the full stellar range,
and then discovered that **most of this had already been published — by this group**. §24.6 is the
most important part of this section; read it before planning further crust-pipeline work.

### 24.1 pMELTS out, MAGEMin in

pMELTS cannot cover the stellar Mg/Si range:

- **Fails outright above molar Mg/Si ≈ 1.6.** Mg/Si 1.7, 1.8, 1.9 and 2.0 all fail at the start
  state at every temperature nudge, because its solution models do not span the
  ferropericlase-bearing assemblages that become stable there.
- **Below ≈ 0.8 it extrapolates badly**, returning 69 wt% SiO₂ rhyolites from what it treats as a
  peridotite.

**MAGEMin** (Riel et al. 2022, Holland–Green–Powell 2018 igneous dataset) converges across
**0.5–2.0 with no failures**, and stabilises `ne` (nepheline) and `fper` (ferropericlase) at high
Mg/Si on its own — independent confirmation, from a different thermodynamic dataset, that the CIPW
desilication cascade of §23.7 was inferring the right phase.

`src/kamino/data/make_crust_compositions.jl` replaces the pMELTS Python script. Three facts that
cost time:

1. **Julia was not actually installed.** juliaup was unpacking the runtime into `~/.julia/juliaup`
   on the **full home partition** and failing with ENOSPC. The fix is
   `JULIAUP_DEPOT_PATH=/home/pt426/data/julia_depot` (the data partition); the existing juliaup
   block in `~/.profile` sets only `PATH`, which is why it broke.
2. **MAGEMin has no isentropic mode.** It is a fixed-(P,T) Gibbs minimiser, so the isentrope is
   tracked by root-finding T at each pressure step to hold entropy constant. This matters:
   prescribing the solid adiabat ignores latent heat and over-melts. At T_p = 1350 the isentrope
   arrives at 1309 °C where the solid adiabat gives ~1362 °C.
3. **Julia buffers stdout to a file**, and the run emits fewer lines than the buffer holds, so the
   log stayed empty for 25 minutes and looked hung. `flush(stdout)` per composition; also narrowed
   the T_p bisection bracket from 1100–2200 to 1150–1650 (every solution lands in 1167–1523),
   roughly halving runtime.

### 24.2 Earth's basalt sets T_p = 1325 °C

`--calibrate` scans T_p at pyrolite against the PRIMELT primary melt:

| T_p | T_end | F | SiO₂ | MgO | Al₂O₃ | CaO | misfit |
|---|---|---|---|---|---|---|---|
| 1300 | 1278.6 | 0.084 | 48.26 | 10.34 | 18.29 | 11.62 | 2.96 |
| **1325** | **1294.5** | **0.117** | **47.77** | **11.06** | **17.67** | **12.31** | **2.22** |
| 1350 | 1309.4 | 0.153 | 47.63 | 11.70 | 16.91 | 12.83 | 2.63 |
| 1400 | 1339.7 | 0.219 | 48.13 | 12.88 | 15.16 | 13.21 | 5.43 |

1325 minimises the misfit *and* is the only value whose F falls in the 0.08–0.12 range of real MORB
primary melts. It sits **25 °C below** the project's `EARTH_MANTLE_POTENTIAL_TEMPERATURE = 1350` —
a model offset, not a claim about Earth: MAGEMin/HGP18 is slightly more melt-productive than
PRIMELT at the same T_p. **Anchor on the melt, not the label.**

### 24.3 T_p is not free — the constant-F closure

The user's argument, and it is right: a mantle that cannot melt also cannot transport heat by
magmatism, so it warms until melting carries the heat out. Holding T_p fixed makes the Mg-rich end
look like planets that "barely melt" (F ≈ 0.003 at Mg/Si 2.0) — an artifact of the closure. At
constant F those planets simply run hotter:

```
Mg/Si   T_p    dT vs Earth   melt
0.50    1167   -158          72 wt% SiO2 (granitic; quartz in the assemblage)
1.25    1325     0           Earth's basalt, by construction
1.60    1317    -8
2.00    1523   +198          MgO 18 wt% (picritic)
```

Supported by **Nature Comms Earth Environ (2022) 3, 261**, which finds thermal and water-cycling
feedbacks buffer the mantle at near-constant homologous temperature. Counterweight worth carrying:
**Korenaga (2016)** argues the adjustment is too slow for true self-regulation, so the buffer is a
tendency, not a constraint.

> ⚠️ **REJECTED: constant homologous temperature T_p/T_solidus**, which is the quantity the
> buffering literature actually identifies. It fails here for an instructive reason. The true
> multi-component solidus is set by the first infinitesimal melt, which **minor alkalis control**.
> Na₂O stays at pyrolite's 0.36 wt% while SiO₂ falls, so silica-poor bulks become nepheline-normative
> and their alkaline eutectic melts LOW — the solidus **collapses 1509 → 1049 °C** from Mg/Si 1.25 to
> 2.0, driving T_p **down 400 °C** and extinguishing melting entirely. The solidus measures a
> trace-driven eutectic, not the bulk refractoriness that governs heat transport. (The literature
> definition uses a *dry peridotite solidus parameterization*, which coincides with the true solidus
> for Earth-like compositions but not across this range.) Do not retry it against the true solidus.

### 24.4 The ultracalcic-melt problem

The user asked whether a Mg-rich crust should contain more olivine. It should, and the model was
wrong. At Mg/Si ≥ 1.6 the melts go **ultracalcic** — CaO/Al₂O₃ reaching **1.77** against MORB's
~0.78 — so the norm assigns Mg to diopside before olivine and returns ~52% diopside with only ~24%
olivine: a clinopyroxenite where a picrite belongs.

The mechanism is documented (**Médard et al. 2004**): once orthopyroxene leaves the residue, melts
turn nepheline-normative with CaO/Al₂O₃ > 1, exactly as observed here (`ol, spl, fper, liq` — both
pyroxenes gone, so all Ca enters the melt while Al is retained in spinel). But the same work states
that CaO up to 19 wt% and CaO/Al₂O₃ up to 1.8 **"exclude an origin from fertile lherzolites at
volatile-absent conditions"** and require volatiles or a cpx-rich source. This run is volatile-free
peridotite — the source the literature rules out.

Tested and **excluded** as the cause: the mantle construction. Building high Mg/Si by enriching MgO
only (holding Ca/Si, Al/Si at pyrolite) instead of trading MgO↔SiO₂ still gives CaO/Al₂O₃ ≈ 1.7.

### 24.5 The fix — stop melting at clinopyroxene exhaustion

Ultracalcic melts appear only *after* melting past cpx-out, so make cpx-out the closure. Melt
fraction is then solved per composition rather than assumed:

```
Mg/Si  fcore  mantle FeO   F      SiO2   MgO    CaO    CaO/Al2O3
1.00   0.70   8.73         0.200  48.56  12.12  12.75  0.90
1.30   0.70   7.95         0.212  47.87  12.66  13.47  0.86
1.60   0.70   7.30         0.152  44.79  11.72  16.94  1.01
1.60   0.90   2.56         0.206  47.93  14.22  14.85  0.84
```

**CaO/Al₂O₃ falls to 0.84–1.01** across almost the whole grid — basalts, not melilitites — and F
comes out at 0.15–0.22, bracketing the 20% Guimond et al. adopt for Earth but derived rather than
assumed. That review explicitly names this as the better-but-undone approach: *"Such a
petrologically informed melt fraction could be identified for all exoplanetary mantles, but given
the assumptions already inherent to these calculations would be beyond the scope of this review."*

### 24.6 ⚠️ The literature refocus — most of this was already published

**[Guimond, Wang, Seidler, Sossi, Mahajan & Shorttle (2024)](https://arxiv.org/abs/2404.15427),
*Rev. Mineral. Geochem.* 90, 259** — "From stars to diverse mantles, melts, crusts and atmospheres
of rocky exoplanets". **Section 4 is "Melts and crusts of rocky exoplanets". Shorttle is a
co-author.**

Their **Figure 5 is what §24.1–24.5 re-derived**: same construction (*"Except for Mg and Si, the
compositions of other elements are kept the same as those of bulk silicate Earth"*), same range
(Mg/Si 0.5, 1.0, 1.5), using Perple_X for mantle mineralogy and **pMELTS for melting** (their
Fig. 8). Established results this session rediscovered the hard way:

- **The boundaries.** *"silicate mantles becoming olivine-free at Mg/Si ≲ 0.8, or orthopyroxene-free
  at Mg/Si ≳ 1.6"* — exactly where the assemblages here lost opx, gained ferropericlase, and the
  melts turned ultracalcic.
- **The olivine answer.** *"Increasing Mg/Si above unity produces more forsterite olivine at the
  expense of enstatite orthopyroxene"* in the **mantle**, while *"melts are less magnesian than
  their corresponding mantles and have broadly similar silica contents"*. So a Mg-rich planet has an
  olivine-rich mantle but a broadly basaltic **melt**. The picritic MgO-18 melts were the anomaly.
- **CaO is already flagged as the quantity that matters, for kamino's exact reason:** *"CaO is also
  higher in melts than in mantles (due to its incorporation into clinopyroxene, which readily melts
  out), which is significant given the key role crustal Ca has in planetary carbon cycles via
  carbonate formation."* That is the Ca-starvation thread (§5, §13, §22.6), already named.
- **Mg/Si is not the dominant control.** [Putirka & Rarick (2019)](https://arxiv.org/pdf/1907.05506),
  surveying >4,000 Hypatia stars, find *"half or more of the range of exoplanet mantle mineralogy is
  controlled by core formation"* via Fe partitioning — held **fixed** at FeO = 8.05 in every sweep
  here. The review agrees **(Mg+Fe)/Si beats Mg/Si** as a predictor of olivine/opx.
- **The exotic ends are rare.** **89% of Hypatia stars** fall in the ordinary olivine + opx field.

**Implication.** Do not chase Mg/Si 1.7–2.0; the models are published as unreliable there and ~11%
of stars are outside the normal field. The additive question for this paper is not "what crust does
Mg/Si make" — that is done — but **what that crust does to ocean chemistry and the carbon cycle**,
which is the gap the review itself names (*"no studies explicitly consider variable carbon mineral
speciation..."*). kamino is the weathering half the review stops short of.

**Ask Shorttle before rebuilding anything** — the melt compositions behind their Figure 9 may
already exist.

### 24.7 Proposed parameterisation, and what remains

A **2-D grid: Mg/Si 0.8–1.6 × Fe_core/Fe_bulk 0.55–0.95**, melt stopped at cpx-out per composition.
The Fe axis earns its place — at Mg/Si 1.6, varying `fcore` 0.55 → 0.90 (mantle FeO 10.6 → 2.6 wt%)
moves CaO/Al₂O₃ **1.34 → 0.84**, i.e. **the Fe parameter controls whether the Mg-rich corner goes
ultracalcic at all**. Invisible on a 1-D sweep with FeO pinned.

Bounds: the review's olivine-free/opx-free limits on Mg/Si; its warning that thermodynamic models
*"become unreliable at very FeO-rich compositions (≳25 wt%)"* on the Fe axis (this grid tops out at
12.5 wt%).

- **The `fcore` → mantle FeO mapping is ours, not theirs** — bulk Fe was back-calculated by
  requiring Earth's BSE FeO = 8.05 at `fcore` = 0.7. Putirka & Rarick define α_Fe = Fe_BSP/Fe_BP on
  a cation weight basis with explicit core Ni and Si; adopt their formulation before citing them.
- **Mg/Si 1.6 with fcore 0.55 is still marginal** (F = 0.099, CaO/Al₂O₃ = 1.34).
- The current `crust_compositions.csv` is the **1-D constant-F** table; its Mg/Si ≥ 1.7 rows should
  not be used as they stand (§24.4).
- `make_crust_compositions.py` (pMELTS) is superseded by the `.jl` version and should be deleted.
- Nothing from §24 is wired into `oxide_composition` yet — the model still reads the PRIMELT
  spreadsheet.

### 24.8 Cross-cutting lessons

- **Search the literature before deriving, not after.** Two days of pipeline work re-derived a
  published figure from this group, including its boundary values. The user's "this seems like it
  has been done before" was worth more than any calculation in this section.
- **A closure can manufacture the anomaly it appears to reveal.** "Mg-rich planets barely melt" was
  an artifact of fixed T_p; "Mg-rich crusts are clinopyroxenites" was an artifact of melting past
  cpx-out. Both looked like results.
- **Check what an end-member actually varies.** Also §23.4 — a recurring failure mode this session.
- **Silence is not progress.** A buffered log made a 25-minute run indistinguishable from a hang.
- **A phase appearing in the residue can mark the edge of a model's validity, not just a mineral
  change.** Ferropericlase coexisting with melt at 1–3 GPa is a lower-mantle phase; its appearance
  was the signal that the calculation had left the range where peridotite melting means anything.
