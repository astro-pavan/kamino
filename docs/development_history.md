# Kamino — Consolidated Development Context

**Compiled 2026-08-17, extended through 2026-08-24** from all Claude Code sessions on this project
(2026-07-14 → 2026-08-24). Later sessions are weighted more heavily; where an early decision was later
reversed, the reversal is what's recorded, with the original noted as history.

Verified against the working copy on branch `NaCl-chemistry`. **§11** lists what is actually in the
code versus what was only proposed. **§18–§25 are the most recent work** and supersede earlier framing
where they conflict; §20 records the iron charge-leak fix, §21 removes reverse weathering from the pore
space (which invalidates the transition evidence in §20.4, see below). **§25 is the current authority
on the CRUST PIPELINE** — it replaces the (T_p, Mg/Si) axes with (Mg/Si, core-formation ΔIW), closes
the CIPW norm with akermanite, and records validation against Guimond et al., Brugman et al. and Katz
et al.; it supersedes §5, §23 and most of §24, though §24.6's literature refocus stands. §22 remains the authority on
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
      + F_prec           # fast ocean precipitation, tau_prec = 100 kyr x (depth/3 km), see §26
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
| `tau_prec` | **100 kyr x (ocean_depth / 3 km)** | fast precipitation; depth-scaled 2026-08-25, see §26 |
| `tau_rw` | 5 Myr | reverse weathering; deliberately **not** depth-scaled, see §26.2 |
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

Databases: **`lt_weathering_pitzer.dat` is the runtime LT database** (`chemistry.py`, overridable
via `KAMINO_LT_DATABASE`); `lt_weathering_sit.dat` is the SIT-based alternative, kept for
sensitivity tests; `hydrothermal.dat` is the HT database. Both LT databases are generated by
`make_database.py`. `hybrid_ocean.dat` and `ocean_chem.dat` are the superseded hand-built files,
retained only as historical reference.

> ⚠️ **CORRECTED 2026-08-24 (§25.12), then RESOLVED (§25.14).** This section previously said
> "`make_database.py` builds these". It did not: it wrote only `lt_weathering_sit.dat`, while the
> model loaded `hybrid_ocean.dat`, generated by a `make_hybrid.py` that was never committed. That
> is fixed. `make_database.py` now builds **either** activity model — `base="pitzer"` (default) or
> `base="sit"` — entirely from databases bundled with the `phreeqc` package, and `chemistry.py`
> loads `lt_weathering_pitzer.dat`. **`hybrid_ocean.dat` is no longer referenced by anything.**
> Read §25.14 before touching the databases; it records what the pitzer base has to graft
> (aluminium, ferric iron, five phases) and the three silent traps involved.

`make_database()` imports the base model,
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
   >
   > ⚠️ **Scope limit, added 2026-08-25 (§26.4).** This holds for the configuration measured here
   > — shallow, `crust_production_rate = 0.01` — and **not in general**. At 20 km and Mg/Si 1.25,
   > Sepiolite(d) precipitates at the 100 kyr baseline (SI +3.778), and lengthening `tau_prec` to
   > 700 kyr changes its flux by 3.6%. The **50× ratio is not the load-bearing quantity**; what
   > matters is which side of saturation each bucket sits on. Do not carry the ratio forward as a
   > design constraint.

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

> 🔴 **Added 2026-08-25 (§27.5): the §22 Earth calibration is now stale and blocks the paper
> figures.** Adopting hedenbergite moved the Earth anchor **−8.2 K at S = 1.0**, so `alpha` and
> `kd_mg_ht` are fitted to a crust the model no longer produces. §25.12 (the reproducible-database
> switch) forces the same exercise, so both should be redone together — one calibration, not two.
> Everything already on disk, including the pilot and the §26 tau runs, predates both changes.

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

---

## 25. The two-parameter crust pipeline, and its validation (2026-08-24)

**Supersedes §5, §23 and most of §24 for the crust pipeline.** The (T_p, Mg/Si) axes are gone; the
PRIMELT spreadsheet is gone; the 1-D constant-F table is replaced by a 2-D grid that is generated,
merged, validated and wired into the model. §24.6's literature refocus still stands and is
strengthened: §25.6 below records what now reproduces from Guimond et al. and from two sources
outside that group.

Full methods reference: **`docs/crust_composition.md`**. `docs/tectonics.md` had a stale
"Crust Composition" section describing the PRIMELT pipeline; it has been rewritten.

### 25.0 The state it was found in

`Planet()` could not be constructed at all on the current machine: `import_primelt_spreadsheet`
needs `xlrd`, which is absent from `.venv`, and `planet.py:130` calls `mineral_composition` at
construction. This was masked on the old box. The baseline was recovered by installing `xlrd`, and
`mineral_composition(1350)` was captured before anything was touched — it matches §5's recorded
assemblage exactly, and is bit-identical after all the changes below.

### 25.1 The two axes, and what was rejected

| axis | range | Earth | controls |
|---|---|---|---|
| mantle molar Mg/Si | 0.5–2.0 | 1.25 | olivine/orthopyroxene, feldspar/feldspathoid |
| core-formation ΔIW | −5 to −1 | −2 | mantle FeO, and hence crust Fe |

Rejected, with reasons recorded in the module docstring so they are not re-proposed:

- **C/O** — rejected on three independent grounds; see §25.9, added after the initial write-up when
  the FGK-only justification was challenged. Short version: the stellar route is rarer for M dwarfs
  than for FGK, the non-stellar (soot) route makes a different population with a kinetically inert
  crust, and the real effect is on the carbon SOURCE, which is already `outgassing`.
- **T_p** — not observable and not free; now *solved* per composition.
- **Ca/Al** — real, but both are refractory so stellar Ca/Al varies least of the candidates.
- **Na₂O** — the strongest **un-swept** lever and should be reported as such (alkalis move the
  solidus more than MgO/FeO does; albite-vs-nepheline sets ocean Na). Excluded because Na is
  moderately volatile and devolatilisation is stochastic, so it does not belong on a grid indexed
  by observables.

### 25.2 The ΔIW axis

Metal–silicate equilibrium, ΔIW = 2 log₁₀(a_FeO/a_Fe) under ideal mixing, inverted for mantle FeO
against the pyrolite non-Fe budget. The activity constant (C = 0.567957) is **calibrated so ΔIW = −2
reproduces BSE FeO = 8.05 wt% exactly**; k = 0.020237 cation mol per wt%.

> **The axis is LOGARITHMIC in FeO.** −5 → 0.26, −4 → 0.82, −3 → 2.59, **−2 → 8.05**, −1 → 24.1 wt%.
> This is the single least intuitive fact about the parameterisation. The grid stops at −1 because
> mantle FeO then hits the ~25 wt% ceiling above which Guimond et al. state the models are unreliable.

**Two fO2 values, not one.** The ΔIW that sets mantle FeO is the *core-formation* value (Earth ≈ −2);
the one `m-class/outgassing/outgassing_model.IW_offset` wants is the *modern melt* value (Earth ≈ FMQ
≈ +3.5). They are separated by post-core-formation self-oxidation (Fe disproportionation, Guimond
§2.4). Resolution: sweep the core-formation value, and carry `DELTA_IW_SELF_OXIDATION = 5.5` in
`constants.py` as an **Earth anchor, not a law**, with `melt_delta_iw()` as the documented handoff.

**Young et al. (2023) validates the formulation and exposes a calibration flaw.** They write ΔIW
identically, ideal-mixing convention and all, and put Earth at −2.2 (mine: −2.00 by construction).
But their reduced anchor disagrees with mine: they give ΔIW ≈ −5 ↔ FeO 0.07 wt%, where this mapping
gives −6.14 — a **~1.1 log unit offset at the reduced extreme**, larger than the ±0.4 slop from
γ_FeO. Cause: single-point calibration at Earth extrapolated three log units down. It is a *labelling*
error, not a physics one (the melt moves 0.08 wt% SiO₂ between ΔIW −5 and −4.5), so the grid was not
regenerated — but **do not quote the reduced-end ΔIW labels against Young's** without a two-point
recalibration (Earth −2.2/8.05 plus E-chondrite −5/0.07).

Young et al. also supply the physical content of the reduced end: E chondrites and aubrites sit at
ΔIW −4.5 to −6.5, and Mercury ≈ −5. The low-ΔIW columns are the enstatite-chondrite regime, not an
arbitrary bound.

> ⚠️ **A sign error made during this session, corrected.** It was asserted that Si's siderophile
> behaviour makes the (low Mg/Si, reduced) corner depopulated, because reducing conditions drive Si
> to the core and raise mantle Mg/Si. Young et al. run the opposite way: their oxidation mechanism is
> Si⁴⁺ + 2Fe⁰ = Si⁰ + 2Fe²⁺, and running it *forward* is what takes an embryo from ΔIW −5.8 to −2.1
> **while** moving Si into the core. So mantle Mg/Si and ΔIW rise *together*, and the depopulated
> corner is (low Mg/Si, **oxidised**). The robust claim — that the two axes are physically coupled and
> the grid does not sample nature uniformly — survives; the sign did not, and the sign is
> model-dependent (classical metal–silicate partitioning favours Si-in-metal at *low* fO2).

### 25.3 Generating the grid

`make_crust_compositions.jl` rewritten for two axes. `MAGEMin_C` was not in the Julia depot on this
machine; `Pkg.add` plus a probe confirmed the API fields (`ph`, `entropy`, `frac_M_wt`, `bulk_M_wt`)
are unchanged.

- **Construction order matters and is the point:** iron first (it changes the size of the silicate
  budget), then Mg/Si within what remains. The other order makes the axes non-orthogonal — the §23.4
  failure mode.
- **Ferric iron stays off** (`O` = 0.0). Enabling it would perturb the Earth anchor; Guimond run
  their grids the same way. See §25.6 for the one place this costs.
- **Closure: F = 0.20 fixed**, Guimond's value. T_p solved per composition by bisection.
- **Sharding by process, not thread** — MAGEMin keeps mutable per-point workspaces in the `data`
  handle, and §23.5 records what a poisoned solver looks like. 7 processes × ~96 min wall for
  153 points; a single process would have been ~5.5 h.

Grid: 17 Mg/Si × 9 ΔIW (0.5 spacing = uniform in log FeO) = **153 points, 0 failures**.
Earth: T_p 1383 °C, F 0.201, SiO₂ 47.88, CaO/Al₂O₃ 0.85. **Misfit to PRIMELT 4.85 wt%**, against
2.22 at F = 0.117 (§24.2) — the quantified cost of adopting Guimond's melt fraction.

New tooling: `merge_crust_slices.py` (refuses a table with holes), `check_crust_table.py` (structure,
anchor, per-oxide mass balance, petrology flags, end-to-end PHREEQC), and `--probe` / `--points` /
`--validate` / `--slice` modes.

### 25.4 Both end-members misbehave, and it is one criterion

**Silica-rich end.** `emit_quartz=False` was discarding excess silica: at Mg/Si 0.5 it **relabelled a
72.7 wt% SiO₂ rhyolite as a 51.7 wt% basalt**, throwing away 0.71 mol SiO₂. Fixed by
`emit_quartz=True` plus `Quartz` in `primary_minerals` (dissolve-only, so it cannot act as a
low-temperature silica buffer). Measured to be a **no-op at Earth** — the assemblage is bit-identical
to the captured baseline.

**Silica-poor end.** 46/153 cells retained a silica deficit after the albite→nepheline cascade,
~5% of their SiO₂ — a genuine mass-balance violation, silently absorbed. All 41 ultracalcic cells
were inside the 46.

**They are the same criterion.** Usability collapses onto mantle **(Mg+Fe)/Si**, the ratio Guimond
name as a better olivine/opx predictor than Mg/Si: usable cells span 0.503–1.694, unusable 1.692–2.691
— a 0.002 overlap. Below ~0.8 the crust is quartz-normative; above ~1.69 the phase set runs out of
desilication capacity. Earth sits at 1.399.

> **The Fe axis controls whether the low-Mg/Si corner is granitic at all.** At Mg/Si 0.5, melt SiO₂
> is flat at 71–75 wt% from ΔIW −5 to −1.5, then collapses to **49.9 wt%** at ΔIW −1. (Mg+Fe)/Si
> crosses 0.8 between exactly those two cells. Invisible on any 1-D Mg/Si sweep with FeO pinned,
> which is how every earlier version of this pipeline ran.

**cpx-out does NOT fix the Mg-rich corner** — tested, contra §24.5. At Mg/Si 1.6 it still gives
CaO/Al₂O₃ = 1.28 (§24.5 predicted 1.01), *with cpx still in the residue*, because at Mg/Si ≥ 1.5 the
source is already orthopyroxene-free, so Médard's mechanism operates regardless of where melting
stops. It also breaks the low end (F = 0.382 at Mg/Si 0.5 — not a crust). cpx-out remains available
as `--closure cpx-out` for sensitivity work, but F = 0.20 is retained.

### 25.5 Akermanite closes the norm

**Akermanite (Ca₂MgSi₂O₇) was already in `hybrid_ocean.dat`** with complete llnl-lineage
thermodynamics on the H₄SiO₄ basis — but was **not** in `make_database.ADDED_PHASES`, so a rebuild
would have silently dropped it. Latent bug, now fixed.

Norm **step 5c**: `2 Diopside → Akermanite + ½ Forsterite + 1.5 SiO₂`, releasing 0.75 mol SiO₂ per
mol diopside. **Clears all 46 cells** at 2.1–7.8 wt% akermanite. The residual-deficit case is now a
`ValueError`, not a warning — no cell triggers it. `check_crust_table.py`: **ALL CHECKS PASSED**,
153/153 mass-balancing every oxide to <0.05 wt%.

> ⚠️ **Larnite was rejected, and this is worth not re-litigating.** Textbook CIPW desilicates to
> larnite (Ca₂SiO₄), and `Kinec_v3.dat` has it **complete — PHASES and RATES** — so it was the
> zero-proxy option, and §23.7's stated blocker ("not in the database") is obsolete. Rejected because
> (i) k_eff at 300 K / pH 6 is **606× wollastonite and ~10⁵× diopside**, so at 5–8 wt% it would supply
> the entire dissolution flux (§5: the Wollastonite→Diopside swap was "one of the most consequential
> changes in the project"; larnite is three orders beyond wollastonite); and (ii) it is a cement
> clinker phase, rare in nature, whereas these melts are melilititic and melilitites crystallise
> **melilite**. Giving larnite a real rate would model a mineral that is not in the rock.

**Akermanite's kinetics are a proxy** — PHASES-only in Kinec_v3, Kinec.v2, llnl, core10 and
Thermoddem. Melilite is a sorosilicate, between orthosilicates and chain silicates, so the rate is
bracketed rather than pinned (`set_akermanite_proxy`): wollastonite / **forsterite (default)** /
diopside. Measured spread over the 46 cells:

| quantity | slow | mid | fast | spread |
|---|---|---|---|---|
| **Ca** | 8.4e−12 | 7.2e−11 | 4.0e−10 | **48×** |
| Mg | 2.5e−10 | 2.8e−10 | 4.6e−10 | 1.8× |
| Alkalinity | 2.0e−09 | 2.2e−09 | 3.2e−09 | 1.6× |

At the fast bound akermanite supplies **98% of all Ca** at 2–8 wt% of the crust. So those cells are
sound for salinity/alkalinity/Mg work and **must sweep the bracket for any Ca or carbonate-burial
result**. A measured åkermanite rate (the CO₂-mineralisation and slag literature is the place to
look) would collapse this entirely. Note also that the *broken* cells delivered Ca = 0.019 mM — the
mass violation was not erring safe, it was starving the carbon cycle's key ion by three orders.

### 25.6 Validation — three independent sources

Figure: `output/crust_validation.png`, from `experiments/plot_crust_validation.py`.

**Katz, Spiegelman & Langmuir (2003)** — the standard anhydrous-peridotite parameterisation,
independent of both MAGEMin and Guimond. Over 1.0–2.0 GPa and F = 0.05–0.20: **mean offset +16 °C,
sd 10 °C, max 29 °C**, MAGEMin consistently marginally hotter for a given F. Their cpx-out melt
fraction at these pressures is F = 0.23–0.26; so three independent estimates of cpx-out — Guimond's
assumed 0.20, our MAGEMin solve 0.213, Katz 0.23–0.26 — agree within ~20%. **The closure is not an
arbitrary inheritance**, which is the most useful single sentence for a methods section.

**Guimond et al. (2024)** — seven claims tested:

| claim | theirs | ours |
|---|---|---|
| olivine-free below Mg/Si | ≲ 0.8 | **0.70** ✓ |
| orthopyroxene-free above Mg/Si | ≳ 1.6 | **1.50** ✓ |
| excess oxides form their own phases | predicted | quartz ≤ 0.6, ferropericlase ≥ 1.7 ✓ |
| melt-vs-mantle relations (§4.2) | 5 relations | **20/20** across 4 compositions ✓ |
| melt vs mantle variance (Fig. 10) | SiO₂/MgO less, others greater | **6/7** ✓ |
| mantle FeO → melting temperature (§4.1) | up to ~100 °C | +53 °C mean ~ |
| F = 0.20 is where cpx is lost | assumed | **0.213** solved ✓ |

Both boundary offsets are *inward*, which is what 20% melt depletion predicts (the residue is more
olivine-rich than their subsolidus mantle) — agreement and its explanation in the same number. The
Fig. 10 SiO₂ failure (1.11) is our uniform grid oversampling the exotic corner: restricting to the
ordinary ol+opx field gives 0.76. FeO stays marginally below 1 in every subset — a genuine small
disagreement, because our ΔIW axis spans 100× in mantle FeO where their Hypatia population barely
varies it.

The (Mg+Fe)/Si ≈ 0.8 boundary **reproduces without being encoded** — see the §25.4 blockquote.

**Brugman, Phillips & Till (2021)** — the only experimental test. Piston-cylinder melting of HEX1, a
hypothetical exoplanet mantle at molar Mg/Si = 1.42, directly on our axis.

- *Bulk:* our construction at Mg/Si 1.42, FeO 8.23 reproduces their independently-designed starting
  material — SiO₂ 42.40 vs 42.00, MgO 40.39 vs 40.04, FeO exact, **molar Mg/Si 1.420 vs 1.421**. The
  Mg–Si–Fe backbone was never fitted to it. (Al₂O₃ −8%, CaO −6%, Na₂O +71% differ by construction.)
- *Melt, matched conditions* (their exact bulk, isobaric 1.5 GPa, F = 0.05): Al₂O₃ **+4.7 wt%**,
  FeO **−3.9 wt%** against their Table 5 average. **This is the weakest link in the chain and should
  be reported as such.** Their table averages 1.0–2.0 GPa and F = 0.004–0.056, so it is not a
  single-condition comparison; but the FeO deficit is the direction the all-ferrous assumption
  predicts, since Fe³⁺ is strongly incompatible and would raise melt FeOtot. This is the concrete
  argument for eventually enabling the ferric component.
- Their solidus claim (HEX1's ≈ Earth's anhydrous peridotite solidus) holds: our T_p at Mg/Si 1.4 and
  1.25 differ by 12 °C.

### 25.7 What is wired, and what is not

Wired: `crust_composition.py` (PRIMELT path deleted, so the `xlrd` dependency is gone),
`constants.py`, `mineral_info.py`, `mineral_rates.py`, `make_database.py`, `planet.py`,
`diagnostics.py`, `experiments/parameter_sweep.py`, `experiments/plot_results.py`.

> ⚠️ **`mg_si_ratio` was RENAMED to `mantle_mg_si`, deliberately.** The old parameter was a
> *multiplier on Earth's 1.23*, so `mg_si_ratio=1` meant Mg/Si = 1.23. Reusing the name would have
> silently changed the meaning of every existing sweep config. An `AttributeError` is the correct
> outcome. `mantle_potential_temperature` is gone; `delta_iw` is new.

Outstanding:

- **The Mg/Si ≤ 1.6 cap is undecided.** It is now a scientific choice, not a forced one: the norm
  closes everywhere. Three independent literature lines put the boundary at 1.6 (Guimond's opx-free
  limit, §24.7's warning, and where MAGEMin stabilises ferropericlase).
- **Two-point ΔIW recalibration** against Young et al.'s E-chondrite anchor (§25.2).
- **A measured åkermanite dissolution rate** (§25.5).
- **Ferric iron** — currently off; §25.6 gives the first concrete evidence it costs something.
- `ggge967-sup-0002-primelt1.xls` was **kept**, against the plan, because it is the provenance for the
  PRIMELT reference melt that `--calibrate` still scores against. Nothing in the model reads it.
- §24.6 stands: **ask Shorttle** whether melilite-normative crusts are something the paper should
  claim at all, and whether the melts behind their Figure 9 already exist.

### 25.8 Cross-cutting lessons

- **A mass-balance violation is not conservative.** The 46 broken cells looked like they were erring
  toward "less reactive crust". They were delivering Ca = 0.019 mM against Earth's 7.3 — a
  three-order starvation of the ion the carbon cycle runs on, arrived at by silently discarding mass.
- **Check whether the phase you need is already in the database.** Akermanite had full thermodynamics
  sitting in `hybrid_ocean.dat` while §23.7 recorded the fix as blocked. Larnite had *complete*
  thermodynamics and kinetics in `Kinec_v3.dat` while the same section recorded it as "not in the
  database". Both statements were true of the wrong file.
- **Complete data is not a reason to use a phase.** Larnite was the zero-proxy option and is still
  the wrong answer. Availability is not petrology.
- **A bracket is a result.** Proxying akermanite's rate three ways converted "unusable, unknown why"
  into "usable, with a 48× uncertainty localised to one ion" — which is actionable, where a single
  guessed rate would have been invisible.
- **Two failure modes can be one criterion.** The granite corner and the melilitite corner looked
  unrelated until both fell out of mantle (Mg+Fe)/Si, with a 0.002 overlap between the usable and
  unusable ranges.
- **Render the figure and look at it.** Panel (c) of the validation figure had `ylim` at 24 while
  SiO₂ is ~48 wt%: the bars were silently clipped *and* their labels drew into a neighbouring panel.
  The validator checks colour, not geometry.
- **A closure can fix one end and break the other.** cpx-out is right for peridotite-like
  compositions, does nothing for the opx-free corner it was proposed to fix, and returns F = 0.38 at
  the quartz-saturated corner.
- **Correct the sign, keep the claim.** The Mg/Si–ΔIW coupling is real; the direction asserted for it
  mid-session was not established and Young et al. reverse it. See the §25.2 blockquote.
- **Check the databases before asserting what is not in them.** Three successive scoping claims about
  carbon phases were made confidently and each was wrong or beside the point — graphite and kerogen
  (the latter with kinetics) turned out to be sitting in files the project already ships. §25.9.
- **Check that the build script builds the file you actually load.** `make_database.py` writes
  `lt_weathering_sit.dat`; the model loads `hybrid_ocean.dat`, which is built by a script that is not
  in the repository. Four separate symptoms were each explained away locally over two months before
  the question "has it always run on hybrid_ocean.dat?" connected them. §25.12.
- **A file that looks broken may only be stale.** `lt_weathering_sit.dat` was three phases short and
  looked unusable; every one of those phases was already in `ADDED_PHASES`, and a rebuild produced a
  complete, working database. Regenerate before diagnosing. §25.12.

### 25.9 Why C/O is not an axis — the full argument (added 2026-08-24, later the same day)

§25.1 originally justified dropping C/O with the FGK statistic alone. That is the wrong stellar
population for this paper — kamino's targets are temperate rocky planets, which are overwhelmingly
M dwarf hosts — and the argument was challenged on exactly that point. It now rests on three
independent legs, written up in full in `docs/crust_composition.md` §1.4.

**(i) The stellar route is rarer for M dwarfs, not commoner.** Nakajima & Sorahana (2016) measured
C and O directly in 46 nearby M dwarfs from K-band CO and H₂O: **none has C/O > 0.8**. Population
limits: Gaidos (2015) puts C/O ≈ 1 below 1.2×10⁻³ (95% conf.), and 6×10⁻⁴ from SDSS DR7 (99%);
Gizis et al. (2016) put 0.8 < C/O < 1 below 1%. Against ~1% for FGK.

> **Counterintuitive and worth remembering: M dwarf C/O is measured BETTER than FGK C/O.** The
> CO + H₂O method is non-differential and needs no assumed solar C/O; the FGK optical C I / [O I]
> route is differential and its literature disagrees violently (some find C/O > 0.8 in >20% of
> stars, others in none). The M dwarf data breaks that tie toward the low-frequency answer, so
> Guimond's ~1% is on the well-supported side. The historical "carbon-rich M dwarfs" were
> artefacts — weak TiO reads as high C/O and is also what low metallicity looks like.

**(ii) The soot route is real, and it is not a C/O problem.** Li et al. (2026), *ApJ Lett.* 997,
L29 — the soot line (~500 K) lies INSIDE the water-snow line (~160 K), so anything forming beyond
the snow line also formed beyond the soot line and should carry refractory organic carbon. Their
soot planet (74% rock / 26% soot) reaches **bulk planetary C/O ≈ 1 at solar stellar C/O ≈ 0.55**.
Location, not composition — which undercuts the "only ~1% of stars can do this" framing entirely.

It still does not become a kamino axis: the carbon is refractory organic CHON (C:H:O ≈ 100:77:14),
not carbide; the archetypes are dry or carry 25–50% H₂O by mass, so neither is a rocky ocean
planet; and reduced carbon is kinetically inert in a cold abiotic ocean, so it acts as a **diluent**
exactly as Quartz does at low Mg/Si.

> ⚠️ **Three wrong scoping statements were made before the right one.** In order: (a) "the PHREEQC
> databases have no graphite, diamond, carbide or organic carbon" — **false**, checked and
> retracted: `C(cr)` (graphite) is in `sit.dat` AND in the project's own `lt_weathering_sit.dat`,
> and `core10.dat` / `Kinec_v3.dat` / `Kinec.v2.dat` carry three kerogen phases **with RATES
> blocks**, of which KerogenC515 (C₅₁₅H₅₉₆O₇₂) matches Li et al.'s soot O/C exactly (14.0 vs 14±3)
> with soot's H/C bracketed by C128 and C292. (b) "adding them makes carbon redox-active, which is
> a big model change" — true but beside the point. (c) The correct statement: **reduced carbon is
> kinetically inert at seafloor conditions**, so the phases being available changes nothing. Only
> diamond is genuinely absent from every database, and that absence costs nothing — see below.
>
> The lesson: three successive scoping claims were each made confidently and each was wrong or
> irrelevant. Check the databases before asserting what is in them.

*Diamond is excluded by physics, not by data.* It needs ~1.76 GPa at 25 °C (Berman–Simon), against
~0.03 GPa at a 3 km ocean floor and ~0.49 GPa at the deepest ocean swept — and liquid water gives
way to ice VI near 1.0 GPa, **below** diamond's field. There is no habitable pressure at which a
liquid-water seafloor sits on diamond-stable rock.

**(iii) The effect is on the SOURCE, and lands on a parameter that already exists.** Crust
composition sets the weathering sink — a feedback bounded by crust production and land area. C/O
sets the carbon source — a linear driver. Kamino carries that as `outgassing`. Two consequences:

- **Magnitude.** Li et al.'s soot is 79.9 wt% carbon, so 26 wt% soot is ~208,000 ppm bulk carbon
  against Earth's mantle 10–100 ppm: a factor of **2,000–21,000**. The sweep tops out at 10×, so a
  soot planet is **200–2,000× beyond the swept range** — and §8 already records the model saturating
  (50 runs pinned at exactly 5.00 bar, `acid_ocean`).
- **Speciation, which fails first.** `outgassing_flux[c_idx]` is inorganic carbon tied to
  alkalinity. A soot planet is reduced enough that volcanogenic carbon comes out as **CH₄ and CO**
  (Li et al.: soot/oxidised-iron ≈ 1.08 against their 0.3 threshold). CH₄ forms no carbonic acid, so
  the ocean chemistry barely registers it — but `get_T_surface(S, P_CO2, albedo)` has no methane
  axis. **A carbon-rich planet breaks this model at the CLIMATE interface, not the crust interface.**

**What C/O does leave behind.** A soot-rich planet is deeply reduced, which places it at the low-ΔIW
end of the axis already built. So that corner now has two independent formation routes motivating
it — enstatite-chondrite-like accretion (Young et al. 2023, §25.2) and soot accretion — where §25.2
gave only one.

### 25.10 The norm switched to pyrolite (2026-08-24, later still)

§25.1–25.9 describe a hand-rolled CIPW norm. It has been replaced by **pyrolite's implementation**
(Williams et al. 2020) plus a documented correction step, for a write-up reason the user gave and
which is sound: a cited implementation with a short list of modifications stays in the main body of
a paper, whereas a bespoke norm has to be described in an appendix. The hand-rolled version is
retained as `_cipw_norm_native` and is now the cross-check.

**The user's construction was the one that worked.** Restricting the OXIDES before the norm runs
(dropping TiO₂, K₂O, MnO, P₂O₅, Cr₂O₃ and renormalising) suppresses magnetite, ilmenite, orthoclase
and leucite at source. An earlier attempt to delete those phases from the OUTPUT was measured and
is wrong: standard CIPW allocates Fe to magnetite/ilmenite and K to feldspar ahead of the
ferromagnesian phases, so deleting the products strands their cations and corrupts the silica
balance — normative olivine came out 60% low and orthopyroxene 4× high at the Earth anchor. With
pre-filtering, feldspar agrees to <0.01 wt%.

> ⚠️ **pyrolite invents ferric iron, and `Fe_correction=None` does NOT disable it.** `normative.py`
> reads `if Fe_correction is None: Fe_correction = "LeMaitre"`, and that default assigns Fe₂O₃/FeO
> from a **TAS rock-type classification** (`_MiddlemostTASRatios`, ratios 0.1–0.5). On an
> all-ferrous melt it produced 3.6 wt% normative magnetite — silently contradicting the entire ΔIW
> axis, which is built on `O` = 0.0 in MAGEMin. The correction is skipped only where **both** FeO
> and Fe₂O₃ are > 0, so `Fe2O3` must be a tiny POSITIVE sentinel (1e-9), never 0.0. Undocumented,
> version-fragile, and the single most likely thing to break on a pyrolite upgrade. `cipw_norm`
> guards it by raising on any unexpected phase and naming magnetite in the message.

**Six correction reactions**, all balanced, all using only database-available phases:

```
hedenbergite + 1/2 forsterite -> diopside     + 1/2 fayalite     silica-neutral
hedenbergite + enstatite      -> diopside     + ferrosilite      silica-neutral  [olivine-free]
ferrosilite  + 1/2 forsterite -> enstatite    + 1/2 fayalite     silica-neutral
2 ferrosilite                 -> fayalite     + SiO2             releases        [olivine-free]
larnite      + 2 diopside     -> 2 akermanite + SiO2             releases
nepheline    + 2 SiO2         -> albite                          reabsorbs
```

The olivine-free fallbacks were not anticipated: 35 cells failed on the first run because
silica-oversaturated melts (Mg/Si ≲ 0.7) carry **no normative olivine at all**, so there is nothing
to exchange iron with. The guard raised rather than emitting a wrong crust. Larnite is routed
through diopside because **zero of the 46 silica-deficient cells carry free enstatite** — the clean
`larnite + enstatite -> akermanite` reaction is unavailable in every case where it is needed, and
the only other balanced route makes wollastonite.

**Agreement with the hand-rolled norm.** 153/153 cells, no failures. Quartz, anorthite, enstatite
and fayalite match to ≤0.01 wt%; **107/153 agree to <0.5 wt% on every phase — exactly the 107 cells
with no silica deficit.** The 46 that differ do so only in how the deficit is absorbed, up to
6.2 wt% on diopside in the Mg/Si ≥ 1.7 corner already flagged ultracalcic. The Earth crust is
unchanged. `check_crust_table.py`: ALL CHECKS PASSED.

**Cost.** pyrolite's CIPW runs at 1.64 s/call against 0.033 ms hand-rolled — a factor of ~50,000.
`mineral_composition` is now `lru_cache`d on (Mg/Si, ΔIW), giving 4.4 µs/call after the first;
without it a sweep would pay the full cost for every run at an identical composition.

**Independent corroboration worth recording:** pyrolite puts normative larnite in **exactly the
same 46 cells** the hand-rolled cascade needed åkermanite for. Two independent implementations
identify the same compositions as silica-deficient.

### 25.11 Hedenbergite — evaluated, not adopted

Asked whether the Fe-pyroxene corrections could be removed by adding the phase instead. The answer
is yes, and it is worth doing, but as a physics decision rather than a tidying step:

- **Thermodynamics exist** for both `Hedenbergite` (CaFe(SiO₃)₂, log_k 19.606) and `Ferrosilite`
  (FeSiO₃, log_k 7.4471) in `llnl.dat` — same lineage as the other primary silicates — and
  `make_database.py`'s `SPECIES_ALIAS` shim already handles the SiO₂→H₄SiO₄ basis conversion. This
  is the §23.7 nepheline workflow.
- **Kinetics would still be proxied**, but far better than åkermanite's: `Augite_ss` has a MEASURED
  rate in Kinec_v3 (Aa 1.52e6, Ea 81834, na 0.7, plus a neutral term) and its composition
  Mg₀.₄₅Fe₀.₂₇₅Ca₀.₂₇₅SiO₃ is a real Fe-bearing clinopyroxene between diopside and hedenbergite.
- **It is not a bookkeeping change.** k_eff at 300 K / pH 6: Augite_ss 1.78e-12 against Fayalite
  3.79e-09 — **2,135×**. At the Earth anchor ~half the crust's iron (clinoferrosilite 7.17 wt% vs
  fayalite 7.78) would change host and drop three orders in dissolution rate. That is arguably MORE
  correct — real basaltic cpx hosts iron and weathers slowly, so the all-Fe-to-fayalite convention
  over-delivers — but §20 was the iron charge-leak fix and §22 is the calibration authority, so both
  need re-checking before it is adopted.

Deferred deliberately. The write-up does not need it: the correction list is already four reactions
plus two fallbacks, and two of them are silica-neutral Fe–Mg exchanges describable in a clause.

### 25.12 The runtime database is not the one the build script makes (2026-08-24)

Prompted by the user asking directly whether the model had always run on `hybrid_ocean.dat` rather
than on `make_database.py`'s output. It had. **§11 has been corrected in place.**

**The finding.** `make_database()` writes `lt_weathering_sit.dat`. `chemistry.py:33` loads
`hybrid_ocean.dat`. These are not two versions of one database:

| | `hybrid_ocean.dat` (runtime) | `lt_weathering_sit.dat` |
|---|---|---|
| activity model | **Pitzer** | **SIT** (ThermoChimie) |
| base | `pitzer.dat` + `ocean_chem.dat` | `sit.dat` |
| PHASES | 85 | 1772 (after rebuild) |
| RATES blocks | 0 (kinetics live in `mineral_rates.py`) | 8 |
| built by | `make_hybrid.py` — **NOT IN THE REPOSITORY** | `make_database.py` |
| carries the §22 calibration | **yes** | no |

> ⚠️ **The runtime database is unreproducible.** `hybrid_ocean.dat`'s own header says "Generated by
> make_hybrid.py"; that script does not exist anywhere in the repo or its history. If the file is
> lost or corrupted, every calibration in §22 goes with it. Reconstructing `make_hybrid.py` from
> the file's header while the provenance is still legible is the cheapest insurance available.

**How it went unnoticed.** Three symptoms were each explained away locally instead of together:

- §23.7 added Nepheline to `hybrid_ocean.dat` **and** to `make_database.ADDED_PHASES` "so the
  documented build path stays complete" — but that build path feeds the file nothing loads.
- §25.5 found Akermanite already in `hybrid_ocean.dat` and absent from `ADDED_PHASES`, and recorded
  it as a "latent bug" in the bookkeeping. It was actually a symptom of this.
- §25.9 found `C(cr)` (graphite) in `lt_weathering_sit.dat` but **not** in `hybrid_ocean.dat` —
  the same divergence seen from the other direction, and again not connected.
- §19 noted zeolites "live only in `lt_weathering_sit.dat`, not the runtime `hybrid_ocean.dat`.
  Not worth pursuing." The divergence was seen and dismissed.

**`lt_weathering_sit.dat` is STALE, not broken.** It loads in PHREEQC and was missing exactly three
of the phases the model needs — Nepheline, Akermanite, Hedenbergite — all of which `ADDED_PHASES`
already knows about. Rebuilding it (`make_database(name=...)`) produces **1772 phases, nothing
missing, loads cleanly, and solves end-to-end** with all nine primary phases equilibrating. The
reproducible database has been one command away from complete the whole time.

**Switching is not a swap, and should not be done casually.** Pitzer and SIT give different
activity coefficients at the ionic strengths these oceans reach (I ≈ 3–4 mol/kg), which is the
regime the model lives in and presumably why Pitzer was chosen. Every number in §22 — `kd_mg_ht`,
`k_na`, the LT/HT flux targets against Coogan & Dosso — was derived on the Pitzer database and would
need re-deriving. The 1772-phase list also means `available_mineral_string` and the precipitating
sets need re-auditing: §19.7's exclusions were reasoned against 85 phases. Note also the §11 naming
trap still applies — silica is `H4(SiO4)` in sit.dat, `SiO2` in llnl/Kinec, `H4SiO4` in
`hybrid_ocean.dat`.

### 25.13 Hedenbergite and the Augite rate — tested, wired, left OFF

§25.11 deferred this; the user asked for it to be tested. It works, and the effect is much larger
than the Earth anchor alone suggested.

**What was added.** `Hedenbergite` (CaFeSi₂O₆) and `Ferrosilite` (FeSiO₃) PHASES blocks, taken from
`llnl.dat` and rewritten onto the H₄SiO₄ basis exactly as Diopside is (`CaFeSi2O6 + 4 H+ = Ca+2 +
Fe+2 - 2 H2O + 2 H4SiO4`, log_k 19.606), inserted **inside** the PHASES block before `PITZER` per
§23.9. Plus `augite_k` in `mineral_rates.py` from Kinec_v3's measured `Augite_ss` RATES block
(Aa 1.52e6, Ea 81834, na 0.7, An 350, En 83000), used as the proxy rate for **both** phases.

**Validated with a control**, which is the §23.9 lesson applied: the test database solves and
returns saturation indices for both phases, while the unmodified runtime database returns
`ERROR: Phase not found in database, Hedenbergite` on the identical input. The phases are genuinely
new and genuinely visible to PHREEQC, not merely to the Python parser.

**The effect, measured across the grid** (`EMIT_FE_PYROXENE` on vs off, k at 300 K / pH 6):

| Mg/Si | Hd wt% | Fs wt% | Δ k[Alk] | Δ k[Fe] |
|---|---|---|---|---|
| 0.50 | 2.96 | 6.00 | **−97.4%** | −100% |
| 0.70 | 6.47 | 9.73 | **−97.9%** | −100% |
| 1.00 | 8.13 | 4.32 | −43.5% | −51% |
| **1.25 (Earth)** | 7.21 | 0.07 | **−21.4%** | −27% |
| 1.60 | 8.64 | 0 | −15.7% | −24% |
| 2.00 | 5.83 | 0 | −12.5% | −22% |

Zero failures anywhere. But at the silicic end essentially **all** iron moves from fast fayalite
into slow ferrosilite and the crust goes nearly inert. That is arguably the right petrology —
granites do weather slowly — and it would greatly increase the compositional lever §23.1 found too
weak. It also introduces a ~40× swing in weathering sensitivity across the Mg/Si axis from a single
modelling choice.

> ⚠️ **The silicic-end collapse rests on a stretched proxy.** `Augite_ss` is a CLINOpyroxene rate
> being applied to an ORTHOpyroxene, and it is 4× slower than enstatite (1.78e-12 vs 7.04e-12 at
> 300 K / pH 6). The −97% figures at Mg/Si ≤ 0.7 are therefore the least trustworthy numbers in the
> table. A measured ferrosilite rate would settle it; none exists in any available database
> (`mineral_rates.py`'s header already recorded "Ferrosilite — PHASES in basic_v2.dat but no
> kinetic data found").

**Left OFF at the time.** `crust_composition.EMIT_FE_PYROXENE = False`, so the norm still
exchanged iron into olivine and nothing changed: `check_crust_table.py` passed unaltered with the
phases present but unused. Adoption was judged to require the §22 recalibration first.

> ✅ **Superseded 2026-08-25 (§27).** Both phases are now emitted unconditionally and the flag has
> been **removed entirely**, along with the two Fe correction reactions it guarded. Adoption ran
> ahead of the §22 recalibration rather than behind it, so that recalibration is now outstanding.
> The measured climate effect is −1.0 K at S = 0.6–0.8 and −8.2 K at S = 1.0 (3 km, Mg/Si 1.25);
> the Mg/Si 0.5 case survives its 2135× rate cut and converges normally. Note the Δk[Alk] column
> above is a *rate-constant* ratio and badly overstates the climate response — see §27.6.

**Note the convergence with §25.12.** Adopting hedenbergite forces a §22 recalibration; switching to
the reproducible SIT database also forces one. Doing both together costs one calibration exercise
instead of two, and ends with the model on a database that can be rebuilt.

**State left on disk (all uncommitted, nothing decided):**

- `hybrid_ocean.dat` — Hedenbergite + Ferrosilite added. This is an edit to the *unreproducible*
  database, i.e. the practice that caused §25.12. Original backed up outside the repo.
- `lt_weathering_TEST.dat` — scratch output of the rebuild test, sitting in `src/kamino/data/`.
  Delete it.
- `lt_weathering_sit.dat` — **still stale**; deliberately not regenerated pending a decision.
- `mineral_rates.py` — `augite_k` added and wired for both phases in `K_FUNCTIONS`.
- `mineral_info.py` — molar masses and `primary_minerals` entries.
- `make_database.py` — `ADDED_PHASES` entries (which, per §25.12, only affect the unused database).

### 25.14 The activity model: measured, and the runtime database made reproducible (2026-08-24)

§25.12 established that the runtime database was unreproducible. This section resolves it, and
answers the question that should have been asked first: **does the activity model actually change
the model's output?**

**It barely does.** Identical inputs, Earth crust, w/r = 3, across T = 275-315 K and
pCO₂ = 10⁻³-10⁻¹ bar:

| quantity | SIT / Pitzer |
|---|---|
| alkalinity weathering flux | **median 1.09x**, range 0.82-1.41 |
| feedback strength `d ln F / dT` | agrees to **~15%** (0.034 vs 0.034 at 10⁻¹ bar) |
| Walker exponent `d ln F / d ln pCO₂` | both ~0 (Pitzer +0.02..+0.08, SIT −0.06..−0.01) |

The pCO₂ = 10⁻⁴ column disagrees by up to 9.6x, but that is a **Pitzer convergence artefact**, not
a real difference: Pitzer's own values there are non-monotonic in temperature (−1.03e-10 at 288 K,
−2.25e-11 at 300 K). It is far below any habitable pCO₂.

> **Why the 65x equilibrium difference does not propagate — the important part.** §25.12 flagged
> b_eq[Ca] = 400 mM under SIT against 6.2 mM under Pitzer as alarming. It is irrelevant, and the
> Maher & Chamberlain form says why:
>
> ```
> F = A_r (b_eq − b_in) / (b_eq/k + A_r/J)
> ```
>
> When b_eq is large the `b_eq/k` term dominates the denominator and **F → A_r·k** — the flux
> becomes independent of b_eq entirely. The model sits in the KINETICALLY-limited (low Damköhler)
> regime, governed by `mineral_rates.py`, which is database-independent. The activity model only
> sets b_eq, and b_eq divides out.
>
> This also predicts exactly where the two WOULD diverge: the thermodynamically-limited, high-Da
> regime. If the model is ever pushed there — see §8's Da investigation — the choice matters again
> and must be re-tested.

**Cost, measured.** Identical PHREEQC solve: Pitzer 2.2 ms, SIT 20.2 ms — **9.1x**. That maps
straight onto the model: 0.57 s per `dY_dt` (~7 solves) and ~8 s per LSODA Jacobian, which is why
a 10 Myr integration would not finish in 500 s when fast_18 did the full 2 Gyr in 6-7 s.

> ⚠️ **The 9x is NOT overhead and cannot be trimmed.** Tested: a trivial pure-water solve costs
> 1.17 ms on SIT against 0.05 ms on Pitzer, so the fixed database-size cost is only ~1.1 ms of the
> 20.2. PHREEQC speciates only the elements actually present, so sit.dat's actinides and organics
> were never being computed. The cost is the complexation network for the model's OWN elements:
> **87 aqueous species against Pitzer's 28**, and 87/28 = 3.1, squared = 9.6, matching the observed
> 8.8x for a Newton solve on a larger system. A filter that drops other elements recovers ~5% and
> breaks the database's dependency graph ("Elements in species have not been tabulated"). The
> attempt and its measurements are recorded in `make_database.py`. **Do not retry it.**

**Conclusion: Pitzer is the default.** 9x faster, for a ~10% difference in the quantity that
matters. The aluminium argument for SIT is real but weaker than it looks — Al is thermodynamically
load-bearing, but the model is kinetically limited, so it barely reaches the output.

#### 25.14.1 `make_database.py` now builds either, and both are reproducible

`make_database(base="pitzer")` (default) or `base="sit"`. The activity model is a documented
parameter rather than an accident of which file was on disk. **`hybrid_ocean.dat` is no longer
referenced by anything**; `chemistry.py` loads `lt_weathering_pitzer.dat`.

Validation against the old hand-built database it replaces:

| database | ms | Albite SI | CO₂(g) SI |
|---|---|---|---|
| `hybrid_ocean.dat` (old, hand-built) | 3.64 | +1.6784 | −9.9539 |
| **`lt_weathering_pitzer.dat` (new)** | 4.18 | **+1.6181** | **−9.8448** |
| `lt_weathering_sit.dat` | 33.7 | +1.3775 | −10.1529 |

Albite within 0.06 log units and CO₂ within 0.11 of the database every §22 calibration was derived
on — close enough that the recalibration should be an adjustment, not a rebuild.

**What the pitzer base needs, and where it comes from** (all bundled with the `phreeqc` package, so
the build needs no file outside the installed environment):

- **Aluminium.** `pitzer.dat` defines 24 elements and **Al is not one of them**. Master species and
  five hydrolysis species are grafted from `Kinec_v3.dat` — the same source `ocean_chem.dat` cites
  for its thermodynamics, and the master-species line is byte-identical across `Kinec_v3.dat`,
  `llnl.dat` and the old `hybrid_ocean.dat`. Kinec_v3's ORGANIC Al complexes (Al(CH3COO)2+ etc.)
  are deliberately excluded: acetate is not a tracked species.
  > ⚠️ This supplies Al **thermodynamics only**. The Pitzer framework has no Al interaction
  > parameters, so Al³⁺ receives the long-range electrostatic term alone. Say so wherever Al
  > speciation is quoted from a Pitzer-based run. This is structural, not an oversight —
  > `pitzer.dat` is not an Al database.
- **Ferric iron.** `pitzer.dat` has only Fe⁺². Goethite (`FeOOH + 3 H+ = Fe+3 + 2 H2O`) needs the
  redox pair, grafted from **`phreeqc.dat`** — not Kinec_v3, which writes the couple against O₂
  (`H+ + Fe+2 + 0.25 O2 = Fe+3 + 0.5 H2O`) where the model's bookkeeping uses the electron form.
  sit.dat already defines `Fe(+3)`, so the graft is skipped there automatically.
- **Five phases** sit.dat has and pitzer.dat lacks: Siderite, Kaolinite, Goethite, Saponite-Na,
  Greenalite. Added base-conditionally — adding them unconditionally would silently override
  ThermoChimie's own versions in the SIT build.

#### 25.14.2 Three silent traps, now handled in code

Each of these produced a *successful-looking* build that failed later, which is what made them
expensive:

1. > ⚠️ **`END` truncates a database.** `pitzer.dat` closes its PITZER block with `END`, and
   > PHREEQC stops reading there. Every appended phase, rate and graft was invisible: the file
   > loaded with **rc = 0 and no error**, and the first symptom was "Phase not found in database"
   > at run time. This is §23.9's trap in a new place, and it is almost certainly why the original
   > `make_hybrid.py` interleaved phases into the base (PHASES at line 259, END at 1145) rather
   > than appending them. The builder now strips standalone `END` lines. **Stripping it also
   > revealed two real errors that had been hidden behind it** — the Fe⁺³ and silica problems
   > below were present all along and simply never read.
2. **Silica naming, in both directions.** `H4(SiO4)` in sit.dat, `H4SiO4` in pitzer.dat and the old
   hybrid_ocean.dat, `SiO2` in llnl/Kinec/Thermoddem phases. §11 records the trap; what it does not
   record is that the alias must point the right way. PHREEQC defines a species on the **right** of
   a reaction whose left side is already known, so the base's master goes on the left. Writing
   `H4(SiO4) = H4SiO4` against a Pitzer base yields "Reaction for species has not been defined".
   The alias block is now generated from whichever base is selected.
3. **Gas-phase naming.** SIT/ThermoChimie uses `O2(g)`, `H2(g)`, `N2(g)`, `CH4(g)`; the Pitzer
   database uses `Oxg(g)`, `Hdg(g)`, `Ntg(g)`, `Mtg(g)`. Requesting an absent gas makes PHREEQC
   reject the **whole `SELECTED_OUTPUT` block**, so the first symptom was a `KeyError: 'si_CO2(g)'`
   three steps removed from the cause. `chemistry.py` now intersects its gas request with what the
   loaded database defines, and points `OXYGEN_GAS_PHASE` at whichever name exists.

#### 25.14.3 Other changes to `chemistry.py`

- **`KAMINO_LT_DATABASE`** environment override, so the two activity models can be A/B tested
  without editing code. This is how §25.14's flux comparison was run.
- **`parse_stoichiometry(keep=...)`.** The SIT base is the full ThermoChimie set; parsing all 1772
  phases demanded a species map for `Am+3`, `PuO2+2`, `Edta-4` and 87 others. Parsing is now
  restricted to the ~24 phases the model can use.
- **`_FUSED_RE`.** Thermoddem writes coefficients fused to species (`2.000H+`, `1.000H4SiO4`);
  PHREEQC accepts this, `_iter_terms` did not. Wollastonite was the phase that exposed it.
- **`-saturation_indices` filtered.** §19 flagged this list as unfiltered at 84 phases; at 1772 it
  would have had PHREEQC computing 1772 saturation indices on every solve.
- **`_OPTIONAL_MINERALS`.** Hedenbergite/Ferrosilite were only reachable behind
  `EMIT_FE_PYROXENE`, so their absence from a database was not an error.
  **Removed 2026-08-25 (§27.1):** both phases are now always emitted, so the exemption would have
  let a database regression silently drop the crust's iron. They are required.

#### 25.14.4 Cross-cutting lessons

- **Measure the thing the model outputs, not the thing the model computes.** Two days of argument
  about activity models and aluminium turned on b_eq — which divides out of the flux entirely in
  the regime this model runs in. The 10-20% flux agreement settled a question that theory had been
  making look decisive.
- **A silent success is worse than a loud failure.** The `END` trap loads with rc = 0. It hid two
  genuine errors for the whole session, and only surfaced them when the truncation was removed.
- **"Reproducible" and "correct activity model" were tangled and are separable.** The original
  motivation for switching to SIT was reproducibility. Once `make_database.py` could build either,
  reproducibility stopped being an argument for SIT at all, and the choice could be made on the
  measurement.

---

## 26. Precipitation timescales: why `tau_prec` scales with ocean depth and `tau_rw` does not (2026-08-25)

Context: the sweep needed for the Mg/Si and ΔIW figures is dominated by water worlds. The 20-run
pilot measured **27.2 min for 10 shallow (3 km) runs against 154.6 min for 10 deep (20 km) runs** —
deep is **5.7× shallow**, worse than the 3–4× estimated beforehand. Extrapolated to the 252-run
cross sweep that is ~38 CPU-hours, of which the water worlds are 85%. The deep runs are not doing
more chemistry; they are **stiffer**, taking 1232–1477 integrator steps against 561 for 3 km. This
section is the investigation into whether that stiffness can be removed without changing the
answer, and what it revealed about the two-bucket precipitation model.

All runs below: `S = 1.0`, `outgassing = 0.1`, `crust_production_rate = 1.0`, `alpha = 2`,
`kd_mg_ht = 0.02`, `k_na_cont_removal = 0.004`, `t_end = 2 Gyr`, all reaching `t_end`
(`termination = 'timeout'`, i.e. not converged and not wall-capped). Concentrations in mM.

### 26.1 The measurement: `tau_prec` alone

| depth | Mg/Si | `tau_prec` | steps | wall | T | pH | pCO₂ | Ca | Mg | Alk |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 km | 1.25 | 100 kyr | 1232 | 587 s | 344.09 | 6.199 | 0.5525 | 3.112 | 17.29 | 37.89 |
| 20 km | 1.25 | 700 kyr | 384 | 140 s | 344.38 | 6.175 | 0.5867 | 3.311 | 17.26 | 38.27 |
| 20 km | 1.25 | 3 Myr | 328 | 98 s | 344.99 | 6.116 | 0.6753 | 3.874 | 17.20 | 39.28 |
| 20 km | 0.50 | 100 kyr | 1477 | 614 s | 355.27 | 5.459 | 4.2736 | 10.015 | 12.22 | 41.56 |
| 20 km | 0.50 | 700 kyr | 243 | 52 s | 355.34 | 5.456 | 4.3084 | 10.128 | 12.12 | 41.64 |
| 20 km | 0.50 | 3 Myr | 250 | 48 s | 355.58 | 5.444 | 4.4275 | 10.509 | 11.88 | 41.93 |

At 700 kyr — the mass-proportional value for a 20 km ocean — the cost falls **3.2× (Mg/Si 1.25)**
and **6.1× (Mg/Si 0.5)** in steps, for **ΔT = +0.29 K** and **+0.07 K** respectively.

A prediction made before the second test was **wrong and is recorded as such**: the acidic Mg/Si
0.5 corner was expected to be the *more* sensitive, on the reasoning that it sits nearer saturation
boundaries. It is the *less* sensitive by a factor of ~6 on Ca (+1.1% against +6.4%), and it is
also the cheaper one. Proximity-to-saturation was the wrong intuition; see §26.3 for the right one.

### 26.2 Why scaling `tau_rw` with it is wrong

The obvious tidy move — scale both timescales with ocean mass, preserving their ratio — was tested
at `tau_prec = 667 kyr`, `tau_rw = 33.3 Myr` (ratio held at 50×):

| depth | Mg/Si | config | steps | T | pCO₂ | Ca | Alk |
|---|---|---|---|---|---|---|---|
| 20 km | 1.25 | 100 kyr / 5 Myr | 1232 | 344.09 | 0.5525 | 3.112 | 37.89 |
| 20 km | 1.25 | 667 kyr / 33.3 Myr | 310 | **319.26** | **0.0338** | 1.110 | 35.44 |
| 20 km | 0.50 | 100 kyr / 5 Myr | 1477 | 355.27 | 4.2736 | 10.015 | 41.56 |
| 20 km | 0.50 | 667 kyr / 33.3 Myr | 269 | 355.34 | 4.3062 | 10.122 | 41.63 |

**Preserving the ratio is what breaks it.** At Mg/Si 1.25 the planet cools **25 K** and pCO₂ falls
**17×**. The mechanism is the carbon budget re-closing: slowing reverse weathering leaves more
alkalinity available to calcite (ocean Ca falls 3.112 → 1.110 as it is consumed), so steady state
moves to a far lower pCO₂. At Mg/Si 0.5 it changes nothing at all, because reverse weathering is
inactive there — Sepiolite(d) is undersaturated at every timescale tested.

### 26.3 The actual distinction: which phases reach equilibrium

Saturation indices recomputed at each attractor's final ocean state (`experiments/probe_saturation.py`;
the two precipitation calls in `dY_dt` are independent PHREEQC equilibrations on the same
`b_ocean`, so the timescales never compete inside one solve — they act through the shared ocean
state over time):

| `tau_prec` / `tau_rw`, 20 km Mg/Si 1.25 | Calcite | SiO2(am) | Kaolinite | Sepiolite(d) | Sep. flux |
|---|---|---|---|---|---|
| 100 kyr / 5 Myr | **+0.004** | **+0.005** | +1.154 | **+3.778** | 1.565e−6 |
| 700 kyr / 5 Myr | +0.021 | +0.026 | +2.508 | +3.776 | 1.621e−6 |
| 3 Myr / 5 Myr | +0.065 | +0.077 | +3.386 | +3.780 | 1.764e−6 |
| 667 kyr / 33.3 Myr | +0.050 | +0.014 | +1.235 | **+8.799** | **1.946e−7** |

This is the whole result in one table.

- The **fast bucket sits at SI ≈ 0**. Those phases genuinely reach saturation, so their flux is set
  by solute **supply**, not by the timescale. `tau_prec` is a numerical relaxation constant and its
  value barely enters the answer — which is why it can be scaled.
- **Sepiolite(d) sits at SI ≈ +3.8**, four log units supersaturated, at *every* `tau_prec`. It never
  reaches equilibrium, so `tau_rw` **is** the reverse-weathering flux, not a relaxation constant.
  Scale it and the flux scales with it (SI +8.8, flux down 8.3×), and the climate follows.

That asymmetry — not the ratio between them — is what the two-bucket model encodes.

### 26.4 Correction to §21

§21 states that *"Amorphous silica precipitates on `tau_prec` = 100 kyr … the reverse-weathering
list runs on `tau_rw` = 5 Myr, **50× slower**, so the silica is gone before it gets a look.
Sepiolite(d) doesn't precipitate at all."*

That is **true for the configuration §21 was written about** (shallow, `crust_production_rate =
0.01`, where the Mg budget showed reverse weathering at −2.7e−06 Tmol/yr) and **false as a general
statement**. At 20 km and Mg/Si 1.25 Sepiolite(d) is already precipitating at the 100 kyr baseline
(SI +3.778, flux 1.565e−6), and moving to 700 kyr changes its flux by **3.6%**. The "50×
separation" was an observation about a symptom in one regime, and was mistakenly carried forward in
this session as a design constraint before being tested. **The ratio is not the load-bearing
quantity.** What matters is that each bucket stays on the correct side of saturation.

### 26.5 Unifying the two timescales, and why it fails

Given the above, the natural simplification is a single timescale — easier to justify in the
methods. It cannot be unified at 100 kyr (that speeds reverse weathering 50×, and asserts that
authigenic clays form as readily as calcite, contrary to observation), so the candidate is
`tau_prec = tau_rw = 5 Myr`:

| depth | Mg/Si | config | steps | T | ΔT | pH | Ca | Calcite SI |
|---|---|---|---|---|---|---|---|---|
| 20 km | 0.50 | baseline | 1477 | 355.27 | — | 5.459 | 10.015 | +0.004 |
| 20 km | 0.50 | 5 / 5 Myr | 239 | 355.74 | +0.47 | 5.436 | 10.783 | +0.035 |
| 20 km | 1.25 | baseline | 1232 | 344.09 | — | 6.199 | 3.112 | +0.004 |
| 20 km | 1.25 | 5 / 5 Myr | 362 | 345.29 | +1.20 | 6.083 | 4.248 | +0.092 |
| **3 km** | 1.25 | baseline | 561 | 316.66 | — | 7.666 | **0.300** | — |
| **3 km** | 1.25 | 5 / 5 Myr | 262 | 323.88 | **+7.22** | **6.946** | **3.326** | **+1.219** |

**The shallow ocean breaks it.** At 3 km the ocean ends up **16× supersaturated in calcite**
(SI +1.219; SiO2(am) +0.388), i.e. the fast bucket stops buffering entirely. Ca then rises
**11-fold**, pH falls 0.72, and T rises 7.2 K. In a small reservoir, 5 Myr is simply not fast
enough to keep up with solute supply.

So the fast bucket is "numerically irrelevant" *only while it stays at SI ≈ 0*, and shallow oceans
are where it stops. **The two timescales stay.** The justification for the methods is one sentence
with two citations — carbonate and amorphous silica precipitate rapidly and hold the ocean at
saturation, while authigenic clay formation is kinetically inhibited and proceeds on Myr timescales
(Michalopoulos & Aller 1995; Isson & Planavsky 2018) — and the model's own SI values (calcite
+0.004 against Sepiolite +3.8) are the demonstration that the two buckets behave differently.

### 26.6 What was implemented

**`planet.py`.** `tau_prec` now defaults to `None` and is resolved in `__init__`:

```python
TAU_PREC_REF = 100e3 * YR
OCEAN_DEPTH_REF = 3000.0   # m; the depth at which TAU_PREC_REF applies
...
if tau_prec is None:
    tau_prec = TAU_PREC_REF * (ocean_depth / OCEAN_DEPTH_REF)
```

| depth | resolved `tau_prec` | `tau_rw` |
|---|---|---|
| 3 km | 100,000 yr (unchanged) | 5 Myr |
| 10 km | 333,333 yr | 5 Myr |
| 20 km | 666,667 yr | 5 Myr |

Resolution happens **before** `planet_config` is built, so saved JSONs record the value actually
used rather than `None`. An explicit `tau_prec=` still overrides, so pinned reproductions and the
tau tests above still work. `tau_rw` is deliberately untouched (§26.2).

**`parameter_sweep.py`.** `WALL_SECONDS_DEEP` raised **1800 → 2700 s**. The pilot's single
`wall_timeout` was `S = 1.0, depth = 20 km, Mg/Si = 1.6, ΔIW = −2.0`, which reached
**t = 1.72 of 2.0 Gyr (86%)** and was still converging — a truncation of real physics, which is
what these budgets exist to prevent. The headroom is cheap now that deep runs take 3–6× fewer
steps.

Expected sweep cost: **~38 → ~17 CPU-hours** for the 252-run cross, i.e. ~2.5 h on the EPYC 7251 at
7 workers rather than ~5.5 h.

### 26.7 Consequences and open items

- **Shallow (3 km) results are bit-identical**; nothing about them changes.
- **Deep-ocean results already on disk — including the pilot — were produced at
  `tau_prec = 100 kyr` and will not reproduce.** They must be regenerated, not merged.
- **Name collision introduced.** `planet.py:40–41` carry `_ABIOTIC_CA_3MYR` and
  `_TAU_PREC_REF_K = 3e6 * YR`, both provably unreferenced anywhere in the repo, and the latter now
  sits 27 lines above `TAU_PREC_REF` with a similar name and a different value. Worth deleting.
- **Still open for the ΔIW figure:** ΔIW = −1 at Mg/Si 1.25 contains åkermanite, whose proxy
  bracket carries a 48× uncertainty on Ca (§25) — on the very axis that figure is about.

### 26.8 Cross-cutting lessons

- **A ratio between two parameters is not automatically meaningful.** The 50× separation looked
  load-bearing because §21 mentioned it in the same breath as a real mechanism. Preserving it was
  the single most damaging thing tested here (25 K); ignoring it cost 0.29 K.
- **Saturation index tells you whether a timescale is physics or numerics.** A phase at SI ≈ 0 has
  reached equilibrium and its rate constant is free; a phase held at SI ≫ 0 is kinetically limited
  and its rate constant *is* the answer. This is a cheap diagnostic and should be the first thing
  checked before tuning any relaxation timescale in this model.
- **Test the regime you think is safe, not just the one you think is risky.** The unification was
  ruled out by the 3 km case, which was included only as a control — the deep cases it was aimed at
  both passed.
- **Notes written about one configuration need their scope recorded.** §21's claim was accurate
  where it was made and misleading everywhere else, and it cost time in this session before being
  tested.

---

## 27. Adopting hedenbergite: iron leaves fayalite for good (2026-08-25)

§25.13 added Hedenbergite and Ferrosilite to the database, measured their effect, and left them
**off** behind `EMIT_FE_PYROXENE = False` pending the §22 recalibration. This section adopts them,
removes the flag, and records what changed. The trigger was a review of whether the phases were
actually wired up — they were, but dormant, and dormant behind a flag that did not work correctly.

### 27.1 What the audit found

Everything downstream of the norm was correct: the database entry
(`CaFe(SiO3)2 + 4 H+ = Ca+2 + Fe+2 + 2 H2O + 2 SiO2`, log_k 19.606, grafted from llnl.dat),
the parsed stoichiometry (`4 alk, 2 Si, 1 Fe, 1 Ca`), the molar mass, `primary_minerals`
membership (dissolve-only, correct for an igneous phase), and the `augite_k` rate. A 50 Myr run
with the flag forced on completed with **zero chemistry fallbacks**.

Two defects sat above it.

**The `lru_cache` did not key on the flag.** `_mineral_composition_cached` keys on
`(mantle_mg_si, delta_iw, cipw_items)`, but `EMIT_FE_PYROXENE` was read as a module global inside
`cipw_norm`. Flipping it after any composition had been computed returned the **stale assemblage,
silently**:

```
flag flipped to True, no cache_clear -> Hedenbergite = 0.0     (correct answer: 0.1573)
identical to the False result? True
```

Any script written to compare the two settings — precisely the script needed to decide whether to
adopt them — would have produced a null result with no error. This is the §25.14 "silent success is
worse than a loud failure" lesson recurring in a different disguise.

**`_OPTIONAL_MINERALS` was load-bearing in the wrong direction.** `chemistry.py` exempted
Hedenbergite and Ferrosilite from the missing-mineral check *because* they sat behind an
off-by-default flag. Had they been switched on while a database lacked them, the crust's iron would
have been dropped from the reactive assemblage silently rather than raising.

### 27.2 The effect on the crust

Iron is **conserved exactly** — 1.085 mmol/g at the Earth anchor either way. Only its *host*
changes, and therefore its dissolution rate. `augite_k` is 1.78e−12 at 300 K / pH 6 against
fayalite's 3.79e−09, a factor of **2135**.

The effect is strongly **non-uniform**, because it depends on how much normative olivine the melt
carries to exchange iron with:

| Mg/Si | ΔIW | Fe in pyroxene | Fe-weighted mean k falls by |
|---|---|---|---|
| 1.25 | −2 | 27% | 1.4× |
| 1.60 | −2 | 24% | 1.3× |
| 0.80 | −2 | 68% | 3.2× |
| **0.50** | **−2** | **100%** | **2135×** |

Below Mg/Si ~0.8 the melts are silica-oversaturated and carry no normative olivine at all, so the
old fallback route (`2 ferrosilite -> fayalite + SiO2`) put **every** ferrous atom in the
fastest-dissolving host in the database. Now none of it is there.

Assemblage at the Earth anchor (Mg/Si 1.25, ΔIW −2), weight fractions:

| phase | before | after | k (300 K, pH 6) |
|---|---|---|---|
| Anorthite | 0.3543 | 0.3543 | 9.44e−12 |
| Diopside | 0.2447 | 0.1817 | 1.55e−11 |
| Albite | 0.1490 | 0.1490 | 1.01e−11 |
| Forsterite | 0.1395 | 0.1603 | 4.49e−10 |
| Fayalite | 0.1106 | 0.0804 | 3.79e−09 |
| **Hedenbergite** | 0 | **0.0721** | 1.78e−12 |
| **Ferrosilite** | 0 | 0.0007 | 1.78e−12 |

### 27.3 The effect on climate

Matched pairs, 3 km ocean, `alpha = 2`, `kd_mg_ht = 0.02`, `k_na = 0.004`, `t_end = 2 Gyr`,
all reaching `t_end`. Concentrations in mM.

| Mg/Si | S | | T | pH | pCO₂ | Ca | Mg | Alk |
|---|---|---|---|---|---|---|---|---|
| 1.25 | 0.6 | before | 294.96 | 6.290 | 0.8021 | 8.126 | — | 37.53 |
| 1.25 | 0.6 | after | **293.94** | 6.305 | 0.7702 | 8.185 | 18.87 | 37.66 |
| 1.25 | 0.8 | before | 298.02 | 6.844 | 0.1630 | 2.757 | — | 26.81 |
| 1.25 | 0.8 | after | **296.93** | 6.870 | 0.1515 | 2.728 | 18.93 | 26.76 |
| 1.25 | 1.0 | before | 316.66 | 7.666 | 0.0255 | 0.300 | 18.77 | 21.62 |
| 1.25 | 1.0 | after | **308.46** | 7.813 | 0.0090 | 0.462 | 13.72 | 11.84 |
| 0.50 | 1.0 | before | 350.57 | 5.929 | 2.4385 | 3.392 | 16.58 | 23.42 |
| 0.50 | 1.0 | after | **339.87** | 6.278 | 0.2869 | 6.331 | 5.44 | 7.02 |

Three things worth recording.

1. **The cold end barely moves** (−1.0 K at S = 0.6 and 0.8). Weathering there is not limited by
   the iron-bearing phases.
2. **The warm end moves ~8–11 K colder**, via a genuine drop in crust reactivity: alkalinity falls
   45% at Mg/Si 1.25 and **70%** at Mg/Si 0.5, and ocean Mg falls 27% and 67%. Ocean Fe is ~0 in
   every run, before and after — so this is **not** an iron-chemistry effect. Removing fayalite,
   the fastest phase in the assemblage, simply makes the crust a poorer weathering substrate.
3. **S = 1.0 moved toward the historical value, not away from it.** fast_18 reported ~308–313 K for
   this configuration; the fayalite-exchange crust gave 316.66 and the Fe-pyroxene crust gives
   308.46. The §25.12 database validation had flagged that overshoot as the one point of
   disagreement with fast_18. It is now smaller.

**The Mg/Si 0.5 case does not break**, which was the open worry: despite a 2135× cut to the
iron-weighted rate, the run converges normally to a sensible 339.87 K with pH 6.28. The 2135× is
real but iron was never carrying the carbon budget.

### 27.4 What was implemented

The flag is **removed entirely**, not merely defaulted on — there is no configuration in which
exchanging iron into fayalite is preferable now that both endmembers exist with a measured proxy
rate.

- **`crust_composition.py`.** `EMIT_FE_PYROXENE` and the `emit_fe_pyroxene` parameter deleted.
  `clinoferrosilite -> Hedenbergite` and `ferrosilite -> Ferrosilite` moved into `_PYROLITE_DIRECT`
  with the other endmembers that map straight onto a database phase. The two Fe correction
  reactions are **deleted as dead code** — they only ever worked around the database lacking
  Fe-pyroxene. `cipw_norm` now converts **one** normative phase (larnite -> åkermanite) rather than
  three, and its docstring is correspondingly shorter. The cache-key defect disappears with the
  global that caused it.
- **`chemistry.py`.** `_OPTIONAL_MINERALS` removed; both phases are now **required**, so a database
  regression fails loudly instead of silently dropping the crust's iron.
- **`planet.py`.** No change — an `emit_fe_pyroxene` field was added to `planet_config` while the
  option still existed, then reverted, because recording a constant is noise.

Verified: no references remain anywhere in the code; all three modules compile; `chemistry.py`
imports clean under the strict check (so both phases really are in the runtime Pitzer database);
and across the full Mg/Si 0.5–2.0 × ΔIW −5…−1 box every assemblage sums to 1.0 with no negative
phases and conserved iron.

### 27.5 Open items this creates

- **§22 recalibration is now due, not optional.** §25.13 said adoption should follow it; adoption
  has happened first. The Earth anchor moved −8.2 K at S = 1.0, so `alpha` and `kd_mg_ht` are
  anchored to a crust the model no longer produces. §25.12 notes the database switch forces the
  same exercise — doing both together still costs one calibration, not two.
- **The ferrosilite rate is the weakest link and is now load-bearing.** `Augite_ss` is a
  CLINOpyroxene rate applied to an ORTHOpyroxene. At Mg/Si 0.5, ΔIW −1 the crust is **39 wt%
  ferrosilite**, so the low-Mg/Si end of the Mg/Si figure rests on it. No measured ferrosilite rate
  exists in any available database. This belongs in the methods as a stated caveat.
- **Every result on disk predates this.** All sweeps, the pilot, and the §26 tau runs used the
  fayalite-exchange crust and will not reproduce.
- **The compositional lever got stronger.** §23.1 found the Mg/Si lever too weak. At S = 1.0 the
  Mg/Si 0.5 → 1.25 temperature difference is now 339.87 − 308.46 = **31.4 K**, against
  350.57 − 316.66 = 33.9 K before — similar in span but reached with far lower pCO₂ and alkalinity
  at both ends.

### 27.6 Cross-cutting lessons

- **"Implemented" and "in use" are different claims, and only one of them was true.** Every
  component was correct and the path was dead. The audit that found this was prompted by a question
  about whether the phase worked, not by any failure.
- **A flag with no valid off-setting is a liability.** It cost a silent cache bug and a
  database check exempted in the wrong direction. Deleting it removed more code than it added.
- **A 2135× change in a rate constant produced an 11 K change in climate.** Iron was never carrying
  the carbon budget; what mattered was that fayalite is the fastest-dissolving phase in the
  assemblage, so moving mass out of it lowers total crust reactivity. Rate-constant ratios are not
  a proxy for model sensitivity.
