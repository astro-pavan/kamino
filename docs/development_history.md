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
et al.; it supersedes §5, §23 and most of §24, though §24.6's literature refocus stands. **§32 now
supersedes §25 on how the melt is COMPUTED** — isobaric batch melting at 1 GPa replaces the
isentrope, the grid is 26 × 25, and `T_p` becomes `T_melt`; §25's norm, akermanite closure and
validation against the literature all stand unchanged. §22 remains the authority on
CALIBRATION — it records the Earth calibration, shows `alpha` is not identifiable from Earth, and
establishes that the LT seafloor flux carries the wrong ion relative to Coogan & Dosso. §22.0 also
corrects two stale entries in §11/§15.

**§33 is the most recent work and supersedes §22/§27 on CALIBRATION**: the seafloor-area fix (§33.3) changes the land-bearing physics, so `KD_MG_CALIB` and `K_NA_CALIB` must be re-fitted before any ocean chemistry at `land_fraction > 0` is quoted. Land-free results are unaffected and verified bit-identical.

**§37 (2026-09-24) is the most recent work.** It records the review of the paper draft (issues tracked in `paper_issues.md` at the repo root) and the model changes it led to: area-scaled shelf precipitation, explicit SO₄ and K⁺ backgrounds in the calibration, crust production 1/130 Myr, a recalibration with the constants moved into `constants.py`, and **no Cl in the sweeps** (Cl kept for the Earth calibration only). It supersedes §3's constants table, the shelf term as described in §3/§33, §35's SO₄ background and §35's chlorine seeding for the sweeps, and corrects §33.4's attribution of the OLR fit. §37.19 replaces the `r_avg` convergence check with a windowed drift test (~15× faster sweeps, same results), and §37.18 sets the default outgassing to 1×. §37.20 replaces the failed joint calibration with an alternating fit (α ≈ 356), and §37.21 calibrates τ_rw to the modern authigenic-clay Mg sink. §37.22 reruns the ocean pilot on the refitted constants: the RWR transition moves to S 0.9–0.95. §37.23 sets the shelf depth to 140 m. §37.24 records the full rerun and the plotting fixes it needed, §37.25 the 'earth' sweep for the calibration figures, and §37.26 the comparison with the old seeded sweep.

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
| **Aug 24** | **The two-parameter crust pipeline** (§25): Mg/Si × ΔIW adopted as the axes, åkermanite closes the norm, three-source validation; the norm switched to pyrolite; the runtime database found unreproducible and made reproducible (§25.12–25.14). |
| **Aug 25** | **Precipitation timescales** (§26): `tau_prec` depth-scaled, `tau_rw` deliberately not. **Hedenbergite adopted** (§27) — iron leaves fayalite, moving the Earth anchor −8.2 K. |
| **Aug 26–27** | **Earth recalibration** (§28) after hedenbergite; `alpha` shown unidentifiable from Earth a second time; **the ocean was silently oxidising** and `pe = −3.0` adopted as the default (§28.3). First full sweep read (§29). Analysis tooling made reproducible and ~30× faster (§30). **Every sweep now runs in both redox states** (§31). |
| **Sep 1** | **`alpha` decision and a second refit** (§34) — production moved to `ALPHA_REF` itself so the sweep and module default cannot drift. ⚠️ *Reconstructed from code on 09-09; no session record exists.* |
| **Sep 3** | **Crust pipeline audited** (§32): the Mg/Si 0.5 / ΔIW −1 "vanishing quartz" traced to a real phase boundary aliased by the ΔIW axis; the isentrope shown redundant under batch melting and **replaced by an isobaric solve** (577 → 219 lines, 130 s → 1 s per point); grid 17 × 9 → **26 × 25**; `check_crust_table.py` found broken since §25.13 and fixed. |
| **Sep 7–8** | **Continental weathering and the crossover** (§33): `continental_baseline.py` rewritten as an Earth-like instellation sweep; a post-runaway hot-branch state found being counted as habitable; the **seafloor-area fix** (§33.3) which invalidates the Earth calibration; land-fraction series, coarse grid and alpha sweep; a TTG second-melt diagnostic showing high Mg/Si cannot make felsic continents; **alpha measured as the largest control** (f\* ∝ α^0.80–0.86). |
| **Sep 9** | **Charge balance and the Cl root cause** (§35): the carbon excess traced to the conservative-ion charge residual, τ_Cl = 5.6 Gyr against a 2 Gyr integration, `K_CL_ANALYTIC` found wrong since §33.3, SO₄ pinned at zero in every sweep; **sedimentation rate extended to all precipitating phases**; **MORB verified against Gale et al.** |

| **Sep 17** | **The sink audit** (§36): reverse weathering shown to be a bimodal switch and a CO₂ *source*, negligible at Earth but the only backstop once `kd_mg_ht` is removed; Greenalite shown inert (pore Goethite starves the ocean of Fe); no Mg carbonate can fill the gap; **crustal Albite shown to contribute exactly zero Na** (supersaturated above ~1 mM Na + `dissolve_only`), quantifying §6.2's switch; nahcolite ruled out as a Na sink. `basic_no_rw` added to `parameter_sweep.py`, unrun. |

| **Sep 24** | **Paper review and the changes it forced** (§37): draft checked line by line against the code (`paper_issues.md`), equations corrected in the draft; **shelf precipitation made area-scaled** (it was a second whole-ocean sink burying 56% of Earth's carbonate); **SO₄ 28.2 / K⁺ 10.2 mM explicit** in the calibration (K added as a fixed element); crust production 1/130 Myr; **recalibration** (η_HT 2.361e-2, η_Na 5.775e-3, α 14.57; Mg +15%); constants consolidated in `constants.py`; 30-run Cl test → **sweeps run with no Cl**. |

*Rows Aug 24 – Sep 1 were missing entirely; Sep 3 and Sep 7–8 were inverted. Both fixed 2026-09-09.*

---

## 3. Model architecture as it now stands

### State vector

> **Superseded 2026-09-24 (§37.6, §37.19).** K is now the last element, and `r_avg` is gone: the state is
> `[P_CO2, P_H2O, *elements]`.

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

> **Superseded 2026-09-24 (§37.7–37.8).** The fitted constants now live only in
> `src/kamino/constants.py`: `KD_MG_HT = 2.361032e-02`, `K_NA_CONT_REMOVAL = 5.775040e-03`,
> `K_CL_SUBDUCTION = 1.961786e-04`, `ALPHA_REF = 14.57`. `parameter_sweep.py` takes them from there, so
> the sweep and module defaults can no longer drift. The table below is kept as history.

> **Updated 2026-09-09.** This table had drifted three refits behind the code (it recorded
> `KD_MG_HT = 0.07`, `K_NA = 2.194806e-03`, `alpha = 1.43` — the §18/§22-era values). The values
> below are read from the working copy. **All three of the fitted constants are stale in a
> different sense: they must be re-fitted after §33.3's seafloor-area fix (§35).**

| Parameter | Value | Role |
|---|---|---|
| `KD_MG_HT` | `1.394755e-02` | HT Mg→Ca exchange; §28.1 gave `1.394362e-02`, refit 2026-09-01 (§34) |
| `K_CL_SUBDUCTION` | `1.373251e-04` | Cl subduction. ⚠️ Its analytic derivation is **wrong** post-§33.3 (§35.2) |
| `K_NA_CONT_REMOVAL` | `4.234317e-03` | Na sink (always-on, `J_total`-scaled); §28.1 gave `4.272026e-03` |
| `alpha` (`weathering.ALPHA_REF`) | `1.100155` | base reactive area per unit crust area; §28.1 gave `0.487612`, refit 2026-09-01 (§34) |
| `f_HT` | `0.0` | vestigial |
| `pe` | **−3.0** | ocean/pore redox; added 2026-08-27, see §28.3. Was silently PHREEQC's default of +4 |
| `tau_prec` | **100 kyr x (ocean_depth / 3 km)** | fast precipitation; depth-scaled 2026-08-25, see §26 |
| `tau_rw` | 5 Myr | reverse weathering; deliberately **not** depth-scaled, see §26.2 |
| `TAU_ATM` | 10 kyr | atmosphere relaxation |
| `P_CO2_CLIMATE_FLOOR` | 1.0 Pa | climate-input clamp (§9.2) |
| `convergence_threshold` | 0.05 /Gyr | over a 50 Myr `convergence_window` since §37.19; was `event_converged` |

> ⚠️ **Memory conflict.** `memory/project_calibration_state.md` records `KD_MG_HT=0.233`,
> `K_NA=3.33e-3`, `f_HT=0.039`, and a convergence threshold of 0.5. **None of these match the
> current code** (values above). That memory predates the Aug 3 switch-back and should be treated
> as stale history, not current calibration.

> ⚠️ **`parameter_sweep.py` pins its own copies and they do not match.** `KD_MG_CALIB` and
> `K_NA_CALIB` (`parameter_sweep.py:34-35`) still hold the §28.1 values, so
> `_warn_constant_drift()` fires on both today and every run is filename-tagged `_kmg…_kna…`
> accordingly. `ALPHA_CALIB = ALPHA_REF` by construction since §34, so alpha cannot drift.
> **After the recalibration, update `planet.py` and `parameter_sweep.py` together** or the new
> sweep is again non-comparable to the module default.

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

> **2026-09-24 additions (§37), all uncommitted:** area-scaled shelf precipitation
> (`Planet.shelf_area_fraction`); K added to `elements` as a fixed background, with SO₄ = 28.2 and
> K⁺ = 10.2 mM in the calibration seed; `EARTH_CRUST_PRODUCTION_RATE_PER_AREA = 1/(130 Myr)`; calibrated
> constants moved to `constants.py`; sweeps run with no Cl (`parameter_sweep.CL_OUTGASSING_RATIO = 0`, run
> names tagged `_cl0`); `get_ocean_state` receives `pe`; sweeps start from a blank ocean by default
> (`KAMINO_SEED_OCEAN`, 2026-09-23); seafloor temperature floor `SEAFLOOR_T_FLOOR = 273.15` K in
> `constants.py` (was 274 K); calibration α anchor = net seafloor flux at 0.9 Teq/yr (§37.16, refit pending);
> sedimentation from cited dust, cosmic-dust floor and porosity 0.7 (§37.17); windowed-drift convergence on a
> stepped LSODA loop, no `r_avg`, pinned SO₄/K Jacobian columns skipped (§37.19); alternating calibration
> (ions at fixed α, then α and τ_rw rescaled to their flux targets), `TAU_RW_REF` in `constants.py` and
> `rw_mg_flux` in the run diagnostics (§37.20–37.21); shelf depth 140 m (§37.23); `t_end_yr` in the output,
> Cl ratio / run length pinned in `plot_results`, per-run `seed` in `run_simulation`, and the 'earth' sweep
> in `continental_baseline.py` (§37.24–37.25). Entries below that conflict
> with these are history.

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
- ✅ Dead Jacobian columns — **applied 2026-09-24** (§37.19): SO₄ and K skipped, `r_avg` removed
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
- ✅ **Sedimentation rate counts every precipitating phase** (added 2026-09-09, §35.4). Was
  carbon-as-calcite plus silicon-as-**quartz**-density only; now each ocean-precipitating mineral
  contributes its own volume via `mineral_info.PRECIPITATE_MOLAR_MASS` / `PRECIPITATE_DENSITY` and
  `precipitation.sediment_volume_rate`. `get_precipitation_by_mineral` returns a fourth item
  (per-mineral molar rates); `get_precipitation`'s 3-tuple signature is unchanged.
- ❌ **The ocean is not seeded.** `time_evolve` zero-fills `Y0` and no sweep passes `b0`, so every
  run starts from an empty ocean and SO₄ — which is *pinned* (`planet.py:466`) and has no source
  term — stays at 0 forever. `calibrate_earth.py` seeds 23.45 mM. **The calibration and the sweeps
  are running different oceans** (§35.1, §35.3).
- ❌ **`K_CL_ANALYTIC` is wrong post-§33.3.** It assumes the Cl source and sink areas cancel, which
  stopped being true at the seafloor-area fix (§35.2).

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

> ✅ **RESOLVED — the redox half, checked 2026-09-09.** The concern was that §28.1 fitted at the
> implicit `pe = 4` while §28.3 moved the default to −3.0. `calibrate_earth.py` does **not** pass
> `pe` at all, so it inherits `planet.PE_DEFAULT = -3.0` (`planet.py:96`), which is exactly
> `PE_REDUCING == PE_DEFAULT_SWEEP` (`parameter_sweep.py:85-87`). Calibration and the production
> sweep arm anchor at the same redox. No action; any refit inherits this automatically.
>
> ✅ **RESOLVED — the `alpha` half, 2026-09-01 (§34).** Production runs at `ALPHA_REF` itself
> rather than a separately pinned round number, so the sweep and the module default cannot drift.
> The *finding* below still stands and is not weakened by the decision: **`alpha` is not
> identifiable from Earth** (concentrations move < 6% across a 41× change) while the land-free
> sweeps are kinetically limited (Da ~ 0.005), where `F ∝ alpha` linearly. The alpha arm is now
> (`ALPHA_REF`, 10, 50). The feedback STRENGTH is alpha-invariant to 7% over 40×, which is what
> the composition figures report.
>
> 🔴 **Superseding both: the whole anchor is stale for a third reason** — §33.3's seafloor-area
> fix, which `KD_MG_CALIB` and `K_NA_CALIB` absorbed as a 1.43×. See §35.

1. **Raise or taper the CO₂ ceiling** to the maximum-greenhouse peak (§9.2/§13), instellation-
   dependent rather than a flat 10 bar. The thermodynamic branch *is* the rising-pCO₂ branch, so it
   hits the wall almost immediately and only ~1 point past each transition survives. This is now the
   main limit on settled fraction (~26%) and on seeing each transition line in full (§20.5).
2. **Run a full-res sweep with the iron fix** (§20.3–20.4) to see the transition across many lines,
   not just the 4 that cleared ≥5 settled points in low-res fast_15.
3. **Seed Cl analytically in the sweep** (§20.4). ⬆️ **PROMOTED to the top of this list and
   root-caused, 2026-09-09 (§35.1).** This is not merely a settling nuisance: it is the single
   largest error in the model's ocean chemistry. Measured τ_Cl = **5571 Myr** against a 2 Gyr
   integration, so a blank start reaches `1 − e^(−2/5.57) = 0.302` of steady state — 0.302 × 780 =
   **235 mM**, the sweep's value to three digits. Because alkalinity is charge-derived (§7), that
   315 mEq Cl deficit *becomes* carbonate alkalinity, and it is 54% of why DIC is 180× Earth's.
   `calibrate_earth.py` seeds; **no sweep does** (`time_evolve` zero-fills, `parameter_sweep.py:241`
   passes no `b0`).
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
   which is what the calcite-saturation bistability requires. ~~**New top-priority item in its place:
   decide `alpha` (§22.9)**~~ — **decided 2026-09-01 (§34): production runs at `ALPHA_REF`.** It
   remains unidentified *by Earth*; the decision is which convention to run, not a measurement.
   **The top-priority item is now item 3 (seed the ocean) plus the §35 recalibration.**
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

> ✅ **DECIDED 2026-09-01 (§34): production runs at `ALPHA_REF`, currently 1.100155.** That is
> option 3 below — the best Earth fit — so read option 3's warning as a live caveat on the
> production configuration, not as a rejected branch. The rest of this section is the reasoning
> as it stood, kept because the argument that Earth cannot identify `alpha` is unchanged.

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

---

## 28. The Earth recalibration, the `alpha` problem, and the redox assumption nobody set (2026-08-26 → 08-27)

Three things, in the order they happened, because each one forced the next. §27 invalidated the §22
calibration, so it was re-run; the re-run exposed that `alpha` is unidentifiable in a way that
matters for the paper; and a question about whether Goethite can form on an anoxic planet revealed
that the model had been silently assuming an oxygenated ocean since it was written.

### 28.1 The recalibration

`experiments/calibrate_earth.py`, re-run after §27. Converged on its own tolerance at **25 of 60**
evaluations, best cost 0.1077.

```
K_CL_SUBDUCTION   = 1.373251e-04   (analytic, unchanged)
K_NA_CONT_REMOVAL = 4.272026e-03   (was 2.194806e-03;  §22 got 3.904380e-03)
KD_MG_HT          = 1.394362e-02   (was 7.0e-02;       §22 got 1.898657e-02)
ALPHA_REF         = 0.487612       (was 1.43;          §22 got 0.908383)
```

> ⚠️ **SUPERSEDED by the 2026-09-01 refit (§34).** The shipped values are now `K_NA =
> 4.234317e-03`, `KD_MG_HT = 1.394755e-02`, `ALPHA_REF = 1.100155`. `K_CL_SUBDUCTION` is
> unchanged — and §35.2 shows its analytic derivation has been *wrong* since §33.3. Note
> `parameter_sweep.py` still pins the §28.1 `K_NA`/`KD_MG` above, so the drift check fires
> today (§3). Both are stale again regardless: see §35.

Earth: `converged`, **T = 294.4 K, pH 7.76, pCO₂ 694 ppm**, `|dlnP/dlnt| = 0.002`, zero fabricated
derivatives.

| | Na | Ca | **Mg** | Alk | C |
|---|---|---|---|---|---|
| this run | −8.3% | +3.0% | **+37.0%** | +33.1% | +32.1% |
| §22 | −0.9% | +2.9% | +4.0% | +24% | +25% |

**The Mg residual is a consequence of §27, not a solver failure.** The deleted correction reaction
was `hedenbergite + ½ forsterite → diopside + ½ fayalite`, which consumed forsterite (Mg-only) and
produced **diopside** (Ca and Mg, 1:1) — it was manufacturing a Ca source out of an Mg-only mineral.
Removing it changes the Earth-anchor assemblage:

| phase | before §27 | after | |
|---|---|---|---|
| Diopside | 0.2447 | 0.1817 | **−25.7%** |
| Forsterite | 0.1395 | 0.1603 | **+14.9%** |

Rate-weighted at Earth conditions: **Ca supply −9.5%, Mg supply +12.3%, Ca/Mg ratio −19.4%**. The
calibration's own split metric moved **−24%** (0.990 → 0.752). Those two independent numbers agree,
so the residual is fully explained by the crust change. `kd_mg_ht` can only trade Mg *for* Ca
mole-for-mole, so with Ca on target there is no way to pull Mg down; the solver stopped at the
compromise where Ca wins.

One fix was needed in the script itself: `TAU_PREC_INIT` was pinned at 100 kyr while §26 made
`tau_prec` depth-scaled, so at Earth's 3700 m the calibration was fitting at a timescale the model
never uses there. It now resolves through `planet.TAU_PREC_REF` (123 kyr at that depth).

### 28.2 `alpha` cannot be identified from Earth, and the sweeps run where it matters most

> **Status 2026-09-01 (§34):** the *finding* here is unchanged and was re-measured, not overturned.
> What changed is the operational choice — production moved from a pinned `alpha = 2` to `ALPHA_REF`
> itself. Everything below about identifiability still holds.

The fit reports `alpha = 0.4876`, but that number is not a measurement. Paired evaluations differing
only in `alpha` return identical oceans, exactly as §22.2 found. Measured directly: across a **41×
change in `alpha`** (0.4876 → 20, with `K_na`/`kd_mg` fixed at the calibrated values), Earth
concentrations move **< 6%**:

| alpha | Na | Ca | Mg | seafloor Alk |
|---|---|---|---|---|
| 0.4876 | −8.3% | +3.0% | +37.0% | 0.11 Teq/yr |
| 20 | −9.7% | +5.7% | +42.6% | 7.57 Tmol/yr primary |

The reason is structural, and it is the worst possible arrangement for this paper: **Earth is
transport-limited, the water worlds are not.** Measured over 19 land-free pilot states, the
Damköhler number has **median 0.005 and 0/19 above 1**, so `F → A_r·k ∝ alpha` linearly, with
elasticity `d ln F/d ln alpha` = **0.998 median**. Earth cannot see `alpha` because continental
weathering dominates at `land_fraction = 0.3`; the sweeps run at land_fraction 0, where `alpha`
carries the entire thermostat.

**What survives: the FEEDBACK is alpha-invariant.** In the kinetic limit `alpha` is a multiplicative
constant on the flux, and a constant cancels out of `d ln F/dT`:

| alpha | F(300 K) | **d ln F/dT** |
|---|---|---|
| 0.4876 | 1.31e−11 | **0.0781** |
| 2.0 | 4.85e−11 | **0.0806** |
| 20.0 | 4.18e−10 | **0.0752** |
| 200.0 | 3.03e−09 | 0.0402 |

Flat to **7% over a 40× range**, degrading only near 200 where the system starts leaving the kinetic
limit. So `alpha` sets the absolute CO₂ offset, not the strength or the composition-dependence of
the feedback — which is what the Mg/Si and ΔIW figures report.

**A degeneracy that only half holds.** At steady state `alpha·g(T, pCO₂) = F_out`, so the solution
should depend on `F_out/alpha` alone. Tested at fixed ratio 0.05:

| alpha | outgassing | T |
|---|---|---|
| 0.4876 | 0.02438 | 306.73 |
| 2.0 | 0.100 | 308.46 |
| 20.0 | 1.000 | **327.83** |

4.1× in `alpha` costs 1.7 K, but the full 41× costs **21 K**. It breaks at the top end, where the
alpha-independent sinks (HT exchange and the Na sink, both ∝ `J_total`) stop rescaling. The
degeneracy is a low-alpha approximation, not a symmetry, and should not be leaned on.

**Production value: `alpha = 2`, chosen not fitted.** `calibrate_earth.py`'s own diagnostic puts the
~1 Tmol/yr primary seafloor anchor at `alpha = 10.4`, but adopting it censors the data:

| S | Mg/Si | alpha = 2 | alpha = 10 |
|---|---|---|---|
| 0.8 | 1.25 | 296.93 | **279.99** |
| 1.0 | 1.25 | 308.46 | 298.62 |
| 1.0 | 0.50 | 339.87 | 316.75 |
| 1.2 | 1.25 | out_of_domain | out_of_domain |

`alpha = 10` does not rescue the hot end and pushes the cold end under: `CROSS_INSTELLATION` starts
at 0.50, and a uniform −17 K puts S ≤ 0.7 at or below freezing, removing the cold half of the
feedback curve. The **alpha arm** (2, 10, 50 — all still in the kinetic limit, Da ≤ 0.13) carries the
argument instead, which is the stronger claim anyway.

> ⚠️ **The `ALPHA_REF` primary-dissolution anchor is 97% dissolved Fe²⁺.** §22.3 measured 88% and
> rejected the anchor on those grounds. §27 made it **worse**, not better: fayalite fell 27% but
> diopside — the main fast Ca source — fell 25.7% while slow anorthite was unchanged, so Ca's share
> collapsed faster than Fe's. At the calibrated state, `Fe 3.682 Tmol/yr (97.3% of charge)` against
> `Ca 0.001`. Anchoring `alpha` to Coogan's measured seafloor flux is therefore still unavailable:
> the model and the observation agree on a magnitude while disagreeing about which ion carries it.

### 28.3 The ocean was oxidising, and nobody had said so

The user asked whether an anoxic atmosphere would affect Goethite formation. It does, and the
question exposed a genuine defect.

**`chemistry._solution_block` set no `pe` and no `redox`**, so PHREEQC fell back to its own default
of **pe = 4.0** — firmly oxidising at seawater pH. Goethite (FeOOH, **ferric**) sat unconditionally
in `clay_minerals`, used for both pore and ocean precipitation. The model is **abiotic by
construction**: there is no oxygenic photosynthesis, so there is no source of free O₂, and its
oceans should be reducing.

The code already contained the contradiction. `chemistry.py:92` charges Fe as **+2** in the
alkalinity balance, justified as *"soluble/conservative under the anoxic conditions where it is
mobile; oxic Fe precipitates as Goethite"* — the anoxic charge convention and the oxic mineralogy at
once.

**`fO2` was plumbed through the whole module and did nothing.** Measured from 1 PAL down to 1e−14
PAL: identical fluxes, Goethite SI 2.939 at every value. The cause is structural — the Pitzer
database's only redox couple is `Fe+2 = Fe+3 + e-`, and its oxygen master `Oxg` is a decoupled inert
gas (`Oxg = Oxg`), so imposing `Oxg(g)` fugacity sets dissolved O₂ and touches nothing else.
`lt_weathering_sit.dat` does carry a proper `O(0)` master and would respond, at ~9× the runtime
(§25.12).

**`pe` works, and the scan is diagnostic** (Earth pore conditions, 50 µM Fe):

| pe | Goethite SI | Siderite SI | |
|---|---|---|---|
| +12 | +7.65 | −9.33 | modern oxic seawater |
| +4 | +7.64 | −1.34 | **PHREEQC's default = the old behaviour** |
| 0 | +5.65 | **+0.43** | siderite (ferrous) takes over |
| −3 | +3.66 | +0.43 | anoxic marine pore water |
| −6 | **−0.34** | +0.43 | goethite finally undersaturated |
| unset | **+7.64** | **−1.34** | identical to pe = 4 |

### 28.4 What anoxia does to the climate

`pe` is now a `Planet` parameter, **default −3.0**, plumbed through `_solution_block` →
`solve_solution` → `get_b_eq` / `get_precipitation` / `get_precipitation_by_mineral` →
`get_weathering_flux`, and into all four `planet.py` call sites (pore weathering, fast ocean,
reverse weathering, shelf). Recorded in `planet_config`.

End-to-end at S = 1.0, 3 km, Mg/Si 1.25, dIW −2, alpha 2:

| pe | T | pH | pCO₂ | ocean Fe | Mg | Alk | fallbacks |
|---|---|---|---|---|---|---|---|
| **None** | **309.79** | 7.868 | 0.0108 | 0 | 15.76 | 15.65 | 0 |
| **4.0** | **309.79** | 7.868 | 0.0108 | 0 | 15.76 | 15.65 | 0 |
| 0.0 | 320.61 | 7.671 | 0.0390 | 0 | 24.13 | 32.16 | 2 |
| −3.0 | **321.31** | 7.656 | 0.0418 | 0 | 24.62 | 33.18 | 0 |
| −6.0 | 321.31 | 7.656 | 0.0419 | 0.08 µM | 24.62 | 33.19 | 9 |

`None` reproduces `pe = 4.0` **to every digit** — the hidden assumption, demonstrated at full-model
level. Anoxia is worth **+11.5 K**.

Two properties make this a well-behaved parameter, unlike `alpha`:

- **It saturates below pe ≈ 0.** pe of 0, −3 and −6 all give 321.3 K, so the value does not need
  tuning; what matters is the binary oxic/anoxic distinction.
- **pe = −3 is the numerically cleanest** (0 fallbacks, against 2 at pe = 0 and 9 at pe = −6) and is
  the middle of the measured range for anoxic marine pore water.

Reference values: modern oxic seawater ~ **+12.5**; anoxic marine pore water ~ **−3 to −5**; Archean
ocean reconstructions ~ **−3 to 0**.

> **A correction recorded because the first answer was wrong.** An initial test removed Goethite by
> hand while leaving `pe` at PHREEQC's default of 4, and reported anoxia as **−19 to −30 K** with
> ocean Fe reaching ~1 mM. That configuration is inconsistent: at pe = 4 **Siderite is also
> undersaturated**, so iron had no sink at all and piled up unphysically. Done properly, Siderite
> (FeCO₃) takes over below pe ≈ 0 — the real Archean iron sink, already in `carbonate_minerals` —
> iron is still removed, and the answer is **+11.5 K, the opposite sign and a third the magnitude**.
> The mechanism is a stoichiometry swap: Goethite removes 2 eq alkalinity per Fe and no carbon;
> Siderite removes 2 eq alkalinity **and a mole of carbon**, so the carbon balance closes elsewhere.

**`pe` is a distinct quantity from the two dIW values.** Those set the oxygen fugacity of the
MANTLE (at core formation, and of the melt); `pe` is the ambient redox of the WATER–ROCK system.
They must not be conflated.

### 28.5 Consequences

- **Every run on disk predates this**, including the §28.1 calibration itself, which was performed at
  the implicit `pe = 4`. The Earth anchor should be re-fitted at the production `pe`.
- **§28.4 roughly doubles the influence of the crust dIW axis.** Iron only matters to the carbon
  cycle when it stays dissolved, and under anoxia ocean Fe spans 273 → 3022 µM across the dIW range.
  The §29 measurement that dIW is a weak, phase-boundary-gated control was made at `pe = 4` and is
  conditional on it.
- **The `fO2` argument threaded through `chemistry.py` remains inert on the Pitzer database.** It is
  left in place because it is correct on SIT, but it must not be mistaken for a working control.

### 28.6 Cross-cutting lessons

- **"No setting" is not "no assumption".** Leaving `pe` unset did not make the model agnostic about
  redox; it made it silently oxidising, because PHREEQC has a default. Every library default that is
  not written down is a modelling choice made by someone else.
- **A plumbed parameter is not a working parameter.** `fO2` was threaded through five functions and
  verified inert only when someone finally scanned it. The database, not the call signature, decides
  whether a control does anything.
- **An inconsistent test can invert a result.** Removing Goethite without setting `pe` produced a
  confident answer of the wrong sign. The check that caught it was asking which sink replaces the
  one being removed.
- **A parameter that saturates is worth more than one that is merely fitted.** `pe` has a
  well-measured natural range and its effect plateaus inside it; `alpha` has neither property, and
  the difference is why one is a default and the other needs a sensitivity arm.

---

## 29. What the first full sweep shows: Mg/Si dominates, ΔIW is a switch, depth inverts (2026-08-26)

2045 runs in `/home/pavan/PhD/sweep_output` — a basic (outgassing × crust) sweep, a 10-point depth
sweep, a near-factorial 8 × 7 composition sweep, and a 3-point alpha arm at the reference crust. All
at `alpha = 2`, the §28.1 constants, and **`pe = 4` (the pre-§28.3 implicit default)**, so every
number here is conditional on an oxidising ocean.

### 29.1 Mg/Si is ~5× stronger than ΔIW

Across 722 in-domain composition runs (3 km, out = 0.1, crust = 1):

| axis | median T range | mean |
|---|---|---|
| **Mg/Si**, at fixed (ΔIW, S) | **45.9 K** | 45.7 K |
| **ΔIW**, at fixed (Mg/Si, S) | **6.8 K** | 8.9 K |

The Mg/Si ordering is also **stable across redox**: Mg/Si 0.5 is the hottest line in every ΔIW
column and Mg/Si 2.0 the coldest, so the Mg/Si result does not depend on which ΔIW is chosen.

### 29.2 ΔIW is a phase-boundary switch, not a lever

The ΔIW effect is not a weak continuous trend — it is flat, punctuated by cliffs:

| Mg/Si | ΔIW T-range | assemblage across the axis |
|---|---|---|
| 1.00 | 3.1 K | unchanged |
| 1.75 | 2.6 K | unchanged |
| 2.00 | 1.6 K | unchanged |
| **0.80** | **16.6 K** | **quartz out, olivine in** |
| **1.50** | **13.3 K** | **åkermanite appears** |
| **1.25** | **12.4 K** | **åkermanite appears** |

**Mean 1.8 K where the mineral assemblage is unchanged; 11.9 K where a phase boundary is crossed.**

**One variable explains the whole grid.** Temperature tracks the rate-weighted Ca+Mg supply from the
crust with **r = −0.992** across all 44 populated cells, slope **−30.6 K per decade of supply**.
ΔIW controls mantle FeO, which changes fayalite and hedenbergite — but iron is not a carbon-cycle
cation (ocean Fe is ~0 in every run at `pe = 4`; it precipitates as Goethite). What matters is when
the extra FeO shifts the norm across silica saturation or the desilication threshold, because that
changes how much **Ca and Mg** the crust delivers.

> ⚠️ **Two of the three steps rest on the åkermanite proxy.** Swapping the §25 bracket:
>
> | step | slow | mid (default) | fast |
> |---|---|---|---|
> | Mg/Si 0.8, silica saturation | −16.8 K | −16.8 K | −16.8 K |
> | Mg/Si 1.25, åkermanite | **−1.2 K** | −4.6 K | **−21.6 K** |
> | Mg/Si 1.5, åkermanite | **−5.1 K** | −8.1 K | **−24.0 K** |
>
> The **silica-saturation step is proxy-independent** — quartz-out/olivine-in is a real
> petrological transition and the supply jumps 3.6× regardless. The two åkermanite steps span
> **−1 to −24 K** on an unmeasured rate and must be reported with the bracket, not as point values.

### 29.3 Ocean depth: three competing effects, and the sign of the temperature response flips

| depth | S = 0.5 | S = 0.8 | S = 1.0 |
|---|---|---|---|
| 3 km | 292.2 | 298.3 | 309.8 |
| 20 km | **266.7** | 308.2 | **343.5** |
| 30 km | **259.8** | 309.2 | **345.6** |

At **S = 1.0 deep oceans run 34 K hotter**; at **S = 0.5 they run 32 K colder**. Three mechanisms,
each verified independently:

1. **Weathering is areal, the ocean is volumetric.** The sink per unit ocean mass scales as
   area/mass = 1/depth, so a 50 km ocean has 1/17th the drawdown per kg of a 3 km ocean. This
   dominates when warm: pCO₂ goes 0.011 → 0.79 bar at S = 1.0.
2. **Pressure suppresses carbonate precipitation.** At identical composition and T, calcite SI falls
   monotonically with pore pressure — **+2.05 at 300 m, +0.02 at 30 km, −1.02 at 50 km**. Below
   ~30 km calcite stops precipitating, so Ca and alkalinity accumulate instead of being buried. This
   dominates when cold: at S = 0.5, 20 km, alkalinity reaches 102 mM against 57 at 3 km, holding
   carbon as DIC rather than atmospheric CO₂.
3. **Chloride is a clock, not an equilibrium.** `Cl × depth` converges to ~58,000 for every deep
   ocean but is only 15,900 at 300 m: deep oceans still hold *all* the Cl ever outgassed, while
   shallow ones have reached the source–sink balance. The relaxation time scales with ocean mass, so
   at 2 Gyr the deep ones are still filling.

So the deep ocean is a **CO₂ buffer whose sign depends on temperature**: warm → it releases (weak
sink), cold → it absorbs (no carbonate burial).

**Salinity changes composition, not just magnitude.** At S = 1.0, 3 km → 30 km: Cl falls
0.58 → 0.07 g/kg while dissolved **Si rises 0.40 → 8.75 g/kg**, becoming 68% of the salt. Shallow
oceans are chloride brines; deep oceans are silica solutions, because silica has no
pressure-suppressed sink the way calcite does.

**Deep-ocean points are the least trustworthy in the sweep**: 115 runs terminated `wall_timeout`,
concentrated at depth, and 50 km is non-monotonic against 30 km at S = 1.0. Treat ≥ 30 km as
indicative.

---

## 30. Analysis tooling: making the figures reproducible and fast (2026-08-26 → 08-27)

Not physics, but it was blocking the physics.

### 30.1 The plotting code was 30 minutes per invocation, and it was PHREEQC

`_add_diag_columns` re-ran the weathering chemistry for every plotted run to recover the Damköhler
number, which drives the solid/dashed line styling and is not in the run JSON.

| | |
|---|---|
| `_diag_from_json`, per run | **921 ms** |
| runs in the sweep | 2045 |
| **total** | **~31 min**, on every invocation |
| load_data + import + style | 3.4 s |

Over 99% of wall time; matplotlib was never involved. Two fixes, in order of correctness:

**`Planet.time_evolve` now writes the diagnostics it already computed.** `dY_dt` produces `Da`,
`pH_seafloor`, the ocean SI dict and `b_pore` on every step and was discarding them. The final-state
re-evaluation that `time_evolve` already performs (§27-era fix for `self._T`) now records:

```json
"diagnostics": {"da": …, "calcite_si": …, "ocean_si": …, "alk_flux": …, "pH_seafloor": …}
```

Cost: **one extra PHREEQC solve per run** — the pore-fluid calcite SI, the only one not part of an
equilibrium `dY_dt` already performs. Verified against the plotting recompute on the same file:
`da`, `ocean_si` and `pH` agree to 15 significant figures, `calcite_si` and `alk_flux` to 8 (solver
tolerance from reconstructing the state). `alk_flux` uses the same `A_SEAFLOOR_EARTH` normalisation
the plots use. Reading the block is **~300,000× faster** than recomputing it.

**A sidecar cache for runs that predate that.** `.plot_diag_cache.json` beside the runs, keyed on
each file's size and mtime so re-running a sweep invalidates its own entries. **31 min → 28 s.**
Verified exact: computing 40 runs with and without the cache gives `maxdiff = 0.000e+00` on all five
columns with matching NaN patterns.

Two bugs surfaced while doing this:

- **`_add_diag_columns` used the FIGURE output directory to locate the run JSONs.** Identical in
  `__main__`, so it worked there, but any caller rendering elsewhere got *every diagnostic silently
  NaN* — the bare `except Exception` swallowed the missing-file errors. Fixed with `RUN_PATH`,
  recorded by `load_data`.
- **`_save_diag_cache` overwrote instead of merging**, so a partial render shrank a complete cache to
  the subset that process touched (observed: 2045 → 40 entries).

### 30.2 `_recompute_T` deleted; the source-side fix works

`plot_results` had been recomputing surface temperature from (instellation, P_CO2) because
`self._T` was a `dY_dt` side effect that could be left on a Jacobian probe's value. `time_evolve`
now re-evaluates the final state, so this is redundant — measured over all 2045 runs: **mean |ΔT|
0.0013 K, max 0.30 K, zero runs above 0.5 K**, against the 11 K error it was written for.

It was also becoming harmful: it always called `get_T_surface_analytic` regardless of the run's
`climate_model`, so it would have silently substituted the analytic model for a clima-interpolator
run. Moved to `plot_legacy.recompute_T`, applied only under `--legacy` and only to analytic runs.

**Every run above 0.01 K is a `wall_timeout`**, and that is mechanistic: the deadline check sits at
the top of `dY_dt`, so when `time_evolve` makes its final-state re-call the deadline has already
passed, the call raises, and `self._T` keeps the aborted probe's value. A real (0.3 K) gap in the
source-side fix; clearing `self._wall_deadline` before that re-call would close it.

### 30.3 Legacy handling split out

`plot_results.py` now assumes the current run schema and nothing else. `experiments/plot_legacy.py`
holds a one-way **schema upgrade** (`upgrade(df)`) rather than a second copy of the plotting code —
old terminations (`snowball`, `hothouse`, `co2_ceiling`, `acid_ocean`, `co2_floor`) map onto
`out_of_domain` + a `domain_wall`, named crust compositions (`_comp_basalt_49`) map onto the
reference crust or are given NaN axes, and discontinued `crust_carbonate_content` runs are dropped.
After `upgrade()`, every current figure renders old data unchanged. Invoked with `--legacy`.

Two of those mappings are **deliberate reclassifications, not translations**, and are documented as
such: `acid_ocean` was not a CO₂ wall but did stop mid-evolution, and `co2_floor` used to count as a
snowball outcome but now becomes "unknown fate".

The current 2045-run sweep exercises **none** of the legacy paths.

### 30.4 Page-based figure sizing

Figures are now sized to the page rather than to whatever looked right on screen:

```python
COLUMN_WIDTH_IN = 240 / 72.27   # 3.32 in -- one MNRAS column
TEXT_WIDTH_IN   = 504 / 72.27   # 6.97 in -- both columns
figure_size(width='single'|'double', height=<inches>)
```

Every publication figure calls it; `diagnostic_size()` covers the wide on-screen grids, which are
explicitly exempt. `_add_figure_legend` measures the rendered legend and drops a column at a time
until it fits inside the panel block, so a single-column figure can never be widened by its own
legend.

**`bbox_inches='tight'` had to go** for publication figures: the style already sets
`constrained_layout`, and `tight` then re-cropped to the content box and came out **larger** than
requested (3.43 against 3.32 in). A figure wider than the column gets scaled down by
`\includegraphics`, shrinking the type below what the style chose.

The style file path was also resolved **relative to the working directory**, so the module raised
`OSError` at import from anywhere but the repository root. Now resolved against `__file__`, with
`KAMINO_PRESENTATION=1` selecting the presentation style.

### 30.5 Crust figures were missing the phases the model emits

Enumerated over all 153 table cells:

| mineral | cells | max wt | in `plot_crust_grid` | in `plot_crust_composition` |
|---|---|---|---|---|
| Hedenbergite | **153** | 0.157 | ✗ | ✗ |
| Ferrosilite | 83 | **0.394** | ✗ | ✗ |
| Akermanite | 46 | 0.117 | ✓ | ✗ |
| Quartz | 29 | 0.462 | ✓ | ✗ |

Hedenbergite is in **every cell** and both scripts dropped it; `plot_crust_composition` was missing
four phases, so its stacked bands did not sum to 1 and the shortfall was invisible. Derived
endmembers now take their parent's hue with a hatch (Akermanite ← Diopside, Hedenbergite ←
Diopside, Ferrosilite ← Enstatite), which keeps `plot_crust_grid`'s CVD-validated hue ring intact.

`plot_crust_composition.py` was additionally **broken outright** — it called
`oxide_composition(T_p, mg_si)` with the pre-ΔIW signature, filtered CSV columns such that
`float("fixed-F")` was attempted, and ignored the `delta_iw` axis entirely (collapsing 153 rows onto
the Mg/Si axis with 9 stacked at each point).

**0 of 153 cells are now mass-violating** — the åkermanite reroute (§25.5) closed that hole
completely, and the "beyond norm" exclusion path no longer triggers anywhere.

---

## 31. The sweep design after `pe`: every sweep in both redox states (2026-08-27)

§28.3 made ocean redox a parameter with no defensible single value — the model is abiotic and so
should be reducing, but every result on disk was produced at the implicit oxidising default, and the
difference is +11.5 K. `parameter_sweep.py` therefore runs **every sweep under both states** rather
than choosing:

```python
PE_REDUCING  = -3.0   # abiotic planet; Siderite (FeCO3) is the iron sink
PE_OXIDISING = +4.0   # oxygenated ocean; ferric Goethite strips dissolved Fe
PE_STATES = [PE_REDUCING, PE_OXIDISING]
```

`+4.0` rather than modern seawater's `+12.5` because it is the value the pre-2026-08-27 sweeps
implicitly ran at, so the oxidising arm reproduces them exactly. The iron system is saturated by
then anyway (Goethite SI +7.64 at pe 4 against +7.65 at pe 12).

**Three new sweeps resolve the axis rather than bracketing it**, on `pe_arm = [12, 4, 0, −3, −6]`,
which spans oxic seawater to below the Goethite saturation boundary (~ −5.8) and straddles the
Goethite→Siderite switch at pe ≈ 0:

| sweep | runs | what it answers |
|---|---|---|
| `pe` | 95 | the redox response curve at 3 km |
| `pe_deep` | 95 | the same at 20 km — worth separating, because Siderite is a CARBONATE and the deep ocean is where carbonate precipitation is pressure-suppressed (§29.3), so the redox switch and the depth effect may not be independent |
| `pe_composition` | 630 | **pe × Mg/Si and pe × ΔIW** — whether the composition signal survives the redox choice |

`pe_composition` is the one that matters most. §29.2 measured ΔIW as a weak, phase-boundary-gated
control, but that was at `pe = 4`, where Goethite strips dissolved iron so it cannot affect the
carbon cycle at all. Under anoxia iron stays in solution and ΔIW is the parameter that sets how much
there is, so **the ΔIW result may be substantially redox-dependent** and §29.2 should not be quoted
without that caveat until this sweep has run.

### 31.1 The resume trap this created, and the guard

A run at the model's own default (`pe = −3`) is deliberately **untagged**, so its filename is
identical to one from before `pe` existed — which was produced at the implicit `pe = 4`. With
`RERUN = False`, a reducing sweep pointed at the existing 2045-run directory would have found those
files and returned them as reducing results: **2045 oxidising runs silently relabelled**. This is
the fast_13 resume trap in a new costume.

`run_simulation` now checks the `pe` recorded in the JSON against what was requested and re-runs on
any mismatch, including the `ABSENT` case that identifies pre-`pe` output:

```
re-running planet_s_1.0_..._kna0.00427203: stored pe='ABSENT' != requested pe=-3.0
  (output predates the pe parameter?)
```

Verified end to end: the same configuration under both states produces **309.79 K (oxidising)** and
**321.31 K (reducing)** with distinct filenames, reproducing §28.4 exactly; and stripping the `pe`
field from a stored run triggers the guard while the correctly-tagged sibling is still reused.

### 31.2 Cost

Doubling every sweep doubles its cost. Measured per-run times from the pilot (2.72 min shallow,
15.46 min deep — the deep figure predates §26's `tau_prec` speedup and is an upper bound):

| sweep | runs | CPU-h | wall-h @ 3 |
|---|---|---|---|
| `cross` | 252 | 11.4 | 3.8 |
| `pe` | 95 | 4.3 | 1.4 |
| `pe_deep` | 95 | 24.5 | 8.2 |
| `pe_composition` | 630 | 28.6 | 9.5 |
| `cross_deep` | 252 | 64.9 | 21.6 |
| `alpha_composition` | 756 | 34.3 | 11.4 |
| `basic` | 1862 | 84.4 | 28.1 |
| `basic_deep` | 1862 | **479.8** | 159.9 |

`basic_deep` and `composition_deep` (411 CPU-h) are now firmly overnight-cluster jobs rather than
laptop jobs. The cheap, high-value order for the paper is **`cross` → `pe` → `pe_composition`**,
which is ~15 wall-hours at 3 workers and answers the three composition figures plus the redox
question.

---

## 32. The crust pipeline audited: the isentrope was redundant and the grid was aliasing (2026-09-03)

Triggered by distrust of the crust-composition generator, and specifically by one pie in
`output/crust_grid.png`: at Mg/Si = 0.5, ΔIW = −1, all normative quartz had vanished and been
replaced by ferrosilite, while every other cell in that column was quartz-rich. It looked like a
broken norm. It was not — but three real problems came out of chasing it, and the pipeline is now
half the size and 35× faster.

### 32.1 The anomaly is a real phase boundary, aliased by the ΔIW axis

The melt at that cell genuinely is what the norm says: SiO₂ 49.9 wt%, FeO 25.2, no free silica, so
pyrolite's own Fe/Mg pyroxene split puts nearly all the pyroxene in the Fe endmember. The jump was
upstream, in the table: between two *adjacent* grid cells the melt moved 21.5 wt% in SiO₂ and the
solved temperature 150 °C.

The mechanism is stoichiometric and needs no thermodynamics to see. Adding FeO at fixed molar Mg/Si
adds a **silica-consuming** component — Fe needs Si for pyroxene — without adding SiO₂, so the
source's free normative quartz is progressively eaten. Solving for where it reaches zero:

| Mg/Si | 0.50 | 0.60 | 0.70 | 0.80 | ≥ 0.90 |
|---|---|---|---|---|---|
| ΔIW at quartz-out | −1.05 | −1.34 | −1.80 | −2.82 | never |

A steep diagonal confined to Mg/Si ≤ 0.8 and ΔIW ∈ [−2.9, −1.0]. Once quartz and feldspar are gone
the rock no longer starts melting at a low-T near-eutectic, so the temperature needed for F = 0.20
climbs steeply. On a 0.5-spaced ΔIW axis that whole transition fell inside **one grid step** at
Mg/Si = 0.5, which is why it read as a discontinuity.

**This is the boundary §6.2 of `crust_composition.md` already claimed as a validation success** —
Guimond's (Mg+Fe)/Si ≈ 0.8 predictor, reproduced without being encoded. The suspicious corner was
the headline result, landing between two cells.

### 32.2 `check_crust_table.py` had not run since Hedenbergite went direct

The validator that should have caught a per-cell problem crashed on the first cell containing
Hedenbergite:

```
File "check_crust_table.py", line 65, in reconstruct
    for ox, k in FORMULA[m].items():
KeyError: 'Hedenbergite'
```

Its `FORMULA` table (used for the every-oxide mass-closure check) was never updated when §25.13's
`EMIT_FE_PYROXENE` work moved `clinoferrosilite → Hedenbergite` and `ferrosilite → Ferrosilite`
into `_PYROLITE_DIRECT`. Hedenbergite is present in **every** cell, so the script could not complete
a single run. Added `Hedenbergite` = CaFeSi₂O₆ and `Ferrosilite` = FeSiO₃ (both verified against
`MINERAL_MOLAR_MASS` to ≤0.002 g/mol); the §6.1 claims are reproducible again.

> **The lesson is the one §30.1 already taught in a different costume:** a validator that fails
> loudly still fails silently if nobody runs it. Two separate documented results — "153/153 cells
> mass-balance" and "0 of 153 cells are mass-violating" — were being quoted from a script that
> could not execute.

Also removed: `cipw_norm`'s `verbose=True` branch referenced undefined `hd`/`fs`, left behind by the
same change, so the debug path raised `NameError`.

### 32.3 The isentrope was doing nothing but manufacturing a label

`isentropic_melt` tracked an isentrope from 3.0 → 1.0 GPa in 2 kbar steps, secant-solving for the
temperature holding entropy constant at each one, inside a bisection on T_p. That is ~130 s per grid
point and the reason `--slice` sharding and `merge_crust_slices.py` existed at all.

It cannot affect the answer. The melting is **batch** — `minim(data, X, P, T)` is called with the
same bulk `X` at every step, melt is never extracted — so every state is a full equilibrium
minimisation of one fixed composition, and equilibrium is path-independent: a function of
(bulk, P, T) alone. Both closures end at 1 GPa with F = 0.20; at fixed bulk and pressure F is
monotonic in T, so that temperature is unique. The two **must** agree.

Measured across all 153 cells of the old grid, isobaric against isentropic:

| quantity | agreement |
|---|---|
| temperature (old `T_end` vs new `T_melt`) | mean **0.60 °C**, max 2.10 |
| melt oxides (SiO₂, FeO, MgO, Al₂O₃, CaO, Na₂O) | mean **0.006–0.026 wt%**, max 0.093 |
| residual assemblage identical | **151/153** (the two are trace cpx/spl within 1 °C of their limits) |

Within the two bisection tolerances. The isentrope's sole product was the potential-temperature
label; the melt composition never depended on it.

**The condition matters and should be quoted with the result:** this holds only for batch melting.
Under fractional melting the residue evolves, history is real, and the path would be required.

### 32.4 `make_crust_compositions.jl`: 577 → 219 lines

Rewritten for the isobaric procedure alone. Removed: `PSTART`/`DP`/`adiabat_T`/`T_at_entropy`/
`isentropic_melt` (the isentrope), `Tp_for_F`/`Tp_for_cpx_out`/`--closure`, `shard`/`--slice`,
`--probe`, `--calibrate`, `--validate`/`--bulk`, `--fixed-p`, `core_mass_fraction`, and the CSV
columns `T_end`, `closure`, `mg_number`, `delta_iw_melt`, `core_mass_fraction`, `warnings` — each
verified to be read by **zero** callers first. `merge_crust_slices.py` deleted with `--slice`.

Two things worth not re-deriving:

- **`T_p` → `T_melt`.** Under an isobaric closure there is no potential temperature, and keeping the
  old column name would have been exactly the mislabelling §31.1 was written about. `T_melt` is the
  melting temperature at 1 GPa (Earth: 1328 °C, against the 1383 °C the isentropic version reported
  as T_p — the difference is the adiabatic gradient plus latent heat). Three readers referenced the
  old name (`crust_composition.py`, `check_crust_table.py`, `plot_crust_composition.py`) and now
  accept either, so old tables still load.
- **`T_for_F` returns NaN when the target is unreachable below 2200 °C**, and `grid_point` warns if
  the converged F is more than 0.05 off target. The old `isobaric_melt` helper had neither guard and
  would have returned a bracket end as though it were a root.

### 32.5 The grid: 17 × 9 → 26 × 25, non-uniform on both axes

At ~1 s per point, density is essentially free, so the axes are now dense where §32.1 says the
assemblage changes and coarse where nothing happens:

| axis | coarse | dense |
|---|---|---|
| Mg/Si | 0.1 over 1.0–1.2, 1.9–2.0 | **0.05** over 0.5–0.9 and 1.3–1.8 |
| ΔIW | 0.5 over −5.0 to −3.0 | **0.1** over −2.9 to −1.0 |

650 cells in **9.0 min serial, 0 failures, 0 warnings** — against 5.5 h for 153. Uniform refinement
would have been waste: below ΔIW −3 the melt moves 0.08 wt% SiO₂ per half-unit.

The transition resolves into a monotonic ramp, which is the whole point:

```
ΔIW    -2.0   -1.5   -1.4   -1.3   -1.2   -1.1   -1.0
SiO₂   72.66  71.34  71.02  70.32  64.21  57.40  49.87     (Mg/Si = 0.5, wt%)
```

| | largest adjacent-cell step along ΔIW |
|---|---|
| old 17 × 9 | **21.5 wt%** SiO₂, 203 °C |
| new 26 × 25 | **7.5 wt%** SiO₂, 114 °C |

Interpolating the old grid across that interval blended a quartz rhyolite with an Fe-pyroxenite — a
mixture that was not itself a MAGEMin solution. Two side benefits: Guimond's (Mg+Fe)/Si ≈ 0.8
threshold is now located at **0.799** rather than bracketed as 0.702 → 0.894, and the
olivine-out/opx-out boundaries are measured at 0.05 resolution (0.55 and 1.55, both in the direction
their depleted-residue argument predicts).

`DIW_SHOW` in `plot_crust_grid.py` gained the 0.1 steps from −1.5 up, so the figure shows the ramp
instead of aliasing it in one row; the grid is now 10 × 10 = 100 pies.

### 32.6 The operational trap that cost the most time: juliaup and a full `$HOME`

Not physics, but it burned an hour and filled a 5 GB home directory twice.

`julia` on this machine is a **juliaup launcher**. juliaup stores toolchains in `<depot>/juliaup`,
and when `JULIAUP_DEPOT_PATH` is unset it falls back to `JULIA_DEPOT_PATH` — which `.bashrc` sets to
`/data/pt426/julia_depot`. But `.bashrc` line 7 is:

```bash
[ -z "$PS1" ] && return          # before line 41's export
```

So in any **non-interactive** shell (a tool, a hook, `ssh host julia …`, an IDE language server)
neither variable is set, juliaup defaults to `~/.julia/juliaup`, finds no toolchain, and downloads
~830 MB of Julia into `$HOME` — where it then fails to extract, because `$HOME` is 5 GB and nearly
full. A working 1.12.7 toolchain was on `/data` the whole time.

Fixed environment-independently, rather than by adding another export that some contexts would also
miss:

```
~/.julia -> /data/pt426/julia_depot
```

Verified with `env -i` (no Julia variables at all): `julia --version` → 1.12.7, no download, home
unchanged. Optional belt-and-braces, not applied: move the `JULIA_DEPOT_PATH` export above
`.bashrc`'s interactivity guard and add `JULIAUP_DEPOT_PATH` explicitly.

### 32.7 Mineral colours re-encoded by family

The pie chart's hues grouped the pyroxenes by STRUCTURE — clinopyroxene pink, orthopyroxene green
— with hatching marking derived endmembers and, per the caption, proxied kinetics. Two hues did
all the work and the texture meant something different in each case.

Now **hue is the mineral family and texture is the cation within it**, so eleven phases read as
six groups, and *dotted means calcium-bearing everywhere it appears*:

| family | hue | plain | dotted (Ca) | lines |
|---|---|---|---|---|
| quartz | amber | Quartz | | |
| feldspar / feldspathoid | green | Albite | **Anorthite** | Nepheline |
| pyroxene, Mg | light blue | Enstatite | **Diopside** | |
| pyroxene, Fe | dark blue | Ferrosilite | **Hedenbergite** | |
| melilite | purple | | **Akermanite** | |
| olivine | red | Fayalite | | Forsterite |

Nepheline is not a plagioclase; it takes the feldspar hue because it *is* the desilicated albite
(norm step 5b converts one into the other), so the cascade shows up as a texture change rather
than a hue change. Akermanite gets its own hue: it is a melilite, and giving it the Mg-pyroxene
hue plus dots would have made it identical to Diopside.

**Three findings from validating this, none of which were guessable:**

1. **Adjacent-pair checking was too weak a test here.** The first search optimised ring adjacency
   and returned worst ΔE 15.3 — but it chose *purple* for Quartz, which co-dominates the low-Mg/Si
   cells with the *dark blue* Fe-pyroxenes without ever touching them in the ring. Re-scored over
   **all pairs**, that assignment is bad. The rule: check all pairs whenever two phases are large
   in the same cells.
2. **The specified green and red do not pass on their own.** At the textbook values the
   green/red pair sits at ΔE **7.2**, inside the 6–8 floor band. Deepening the green to `#0a5c28`
   and lightening the red to `#ef5350` lifts the whole palette to **worst all-pairs ΔE 11.8,
   tritan 9.5, normal-vision 18.8** — better than any previous version of this figure.
3. **Hatching is drawn in the EDGE colour, so it vanishes on dark fills.** Near-black dots on the
   dark blue Fe-pyroxenes (relative luminance 0.032) or the dark green feldspars (0.079) are
   invisible. `hatch_ink` flips the texture to the surface colour below luminance 0.22, so
   Hedenbergite, Anorthite and Akermanite carry *white* dots and the rest carry dark ones.

The caption changed with the semantics: hatching no longer means "proxied kinetics", so that
caveat moved into its own sentence rather than riding on a visual channel that now means Ca.

> ⚠️ **The mineral and oxide palettes now overlap, and that is unavoidable.** The diagnostic
> figure (§32.8) shows oxide pies beside a mineral pie, and the two encodings are independent —
> mineral hue means family, and no oxide palette can mirror that. They collide exactly once:
> Al₂O₃ and Quartz are both `#eda100`. Identity is carried by direct labels on every significant
> slice and by two separate legends, but if the collision proves confusing the fix is to move
> Quartz off amber, not to re-anchor the two palettes to each other.

**`plot_crust_composition.py` now imports this palette** rather than defining a second one. It
had carried its own (Paul Tol 'bright', hue by structure) which was internally fine but meant a
mineral changed colour between the stack plot and the pie charts. It imports `MINERALS`, `COLORS`,
`HATCHED` and `hatch_ink` directly, so order, hue and texture cannot drift apart again, and it now
runs with **no CLI arguments** (`--csv` defaults to `CRUST_TABLE`).

Two defects surfaced while doing it, both of the same shape — a hard-coded value that was right
under the old palette and wrong under the new one:

- **Direct band labels were hard-coded white.** Legible on the old dark hues; nearly invisible on
  amber Quartz and light-blue Enstatite/Diopside. They now take `hatch_ink` of their own band.
- **A `T_p` reference line on a `T_melt` axis.** The panel drew `EARTH_TP = 1325` — a *potential*
  temperature — as "Earth" on an axis that has plotted the melting temperature since §32.4. The
  two agree to 3 °C at the anchor purely by coincidence. The line is now read off the plotted
  slice (`np.interp` at Earth's Mg/Si), the axis is labelled `T_melt`, and the local variable is
  no longer called `T_p`.

### 32.8 `plot_crust_diagnostic.py`: the pipeline as three pies

New figure, because the grid figure shows only the *output* of the pipeline and the question that
started this whole section was about the steps before it. Per composition, three pies side by side
— **bulk mantle oxides → primary melt oxides → normative minerals** — every slice labelled with
its own percentage, under a header carrying `T_melt` and F. One row per composition reads left to
right as one parcel of rock going through the whole calculation.

Two forms: individual figures for a selection of cells (`--points`, default Earth plus the
quartz-out transition and two corners), and the requested **3 × 3 grid** at Mg/Si = [0.5, 1.25,
2.0] × ΔIW = [−1, −2, −5], 27 pies in nine cells. A companion CSV carries every value, including
the slices too small to label.

Three things worth keeping:

- **`mantle_composition` now exists in Python** (`crust_composition.py`), mirroring the `.jl` — the
  bulk mantle was previously computable only inside the generator, so nothing downstream could
  plot the pipeline's input. Verified identical to the Julia construction (Mg/Si 0.5, ΔIW −1:
  SiO₂ 51.25, MgO 17.19, FeOt 24.14).
- **The oxide pies share the mineral hue ring**, and the shared hues mean the same thing on both
  sides: SiO₂/Quartz blue, MgO/Mg-pyroxene pink, FeOt/Fe-pyroxene green, Na₂O/Albite orange. So a
  wedge can be followed by colour from mantle to crust. The four unanchored hues were assigned by
  exhausting all 24 permutations against the CVD checks and keeping the best: worst adjacent
  pair ΔE 9.2, tritan 9.6, normal-vision 20.8.
- **Pie labels need a de-collision pass.** Two thin adjacent slices put their labels at nearly the
  same angle and matplotlib stacks them; `_spread_labels` pushes them apart per side, which is what
  makes 9-phase assemblages legible. Also: do NOT override matplotlib's per-angle horizontal
  alignment on pie labels — with `ha='center'` long names run back over the wedges.

### 32.9 A MORB reference, and a paper-sized grid

**MORB reference pie.** Every grid figure now carries average mid-ocean ridge basalt put through
the *same* CIPW norm as the computed crusts, so the grid is read against something measured rather
than only against itself. `MORB_OXIDES` lives in `crust_composition.py` beside the pyrolite
constants; nothing in the model consumes it.

The norm returns a textbook basalt — **53 wt% plagioclase (Anorthite 28.3, Albite 24.3), 25 wt%
clinopyroxene (Diopside 13.0, Hedenbergite 11.5)**, the rest olivine and orthopyroxene — which
sits visibly close to the Earth anchor cell, as it should.

> ✅ **VERIFIED 2026-09-09 against the primary source.** The paper was obtained and every one of
> the ten oxides in `MORB_OXIDES` matches Gale, Dalton, Langmuir, Su & Schilling (2013), G3 14,
> 489, **Table 1, "ALL MORB", arithmetic mean** — which the table's own caption names as the
> preferred composition ("the arithmetic mean (bold font) is our preferred ALL MORB composition").
>
> | | SiO₂ | TiO₂ | Al₂O₃ | FeOT | MnO | MgO | CaO | Na₂O | K₂O | P₂O₅ |
> |---|---|---|---|---|---|---|---|---|---|---|
> | Gale Table 1 | 50.47 | 1.68 | 14.70 | 10.43 | 0.184 | 7.58 | 11.39 | 2.79 | 0.160 | 0.184 |
> | `MORB_OXIDES` | 50.47 | 1.68 | 14.70 | 10.43 | 0.18 | 7.58 | 11.39 | 2.79 | 0.16 | 0.18 |
>
> Two conventions confirmed rather than assumed: the mean **excludes back-arc spreading centers**
> (that is Gale's "ALL MORB PLUS BAB", a different column), and iron is **FeOT, total iron as
> FeO**, which is what the `FeOt` key means and what the norm needs. The internal-consistency
> checks made when the numbers were unverified all hold at the verified values — CaO/Al₂O₃ =
> 0.7748, Mg# = 0.564, sum 99.568.
>
> Note the `pdftotext` extraction of Table 1 is **column-scrambled**: the oxide labels and the
> log-normal column are offset by a line, so a naive read pairs SiO₂ with 10.43. The assignment
> above is the physically unambiguous one and was cross-checked against the log-normal column
> (SiO₂ 50.41, MgO 7.69, FeOT 10.07, …), which is self-consistent only under this pairing.

**`--paper`.** `plot_crust_grid.py --paper` writes `crust_grid_paper.*` at MNRAS text width
(504 pt) with **5 × 5 = 25 pies** instead of 100: Mg/Si [0.5, 0.9, 1.25, 1.6, 2.0] × ΔIW
[−5, −3, −2, −1.3, −1]. The −1.3 row is kept deliberately — it is the transition §32.1 is about,
and dropping it would put the figure back to aliasing the thing the section exists to show. The
title block is dropped (a figure in a paper is captioned by LaTeX; repeating it above the axes
costs a third of the height the pies need) and all type is stepped down to stay legible at print
size.

**Two traps hit while building this:**

- **Importing `plot_results` for its page-width constant applies a global matplotlib style.**
  `plot_results.py:46` runs `plt.style.use()` at import, and that style sets
  `constrained_layout: True`, which silently disables the `subplots_adjust` this figure's whole
  layout depends on — matplotlib says so in a warning and carries on. It would have broken the
  existing full-size figure too, not just the new one. `TEXT_WIDTH_IN` is now duplicated locally
  with a comment saying why.
- **An axes added to a figure that uses subfigures is painted over by them.** The MORB pie was
  invisible on the diagnostic grid until it was added to the title *subfigure* rather than to the
  parent figure.

Also fixed here: a stale `wedge.set_edgecolor(INK)` left one line below the new
`hatch_ink` call in `plot_crust_grid`, so the pies were drawing dark hatching while the legend
drew light. The two disagreed for exactly as long as §32.7's hatch-contrast fix had been in.

### 32.10 Why the Earth cell is not MORB, and what that costs

Adding the MORB reference (§32.9) made an offset visible, so it was worth pinning down. It is not
a bug, and the causes are separable and quantifiable.

| oxide (wt%) | ours, F = 0.20 | MORB | diff |
|---|---|---|---|
| SiO₂ | 47.94 | 50.88 | −2.94 |
| TiO₂ | 0.93 | 1.69 | −0.77 |
| Al₂O₃ | 15.72 | 14.82 | +0.90 |
| FeOt | 7.71 | 10.51 | −2.80 |
| MgO | 12.48 | 7.64 | **+4.84** |
| CaO | 13.33 | 11.48 | +1.85 |
| Na₂O | 1.74 | 2.81 | **−1.07** |
| **Mg#** | **0.743** | **0.564** | |

**Cause 1 — MORB is not a primary melt.** Ours is a liquid in equilibrium with mantle residue;
erupted MORB has crystallised in crustal chambers. That is the whole MgO/Mg# gap, and our values
are where back-calculated *primary* MORB sits.

**Cause 2 — F = 0.20 is about twice MORB's melt fraction.** In our own melting Na₂O and TiO₂ are
near-perfectly incompatible (implied **D = 0.007 and 0.020**), so C_liq ≈ C_source/F. Inverting
that for MORB's concentrations against *our* source:

```
Na2O  ->  F = 0.128
TiO2  ->  F = 0.119
```

Two independent elements agreeing on **F ≈ 0.12**, inside the 0.08–0.12 range for real MORB that
§24.2 already records. At F = 0.20 every incompatible is diluted ~1.7×, which is precisely the
Na₂O and TiO₂ deficit. §6.1's PRIMELT misfit says the same thing from the other direction: 4.84
wt% at F = 0.20 against 2.22 at F = 0.117.

**A negative result worth keeping.** Fractionation alone does *not* reconcile the two. Removing
12 wt% of equilibrium olivine (Fo86; Fe–Mg Kd = 0.30, Roeder & Emslie 1970 — a one-off diagnostic,
not model code) puts MgO exactly on MORB's value but drives Al₂O₃ to 17.8 (MORB 14.8) and CaO to
15.1 (MORB 11.5), and leaves Na₂O, TiO₂ and FeOt still short. Real MORB also fractionates
plagioclase and cpx, which is what holds Al₂O₃ and CaO down. So cause 1 is real but partial, and
the incompatible deficit belongs to cause 2.

**What neither explains:** FeOt (7.71 vs 10.51) and CaO/Al₂O₃ (0.848 vs 0.775). The latter is the
ultracalcic bias at F = 0.20 past cpx-out already in §7 of the methods doc. The FeOt gap is partly
that ferric iron is off (MAGEMin `O` = 0, while MORB's FeOt includes ~10% Fe³⁺) and partly the
single 1 GPa segregation pressure, where real MORB pools a column extending deeper and
higher-pressure melts are more Fe-rich.

**What it costs the weathering model — this is the part that matters.**

| | ours | MORB |
|---|---|---|
| olivine | **24.1** | 11.3 |
| orthopyroxene | 0.2 | 11.5 |
| Anorthite / Albite | 35.4 / 14.9 | 28.3 / 24.3 |

Total plagioclase is nearly the same (50.4 vs 52.6) but its Ca/Na split is not, and we carry
**twice the olivine**. Against this repository's own `K_FUNCTIONS` at 25 °C, pH 6.5, the two
olivines are 1.3–3.3 orders of magnitude faster than every other phase in the Earth cell:

```
log10 k_eff (mol/m2/s):  Fayalite -8.64  Forsterite -9.57
                         Diopside -10.90  Anorthite -11.02  Albite -11.04  Hedenbergite -11.92
```

Use `K_FUNCTIONS`, not `RATE_FUNCTIONS`, for any such comparison: the latter omits Nepheline,
Hedenbergite, Ferrosilite and Akermanite (which are proxied or carry effective constants) and its
members do not share a call signature, so ranking across it silently mixes conventions. Ranked
over the full assemblage the fast group is wider than olivine alone — **Akermanite ties Forsterite
at −9.57** (it is proxied on forsterite) and **Nepheline is −9.82**.

So the modelled crust weathers faster than real oceanic crust and **over-delivers Mg and Ca while
under-delivering Na** — Na being exactly the ion that needed the K_NA sink to behave (§6). The
Earth calibration was performed with this crust, so much of the bias is absorbed into the tuned
constants; that is the caveat, not the defence. Those constants are not independently meaningful,
and changing F would require recalibrating.

**Verdict: keep F = 0.20.** It is independently corroborated inside this pipeline (the cpx-out
closure returns F = 0.213 at the Earth anchor, §24.5) and it is the right closure for a *generic*
planet, where Earth's particular melting regime cannot be assumed. Earth simply melts less than
the cpx-out limit. The honest framing for the paper is: primary melts at a cpx-out-anchored melt
fraction, validated against PRIMELT primary melts, with **MORB shown for orientation and never as
a validation target**. Reproducing erupted crust needs a different closure — lower F plus an
explicit fractionation step — which is a recalibration, not a parameter tweak.

Trends across the grid are much safer than absolute values: the offset is systematic, so relative
behaviour along Mg/Si and ΔIW carries over even where the absolute fluxes do not.

**Tested directly: what F = 0.12 actually gives.** `--ftarget` was restored to
`make_crust_compositions.jl` for this (it had been cut as unused tooling in §32.4; the F-sensitivity
question is exactly what it is for):

```
julia src/kamino/data/make_crust_compositions.jl --points "1.25,-2.0" --ftarget 0.12
```

At the Earth composition, F = 0.120, **T_melt 1296 °C** (against 1328 at F = 0.20):

| phase (wt%) | F = 0.20 | F = 0.12 | MORB |
|---|---|---|---|
| Albite | 14.9 | 16.9 | 24.3 |
| Anorthite | 35.4 | 36.9 | 28.3 |
| **Nepheline** | 0.0 | **3.0** | 0.0 |
| Diopside | 18.2 | 14.5 | 13.0 |
| Hedenbergite | 7.2 | 6.0 | 11.5 |
| Forsterite | 16.0 | 15.0 | 5.4 |
| Fayalite | 8.1 | 7.8 | 6.0 |
| Enstatite / Ferrosilite | 0.1 / 0.1 | 0.0 / 0.0 | 5.7 / 5.8 |

Grouped: plagioclase 50.4 → **53.8** (MORB 52.6, now matching), clinopyroxene 25.4 → **20.4**
(MORB 24.6, worse), olivine 24.1 → **22.8** (MORB 11.3, essentially unmoved).

> ⚠️ **Lowering F does not make the crust more MORB-like — it makes it silica-UNDERSATURATED.**
> Nepheline appears at 3.0 wt%. A low-degree melt is alkali-rich relative to silica, i.e. an alkali
> basalt; MORB is a tholeiite and is never nepheline-normative. Na₂O improves (1.74 → 2.61 against
> MORB's 2.81) but Al₂O₃ degrades in step (15.7 → 17.6 against 14.8), so the RMS oxide misfit to
> MORB barely moves: **2.40 → 2.29 wt%**.

**Why that would be expensive here specifically.** Nepheline dissolves **1.22 decades faster than
albite** (−9.82 vs −11.04), and albite-versus-nepheline is what sets ocean Na in this model (§6).
Adopting F = 0.12 would therefore raise the Na flux twice over — more Na₂O in the melt, and that Na
sitting in a far more reactive phase. It is a recalibration, not a refinement.

**A third cause of the offset, revealed by the same run.** K₂O goes 0.14 → 0.23 and now *overshoots*
MORB's 0.16. Source K₂O / F = 0.029 / 0.12 = 0.24, i.e. perfectly incompatible — and MORB's K₂O
inverts to F ≈ 0.18, contradicting the F ≈ 0.12 that Na₂O and TiO₂ give. **No single F reconciles
all three**, which is itself the evidence: we melt BSE pyrolite, whereas MORB comes from mantle
already stripped of its most incompatible elements, and K is more incompatible than Na. Part of the
offset is source depletion, not melt fraction, and no choice of F can absorb it.

Neither of the two §7 limitations is relieved either: at F = 0.12 and ΔIW −2, Mg/Si 0.5 is still
strongly quartz-normative (SiO₂ 72.4) and Mg/Si 2.0 is still ultracalcic (CaO/Al₂O₃ 1.78).

All of which **reinforces keeping F = 0.20**: the closure is defensible on its own terms, and the
one thing lowering F buys (Na₂O) costs silica saturation, Al₂O₃, K₂O and a recalibration.

### 32.11 Where this leaves the crust pipeline

- ✅ 650-cell isobaric table in place; `check_crust_table.py` passes every check including PHREEQC
  at all three Mg/Si extremes; Earth anchor unchanged in composition (SiO₂ 47.88, CaO/Al₂O₃ 0.85).
- ✅ `docs/crust_composition.md` rewritten for the isobaric method (§3, §5, §6, §7).
- ⚠️ **Every stored sweep result predates this table.** The composition axis values themselves are
  unchanged at the shared cells (§32.3), so results at those points stand; anything that
  *interpolated* between ΔIW −1.5 and −1.0 at Mg/Si ≤ 0.7 was reading a blend of two rock types and
  should be regenerated.
- ⚠️ The cpx-out closure (`Tp_for_cpx_out`, F = 0.213 at the Earth anchor, §24.5) and the Brugman
  experimental comparison (`--validate --bulk`, §6.3 of the methods doc) are recorded results whose
  code is now deleted. Restore from git if either needs re-running.
- ⚠️ `T_melt` is not a potential temperature. If the paper wants T_p, run the isentrope once for
  that single number rather than per grid point.
- ⚠️ **The Earth cell is a primary melt and does not reproduce erupted MORB** (§32.10): Mg# 0.74
  against 0.56, Na₂O 1.74 against 2.81, twice the normative olivine. Two causes, both quantified —
  MORB is differentiated, and F = 0.20 is ~1.7× the melt fraction MORB's own incompatibles imply
  (F ≈ 0.12). Keep F = 0.20, but never present MORB as a validation target, and expect the crust
  to over-deliver Mg and Ca and under-deliver Na.
- ✅ **The MORB constant is verified against Gale et al.** (§32.9, checked 2026-09-09): all ten
  oxides match Table 1's ALL MORB arithmetic mean exactly. Safe to publish.
- ✅ F = 0.12 tested directly and rejected (§32.10): it fixes Na₂O but turns the Earth melt
  nepheline-normative, and nepheline is 1.2 decades faster-dissolving than albite. `--ftarget` is
  back in the generator if the question needs revisiting.

---

## 33. Continental weathering, the seafloor-area fix, and what actually controls the crossover (2026-09-07 → 09-08)

The question this session set out to answer was narrow — give the model a habitable zone for
Earth-like planets, to compare against the land-free ocean worlds — and it turned into a
sequence of corrections, one of which invalidates the Earth calibration. **§33.3 is the one to
read if you read only one.**

### 33.1 `continental_baseline.py` rewritten

The old script was dead code: it named `basalt_49` (replaced by the Mg/Si–ΔIW axes in §25), wrote
to a `/data/pt426/kamino_experiments_fast_3` path that no longer exists, and set `f_HT = 0.01`,
which `plot_results.plot_continental_baseline` filters out — so it could not have produced the
figure it existed to feed.

It now runs one instellation line (0.30 → 1.45, step 0.05) at `land_fraction = 0.3` with every
other axis at Earth: 1× outgassing, 1× crust production, 3 km, Earth crust, reverse weathering
on, reducing. The land-free arm runs alongside as the comparison; its run names are identical to
the ones `sweep_basic` already wrote (the `_land` tag is suppressed at 0), so those come off disk
for free. `parameter_sweep._run_name` and `run_simulation` gained a `land` argument, both
defaulted, so no existing sweep or run name moves.

`SWEEP` at the top of the file selects `'baseline'`, `'land'`, `'grid'` or `'alpha'`.

**Result at Earth:** T = 295.5 K, pCO₂ 873 ppm (§22's calibration reports 294.4 K). Habitable
band **0.48 ≤ S ≤ 1.125**. The land-free arm at 1× outgassing has **no habitable band at all** —
19 of 24 runs peg the 10 bar CO₂ ceiling — which is why the ocean-world sweeps run at 0.1×
outgassing. Seafloor weathering alone cannot balance Earth's outgassing.

### 33.2 A real bug: post-runaway states counted as habitable

`get_T_surface_analytic` scans T upward and returns the first sign change. The OLR is **not
monotonic in T** — at 1 Pa CO₂ it peaks at 271.5 W/m² near 320 K (the Simpson–Nakajima limit),
dips, then climbs. When absorbed instellation exceeds that peak there is no cool-branch root, and
the solver silently returns a value on the **hot branch beyond the runaway**.

| S | F_in (W/m²) | T returned | branch |
|---|---|---|---|
| 1.10 | 261.8 | 298.3 | cool |
| 1.14 | 271.1 | 317.6 | cool |
| 1.15 | 273.7 | **357.4** | **hot — past runaway** |

357.4 K passed a `T < T_RUNAWAY = 360` test, putting the inner edge one grid point too far out.
`hz_edges` now rejects states above the OLR limit and above the fit's 350 K validity ceiling.
Inner edge 1.150 → **1.125**; the climate model's own runaway threshold at the 1 Pa floor is
S = 1.141, consistent with that bracket.

### 33.3 ⚠️ The seafloor-area fix — THE CALIBRATION IS NOW STALE

`J_total` comes from `EARTH_HYDROTHERMAL_FLUX_PER_AREA`, normalised on `A_SEAFLOOR_EARTH`, and
`flux_LT` inherits that basis — but `F_diss`, the HT Mg→Ca exchange, the Na sink and Cl
subduction all applied it over the **full sphere**. `ocean_water_mass` had the same error, and
since every one of those terms is `(per-area rate × area) / ocean_water_mass` the two cancelled,
so the seafloor terms were right as concentration rates. What was actually wrong is that
**`F_cont` and `F_vol` were diluted into 1/(1 − f) too much water** — 1.43× too weak at Earth's
land fraction, overweighting the seafloor sink relative to continental weathering by that factor.

`ocean_depth` is the true depth over the seafloor (that is how `dY_dt` uses it for the pore
pressure), so the ocean covers (1 − f) of the surface. `Planet.seafloor_area` added and used for
`ocean_water_mass` and the four hydrothermal terms, plus a guard on `land_fraction ≥ 1`.

**`land_fraction = 0` is bit-identical under this** ((1 − 0) = 1), verified — so all 7586
land-free runs, i.e. every main sweep, remain valid.

**But `calibrate_earth.py` fits at `LAND_FRAC = 0.3`**, so `KD_MG_CALIB` and `K_NA_CALIB`
absorbed the 1.43×. Measured at S = 1, land 0.3:

| | old | new | vs Earth |
|---|---|---|---|
| T | 295.5 K | 295.2 K | — |
| pCO₂ | 873 ppm | 824 ppm | — |
| Na | 491 | 679 mM | 1.02× → 1.41× |
| Mg | 46 | 77 mM | 0.88× → 1.45× |
| Cl | 165 | 235 mM | 0.30× → **0.43×** (toward Earth) |

Climate barely moves; the ocean does. The mechanism is clean: the Na and Cl *sinks* are
`J·A/M` and unchanged, while their sources are now 1.43× stronger. **Re-run
`calibrate_earth.py` before quoting ocean chemistry.** The 33 land-bearing runs on disk at the
time were deleted rather than left for the resume path to reuse.

### 33.4 Why the habitable zone does not look like Kopparapu 2013

It is not supposed to. The climate model's *radiative* limits do match:

| | Kopparapu 2013 | this model |
|---|---|---|
| Earth (S = 1, 280 ppm) | 288 K | 290.3 K |
| runaway greenhouse | 1.06 (1.107 Earth-mass) | 1.141 |
| maximum greenhouse | 0.343–0.35 | ~0.41 |

The OLR fit is Haqq-Misra et al. (2016), ApJ 827, 120 — a polynomial fit *to* Kopparapu's 1-D
radiative-convective columns, valid for 10⁻⁵–10 bar CO₂ and **150–350 K**. That is an
OLR(T, pCO₂) parameterisation, not the S_eff boundary polynomials.

> ⚠️ **Corrected 2026-09-24 (§37.10).** The fit in `climate/analytic.py` is **Kadoya & Tajika (2019),
> ApJ 875, 7**, not Haqq-Misra et al.: its form (ξ = 0.01(T − 250), I₀ = −3.1 W m⁻², separate
> polynomials above and below 1 bar) is theirs, fitted to Kopparapu et al.'s model and valid for 150–350 K
> and 10⁻⁵–10 bar. Haqq-Misra et al. (2016) write their OLR fit in log₁₀T; Kadoya & Tajika use Haqq-Misra's
> *albedo* fit, which this code does not (it has its own Rayleigh-scattering albedo). The validity range
> quoted above is unchanged. The paper cites Kadoya & Tajika correctly.

What differs is the **outer** edge, and it is physics not error. Kopparapu's maximum greenhouse
asks what the best possible CO₂ can do; the carbon cycle asks what CO₂ you actually get. At
S = 0.5 the maximum greenhouse would use 4.3 bar (302.9 K) and holding 273 K needs 0.79 bar — the
sweep supplied **0.474 bar**, giving 261 K. The reason is the WHAK parameters: with β = 0.3 the
CO₂ term cancels the temperature term, so the sink never shuts off:

| S | T (K) | pCO₂ (bar) | CO₂ term | T term | **f** |
|---|---|---|---|---|---|
| 0.50 | 261.0 | 0.474 | 9.30 | 0.204 | **1.90** |
| 1.00 | 295.5 | 8.7e-4 | 1.41 | 1.555 | **2.19** |

Weathering sits near 2× modern Earth across the *entire* zone. Sustaining 4.3 bar at S = 0.5
would need it at **43×**, which 1× outgassing cannot feed. So the outer edge is CO₂-supply
limited at 0.48, not radiation limited at ~0.41. There is also no runoff term, no supply limit,
and no ice-albedo feedback (`LAND_ALBEDO == OCEAN_ALBEDO == 0.3`), so weathering never shuts down
through glaciation — which is how CO₂ reaches multi-bar levels in max-greenhouse calculations.

⚠️ At S ≥ 1.1 the carbon cycle drives pCO₂ to 10⁻⁴–10⁻⁹ bar, **below the OLR fit's 10⁻⁵ bar
floor**. The code clamps to 1 Pa so radiation stays in range, but the climate is then evaluated
at a CO₂ up to three orders of magnitude above what the chemistry says. The inner edge is set by
CO₂ exhaustion, so this sits exactly where it matters.

### 33.5 The land-fraction series and the crossover

`LAND_FRACTIONS = [0.3, 0.2, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0]`, log-spaced because at
land 0.3 continental alkalinity is ~21 Tmol eq/yr against a seafloor ~0.018 — a ratio near 1200.

At the Earth reference the crossover is **f\* = 3.2e-4 to 1.7e-3**: seafloor weathering only
dominates below a land fraction of ~10⁻³. Over the whole plausible terrestrial range continental
weathering wins by 10²–10³.

The coarse grid (`SWEEP = 'grid'`, 810 combos: instellation × land × outgassing × crust × Mg/Si)
shows f\* moving **~2 orders of magnitude**, 3e-4 to 0.075. It rises with crust production
(which drives `J_total` directly), falls with outgassing, and at 10× outgassing with Mg/Si 1.25
there is no crossover anywhere sampled. Mg/Si 1.8 shifts it up 3–10×.

Cross-design check: the Mg/Si 1.25, crust 1×, out 1× cell gives 0.00031–0.0011 against the fine
series' 0.00032–0.0017 — two independent designs agreeing.

### 33.6 The Mg/Si trend is climate-mediated, not a reactivity effect

`mantle_mg_si` reaches the model only through `Planet.crust_composition`, and
`crust_production_rate` only through `J_total`. Neither is an argument to
`get_continental_weathering_flux`, which sees T and pCO₂ against an `F_alk_ref` pinned to modern
Earth and a cation split fixed to modern river chemistry. **Continental weathering cannot respond
to crust chemistry or tectonic rate except through the shared climate.**

Measured at land 0.003, S = 0.8, Mg/Si 1.25 → 1.8: the seafloor flux rises only **1.3–1.8×**
while the continental flux **falls to 0.48–0.68×**, because the stronger seafloor sink draws
pCO₂ 1.45 → 0.76 bar and cools the planet 330.5 → 321.6 K, weakening WHAK. In 4 of 5 cells the
climate-mediated continental change is the larger factor.

### 33.7 The crust-reactivity diagnostic: high Mg/Si cannot make felsic continents

To test whether the Mg/Si signal is real reactivity, a **second melting stage** was added:
`src/kamino/data/make_continental_compositions.jl` melts the stage-1 oceanic crust hydrously at
15 kbar / 3 wt% H₂O to F = 0.20, keeping the melt and discarding the residue — the TTG model of
Archean continental crust (Rapp & Watson 1995, J. Petrol. 36, 891; Moyen & Martin 2012, Lithos
148, 312). Not part of the planet model; a diagnostic only.

The pressure was chosen from a probe, not assumed: at 10 kbar the residue is spinel + feldspar
with no garnet and the melt stays basaltic (SiO₂ 48); at 15 kbar the residue is garnet + cpx +
amphibole and SiO₂ goes 47.9 → 56.2, Na₂O 1.74 → 5.87. 650/650 grid points, 0 failed.

`experiments/crust_reactivity.py` norms both crusts through the **same** norm and evaluates
`get_k` at one fixed (T, pH). Two choices matter:

- **Both sides use `_cipw_norm_native`.** The stage-2 melts are peraluminous (361/650) and the
  pyrolite `cipw_norm` refuses a corundum-normative rock. Mixing norms would have put a **1.26×
  artefact** into the ratio (k_oceanic 2.92e-10 pyrolite vs 3.68e-10 native at the Earth point).
- **One fixed (T, pH) for both**, or a chemistry difference is folded into a reactivity ratio.

**Earth point: 2.49** — the oceanic crust releases alkalinity 2.5× faster per unit mass than the
trondhjemite derived from it. Across the grid the ratio spans only 0.25×–2.5× Earth's.

At **Mg/Si 1.8, ΔIW −2 the ratio is 0.40× Earth's**, and 100% of the Mg/Si 1.8 column is below
Earth. Decomposed:

| | Mg/Si 1.25 | Mg/Si 1.8 | change |
|---|---|---|---|
| k_oceanic | 3.68e-10 | 4.95e-10 | ×1.34 |
| k_continental | 1.48e-10 | 4.95e-10 | **×3.34** |
| ratio | 2.49 | 1.00 | ×0.40 |

The oceanic crust *does* get more reactive — ×1.34, matching the ×1.3–1.8 seafloor flux change in
§33.6. But its continental derivative gets **more** reactive still. The mechanism is visible in
the compositions: at Earth's Mg/Si, melting gives a genuinely felsic rock (SiO₂ 47.9 → 56.7, MgO
12.5 → 4.1); at Mg/Si 1.8 the parent is silica-poor and the melt barely differentiates (44.7 →
46.0, MgO 23.1 → 12.7). **A high-Mg/Si planet cannot make felsic continents from this route** —
its "continental" crust stays basaltic and weathers nearly as fast as its seafloor.

Composing the model's f\* trend with the missing continental term (÷3.34 — the model already
contains the seafloor half) flips the sign in 3 of 5 (crust, outgassing) cells:

| crust | out | model | corrected |
|---|---|---|---|
| 0.1 | 0.1 | 4.23× | 1.27× |
| 1 | 1 | 6.88× | 2.06× |
| 1 | 0.1 | 1.85× | **0.55×** |
| 10 | 0.1 | 2.02× | **0.60×** |
| 10 | 1 | 2.74× | **0.82×** |

⚠️ That composition assumes continental weathering is **kinetically** limited. On Earth it is
substantially **supply** limited (West et al. 2005), where a more reactive rock buys little — in
which case the correction collapses toward 1 and the model's original trend stands. The seafloor
law carries transport/supply limitation (sedimentation, Damköhler); the continental law has
neither, so that coupling has no route into the model. Both readings should be quoted as bounds.

### 33.8 Habitable planets are kinetically limited — so `get_k` is the right comparator

Measured across the steady states on disk at 3 km:

| population | median Da | fraction at Da ≥ 1 |
|---|---|---|
| habitable (converged/timeout) | 0.0067 | 10.5% |
| non-habitable | 0.12 | 35.3% |

and Da falls with land fraction: 8.8e-3 at land 0, 3.3e-3 at 0.03, **2.2e-4 at land 0.3**.
Continental weathering keeps the planet cool and low-CO₂, so seafloor weathering never approaches
saturation. Since crossing the thermodynamic limit is what flips the feedback (§21), the habitable
population is essentially all kinetic, and a kinetic rate-constant comparison is the operative
one rather than a partial view.

⚠️ `parameter_sweep.py` states "Earth is transport-limited (Da ≫ 1)" in its `alpha` argument. The
land-free half of that comment matches measurement (8.8e-3), but the land-bearing runs come out
at Da = 2.2e-4 — as kinetic as anything on the grid. The parenthetical appears to conflate
real-Earth *continental* supply limitation with the model's *seafloor* Damköhler. It is
load-bearing for the alpha argument and should be corrected.

### 33.9 `alpha` is the largest control, and it is sub-linear

`SWEEP = 'alpha'`: 405 combos, alpha [1.1, 10, 50] × outgassing [0.1, 1, 10] × land × 9
instellations, at Earth crust production and Mg/Si. 0 failed.

| outgassing | α = 1.1 | α = 10 | α = 50 | exponent |
|---|---|---|---|---|
| 0.1× | 0.0038 | 0.0585 | 0.0725 | **0.80** |
| 1× | 0.00059 | 0.0042 | 0.0154 | **0.86** |
| 10× | no crossover at any α | | | — |

**d log f\* / d log α = 0.80–0.86**, against the α¹ the kinetic limit predicts. Both arms are
sub-linear, so the climate feedback damps alpha's leverage by ~15–20%. Over alpha's unconstrained
1.1–50 range that is still a **19–26× lever on f\***, larger than any planetary property
measured: outgassing 6–21×, crust production 1.4–4.7×, crust composition ~1× once continental
crust is allowed to track the mantle. **And alpha is a model parameter §28.2 records as
unidentifiable from Earth.**

Outgassing sets *whether* a crossover exists (none at 10× at any alpha); alpha sets *where*.

⚠️ These numbers are the **corrected** ones. Measured before the wall-timeout recovery (§33.11)
they were 0.78–1.16, i.e. one arm apparently super-linear, which supported the opposite
conclusion — that the feedback does not damp alpha at all. Several crossover points were being
read off a grid with holes in it. Do not cite the earlier figures.

⚠️ The α = 50, outgassing 0.1× cell reaches f\* ≈ 0.073 against a highest sampled land fraction
of 0.3, so it is saturating against the grid edge; its 0.80 exponent is partly that, not physics.

### 33.10 What the gaps in the ratio maps are

86 missing cells in the alpha grid, categorised rather than assumed:

| cause | n | what it is |
|---|---|---|
| out_of_domain: CO₂ ceiling | 61 | pCO₂ 2.8–10 bar; the sink cannot balance outgassing |
| out_of_domain: runaway | 9 | T → 389 K, the inner edge |
| net alkalinity sink | 7 | steady state, but seafloor flux **negative** |
| wall_timeout | 5 | still slow at 3600 s |
| fallback_limit | 4 | chemistry, see §33.12 |

Only 9 of 86 are computational; 77 are the model reporting a physical outcome. The CO₂-ceiling
cells fill the corner where large carbon supply meets small continental sink — 65 of 86 are at
10× outgassing, 35 at land 3e-4 — which is the same statement as §33.1's land-free arm.

The **net alkalinity sink** cells were previously dropped *silently*, leaving unmarked blanks
indistinguishable from "no run". They are steady states where pore precipitation exceeds basalt
dissolution, occurring where continental weathering runs at ~154 Tmol/yr and floods the ocean
with cations. A negative ratio has no place on a log scale, but it is a real outcome:
`_ratio_cells` now returns them separately and the maps mark them with an open circle against the
`×` used for never-reached-steady-state.

### 33.11 Finishing the wall-timeout runs, and two spawn hazards

`ProcessPoolExecutor` spawns workers that **re-import `parameter_sweep`**, so a value assigned in
the parent never reaches them: monkeypatching `RERUN = True` would have produced a sweep that
reported success having recomputed nothing. `RERUN`, `WALL_SECONDS_SHALLOW` and
`WALL_SECONDS_DEEP` now read from the environment (`KAMINO_RERUN`, `KAMINO_WALL_SHALLOW`,
`KAMINO_WALL_DEEP`), which does cross the spawn boundary. Defaults unchanged.

The same hazard bit the one-off runner: without an `if __name__ == '__main__'` guard every worker
re-executed the launch and the pool died with `BrokenProcessPool`. `experiments/rerun_wall_timeouts.py`
has the guard, and asserts every rebuilt run name **round-trips to an existing file** first —
`_run_name` formats `1` and `1.0` differently, so a rebuilt combo can easily write a *new* file
and leave the original wall_timeout untouched.

100 runs re-run at 3600 s (4× the cap): **67 reached 2 Gyr**, 24 still wall_timeout, 8
fallback_limit, 1 chemistry_void. Alpha-grid gaps 99 → 86, wall_timeout 22 → 5.

The 8 new fallback_limits are informative: given four times the wall clock they burned through
5000 chemistry fallbacks instead. Those runs were never merely slow — the short cap was masking
the real failure. Raising the budget again would convert wall_timeouts into fallback_limits, not
into results.

### 33.12 What `fallback_limit` actually is

When PHREEQC fails to converge, `dY_dt` catches the `ChemistryError`, **reuses the last
derivative that did converge**, and counts a fallback; past 5000 the run is abandoned. So it is
not "ran out of time" but "spent 5000 steps integrating physics it never computed".

The comment at that `except` says *"typically high P_CO2 where PHREEQC cannot converge"*. **That
is not what these are.** Three of them — S = 0.90, out 1×, land 0.3, at all three alphas — sit at
283 K and 3400 ppm CO₂. Failing identically across a 45× range of alpha says the problem is the
chemical state, not flux magnitude.

The converged neighbours show it: Alk **598 mM**, DIC 363–522 mM, **Ca 0.36 mM** (a thirtieth of
seawater), Na 650–679, Cl 235. The charge balance closes exactly — Na 678.6 + 2×Mg 76.8 + 2×Ca
0.36 − Cl 235.1 = 597.8 against the 597.7 recorded. This is the **Cl deficit** (§7): the missing
anion charge is carried by carbonate alkalinity, which drags DIC up and precipitates calcite
until Ca is nearly exhausted. Near-zero Ca against enormous carbonate alkalinity is a nasty
corner for a speciation solver — calcite's SI becomes hypersensitive to tiny Ca changes.

⚠️ A `fallback_limit` run stores an **empty `data.y`** (0 rows against the expected 13), so
salinity and the ion state are unrecoverable; only the abort-state T, pCO₂ and diagnostics
survive. Correctly excluded from the figures either way.

**These are a symptom of the Cl budget, not an independent numerical problem.** Fixing Cl would
likely remove them, and matters far more for the ocean chemistry results than for these 4 cells.

### 33.13 Coupled axes: what should not be swept independently

Two pairs in the grid are not physically independent, and both were treated as such at some point
this session:

- **`crust_production` × `land_fraction`.** A stagnant-lid planet has no subduction, but felsic
  crust can still form by intracrustal melting of a thickened hydrated basaltic pile — drip
  tectonics/sagduction (Sizova et al. 2015; Rudnick 1995 and Jagoutz & Kelemen 2015 for the
  modern arc+delamination alternative). It needs the pile to reach garnet-amphibolite facies.
  Using the model's own 50 Myr resurfacing at 1× and ~7 km per turnover, over 4.5 Gyr with no
  recycling: 1× → 630 km, 0.1× → 63 km, 0.03× → 19 km, 0.01× → 6 km. There is a floor near
  **0.05–0.07×** below which the crust never gets deep enough to melt. At 0.1× only ~23 km ever
  passes below 40 km, giving ~4–5 km of TTG as a global layer against the **12 km** Earth's
  continents represent. So the (crust 0.1×, land 0.3) corner is unpopulated. Venus tesserae (~8%
  of the surface, possibly felsic) and Mars (felsic only in ancient highlands) are the
  observational anchors: some felsic crust, nothing like continents.
- **`ocean_depth` × `land_fraction`.** These are two readings of water inventory and hypsometry,
  so they trace a curve, not a plane. At land 0.3, 3.7 km gives 0.98× Earth's ocean mass (which is
  why `calibrate_earth.py` uses 3700 m); 10 km at land 0.3 demands 2.6× Earth's water while
  keeping 30% dry. Depth is a free axis for `land = 0` only. Note the continental baseline runs at
  3000 m, i.e. **0.79× Earth's ocean mass** — chosen for comparability with the ocean-world
  sweeps, but 3700 m is the Earth-consistent value.

### 33.14 Figures and analysis added

All through `plot_results`, so they share its style, geometry and termination markers.

- `continental_vs_ocean_{tp,chem}` — the two arms at Earth outgassing, Damköhler-styled.
- `continental_habitable_zone` — T vs S plus the zone as bars; edges labelled `crossing`,
  `bracketed` or `open` so a bound is never read as a measurement.
- `weathering_ratio_map` — seafloor/continental flux over instellation × land fraction.
- `weathering_ratio_grid_mgsi{1.25,1.8}` — the same faceted over crust × outgassing.
- `weathering_ratio_alpha_grid` — faceted over alpha × outgassing.
- `alpha_scaling` — f\* against alpha with the α¹ reference.
- `crust_reactivity_ratio{,_absolute}` — the (Mg/Si, ΔIW) reactivity map.
- `plot_results` gained `CONTINENTAL_HZ_OUTER = 0.480`, `CONTINENTAL_HZ_INNER = 1.125` and
  `SHOW_HZ_EDGES`; the lines are drawn on **temperature panels only**, and the legend entry
  appears only where they do. `continental_baseline._report` warns when its measured edges
  disagree with the constants, so re-running cannot silently leave every other figure stale.

All ratio maps are diverging about 1 with a neutral midpoint (`cmr.fusion_r`, chroma 0.000 at its
centre, checked not assumed), **not** forced symmetric — the ratio spans ~1e-3.5 to ~1e0.5, so a
symmetric range would reserve half the ramp for values that never occur. Faceted maps share one
colour scale across all panels, or each panel would renormalise and look alike whatever its
numbers were.

### 33.15 Status

- ✅ Continental baseline, land-fraction series, coarse grid and alpha sweep all run, 0 failures.
- ✅ Seafloor-area fix; land-free runs verified bit-identical.
- ✅ Crust-reactivity diagnostic; 650/650 stage-2 melts.
- ⚠️ **`calibrate_earth.py` must be re-run** (§33.3). `KD_MG_CALIB` / `K_NA_CALIB` absorbed the
  1.43×. Every land-bearing result in this section uses post-fix physics with pre-fix constants:
  Na 1.41× Earth, Mg 1.45×. The *ratios* are robust (10²–10³ effects against a 1.43× shift); the
  absolute ocean chemistry is not.
  > **2026-09-09 (§35):** `K_CL_SUBDUCTION` absorbed the same 1.43× and was **not** on this list —
  > its analytic derivation assumes the source and sink areas cancel, which §33.3 broke. Fix it
  > *before* the refit; it is an input to the fit, not an output. Two further defects found in the
  > same place: the ocean is never seeded (τ_Cl = 5.6 Gyr, so Cl reaches only 30% of steady state
  > in 2 Gyr) and SO₄ is pinned at zero in every sweep while the calibration seeds 23.45 mM.
- ⚠️ `CONTINENTAL_HZ_OUTER/INNER` were measured before the area fix. T moved −0.1%, so they should
  be unchanged, but the drift check will say so on the next run.
- ✅ The Cl deficit (§7) is now implicated in the `fallback_limit` runs as well as the ocean
  chemistry (§33.12) — **root-caused 2026-09-09 (§35.1)**: it is an unseeded initial condition
  against a 5.6 Gyr residence time, not a mis-set parameter. Seeding recovers Ca (0.36 → 12.4 mM),
  which should remove the `fallback_limit` corner.
- ⚠️ Both melting stages use MAGEMin's `"ig"` database. A metabasite set would be more defensible
  for stage 2 if it reaches the paper.
- ⚠️ `parameter_sweep.py`'s "Earth is transport-limited (Da ≫ 1)" comment contradicts measurement
  (§33.8).

---

## 34. The `alpha` decision and the second refit (2026-09-01)

*Numbered after §33 but chronologically before it: this session was only discovered on 09-09.*

> ⚠️ **This section is RECONSTRUCTED FROM THE CODE on 2026-09-09, not written from a session
> record.** No section covered 2026-09-01 — `development_history.md` jumped from §31 (08-27) to
> §32 (09-03) — and the work was found only because `weathering.py`'s `ALPHA_REF` and
> `parameter_sweep.py`'s comment block both carry that date. What follows is what the code
> asserts. **The reasoning behind the numeric values is not recorded anywhere and could not be
> recovered**; if the fit's own output still exists, attach it here.

### 34.1 What changed

`weathering.py:16-19` records a joint least-squares refit "after the `crust_composition.py`
rewrite invalidated the prior fit":

| | §28.1 (08-27) | 09-01 refit | shipped today |
|---|---|---|---|
| `ALPHA_REF` | 0.487612 | **1.100155** | 1.100155 |
| `KD_MG_HT` | 1.394362e-02 | **1.394755e-02** | 1.394755e-02 |
| `K_NA_CONT_REMOVAL` | 4.272026e-03 | **4.234317e-03** | 4.234317e-03 |

`alpha` moved 2.26×; the other two moved by <1%, which is consistent with §22/§28's repeated
finding that `alpha` is nearly free on Earth while `K_na` and `KD_mg` are tightly pinned by the
Na/Ca/Mg targets.

### 34.2 The operational decision, which is the durable part

`parameter_sweep.py:45-47` — **production runs at `ALPHA_REF` itself rather than a separately
pinned round number**, so `ALPHA_CALIB = ALPHA_REF` by construction and the drift check
"becomes a tautology for alpha specifically, which is the point". The sensitivity arm became
`[ALPHA_REF, 10, 50]` (was `[2, 10, 50]`).

This closes the item §15 and §22.9 both called top-priority — but as a **convention, not a
measurement**. `parameter_sweep.py:36-42` is explicit that nothing about identifiability changed:
ocean concentrations still move <6% across a 41× change in `alpha`, because Earth is
transport-limited while the land-free worlds the sweeps target are kinetically limited (Da ~ 0.005
over 19 pilot states, 0/19 with Da > 1), where `F ∝ alpha` linearly.

A domain-coverage check was re-run with `ALPHA_REF` added as a column:

| S | Mg/Si | alpha=2 | alpha=1.10 | alpha=10 |
|---|---|---|---|---|
| 0.8 | 1.25 | 298.33 K | 304.49 K | 280.01 K |
| 1.0 | 1.25 | 321.31 K | 327.22 K | 310.35 K |
| 1.0 | 0.50 | 346.11 K | 348.34 K | 335.99 K |
| 1.2 | 1.25 | \[out of domain, all three — alpha-independent\] |

`ALPHA_REF` is warmer than `alpha = 2` everywhere (lower alpha → weaker weathering → less
cooling), so the move **relaxes** the cold-end domain constraint rather than tightening it, and
the S = 1.2 wall is unrelated to `alpha`.

### 34.3 The caveat nobody wrote down

§22.9 listed four options for `alpha` and warned about option 3 — the best Earth fit — in these
terms: the net seafloor alkalinity flux falls to ~0.005 Teq/yr, **~180× below Coogan**, which
"effectively switches seafloor weathering off, which is fatal on land-free worlds." `ALPHA_REF =
1.100155` **is** that branch (§22.9 quoted 0.908 for it at the time).

Nothing in the code records that this warning was revisited when the decision was made. The
domain-coverage table above checks that the runs stay in the climate model's domain — it does not
check the seafloor alkalinity flux against Coogan. **Measure the net seafloor flux at the refit
value before the paper leans on the land-free thermostat.**

---

## 35. Charge balance, the Cl root cause, and the sedimentation rate (2026-09-09)

Entry point was a plotting question (log y-axis on `continental_baseline_ions`), which exposed the
ion panel: DIC 180× Earth, alkalinity 250×. The session then root-caused it.

### 35.1 The carbon excess is a charge-balance artifact, and Cl is over half of it

Tracked alkalinity equals the conservative-ion charge residual **exactly** — this is §7's design
working, not a bug:

```
Na 678.6 + 2(76.8) + 2(0.36) − 235.1 = 597.8 mEq    vs    tracked Alk = 597.7
```

So DIC is *slaved* to alkalinity, and alkalinity is a small difference of large numbers.
Decomposing the 587 mEq excess against Earth, at S = 1, land 0.3, 3000 m:

| ion | model | Earth | Δcharge | share |
|---|---|---|---|---|
| Cl | 235.1 | 550 | **+315** | 54% |
| Na | 678.6 | 480 | **+199** | 34% |
| SO₄ | 0.0 | 28 | +56 | 9% |
| Mg | 76.8 | 52.8 | +48 | 8% |
| Ca | 0.36 | 10.3 | −20 | −3% |
| K (untracked) | — | 10.2 | −10 | −2% |

**Cl never converges.** Measured τ_Cl = **5571 Myr** against a 2 Gyr integration; from a blank
ocean that reaches 30.2% of steady state, and 0.302 × 780 mM = 235 mM — the run's value to three
digits. This also explains why every baseline run terminates `timeout` rather than `converged`.
It is §15 item 3, promoted from a settling nuisance to the dominant error term.

Seawater Cl is an inherited inventory from early degassing, with modern volcanic Cl being recycled
subducted seawater Cl rather than primordial — [Kendrick et al., PNAS 2021](https://www.pnas.org/doi/10.1073/pnas.2116083118);
[Sharp & Draper 2013, EPSL](https://www.sciencedirect.com/science/article/abs/pii/S0012821X13001192).
The blank-ocean start is the unphysical case, so seeding is the literature-consistent choice, not
a convenience.

### 35.2 `K_CL_ANALYTIC` has been wrong since §33.3

`calibrate_earth.py:126` derives it from a steady-state balance in which the source and sink areas
cancel. **§33.3 broke that cancellation**: the Cl source is `F_vol`, over `surface_area`, while the
sink is over `seafloor_area`. It is missing a factor `A_surf/A_sf = 1/0.7 = 1.43` — the same 1.43×
`KD_MG_CALIB` and `K_NA_CALIB` absorbed, but here in closed form rather than through a fit.

Consequence: the Cl *steady state* is 780 mM, not the 550 targeted. This is an **input** to the
calibration (`calibrate_earth.py:269` passes it into the `Planet` the solver runs against), so it
biases the fit for `K_na`, `alpha` and `KD_mg` — it must be fixed **before** the refit, not after.

### 35.3 SO₄ is pinned at zero in every sweep

`planet.py:466` sets `F_net[so4_idx] = 0.0` by design, so sulfate can only enter through `b0`.
`calibrate_earth.py:244` seeds 23.45 mM; **no sweep seeds anything.** The calibration and the
production sweeps have been running oceans that differ by 56 mEq of charge.

Two follow-ons this raises, both **open decisions** rather than findings:

- **Sulfate is redox-coupled and the model does not couple it.** 28 mM is an oxygenated ocean;
  since §31 every sweep runs both redox arms, and the reducing arm should not carry oxic sulfate.
- **The seed's scaling law is a first-order lever.** Ocean depth spans 300 m – 50 km in the sweep,
  a **167× range in ocean mass**. Fixed-concentration and fixed-inventory seeding differ by that
  factor. This interacts with §33.13: `calibrate_earth.py` runs 3700 m (0.98× Earth's ocean mass)
  while the continental baseline runs 3000 m (0.79×), and `tau_prec` is depth-scaled, so the fit
  is anchored at 123 kyr and the "Earth" figure runs at 100 kyr.

### 35.4 The sedimentation rate now counts every precipitating phase

`S_sed` feeds `weathering.seafloor_reactive_area` through the burial timescale `t_cover`. It was
computed from carbon-as-calcite plus silicon-at-**quartz** density only, which both under-counted
(clays, evaporites and the reverse-weathering phases contributed nothing) and mis-counted (a mole
of Sepiolite(d) carries 6 mol Si but occupies 287 cm³, not 6 × 22.7 cm³ of SiO₂(am)).

Now each ocean-precipitating mineral contributes its own volume. `get_precipitation_by_mineral`
returns a fourth item — per-mineral **molar** rates, which the aqueous flux vectors cannot express
(they carry only tracked elements, so a phase's H and O are invisible and halite's mass is split
across two entries). `get_precipitation`'s 3-tuple signature is unchanged, so its ~15 call sites
were untouched; the clamped sum was extracted to `sum_precipitation` because `Planet` applies it
to the fast and reverse-weathering assemblages separately.

Molar masses are computed from the formula the **runtime database** uses, not from the mineral
name — Sepiolite(d) is the 6H₂O hydrate (647.8 g/mol) and Saponite-Na the non-integer Na₀.₃₄
endmember (386.7 g/mol); an idealised formula would be badly wrong for both. Densities are
Handbook of Mineralogy `D(meas.)`, or `D(calc.)` where none is given.

Measured effect:

| state | new/old | dominant phases |
|---|---|---|
| Earth seawater | 1.004× | Calcite 99% |
| model blank-run ocean | 1.203× | SiO₂(am) 77%, Calcite 16%, Sepiolite(d) 7% |
| Si-rich, low Ca | 1.081× | Calcite 86%, Sepiolite(d) 14% |

Negligible at Earth (calcite dominates and was already counted), ~20% in the model's own high-Si
ocean-world chemistry. End-to-end runs move <0.1%: `t_cover` combines harmonically with `t_clog`,
so a 1.2× change in `S_sed` is heavily damped. **One behaviour change beyond "count everything":
`SiO2(am)` density is now 2200, not quartz's 2650.**

A phase missing from either table now **raises** rather than being silently skipped — a silent
skip is exactly the failure this change exists to fix.

### 35.5 What a corrected charge balance actually does

Cumulative A/B, same planet, S = 1, land 0.3, 3000 m. `K_NA`/`KD_MG` here are hand-scaled by
§33.3's 1.41/1.45, **not** re-fitted, so row 4 is an indication of direction and magnitude only:

| | T | pCO₂ | Alk | C | Ca | Mg | Na | Cl |
|---|---|---|---|---|---|---|---|---|
| as-is | 295.2 | 824 | 597.7 | 363.4 | 0.36 | 76.8 | 678.6 | 235.1 |
| + seawater seed | 295.0 | 797 | 160.1 | 103.1 | 0.36 | 80.6 | 671.5 | 619.3 |
| + `K_CL` area fix | 295.0 | 796 | 226.0 | 143.0 | 0.36 | 79.4 | 669.4 | 549.0 |
| + `K_NA`,`KD_MG` ×1.41/1.45 | 294.3 | **687** | **3.4** | **3.1** | **12.4** | 72.3 | 437.2 | 549.2 |
| Earth | 288 | 280 | 2.3 | 2.0 | 10.3 | 52.8 | 480 | 550 |

DIC goes 363 → 3.1 mM. **Ca recovers 0.36 → 12.4 mM**, confirming the Ca collapse was an
alkalinity artifact, not a missing sink — and confirming §33.12's reading of the `fallback_limit`
runs, which sit in exactly that near-zero-Ca corner and should disappear with the seed.

The seed and the `K_CL` fix are **complementary, not redundant**: once the steady state *is* 550,
seeding at 550 makes Cl stationary and the 5.6 Gyr timescale stops mattering.

> ⚠️ **Two cautions.** (1) Row 3 → row 4 is violent: Alk 226 → 3.4 for a 1.41× change in `k_na`.
> Near Earth's balance alkalinity is a ~0.4% residual of ~600 mEq terms, so it is **ill-conditioned
> as a fit target** — fit the conservative ions and let Alk and DIC fall out. (2) pCO₂ is still
> 687 ppm, 2.5× Earth, after all of this. That is **not** charge balance and needs its own
> treatment.

### 35.6 Mg will not be fixed by the recalibration

The combined-fix run lands at Mg = 72.3 mM against Earth's 52.8 — **+36.9%**. §28.1 independently
reports **+37.0%** and attributes it to §27: removing the hedenbergite correction reaction shifted
the Earth assemblage (Diopside −25.7%, Forsterite +14.9%, Ca/Mg supply ratio −19.4%). Two
independent routes to the same number, so treat it as structural. `kd_mg_ht` trades Mg for Ca
mole-for-mole, so with Ca on target there is no way to pull Mg down. `kd_mg_ht` still has three
disagreeing anchors: Earth fit 0.019, Coogan HT 0.005–0.009, §18 first-principles 0.07, against
the shipped 0.0139.

### 35.7 MORB verified

§32.9's ⚠️ is closed: `MORB_OXIDES` matches Gale et al. (2013) Table 1's ALL MORB arithmetic mean
on all ten oxides. Details and the extraction trap are recorded there.

### 35.8 Status going into the recalibration

- ✅ MORB constant verified against the primary source.
- ✅ Sedimentation rate counts all precipitating phases.
- ✅ `pe` concern from §27.5/§15 checked and closed — calibration and sweeps share `pe = −3.0`.
- 🔴 **`K_CL_ANALYTIC` must be fixed before the refit** (§35.2) — it is an input to the fit.
- 🔴 **`planet.py` and `parameter_sweep.py` constants must be updated together** after the refit,
  or the drift check fires and every run is tagged non-comparable (§3).
- ⚠️ **Open decisions, needed before the run:** the SO₄ background and whether it couples to `pe`;
  fixed-concentration vs fixed-inventory seeding; and whether the continental baseline moves to
  3700 m to match the calibration anchor or the calibration to 3000 m to match the sweeps (§33.13).
- ⚠️ `alpha`'s seafloor alkalinity flux was never checked against Coogan at the 09-01 value (§34.3).
- ⚠️ Old sweep output has been moved aside; the new sweep starts from an empty directory, so the
  §32.11 "every stored result predates this table" caveat no longer applies.

---

## 36. The sink audit: reverse weathering, the Mg sink, and the Na sink (2026-09-17)

Three questions, asked in sequence, each answered by measurement rather than by reading the budget
shares: **is reverse weathering doing anything?**, **what happens if the HT Mg–Ca exchange is
switched off?**, and **could nahcolite replace albitization as the Na sink?** The answers share one
shape — *every secondary sink in this model is a backstop whose setpoint sits outside the ocean
states the model actually produces* — and the Na half of it exposes a structural duplication that
had not been noticed before (§36.9).

**Nothing was changed in `src/kamino`.** The only repository change is a new `basic_no_rw` sweep in
`parameter_sweep.py` (§36.11), which has not been run.

### 36.0 Method, and one trap avoided

Two independent lines of evidence throughout:

1. **Replay.** `Planet` rebuilt from each run's saved config, `dY_dt` evaluated once on the stored
   final state, and `_flux_terms` / `_state` read back. 900 runs sampled from `sweep_output`,
   stratified over (tag, out, crust, mgsi, diw, depth); 895 succeeded. Population-weighted
   estimates (reweighting each stratum back to the full 9 721 runs) agree with the unweighted ones
   to within 0.02 everywhere, so the sample is representative.
2. **Controlled A/B.** Both arms re-run from scratch through `parameter_sweep.run_simulation` at
   the current constants, rather than differencing against output on disk. 113 runs total: 40 RW
   on/off pairs, 39 HT arms, 34 Na arms.

§16's warning was respected: nothing is concluded from `dY_dt` on a fixed `Y` alone. Every causal
claim below comes from a paired `time_evolve`.

⚠️ **`diagnostics.planet_from_config` is stale and `diagnose_run` raises as written.** It passes
`crust_composition` and `f_bio`, neither of which is a `Planet` parameter, and it drops `pe`,
`tau_prec`, `tau_rw`, `kd_mg_ht`, `k_na_cont_removal`, `k_cl_subduction` and `water_rock_ratio`
from the config. `diagnose` also omits `pe=planet.pe` on its per-mineral call, which matters for
the oxidising arm. A private replacement was used for this audit; the module was left alone.

### 36.1 Reverse weathering is a bimodal switch, not a background process

In ~60% of the sweep the RW flux is numerically zero. Where it fires it is frequently the largest
single sink in the ocean. Population-weighted over 9 721 runs:

| RW share of total sink | Alkalinity | Mg | Si |
|---|---|---|---|
| > 1% | 35% of runs | 39% | 34% |
| > 10% | 19% | 24% | 22% |
| > 50% | **7%** | **14%** | **12%** |

Against the alkalinity *source* (`seafloor LT` + continental), `|F_rw|` exceeds 5% in 23.5% of runs,
20% in 14.5%, and 50% in 7.1%.

**Sepiolite(d) carries ~98% of it** (median share of the RW alkalinity flux 0.976). The other two
entries in `reverse_weathering_minerals` do essentially nothing:

- **Saponite-Na** is Al-limited to ~10⁻⁵ Tmol/yr, exactly as §6.1 predicted.
- **Greenalite never saturates anywhere in the sweep** — median SI −16.19, **maximum −1.00**,
  `frac SI>0 = 0.000`. See §36.5.

So RW is an **Mg and Si** sink. It is not a Na sink: Na removal is 95.5% albitization, 4.5% shelf
carbonate, and 0.01% RW.

The active corner is **low crust production (0.01–0.1), low-to-moderate outgassing, and
T_seafloor above ~320 K**. At crust ≥ 3 it is off everywhere; below 300 K it is off. Sepiolite SI
by seafloor-temperature band:

```
T_sf       270–280  280–300  300–320  320–340  340–360
median SI    −3.32    −3.23    −2.47    +2.01    +2.95
frac SI>0     0.37     0.32     0.38     0.58     0.74
```

Where RW is active, the Mg residence time against it alone has median 697 Myr and **p10 = 5 Myr —
exactly `tau_rw`**. That is §26.3's result restated from the other end: Sepiolite is so far
supersaturated that `tau_rw` *is* the flux, and the whole excess inventory is removed on that
timescale.

### 36.2 The mechanism is retrograde solubility, and it makes RW a CO₂ source

The runtime database settles it:

```
Sepiolite(d)  Mg4Si6O15(OH)2:6H2O + 8 H+ = 4 Mg+2 + 6 SiO2 + 11 H2O
              -delta_H -157.339 kJ/mol
```

Dissolution is strongly exothermic, so the mineral is **less** soluble when hot. Hotter ocean →
more Sepiolite → less alkalinity → more CO₂ → hotter still. Precipitation consumes 8 eq of
alkalinity per formula unit, i.e. **Alk : Mg = 2 : 1**, which is exactly the ratio §21.3 measured
(−34.10 Alk against −17.05 Mg).

This is the same retrograde behaviour §21.2 established for `b_eq`, now attributed to a specific
phase and its measured enthalpy rather than inferred from a temperature scan. **RW is a positive
climate feedback and a CO₂ source. It is never a carbon sink.**

### 36.3 The controlled A/B: 40 pairs, `reverse_weathering` True vs False

| | |
|---|---|
| median \|ΔT\| | **0.000 K** |
| pairs with \|ΔT\| > 1 K | **4 / 40** |
| max ΔT | **+23.9 K** |
| pairs with pCO₂ shifted > 2× | 3 / 40 |

The sign is always the same: RW on is hotter. The extreme case is S = 0.85, out = 0.01,
crust = 0.01 (297.8 K / 0.083 bar → 321.7 K / 0.46 bar), and its budget shows the mechanism:

```
                        RW on      RW off
seafloor LT alk       +1.141      +0.172     <- 6.6x, driven by the extra heat
reverse weathering    -0.978       (none)    <- 87% of the total sink
```

RW removes alkalinity until the planet is hot enough for seafloor weathering to outrun it; the
budget re-closes at a far hotter attractor. This is §26.2's `tau_rw` 5 → 33 Myr result (−25 K,
pCO₂ ÷17) reproduced by removing the sink entirely rather than by slowing it.

**At Earth it is calibrated-plausible and climatically irrelevant.** land = 0.3, S = 1, out = 1,
crust = 1, both arms converged: RW removes **0.083 Tmol Mg/yr** against Dunlea et al. (2017)'s
0.02 Tmol/yr observational authigenic Mg sink — the right order, ~4× the lower bound — and it is
0.9% of the alkalinity sink. Switching it off costs **0.09 K** (pCO₂ 686 → 673 ppm).

The two continental pairs at out = 0.01 both hit the wall cap with pCO₂ at 10⁻¹⁰ bar in a snowball
state; their apparent 13× and 0.25× pCO₂ ratios are numerical noise on a vanishing number and
should not be read as signal.

### 36.4 ⚠️ Sepiolite's setpoint is ~100 mM Mg — too high to be a working sink

This is the finding that generalises. RW only engages once ocean Mg exceeds roughly 100 mM, against
Earth's 53 mM. Between Earth-like Mg and that threshold **the model has exactly one Mg sink and no
redundancy**. RW looks negligible at the Earth anchor (0.09 K) not because it is weak but because
it has not switched on yet — and §36.6 shows it is the only thing preventing a 312 K runaway once
the HT exchange is removed.

**It is a backstop, not a contributor.** That distinction matters for how the paper describes it.

### 36.5 Why Greenalite is inert: the ocean never receives any iron

Not a race between Fe phases — both precipitation calls in `dY_dt` (lines 336 and 343) take the
**same** `b_ocean`, so nothing consumes Fe before Greenalite sees it. Three measurements:

| phase | median SI | max SI | frac SI > 0 |
|---|---|---|---|
| Greenalite | −16.19 | **−1.00** | **0.000** |
| Goethite | −3.21 | +3.33 | 0.234 |
| Siderite | −4.92 | +1.45 | 0.113 |

`SI(Goethite) − SI(Greenalite) > 0` in **100.0%** of states (median +11.2); the Siderite comparison
holds in 98.8%. Greenalite is the least favourable Fe phase everywhere.

The cause is **Fe starvation, upstream of the ocean**. Median ocean Fe is 3.4 × 10⁻¹⁰ mol/kgw and
69% of runs sit below 10⁻⁷. The seafloor source delivers Fe ~10⁷× more slowly than Mg:

```
seafloor LT  Fe : median 1.8e-07 Tmol/yr
seafloor LT  Mg : median 6.9e-01
seafloor LT  Si : median 2.7e+00
```

**Pore-space Goethite is what removes it.** In one run checked directly, pore Goethite sits at
**SI +4.66** while ocean Goethite in the same run is at **−4.76**, and the net Fe delivered to the
ocean is −2 × 10⁻⁹ Tmol/yr. Primary dissolution liberates Fe; `pore_precipitating_minerals`
(= `clay_minerals` = Kaolinite + Goethite, per §21.3's fix) captures it before the fluid reaches the
ocean.

Stoichiometry finishes the job. Greenalite is Fe**₃**Si₂O₅(OH)₄, so SI moves 3 decades per decade of
Fe. A counterfactual scan at fixed T, pH and pe needs **10⁴–10⁵× more Fe** to reach saturation.
Redox cannot rescue it either — scanning pe at fixed composition, Greenalite is flat at −9.0 for
pe ≤ −1 and collapses above it (−13.5 at pe 4, −25.5 at pe 8) while Goethite rises monotonically.

In the 314 runs where Fe *is* removed from the ocean, it is Goethite (205) or Siderite (100),
**never both and never Greenalite** — the two split cleanly by redox, as §28.3's `pe` note implies.

✅ **Conclusion: the Fe arm of `reverse_weathering_minerals` is dead code in practice.** Greenalite
is a primary BIF precipitate requiring anoxic, ferruginous, silica-rich seawater (Rasmussen et al.
2019; Johnson et al. 2018), which this model's oceans never are — because pore Goethite oxidises
the Fe out upstream. Note also that Tosca et al. (2021) find a *small* amount of Fe(III) triggers
greenalite nucleation in simulated Archean seawater; the model routes all Fe(III) to Goethite and
structurally cannot represent that pathway.

### 36.6 Switching off the HT Mg–Ca exchange: RW takes over, but Ca is the casualty

13 configurations × 3 arms (`kd_mg_ht` at the calibrated value; `kd_mg_ht = 0`; `kd_mg_ht = 0` with
RW also off), 39 runs. The cleanest ocean-world case, S = 0.85, out = 0.1, crust = 1:

| arm | T | pCO₂ | Mg | Ca | F_HT[Mg] | F_RW[Mg] | SI(Sep) |
|---|---|---|---|---|---|---|---|
| HT on | 288.3 K | 0.034 bar | 17.5 mM | 178 mM | −0.688 | **0** | **−6.40** |
| HT off, RW on | 322.7 K | 0.498 bar | **307 mM** | **0.19 mM** | 0 | **−1.88** | **+9.59** |

RW goes from completely dead to carrying the entire Mg sink — but only after Mg rises 18×, and the
planet warms **34 K** with **15× pCO₂**. The same at S = 1.0 / out = 0.1 / crust = 1 (Mg 17.5 → 209
mM, +36 K, 94× pCO₂), S = 1.0 / out = 0.01 / crust = 1 (Mg 0.087 → 192 mM, a factor of **2200**,
+36 K) and S = 0.85 / out = 0.01 / crust = 1 (Mg 0.085 → 217 mM, +36 K).

**The larger damage is to Ca, not Mg.** On a land-free planet the HT exchange is ~99.9% of the Ca
source (the LT flux delivers 0.0002 Tmol/yr against HT's 0.15). Removing it collapses Ca by two to
three orders of magnitude — 178 → 0.19, 177 → 1.33, 191 → 2.21, 78.9 → 0.33 mM. With no Ca there is
no calcite, the carbon sink dies, and *that* is what drives the warming. **`kd_mg_ht` is carrying
the carbon cycle on ocean worlds, not just the Mg budget.**

**Mg does not run away with both sinks off.** S = 1.0 / out = 0.1 / crust = 1 settles at Mg = 278 mM
with **F_LT[Mg] = −0.0033 Tmol/yr** — the seafloor source has gone *negative*, because the pore fluid
saturates against the primary Mg silicates and dissolution stops. There is a thermodynamic backstop
at roughly 250–500 mM. This **refines §12's "Mg simply accumulates"**: the flux statement is right,
but the endpoint is bounded by source shutdown rather than unbounded.

The both-off arm also degrades numerically, reproducing §12's pathology as a controlled experiment:

```
HT_on        converged 1   timeout 8   out_of_domain 4
HToff_RWon   converged 1   timeout 6   out_of_domain 6
HToff_RWoff  converged 0   timeout 3   out_of_domain 6   wall_timeout 3   fallback_limit 1
```

Only the both-off arm produces zero converged runs and the only `fallback_limit`.

**The continental case separates the three arms cleanly** (land = 0.3, S = 1, out = 1, crust = 1):

| arm | T | pCO₂ | Mg | Ca | Alk | termination |
|---|---|---|---|---|---|---|
| HT on | 294.33 K | 686 ppm | 74.4 mM | 10.1 mM | 3.56 mM | converged |
| HT off, RW on | 295.14 K | 814 ppm | 199 mM | **0.0 mM** | **262 mM** | converged |
| HT off, RW off | **312.12 K** | **14 779 ppm** | — | — | — | **wall_timeout** |

With continents the climate barely notices losing the HT exchange (+0.8 K) because continental
weathering supplies the alkalinity — while the ocean becomes chemically absurd (Ca → 0,
Alk 74× modern). ⚠️ **Temperature alone does not diagnose this model.** Lose RW as well and the
planet goes to 312 K at 21× pCO₂ and never converges.

**Regime caveat:** at out = 1, crust = 0.01 all three arms are identical (Mg 93.3 / 94.1 / 94.0),
because `F_HT` was already only −0.037 and RW already ~zero. None of this is universal.

### 36.7 No magnesium carbonate can fill the gap

Pre-precipitation SI for every Mg phase in the runtime database, across 298 ocean states. (SI is
only reported for phases in `available_mineral_string`, which excludes these — the survey had to
extend it explicitly, which is worth knowing before anyone repeats the measurement.)

| phase | median SI | frac SI > 0 (all) | **frac SI > 0 where RW is off** |
|---|---|---|---|
| Dolomite | +0.19 | 0.537 | 0.360 |
| Huntite | −3.41 | 0.332 | 0.180 |
| Magnesite | −0.81 | 0.265 | 0.143 |
| Nesquehonite | −3.10 | 0.037 | **0.012** |
| Brucite | −5.90 | 0.081 | **0.000** |
| Artinite | −6.14 | 0.040 | **0.000** |

The split lands exactly on the kinetic trap §12 identified:

- The phases that **would** take Mg — Dolomite, Huntite, Magnesite — are the ones that do not
  precipitate abiotically at low temperature. Land (1998) failed to nucleate dolomite from
  supersaturated solution at 25 °C over a **32-year** experiment; the inhibition is attributed to
  Mg²⁺ dehydration kinetics, and magnesite's scarcity in modern surface environments rivals
  dolomite's for the same reason.
- The phases that **are** kinetically defensible — nesquehonite, artinite, brucite — are
  supersaturated in **0.0–1.2%** of the runs where RW is off. Zero in exactly the corner that needs
  a sink.

This confirms §12 at the current code state and sharpens it: §12 measured nesquehonite −0.40 and
artinite −0.91 at one state; across 298 states they sit at −3.10 and −6.14 median.

There is also a **thermodynamic** reason independent of kinetics. Where RW is off the median state
is **pH 5.89, pCO₂ 0.80 bar, T_seafloor 277 K** — acidic and cold. Carbonate solubility rises with
pCO₂, so the regime that needs a Mg sink is structurally the regime where carbonates are least able
to form. Adding one cannot help there.

⚠️ One recent wrinkle if this is revisited: Kim et al. (2023, *Science*) show dolomite *can* grow
near ambient conditions under dissolution–recrystallisation cycling. That weakens the blanket
kinetic prohibition but requires a fluctuating saturation state this steady-state model does not
represent.

### 36.8 The Na source is zero, not small — crustal Albite is inert

Albite is **14.92 wt%** of the Earth-reference crust (Mg/Si 1.25, ΔIW −2), a major phase. The
measured seafloor LT flux nevertheless gives a **Na/Mg ratio of 5.9 × 10⁻⁴** against a crust molar
Na/Mg of ~0.15 — a ~250× suppression.

The weathering law's driving force is `b_eq − b_input` (`weathering.py:93`). Scanning Albite's
saturation index in the pore fluid against ocean Na:

| ocean Na | SI(Albite) |
|---|---|
| 0 | −5.89 |
| 0.01 mM | −1.89 |
| **1 mM** | **+0.11** ← crosses saturation |
| 480 mM (seawater) | **+2.74** |
| 6400 mM | +3.88 |

Above ~1 mM Na, Albite is supersaturated. Primary minerals carry PHREEQC's `dissolve_only`
modifier (`chemistry._equilibrium_block`), so a supersaturated primary phase can neither dissolve
nor precipitate — it is **completely inert**. At seawater Na the driving force is slightly
*negative*.

**It is Albite's own solubility, not competition for Al:**

| assemblage | ocean Na = 0 | ocean Na = 480 mM |
|---|---|---|
| Albite alone (100 wt%) | +0.0017 mM Na | **−0.0006 mM** |
| full crust | +0.0002 mM | +0.05 mM † |
| full crust, no Anorthite | +0.0010 mM | +0.05 mM † |

† water/rock mass bookkeeping, not dissolution — see below.

Even 100 wt% Albite into Na-free water yields **0.0017 mM** equilibrium Na, five orders of magnitude
below seawater; stripping Anorthite changes it by a factor of 5. Albite is simply very insoluble at
the pore fluid's pH ~8.0–8.6, which is the feldspar solubility minimum. For comparison the same
equilibration gives `b_eq[Mg] = 78 mM` from Forsterite + Diopside — the mafic phases deliver
~350 000× more Mg than Albite delivers Na.

⚠️ **A diagnostic trap worth recording:** `b_eq[Na] / b_in[Na] = 1.0001` at *every* ocean Na from
10⁻³ to 6400 mM. That constant 0.01% is the water/rock mass balance, not dissolution. Read as a
chemical signal it looks like a tiny but real Na source; it is not one.

**Consequence: land-free planets have no Na source at all.** With albitization on, it drains the
480 mM seawater seed to **0.0076–0.008 mM** — and, because `F_na_rw` removes Alk 1:1 with Na
(`planet.py:376`), it takes ~480 mM of alkalinity with it. Switching it off leaves Na pinned at the
seed (468–470 mM) for 2 Gyr, because there is nothing to move it.

That also means switching it off **warms** ocean worlds, by removing an alkalinity sink:

| S = 1.0, out = 1, crust = 1 | T | pCO₂ | Na | Alk |
|---|---|---|---|---|
| albitization on | 341.6 K | 0.365 bar | **0.0076 mM** | 5.85 mM |
| albitization off | **352.1 K** | **2.96 bar** | 470 mM | **245 mM** |

⚠️ Ocean-world Na of 0.008 mM is precisely the low-Na state **§6.2** identifies as killing HT Ca
release. The over-drained Na and the HT Ca collapse are the same problem seen from two sides.

### 36.9 This quantifies §6.2's "Albite is the switch" — and §6.3 already tried the obvious fix

§6.2 established that **Albite is the switch** controlling the HT path: *"Low ocean Na drives Albite
dissolution, which floods Na and simultaneously collapses Ca release."* §36.8's saturation scan is
the same curve measured from the other side, and it puts a number on the switch point:

```
SI(Albite) crosses zero at ocean Na ~= 1 mM
   below it  -> undersaturated -> Albite dissolves and floods Na   (the 6.2 failure mode)
   above it  -> supersaturated -> dissolve_only makes it inert     (the 36.8 zero source)
```

Both observations are one saturation boundary. §6.2 saw it at 473 K in the HT path at Na = 27 mM;
§36.8 sees it at ~290 K in the LT path at Na = 480 mM. The HT path is currently off (`f_HT = 0.0` in
every sweep), so only the second regime is live.

**The structural oddity is that albitization is represented twice, in opposite directions, and the
two are not coupled:**

- The crust's Albite sits as an inert `dissolve_only` primary phase, supersaturated at SI +2.74.
- The Na sink is a **separate fitted parameter**, `-k_na . b_Na . J_total` (`planet.py:375`), with
  no thermodynamic link to that supersaturation.

PHREEQC knows Albite should be precipitating; the model never asks it to. Hence `k_na = 0` gives no
sink *at all* rather than degrading to a thermodynamic fallback, and the fitted linear form
over-drains ocean worlds to 0.008 mM instead of shutting off near the saturation boundary.

⚠️ **The obvious fix has already been tried and failed — see §6.3.** "Non-`dissolve_only` + Albite
precipitate-only" was the chemically best PHREEQC result of that arc (dCa/-dMg 0.8-1.1, Ca robust at
low Na), and **the full sweep still hothoused everything**, with *"Na still drained to 0
(Albite-precip-only removed the only Na source on a landless world)"*. §36.8 explains why that
happened: making Albite precipitate-only removes the one regime (Na < ~1 mM) in which it could have
been a source, so Na has a sink and no source and goes to zero — which is exactly what the *current*
parameterisation also does, by a different route.

So this is **not** a recommendation to retry it. The measurement says something narrower and more
useful: **on a land-free world the model has no Na source in either configuration**, and the fitted
`k_na` is therefore not balancing a flux — it is draining a seed. That is worth stating explicitly
in the methods, because it means ocean-world Na is set by the initial condition and the sink
constant, not by any weathering process.

The Earth continental case is the only configuration with a genuine Na source (3.65 Tmol/yr
continental), and it is where the missing thermodynamic fallback bites:

| arm | T | pCO2 | Na | Alk | termination |
|---|---|---|---|---|---|
| albitization on | 294.33 K | 686 ppm | 427 mM | 3.56 mM | **converged** |
| albitization off | 293.3 K | 548 ppm | **6417 mM** | **5887 mM** | timeout |

Na runs away 15x and alkalinity 1650x. §6.1's requirement — *"Na must still reach steady state or
the ocean runs away"* — is exactly what fails here, and because `F_na_rw` debits alkalinity 1:1
with Na (`planet.py:376`), the alkalinity follows it up. Note again that the **climate barely
moves** (-1 K) while the ocean becomes a 6 M brine.

### 36.10 Nahcolite cannot be the Na sink

Nahcolite (NaHCO₃) is already in `carbonate_minerals`, therefore already in
`fast_ocean_precipitating_minerals` — it is one of the Na-carbonates §6.1 added to
`make_database.py` precisely because *"the only Al-free Na sinks are Na-carbonates"*. **No code
change was needed to test it. It has simply never fired.**

Across 895 replayed sweep states its SI has median −4.58 and **maximum −0.00**, with
`frac SI > 0 = 0.0000`; across all 34 albitization runs `F_prec[Na] = 0` in every cell. The closest
approach anywhere is **SI −0.28**, in an ocean holding 6.4 M Na and 5.9 M alkalinity.

It does, however, pin at SI ≈ 0 in the most extreme sweep states (Na 1300–3200 mM, Alk 600–1900 mM),
so it is acting as a ceiling — at a setpoint 3–6× seawater Na. The same pattern as Sepiolite in
§36.4, one notch further out of reach.

This is physically correct rather than a gap. Nahcolite is the stable Na-carbonate above
~1125 ppm CO₂ — the basis of the Green River nahcolite palaeobarometer (Lowenstein & Demicco 2006;
Jagniecki et al. 2015 give 680–1260 ppm) — and it is a **closed-basin lacustrine evaporite**
concentrated far beyond seawater, not a marine phase. Modern seawater at 470 mM Na does not
precipitate nahcolite and a model of open seawater should not either.

Worth recording for the carbon budget: nahcolite removes Na, alkalinity **and carbon** 1 : 1 : 1,
whereas albitization removes only Na and alkalinity. At Earth's 3.65 Teq/yr Na sink, routing that
through nahcolite would add a ~3.65 Tmol/yr carbon sink — comparable to the entire shelf carbonate
sink (4.06). It would matter *if it could fire*. It cannot, at marine Na.

If a second Na sink is wanted, the marine candidates are authigenic Na-clay (Saponite-Na, already
present and Al-limited to ~10⁻⁵ Tmol/yr) or halite in restricted basins (already present as
`evaporite_minerals`, gated on `land_fraction > 0`). Neither is a carbonate.

### 36.11 `basic_no_rw`, and what was left alone

`parameter_sweep.py` gained one sweep and nothing else:

- `reverse_weathering_off = [False]` beside `reverse_weathering_default`
- `sweep_basic_no_rw()` — identical to `sweep_basic` on every axis, RW off
- registry entry and `_sweep_size` sizer

**1862 runs, ≈84 CPU-h.** Verified before shipping: both grids are 1862 combos differing *only* in
the RW flag, and there are **zero filename collisions** with the existing `basic` output, because
`_run_name` appends `_rw` only when reverse weathering is on — these land as the untagged variant
(`..._depth_3000_mgsi1.25_diw-2` against `..._depth_3000_rw_mgsi1.25_diw-2`). Every new run pairs
one-to-one with one already on disk. **Not run.**

⚠️ If it is run for comparison, consider re-running `basic` alongside rather than differencing
against the existing `sweep_output`: those runs were produced across several sessions and predate
the §33.3 seafloor-area fix, so the constants baked into them may not match. `_warn_constant_drift()`
reports this at launch.

### 36.12 Status

- ✅ Reverse weathering characterised: bimodal, Sepiolite-carried, a CO₂ source, negligible at the
  Earth anchor but the only backstop once `kd_mg_ht` is removed.
- ✅ Greenalite shown inert, with the cause traced upstream to pore Goethite.
- ✅ The Mg-carbonate question closed: no phase in the database can fill the gap.
- ✅ The Na source shown to be zero by thermodynamics, not small by kinetics.
- ✅ Nahcolite ruled out as a Na sink, with a literature basis for why that is correct.
- 🔴 **`diagnostics.planet_from_config` is broken** (§36.0) — `diagnose_run` raises as written.
- ⚠️ **Sepiolite's ~100 mM setpoint and Nahcolite's ~1500 mM setpoint** mean both "sinks" are
  inactive across the entire Earth-like part of parameter space. If the paper describes either as
  part of the steady-state budget, that needs qualifying.
- ⚠️ **Land-free worlds have no Na source in any tested configuration** (§36.8-36.9). `k_na` is
  draining the seawater seed, not balancing a flux. The obvious thermodynamic fix was already tried
  and failed in §6.3; this should be stated as a model limitation rather than re-attempted.
- ⚠️ `Greenalite` can be removed from `reverse_weathering_minerals` with no effect on any result,
  or kept as documentation of a pathway the model cannot currently reach.

### 36.13 Cross-cutting lessons

- **A sink that is inactive at the calibration point is not a negligible sink.** RW costs 0.09 K at
  Earth and prevents a 312 K runaway two experiments later. Calibration-point sensitivity is the
  wrong test for whether a term belongs in the model.
- **Temperature is a poor diagnostic for this model.** Two separate experiments produced oceans with
  Ca = 0 and Alk 74–1650× modern while T moved by less than 1 K. Check the ion budget, not the
  climate.
- **Ask what the model's own thermodynamics already says before adding a phase.** Nahcolite was
  already in the precipitating list and had never fired; the answer cost one SI query, not a
  new sink.
- **A supersaturated `dissolve_only` phase is invisible.** Albite is 15 wt% of the crust and
  contributes exactly nothing, and nothing in the output says so. §21.6's "ask which reservoir a
  sink acts on" has a companion: *ask which side of saturation a phase is on*.
- **Check the history before proposing a fix.** "Let Albite precipitate" looked like the obvious
  correction to §36.9 until §6.3 turned out to have tried it and recorded why it failed. The new
  measurement's value was in *explaining* that failure, not in reversing it.
- **Flux shares and causal importance are different quantities.** HT exchange is 89% of the Mg sink
  by flux, but its removal damages the Ca and carbon budgets far more than the Mg budget.

### 36.14 References added this session

- **Dunlea et al. (2017)**, *Nat. Commun.* **8**, 844 — already in §17; used here as the
  observational check on the Earth RW Mg flux (0.02 against the model's 0.083 Tmol/yr).
- **Rasmussen et al. (2019)**, *Precambrian Res.* — widespread greenalite deposition forming BIFs
  before the GOE; greenalite as a primary precipitate from anoxic ferruginous seawater.
- **Johnson et al. (2018)**, *Geophys. Res. Lett.* **45** — low-Fe(III) greenalite as a primary
  Neoarchean ocean mineral.
- **Tosca et al. (2021)**, *Geology* **49**, 905 — ferric iron triggers greenalite formation in
  simulated Archean seawater; the pathway this model cannot represent.
- **Land (1998)** — the 32-year failure to precipitate dolomite at 25 °C; the canonical statement of
  the dolomite problem.
- **Kim et al. (2023)**, *Science* — dissolution enables dolomite growth near ambient conditions;
  the partial counter-argument.
- **Lowenstein & Demicco (2006)**, *Science* **313**, 1928 — nahcolite + halite coprecipitation
  requires pCO₂ > 1125 ppm; Eocene CO₂.
- **Jagniecki et al. (2015)**, *Geology* **43**, 1075 — Eocene atmospheric CO₂ from the nahcolite
  proxy, 680–1260 ppm.

---

## 37. The paper review, and the model changes it forced (2026-09-24)

The session started as a line-by-line check of the paper draft, `ocean_chemistry.tex`, against the code. The
issue list lives in `paper_issues.md` at the repo root and is the authority on what is wrong with the *draft*.
This section records what the review found about the *model*, and what was changed as a result.

Everything is uncommitted. The sweep running when the session began (`basic_low_mgsi` then
`basic_high_mgsi`, blank ocean, files dated 23–24 Sep) used the pre-change code and is superseded.

### 37.1 The paper review

The first pass found ~105 issues. At the user's request, the equations were then corrected directly in the
draft:
- **τ_cover:** the `V_o/d` factor had the wrong units; it is now `ρ_w d`, with the floor and terrigenous terms
  added.
- **f_diss:** `ḃ_prec J` had the wrong units; it is now `J Δb_p`, because pore precipitation is instantaneous.
- **θ_r:** the notation was made consistent, and it is no longer called a "fraction" (it exceeds 1).
- **Da:** now the charge-weighted version the code computes.
- **Eq. (1):** now includes `F_cont` and `F_shelf`.
- **Eq. (3):** the missing inputs were added.
- **Eq. (6):** the pCO₂ function label was wrong.
- **LaTeX:** the pCO₂ double subscript was fixed.
- **Sedimentation rate:** renamed from `S` to `ḣ_sed`, because `S` is the instellation.

Facts established along the way, now recorded in the issue list:
- The Fig. 1 caption had the wrong sign: olivine *increases* with ΔIW.
- 40 ocean masses is 0.9 % of Earth's mass, not 0.5 %.
- Ice VI forms at ~0.63 GPa at 273 K; 1 GPa applies only near 300 K.
- The 274 K seafloor floor engages whenever the surface is below 285 K.
- The blank sweep spans pH 4.75–11.8, not the draft's 5.5–8.

### 37.2 How `k_p` weights minerals: weight fraction is defensible, mole fraction is not

The draft said `k_p` is weighted by mole fraction; the code (`get_k`) uses weight fraction. Physically, the
weight should be each mineral's share of the reactive surface area. Two standard approximations exist
(Beckingham et al. 2016):
- weight fraction, if every mineral has the same surface area per gram;
- volume fraction, if every mineral has the same grain size.

Mole fraction depends on how the formula unit is written, so it has no physical basis. Measured on the
alkalinity-weighted `k` at 290 K and pH 8:

| crust Mg/Si (ΔIW −2) | 0.5 | 0.8 | 1.25 | 1.6 | 2.0 |
|---|---|---|---|---|---|
| volume / weight | 0.98 | 0.78 | 0.76 | 0.82 | 0.85 |
| mole / weight | 0.41 | 0.94 | 1.15 | 1.01 | 0.99 |

The code was kept; the draft's wording was fixed. Most of the volume-weighting offset is uniform and would be
absorbed by α.

### 37.3 Code questions found and left open

- **The weathering law is applied element by element** (`weathering.py`, `F_primary`). In the kinetic limit,
  each element is scaled by its own affinity factor (1 − b_in/b_eq), not by each mineral's saturation state.
  At Earth pore conditions (~300 bar, 286 K, 30 Pa CO₂, pe −3):

  | pore-fluid input | Si | Mg | Ca | Na |
  |---|---|---|---|---|
  | seawater | 0.973 | 0.121 | 0.001 | 0.000 |
  | blank | 1 | 1 | 1 | 1 |

  So dissolution is not stoichiometric in seawater-like oceans. `k` also sums over minerals that their own
  saturation state says shouldn't dissolve. This is probably part of why the seafloor Ca source is ~0
  (§36.6). The standard alternative is per-mineral transition-state kinetics, r_i = k_i A_i (1 − Ω_i^{1/σ})
  (Aagaard & Helgeson 1982; Lasaga 1984). The `*_rate(T, pH, omega)` functions already exist, but combining
  them with the transport limit has no closed form.
- **The pore fluid is held at the atmospheric pCO₂**, via a CO₂ gas phase at ~300 bar in `_equilibrium_block`.
  This is an open system, and it is how pCO₂ enters the weathering law.
- **The sediment volume has no porosity.** Deep-sea sediment in the top ~100 m is 60–80 % water
  (Hamilton 1976), so τ_cover is probably ~3× too long. The terrigenous 5 m/Myr term is presumably a bulk
  thickness, so the two terms may be on different bases.
- **The §7 positivity clamp is inactive.** None of the 2,498 finished runs in the blank sweep ends with
  alkalinity pinned at 0 against a negative charge balance.

### 37.4 Crust production rate: 1/50 → 1/130 Myr (changed by the user)

`EARTH_CRUST_PRODUCTION_RATE_PER_AREA` is now `1/(130e6 YR)`, derived from 2.7 ± 0.2 km²/yr. Cogné &
Humler (2004) give a half-rate of 1.30 ± 0.28 km²/yr, and the draft now cites Cogné et al. (2006). J⊕ is
unaffected, because J = J⊕ R̃; the rate enters only θ_r.

Replaying stored final states confirmed that the running job used 1/50 throughout: the stored Damköhler
numbers are reproduced exactly with 1/50 and are 10–50 % off with 1/130, for runs from both sweeps. **Note:
`parameter_sweep.py` starts a new `ProcessPoolExecutor` for each sweep in a job, so later sweeps in the same
job do import an edited `src/`.**

### 37.5 Shelf precipitation is now area-scaled

**Before.** The shelf term was the whole ocean relaxing towards calcite saturation at 1 km pressure, with the
full τ_prec. It was added on top of the deep precipitation and switched on at full strength for any γ > 0,
whatever the land area. So land planets relaxed calcite twice as fast as ocean worlds, and the land-fraction
series had a step at γ = 0. Replaying the 9 Sep calibrated Earth: outgassing +7.50 Tmol C/yr, deep carbonate
burial −3.32, shelf −4.18 (56 % of the burial).

**After.** Carbonate precipitation is split by area rather than counted twice:
- F_carb = (1 − f_s) F_deep + f_s F_shelf, where f_s is the shelf's share of the seafloor.
- f_s = min(EARTH_SHELF_AREA × land_area / A_LAND_EARTH / seafloor_area, 1).
- `EARTH_SHELF_AREA = 2.7e13` m² (Harris et al. 2014: ~27 million km², ~7 % of the ocean area).

| γ | 0 | 0.001 | 0.03 | 0.3 | 0.9 |
|---|---|---|---|---|---|
| f_s | 0 | 1.8e-4 | 0.0055 | 0.0756 | 1 (capped) |

- **Sedimentation rate unchanged.** `moles_prec` is left unweighted, because the deep share of the carbonate
  lands on the deep share of the seafloor.
- **Land-free runs unchanged,** because the split is guarded on f_s > 0.
- **The old calibration no longer balances.** At the old calibrated state, burial becomes deep −3.07 plus shelf
  −0.32 (9 %) = −3.39 Tmol/yr against +7.50 of outgassing, so the calibration had been relying on the
  double-counted sink.
- **Still open:** the shelf is still 1 km deep and uses the deep-water temperature. The mean shelf break is
  140 m (Harris et al. 2014).

### 37.6 SO₄ and K⁺ as explicit fixed backgrounds

**Before.** The calibration used SO₄ = 23.45 mM, derived from the charge balance to give alkalinity 2.3
meq/kg. It was standing in for K⁺ and every minor ion. A PHREEQC check at the calibration seawater (288 K,
same alkalinity and dissolved carbon) shows the proxy itself is harmless:

| | pH | pCO₂ | calcite SI |
|---|---|---|---|
| 23.45 mM SO₄, no K⁺ | 8.347 | 302 ppm | +0.702 |
| 28.2 mM SO₄ + 10.2 mM K⁺ | 8.352 | 297 ppm | +0.695 |

Using 28.2 mM without K⁺ would force alkalinity to −7.2 meq/kg.

**Implemented (the user's decision):**
- **`chemistry.py`:**
  - K is added to the end of `elements`, so existing indices are unchanged, with `ION_CHARGE` +1.
  - `'K+'` is mapped in the stoichiometry parser and removed from `IGNORED_SPECIES`.
  - `SEAWATER_SO4 = 28.2e-3` and `SEAWATER_K = 10.2e-3` (Millero et al. 2008), and `seawater_seed()` sets both.
- **`planet.py`:** K's net flux is set to zero, like SO₄'s.
- **`calibrate_earth.py`:** the seed, the flux-anchor evaluation and the static α check all use the new values.
- **Plotting:** `diagnostics.PLOT_ELEMENTS` excludes K; `plot_results` salinity includes K, and its Earth
  reference now includes SO₄ and K (35.5 g/kg, previously 32.4).

**Consequence, and an open decision.** The untracked Br⁻ (0.84 mM), F⁻ (0.07) and Sr²⁺ (0.09) now leave
+0.73 meq/kg in the derived alkalinity. The seed's alkalinity is 3.00 rather than 2.30, and at the seed's fixed
dissolved carbon its pH is 9.16 rather than 8.35. Folding Br⁻ into Cl, which is what chlorinity means, brings
alkalinity to 2.16 (−0.14 meq/kg). **Not done yet.**

For the paper: the alkalinity definition is the explicit conservative alkalinity of Wolf-Gladrow et al.
(2007), and SO₄ = 0 on anoxic ocean worlds is supported by Archean sulfate below 2.5 µM (Crowe et al. 2014).
Old output files have one fewer element row, and resume only checks `pe` and the seeding flag, so don't mix
them with new runs.

### 37.7 Recalibration (run by the user, 2026-09-24)

This includes the shelf, K⁺/SO₄ and crust-rate changes, but not the Br⁻ fold.

| constant | before | after |
|---|---|---|
| `KD_MG_HT` | 1.969604e-02 | 2.361032e-02 |
| `K_NA_CONT_REMOVAL` | 6.099720e-03 | 5.775040e-03 |
| `ALPHA_REF` | 4.9 | 14.57 |
| `K_CL_SUBDUCTION` | 1.961786e-04 | unchanged (analytic) |

Best evaluation (`output/calib_ls_027.json`, converged at 0.91 Gyr):
- T 294.4 K, pCO₂ 697 ppm.
- Na 455.8 mM (−3 %), Ca 10.44 (+1 %), Mg 60.51 (+15 %).
- Cl 546 (seeded and stationary), SO₄ 28.2 and K⁺ 10.2 (fixed).
- Alkalinity 5.51 meq/kg (×2.4), dissolved carbon 4.82 mM (×2.3), Si 4.83 mM.
- **Net seafloor alkalinity flux 0.042 Tmol/yr.** The α residual anchors the *primary-dissolution* flux at
  1 Tmol/yr, and the pore clays remove almost all of it.

The Mg offset falls from +36–40 % (§28.1, §35.6) to +15 %.

Before the refit, the same Earth run with the new code but old constants gave T 294.4 K, pCO₂ 701 ppm, Na 433,
Mg 72.0, Ca 10.5, alkalinity 5.43 and a net flux of 0.012 Tmol/yr.

The run behind Fig. 4 (the continental baseline at S = 1: 3 km, blank start, 2 Gyr, still with Cl at the
time) gave:
- Cl 219 mM (40 % of seawater);
- Ca 0.71 mM, alkalinity 366 meq/kg, dissolved carbon 231 mM;
- surface pH 9.67.

**Fig. 4 must be drawn from the calibration run.**

### 37.8 The calibrated constants now live only in `constants.py`

- **`constants.py`:** the only definition of `KD_MG_HT`, `K_NA_CONT_REMOVAL`, `K_CL_SUBDUCTION` and
  `ALPHA_REF`.
- **`planet.py` and `weathering.py`:** import them, so `from kamino.planet import KD_MG_HT` and similar still
  work.
- **`parameter_sweep.py`:** `KD_MG_CALIB = KD_MG_HT` and `K_NA_CALIB = K_NA_CONT_REMOVAL`. This closes the
  drift in §3 and §34, and `_warn_constant_drift` is silent.
- **`calibrate_earth.py`:** starts from the shipped constants (`K_NA_START` and the others are now `None`) and
  prints paste instructions for `constants.py`.

**Stale:** the α arms now centre on 14.57 (`alpha = [ALPHA_REF, 10, 50]`, and `ALPHA_PLANE`), but their
comments still assume α increases from `ALPHA_REF` up to 50.

### 37.9 Chlorine: removed from the sweeps

**The problem.**
- The Cl relaxation time is τ_Cl = ρ_w d / (η_Cl J) = 3.90 Gyr × (d / 3 km) / R̃: 4.81 Gyr at 3.7 km and
  26.0 Gyr at 20 km.
- From a blank start at 3 km, Cl reaches 40 % of steady state by 2 Gyr and 64 % by 4 Gyr (7 % and 14 % at
  20 km). So running to 4 Gyr doesn't fix it.
- The steady state, 0.38 M × F̃/R̃, is 38 mM at the defaults but spans 0.4 mM–380 M across the grid.
- The blank sweep's median Cl of 15 mM is exactly 40 % of 38 mM.
- Physically, a planet's Cl inventory is set during accretion (Sharp & Draper 2013), not by these fluxes.

**The test.** 30 runs: 3 km, R̃ = 1, Mg/Si 1.25, ΔIW −2, pe −3, reverse weathering on, α = 14.57 and the
recalibrated η. Four cases at S = 0.6–1.1:
- the current model (Cl arriving as volcanic HCl);
- no Cl (Cl/C = 0);
- a fixed NaCl background of 55 mM;
- a fixed NaCl background of 546 mM (Na = Cl seeded, with no Cl source or sink and no Na sink).

The first two were run at 0.1× and 1× outgassing, the NaCl cases at 0.1× only.

0.1× outgassing: T (K) / pCO₂ (mbar) / seafloor pH / salinity (g/kg):

| S | current (Cl 15 mM) | no Cl | NaCl 55 mM | NaCl 546 mM |
|---|---|---|---|---|
| 0.6 | 242.6 / 5.3 / 7.90 / 2.3 | 244.1 / 8.5 / 7.97 / 2.6 | 243.9 / 8.0 / 8.01 / 5.8 | 243.1 / 6.2 / 8.05 / 34.5 |
| 0.8 | 269.7 / 5.3 / 7.90 / 2.3 | 271.7 / 8.5 / 7.97 / 2.6 | 271.4 / 8.0 / 8.01 / 5.8 | 270.3 / 6.2 / 8.05 / 34.5 |
| 0.9 | 285.4 / 5.5 / 8.11 / 2.2 | 286.8 / 7.3 / 8.28 / 2.5 | 286.6 / 7.1 / 8.29 / 5.7 | 285.9 / 6.2 / 8.24 / 34.4 |
| 1.0 | 300.6 / 2.4 / 8.48 / 2.1 | 304.9 / 5.2 / 8.45 / 2.4 | 303.6 / 4.2 / 8.53 / 5.6 | 300.7 / 2.4 / 8.60 / 34.2 |
| 1.1 | 360.8 / 150 / 6.88 / 2.3 | 364.7 / 301 / 6.90 / 2.2 | 363.4 / 231 / 7.00 / — ‡ | 359.0 / 113 / 7.17 / — ‡ |

1× outgassing, current model (Cl 153 mM) vs no Cl:

| S | current | no Cl |
|---|---|---|
| 0.6 | 311.3 / 1860 / 6.52 / 22.8 | 312.8 / 2030 / 6.76 / 26.3 |
| 0.8 | 323.3 / 873 / 6.89 / 20.7 | 326.9 / 1140 / 7.07 / 23.9 |
| 0.9 | 333.3 / 753 / 6.97 / 19.3 | 341.2 / 1730 / 6.95 / 22.5 |
| 1.0 | 346.2 / 1190 / 6.77 / 17.8 | 354.1 / 3740 / 6.69 / — ‡ |
| 1.1 | 354.4 / 1920 / 6.60 / — ‡ | 363.6 / 6640 / 6.46 / — ‡ |

‡ Stopped at the 900 s wall-clock cap. The values are the state reached, and salinity is unavailable because
no trajectory is stored.

- **Temperature:** in temperate states the Cl choice moves T by ≤ 2 K at 0.1× outgassing and ≤ 4 K at 1×.
  Near the inner edge the spread grows to 8–9 K and ×3 in pCO₂.
- **pH:** moves by ≤ 0.3 everywhere.
- **Salinity:** without Cl it is weathering-derived Mg bicarbonate (alkalinity ~32 meq/kg at 0.1×, ~300 at 1×;
  Ca < 1 mM). HCl-derived Cl *lowers* salinity by 10–15 %, because each mole of HCl removes a mole of
  alkalinity.
- **Convergence:** every current-model run timed out at 2 Gyr. The no-Cl and NaCl runs converged, except the
  capped hot cells.

**Decision (the user's): the sweeps run with no Cl.** `parameter_sweep.CL_OUTGASSING_RATIO = 0.0` is passed to
every Planet built by `run_simulation`: the ocean worlds, the continental baseline, and the land and α grids.
Run names carry a `_cl0` tag. Cl stays only in `calibrate_earth.py`. The paper will note that the Cl inventory
is a secondary effect and beyond its scope.

Verified: the default S = 0.9 run through the sweep path converged at 286.8 K with Cl ≡ 0. The test scripts and
outputs exist only in the session scratchpad.

**Open:** how far the RWR limit moves with Cl near the inner edge.

### 37.10 Corrections and smaller changes

- **§33.4 corrected:** the OLR fit is Kadoya & Tajika (2019), not Haqq-Misra et al. (2016). The comment at
  `plot_results.py:1310` still repeats the old attribution.
- **`get_ocean_state` receives `pe`,** so the surface pCO₂ and pH are now computed at the model's redox state.
  This change predates the recalibration, which includes it.
- **The user changed sweep defaults** around this session:
  - blank ocean by default (`KAMINO_SEED_OCEAN`, 2026-09-23);
  - `pe=[PE_DEFAULT]` for most sweeps;
  - an edited `DEFAULT_SWEEPS`.
- **Seeded runs still carry Cl.** With `KAMINO_SEED_OCEAN=1`, the seed puts 546 mM of Cl into sweeps, which then
  drains slowly through subduction because there is no Cl source.

### 37.11 Status: before the rerun

This mirrors `paper_issues.md` §10. Settle every item, recalibrate once, then rerun once.
- [ ] Fold Br⁻ into Cl (chlorinity) for the calibration (§37.6).
- [x] Shelf depth: 140 m (§37.23). [ ] Shelf-water temperature: still the deep-water value (§37.5).
- [x] Seafloor temperature floor: lowered to 273.15 K (§37.14).
- [x] Sediment porosity in τ_cover, with cited dust and cosmic-dust floor (§37.17).
- [ ] Element-by-element weathering law and the pore fluid open to atmospheric pCO₂: keep and justify, or
  change (§37.3).
- [x] Calibration flux target: net flux, 0.9 Teq/yr (§37.16). **Still to do:** refit (expect α ≈ 425–440),
  then recentre the α arms. The new-sedimentation pilot is done and the default outgassing is 1× (§37.18).
  The calibration now alternates, and a test converged at α = 356 (§37.20). The real refit with τ_rw is
  running (§37.21).
- [x] τ_rw: calibrated to Dunlea et al. (2017)'s 0.02 Tmol Mg/yr; τ_rw = 39.2 Myr (§37.21).
- [x] Refit: α = 348.3, η_HT 3.055e-2, η_Na 5.254e-3 in `constants.py` (§37.21). Still to do: recentre the α arms.
- [x] Convergence check: windowed drift, no `r_avg` (§37.19).
- [ ] RWR limit with and without Cl at two crust production rates (§37.9).
- [ ] Rerun into a fresh output directory, or with `KAMINO_RERUN=1`.
- [ ] Update the α comments in `parameter_sweep.py` (§37.8).
- [ ] Redraw Fig. 4 from the calibration run (§37.7).

### 37.12 Cross-cutting lessons

- **A background ion that closes the charge balance carries more than its name.** SO₄ = 23.45 mM was standing
  in for K⁺ and every minor ion; making it literal moved +0.73 meq of untracked charge into alkalinity. In a
  model where alkalinity is derived from charge, a change to any background ion must conserve the net
  background charge.
- **A term that only exists on land planets can still dominate the calibration.** The shelf sink buried 56 % of
  Earth's carbonate because it was a second whole-ocean sink, and the calibrated constants had absorbed it.
- **A quantity whose relaxation time exceeds the planet's age across most of the grid is an input, not a
  result.** For Cl, τ ∝ d/R̃. Measure whether it matters before adding a parameter to sweep: here it didn't,
  for temperate climates.
- **Check which code a stored run used by replaying a diagnostic,** rather than reasoning about process
  lifetimes. The stored Damköhler number settled which crust rate the job had used.
- **Check citations of the code against the code.** The paper cited the OLR fit correctly; this history had it
  wrong since §33.4.
- **State the physical basis of a weighting.** Weight fraction for `k_p` is defensible (equal surface area per
  gram); mole fraction would depend on how the formula unit is written.

### 37.13 References added this session

- **Beckingham et al. (2016)**, *GCA* **188**, 310 — reactive surface area estimates (volume-fraction,
  BET/mass-specific, image-based).
- **Harris et al. (2014)**, *Mar. Geol.* **352**, 4 — continental shelves ~27 million km² (~7 % of the ocean),
  mean shelf-break depth 140 m.
- **Cogné & Humler (2004)**, *EPSL* **227**, 427 — seafloor generation half-rate 1.30 ± 0.28 km²/yr.
- **Cogné et al. (2006)**, *G3* **7**, Q03012 — trends and rhythms in seafloor generation; cited in the draft
  for 1/130 Myr.
- **Gerlach (2011)**, *Eos* **92**, 201 — volcanic CO₂ 0.13–0.44 Gt/yr, preferred 0.15–0.26.
- **Millero et al. (2008)**, *DSR I* **55**, 50 — reference composition of seawater (K⁺ 10.2, SO₄ 28.2
  mmol/kg).
- **Wolf-Gladrow et al. (2007)**, *Mar. Chem.* **106**, 287 — explicit conservative total alkalinity.
- **Crowe et al. (2014)**, *Science* **346**, 735 — Archean seawater sulfate < 2.5 µM.
- **Sharp & Draper (2013)**, *EPSL* **369**, 71 — Earth's halogen inventory set during accretion.
- **Hamilton (1976)**, *J. Sediment. Petrol.* **46**, 280 — porosity and density of deep-sea sediments with
  depth.
- **Aagaard & Helgeson (1982)**, *Am. J. Sci.* **282**, 237; **Lasaga (1984)**, *JGR* **89**, 4009 —
  transition-state rate laws with per-mineral affinity.
- **Kadoya & Tajika (2019)**, *ApJ* **875**, 7 — the OLR fit used in `climate/analytic.py`.
- **Jickells et al. (2005)**, *Science* **308**, 67 — aeolian dust deposition to the oceans, ~450 Tg/yr (§37.17).
- **Love & Brownlee (1993)**, *Science* **262**, 550 — cosmic dust accretion, (4 ± 2) × 10⁷ kg/yr (§37.17).
- **Coogan & Dosso (2022)**, *GCA* **329**, 22 — net low-temperature seafloor alkalinity flux, 0.90 Teq/yr
  (1σ 0.16–1.64); now the α anchor (§37.16; see also §22.4).

### 37.14 Seafloor temperature floor: 274 → 273.15 K (later the same day)

The 274 K floor was chosen to keep water below freezing out of PHREEQC. It isn't needed for that: every
PHREEQC call already clamps the temperature to 0.01 °C (`chemistry.py`, `_solution_block`). So the floor's only
role is physical, setting the deep-water temperature that ocean precipitation sees. Weathering runs at the floor
+ 9 K, so it is always at 282 K or above.

A composition-dependent freezing point was rejected as too complex for a second-order effect. Seawater freezes at
−1.9 °C, and the model's fresher oceans (2–26 g/kg) at roughly −0.1 to −1.4 °C.

**The floor is now `SEAFLOOR_T_FLOOR = 273.15` K** in `constants.py`: the freshwater freezing point and the lower
limit of the thermodynamic data, within ~1.5 K of every model ocean's freezing point. It is used everywhere the
floor appeared:
- `planet.py` (the physics);
- `calibrate_earth.py` (the flux-anchor evaluation);
- `plot_results.py` (the "at floor" marker and its legend label);
- `probe_saturation.py`.

**Effect.** Only planets whose surface is below (273.15 + 16.7)/1.02 = 284.2 K are affected. Pore rates fall by
~9 % for the 0.85 K drop (~10 % per K at 60–80 kJ/mol). The Earth calibration is unaffected: its seafloor sits at
283.6 K, above both floors. Checked: an S = 0.6 planet at a 243.8 K surface gets a 273.15 K seafloor and a
282.15 K pore.

**A related limitation, for the paper.** On planets with a frozen surface, the surface pCO₂ is computed from ocean
water clamped to 0.01 °C, because the model has no ice.

### 37.15 Initial conditions don't change the result

32 runs with the current code at sweep settings (no Cl, 273.15 K floor, α = 14.57, 3 km, Mg/Si 1.25, ΔIW −2,
0.1× outgassing): S = 0.4–1.1 at R̃ = 1, and S = 0.8 and 1.0 at R̃ = 0.1, each from four starts:
- blank ocean with pCO₂ 1000 Pa (the default);
- blank with pCO₂ 10 Pa;
- blank with pCO₂ 1 bar;
- a concentrated start: 100 mM Mg(HCO₃)₂ + 10 mM Ca(HCO₃)₂ + 1 mM Si, charge-consistent.

**Every configuration reaches the same final state from all four starts,** within 0.3 K and ~2 % in pCO₂, with
the ocean chemistry within ~1 %. So there is a single steady state, approached from below and from above. The
calcite bistability seen in the calibration came from the old Cl charge problem and doesn't appear in Cl-free
ocean worlds.

- **Only the convergence time differs.** The concentrated start takes up to 1.6 Gyr at R̃ = 1. At R̃ = 0.1 it is
  the only start that meets the convergence criterion within 2 Gyr, but the blank starts end within 0.3 K of it.
  So a 2 Gyr timeout is effectively steady state there.
- **Initial pCO₂ is forgotten within ~10 kyr,** as the atmosphere relaxes to the ocean.
- **S = 0.4–0.8 give identical chemistry and pCO₂ (10.3 mbar);** only the surface T differs. All three sit on the
  seafloor temperature floor.
- **The existing blank-era sweep never hit the cold wall.** Its 186 hot-wall stops before 1 Myr are all at
  S ≥ 1.15, which is past the runaway even at the 1 Pa CO₂ floor. So the domain walls don't depend on the start.

### 37.16 The α flux anchor: from primary to net flux

**The defect.** `calibrate_earth.seafloor_alk_flux_tmol`, which α was fitted against, called
`get_weathering_flux(..., precipitating_minerals=[])` without `pe`. It therefore measured the *primary* flux,
before the pore kaolinite and goethite, at PHREEQC's default pe = +4. It also used a textbook seawater
composition (Si 0.1 mM) rather than the run's own ocean. Decomposed at the recalibrated Earth
(`calib_ls_027`, α = 14.57), in Teq/yr:

| step | alkalinity | Fe | Mg | Al |
|---|---|---|---|---|
| A. anchor as coded (primary, pe +4, textbook seawater) | 1.177 | +1.098 | +0.060 | +0.019 |
| B. A at pe −3 | 0.113 | +0.007 | +0.058 | +0.048 |
| C. B with pore kaolinite and goethite | 0.059 | 0 | +0.058 | +0.001 |
| D. C with the run's actual ocean | 0.042 | 0 | +0.042 | 0 |
| E. what `dY_dt` applies | 0.042 | | | |

So 93 % of the fitted "1 Tmol/yr" was iron that only dissolves in an oxidising pore fluid. The model's actual net
flux, 0.042 Teq/yr, is all Mg with no Ca. It is ~20× below Coogan & Dosso's (2022) net 0.90 Teq/yr (1σ range
0.16–1.64), which is mostly Ca release (§22.4).

**The fix, now in `calibrate_earth.py`:**
- `run_planet` returns the run's own `diagnostics["alk_flux"]`, the net flux `dY_dt` applied at the final state.
- `seafloor_alk_flux_tmol` now simply reads it, so the pore clays, the redox state and the ocean composition are
  all the model's.
- `FLUX_TARGET = FLUX_TARGET_NET = 0.9` Teq/yr (Coogan & Dosso 2022).
- The start-up α diagnostic now uses the pore clays and `pe`.
- The old primary-only function is deleted.

**Earth α scan** (calibration setup, current η, old sedimentation):

| α | net flux (Teq/yr) | Da | T (K) | pCO₂ (ppm) | Mg (mM) | Ca (mM) |
|---|---|---|---|---|---|---|
| 14.57 | 0.042 | 0.004 | 294.4 | 697 | 60.5 | 10.4 |
| 150 | 0.377 | 0.038 | 294.2 | 671 | 64.2 | 11.7 |
| 300 | 0.670 | 0.076 | 294.1 | 648 | 67.4 | 12.8 |
| 450 | 0.914 | 0.112 | 293.9 | 629 | 70.0 | 13.7 |

- **0.9 Teq/yr needs α ≈ 440.** The flux grows more slowly than α.
- **Earth's climate barely moves,** because continental weathering dominates its carbon budget.
- **Mg drifts up**, which the joint refit has to re-balance.

**The transfer problem, quantified.** At the *same* α = 14.57, the ocean worlds already deliver 1.4–2.0 Teq/yr of
net flux (normalised to Earth's seafloor area), against Earth's 0.042. The steady state needs 2 × the outgassed
carbon, i.e. 1.5 Teq/yr at 0.1× outgassing. That's 30–50× more per unit α, for two reasons:
- a ~7× larger reactive fraction of the crust, because Earth's crust is buried by terrigenous sediment (§37.17);
- no element-by-element suppression of Ca and Mg release in dilute oceans (§37.3).

**Ocean-world pilot at α = 440** (old sedimentation, since it ran before §37.17), against α = 14.57, 3 km, R̃ = 1,
no Cl. Each entry is T (K) / pCO₂ (mbar) / Da:

| S | 0.1×, α 14.57 | 0.1×, α 440 | 1×, α 14.57 | 1×, α 440 |
|---|---|---|---|---|
| 0.4 | 216.4 / 10.3 / 0.026 | 211.1 / 0.38 / 4.8 | 268.6 / 2760 (CO₂ ceiling) | 221.9 / 59 / 0.08 |
| 0.6 | 244.7 / 10.3 / 0.026 | 237.3 / 0.38 / 4.8 | 312.8 / 2030 / 0.05 | 253.5 / 60 / 0.08 |
| 0.8 | 272.6 / 10.3 / 0.026 | 261.9 / 0.38 / 4.8 | 326.9 / 1140 / 0.35 | 284.8 / 57 / 0.09 |
| 0.9 | 286.7 / 7.3 / 0.049 | 275.4 / 0.38 / 4.8 | 341.2 / 1730 / 2.4 | 334.0 / 822 / 51 † |
| 1.0 | 304.9 / 5.2 / 1.10 | 297.9 / 1.4 / 20 | 354.1 / 3740 / 8.8 † | 352.6 / 3130 / 259 † |
| 1.1 | 364.7 / 301 / 299 | 364.6 / 296 / 9030 | 363.6 / 6640 / 20 † | 362.7 / 6000 / 595 † |

† stopped at the 900 s wall-clock cap.

- **At 0.1× outgassing, α = 440 leaves no negative-feedback band.** Every S ≤ 0.9 sits on the seafloor floor and is
  thermodynamically limited (Da ≈ 5), and S = 1.0 is already on the positive branch.
- **At 1×, S = 0.4–0.8 are temperate or cold and kinetically limited (Da ≈ 0.08).** The RWR transition falls
  between S = 0.8 and 0.9, at the same instellation as at α = 14.57, but 40–60 K cooler below it.
- **So the default outgassing would move towards 1×,** and the window with a negative feedback (kinetic *and*
  above the floor) is narrow.
- **This needs re-running with the §37.17 sedimentation and the refitted constants** before any conclusion.

### 37.17 Sedimentation: cited dust, a cosmic-dust floor and porosity

**The problem.** The terrigenous rate (5 m/Myr at 30 % land) and the 0.3 m/Myr floor had no source, and there was no
porosity. `h_cover` = 100 m is a bulk thickness, but the precipitates were counted as solid mineral volume. On
Earth the uncited 5 m/Myr was most of the sedimentation (plus 1.39 m/Myr of precipitates), so θ/α = 0.12. Ocean
worlds had only 0.24–0.34 m/Myr, at or near the floor, so θ/α ≈ 0.82. That made the terrigenous number a hidden ~7×
multiplier on ocean-world weathering relative to the calibrated Earth.

**Implemented** (constants in `constants.py`):
- `EARTH_DUST_FLUX_TO_OCEAN = 450 Tg/yr` (Jickells et al. 2005): dust reaching ridge flanks, scaled by land area
  relative to Earth's and spread over the planet's seafloor. That is 0.476 m/Myr of solid at Earth, and zero on
  land-free planets.
- `COSMIC_DUST_FLUX_PER_AREA` from 4×10⁷ kg/yr (Love & Brownlee 1993), about 1×10⁻⁴ m/Myr bulk. It replaces the
  0.3 m/Myr floor. A land-free planet has no continental dust, and volcanic ash remains unmodelled.
- `SEDIMENT_GRAIN_DENSITY = 2650` kg/m³.
- `SEDIMENT_POROSITY = 0.7` (Hamilton 1976; ~0.6–0.8 in the upper 100 m).
- `planet.py` computes the bulk rate as (precipitates + dust solid volume) / (1 − φ).
- `seafloor_reactive_area` uses the cosmic-dust floor, and falls back to Earth's dust alone (1.59 m/Myr bulk)
  only when no rate is passed.
- The old `_S_TERR_EARTH` / `_s_terr` are removed.

**Effect:**

| | θ/α before | θ/α after |
|---|---|---|
| Earth calibration | 0.120 | 0.124 |
| ocean world S = 0.6 / 0.9 / 1.0 | 0.827 / 0.808 / 0.828 | 0.557 / 0.523 / 0.624 |

- **Earth barely moves,** because the cited dust plus porosity happens to reproduce the old total.
- **Ocean worlds lose ~30–35 % of their reactive area** relative to the calibrated Earth.
- **With no sediment at all, θ/α → 1**: the crust is never buried.

**Earth α scan with the new sedimentation** (calibration setup, current η):

| α | net flux (Teq/yr) | Da | T (K) | pCO₂ (ppm) |
|---|---|---|---|---|
| 300 | 0.663 | 0.075 | 294.1 | 649 |
| 400 | 0.827 | 0.099 | 294.0 | 636 |
| 500 | 0.976 | 0.123 | 293.9 | 625 |

So 0.9 Teq/yr now needs **α ≈ 450** (was ≈ 440). This is a scan, not the refit, which is still pending.

### 37.18 The new-sedimentation pilot, and the default outgassing (1×)

> **Superseded by §37.22** for the numbers: that section reruns this pilot on the refitted constants. The
> choice of 1× as the default stands.

The user changed `parameter_sweep.outgassing_default` to `[1]`. Tested at α = 450, 3 km, R̃ = 1, no Cl, blank
ocean, with the §37.17 sedimentation. The numbers below are from the rerun under §37.19's convergence check. In the
first run, 6 of the 24 hit the 900 s wall-clock cap and 5 ran to 2 Gyr. The rerun matches every state that both
runs reached (§37.19).

Seafloor T = 1.02 T_s − 16.7 K, so the 273.15 K floor binds whenever T_s < 284.2 K.

| S | 1× outgassing: T_s (K) / pCO₂ (bar) / Da | 0.3× | 3× |
|---|---|---|---|
| 0.4 | 244.0 / 0.73 / 0.025 (floor) | 214 / 0.0038 / 0.35 (floor) | CO₂ ceiling at 43 Myr |
| 0.6 | 286.9 / 0.58 / 0.028 | 242 / 0.0038 / 0.35 (floor) | 322.1 / 4.18 / 0.30 |
| 0.8 | 297.3 / 0.155 / 0.13 | 269 / 0.0038 / 0.35 (floor) | CO₂ ceiling |
| 0.85 | 304.7 / 0.14 / 0.41 | | |
| 0.9 | 334.2 / 0.85 / 21 | 284 / 0.0039 / 0.34 (floor) | CO₂ ceiling |
| 0.95 | 346.0 / 2.01 / 59 | | |
| 1.0 | 352.7 / 3.15 / 103 | 339.7 / 0.28 / 229 | CO₂ ceiling |
| 1.1 | 362.7 / 6.00 / 231 | 358.9 / 1.03 / 1128 | CO₂ ceiling |

Crust production at 1× outgassing: 0.1× hits the CO₂ ceiling at S 0.8 and 0.9. 10× sits on the floor at
pCO₂ ≈ 1 mbar (S 0.8: 264 K; S 0.9: 279 K).

- **1× is the only one of the three with a negative-feedback band above the floor.** Over S 0.6–0.85, pCO₂
  falls from 0.58 to 0.14 bar as S rises, with Da < 1 and T_s 287–305 K. So 1× is the right default. It is
  also Earth's outgassing, the value used in the calibration.
- **The transition is between S 0.85 and 0.9.** Da jumps from 0.4 to 21 and pCO₂ rises sixfold. With the old
  sedimentation (§37.16) it was between S 0.8 and 0.9, so it has barely moved.
- **0.3× is floor-bound.** Every S ≤ 0.9 has the same pCO₂, because the seafloor never leaves 273.15 K, and
  S = 1.0 jumps straight to the hot branch. This matches 0.1× in §37.16.
- **3× has no weathering steady state, except a hot one.** The kinetic flux tops out below twice the
  outgassed carbon, so pCO₂ runs to the maximum greenhouse at every S except 0.6. S 0.6 settles hot, at
  322 K and 4.2 bar.
- ⚠️ **Past the transition the planet settles hot rather than running away.** At 1×, S 0.9–1.1 are steady at
  334–363 K and 0.85–6 bar, below the maximum greenhouse. Before §37.19 these runs hit the wall-clock cap or
  ran to 2 Gyr, which made them look like runaways. The paper's "retrograde weathering runaway" wording needs
  to match: in this part of parameter space it is a jump to a hot branch.

### 37.19 Convergence check: windowed drift replaces the smoothed rate

**The old check.** `r_avg` was an extra ODE state that relaxed towards max|F_net|/b with τ = 30 Myr, and
`event_converged` fired when it fell below 0.05/Gyr. It dates from commit 608c8df. It replaced an event that
called `dY_dt` itself, which was costly and impure, since solve_ivp events must be pure functions of (t, y).
Its problems:
- **It averaged the magnitude of noisy rates.** PHREEQC jitter in Si, Na and Ca kept it above threshold long
  after the state had settled. Replaying 50 saved runs: one was steady by ~470 Myr but was held until 1883 Myr.
- **It throttled the integrator.** `r_avg` sat in the error norm with a tight atol and a noisy, non-smooth
  right-hand side. Runs that end at the same wall now take 2–50 % fewer steps.
- **It cost a Jacobian column** (2 `dY_dt` calls per Jacobian).
- **Its start value set a minimum run length.** Decaying from 1/Myr to 0.05/Gyr takes ~300 Myr.

**The new check (`planet.time_evolve`).**
- **The rule:** a run is converged once no species above 1e-6 mol/kgw (1e-7 until §37.22) has changed by more than
  `convergence_threshold` (0.05/Gyr) × `convergence_window` (50 Myr) = 0.25 % over the last window. The
  denominator is floored at 1e-6, as before.
- **Why a window:** measuring net displacement over a baseline cancels step-to-step noise instead of averaging
  its magnitude.
- **The loop:** `solve_ivp` is replaced by a loop over `scipy.integrate.LSODA.step()`, the same integrator,
  tolerances and `max_step`, so the check can see the history.
- **The domain event:** checked after each accepted step. The crossing is located with `brentq` on the
  step's dense output, copying solve_ivp's `solve_event_equation`.
- **Clean-up:** `r_avg`, `tau_r_avg`, `event_converged` and `min_time` are gone. The state is now
  `[P_CO2, P_H2O, *elements]`.
- **Jacobian:** `macro_jacobian` also skips the pinned SO₄ and K columns. This is exact: their rows are zero,
  so their Newton updates are zero. It closes §10's "dead Jacobian columns" item. A Jacobian now costs 20
  `dY_dt` calls, down from 26.
- **Readers:** `diagnostics.diagnose` truncates Y to `2 + len(elements)`, so old outputs with the `r_avg` row
  still load. `plot_results` reads `len(elements)` ions.
- ⚠️ **Old and new files can't be told apart by row count.** New files have 13 rows
  (`[P_CO2, P_H2O, 11 elements]`), and so do the pre-K sweep files (10 elements + `r_avg`). In a pre-K file the
  13th row is `r_avg`, ~1e-17 s⁻¹, which would be read as K ≈ 0. Tell them apart by the `_cl0` tag or the date.

**Test: A/B against the old check,** same code otherwise.

| | old check | new check |
|---|---|---|
| ocean pilot (§37.18), total wall time, 24 runs | 7805 s | 505 s |
| runs stopped by the 900 s cap / reaching 2 Gyr | 6 / 5 | 0 / 0 |
| Earth α = 300 / 400 / 500, wall time per run | 15 / 15 / 19 s | 10 / 10 / 10 s |

- **Converged states:** ΔT ≤ 0.04 K, ΔpCO₂ ≤ 0.4 %, major species ≤ 0.8 %. The worst is Na in the slow 1×
  S 0.4 run.
- **Domain-wall runs:** the same wall at the same time, to 1 Myr, with the same state.
- **Earth:** pCO₂ within 0.002 % and species within 0.1 %. The α–flux relation is unchanged.
- **Formerly capped runs:** each converged state matches the capped run's last state. For example, 1× S 1.0
  gives 3.15 bar at 442 Myr, against 3.13 bar at 741 Myr before.

⚠️ **`diagnostics.diagnose` is still broken, independently.** It reads `planet.f_bio` and a
`crust_composition` config key, and neither exists any more. The file carries its own "STALE" note.

### 37.20 The first net-flux refit failed; the calibration now alternates

**What the user's refit returned:** α = 55.43, η_HT = 5.806×10⁻², η_Na = 4.831×10⁻³. It was not a converged
answer. `calibrate()` keeps the lowest-cost evaluation, and that was `calib_ls_045`, a finite-difference
probe, not a trust-region iterate. The iterates ended at α ≈ 29.4. That best point gave:

| | target | calib_ls_045 |
|---|---|---|
| net seafloor flux | 0.9 Teq/yr | 0.40 |
| Mg | 52.8 mM | 26.3 |
| Ca | 10.3 mM | 7.5 |
| Na | 469 mM | 531 |

**Cause 1: the finite-difference steps were far too large.** `diff_step=0.2` on ln(x) means a step of
0.2 × |ln x|, and scipy scales `diff_step` by max(1, |x|). That gave probes of ×0.36 in K_na, ×0.47 in η_HT
and ×1.7 in α. Every K_na and η_HT probe landed on the Ca-collapsed branch (Ca 0.70 mM, Na 1400–1700 mM), so
two of the three Jacobian columns measured a jump between branches, not a slope.

**Cause 2: the flux residual dominated the cost.** ln(0.05/0.9) ≈ −2.9. The fitter found that raising η_HT
lowers Mg, which weakens the element-by-element Mg suppression and raises the flux. So it traded Mg (down to
24 mM) for flux while α crept up.

**The fix (`calibrate_earth.py`).**
- **Inner step:** `least_squares` on (ln K_na, ln η_HT) against (Na, Ca, Mg) at fixed α, with
  `diff_step` ≈ 0.05 in ln(p).
- **Outer step:** α ← α × 0.9/flux. The flux is ~linear in α (exponent 0.99), and α moves the ions by ~1 % over
  a ×1.7 step. The docstring's objection to separating the fit (§36: "the ocean moves the flux") is handled by
  iterating.
- **Probe direction:** a 5 % *drop* in K_na still tips Earth onto the collapsed branch (Na 475 → 551 mM,
  Ca 10 → 0.8 mM). So the parameters are shifted, x = ln p + 20, which makes scipy's forward probes step
  *up*, away from the collapse. Earth's calibrated state sits within ~5 % of K_na of the calcite tipping point.
- **Round cap:** `MAX_RUNS_PER_ROUND` caps the solver steps per round. scipy's `max_nfev` doesn't count
  Jacobian probes, so a round takes ~20–50 runs.
- **`ALPHA_PINNED`:** still available, to fit the ions at a fixed α.

**Scratch test** (before the probe-direction fix). It converged in 4 rounds and 78 runs (~13 min). α went
55 → 259 → 340 → 356; at α = 55 the ions fitted to cost 0.0015 but the flux was 0.19 Teq/yr. The final state:

| Quantity | Value |
|---|---|
| α | 356.4 |
| η_HT | 2.864×10⁻² |
| η_Na | 5.390×10⁻³ |
| net seafloor flux | 0.900 Teq/yr |
| Na / Ca / Mg | 461 mM (−1.7 %) / 10.7 mM (+4.2 %) / 57.6 mM (+9.0 %) |
| Alk / C | ×2.3 / ×2.2 (the usual abiotic offset) |
| T / pCO₂ | 293.9 K / 631 ppm |

Mg was still improving when the round cap hit; every K_na probe was still landing on the collapsed branch.
`constants.py` was set to these values as the starting point for the real refit (§37.21).

**Side effect.** The first scratch attempt wrote `output/calib_ls_000.json` into the repo, because
`Planet.output_path` isn't set from `OUTPUT_DIR`. So the user's original eval 000 was overwritten.

### 37.21 τ_rw calibrated to the modern authigenic-clay Mg sink

**The gap.** τ_rw = 5 Myr had no source (the draft's own comment at L432 calls it a guess). Sepiolite(d) stays
~4 log units supersaturated, so τ_rw *is* the reverse-weathering flux (§26.3). It is not a relaxation
constant, and it matters on hot worlds: 5 → 33 Myr cooled a 20 km world by 25 K (§26.2).

**The constraint.** Dunlea et al. (2017, Nat. Commun. 8, 844) quantified authigenic clay in South Pacific
Gyre sediment:
- typical deep-sea sediment removes ~**0.02 Tmol Mg/yr**;
- if Si-rich (chert-forming) sedimentation covered 50–100 % of the seafloor, it would remove
  **0.4–0.8 Tmol Mg/yr**.

Deltaic clay formation (Michalopoulos & Aller 1995, Science 270, 614) is a river-fed process that ocean worlds
lack. At Earth, with τ_rw = 5 Myr and the current code, reverse weathering removes ~0.16 Tmol Mg/yr, so
τ_rw ≈ 40 Myr is expected.

**Implemented:**
- **`constants.py`:** `TAU_RW_REF` moved here from `planet.py`, into the calibrated block, which now has five
  constants to update together.
- **`planet._final_diagnostics`:** records `rw_mg_flux`, the Tmol Mg/yr removed by reverse-weathering clays
  (positive means a sink).
- **`calibrate_earth.py`:**
  - `RW_MG_TARGET = 0.02` (Dunlea et al. 2017).
  - The outer loop also updates τ_rw ← τ_rw × rw/0.02, since flux ~ excess/τ_rw.
  - It stops when the flux is within 3 % *and* the Mg sink is within 5 %.
  - `RW_MG_TARGET = None` holds τ_rw fixed. Use it for the Si-rich sensitivity case (τ_rw ≈ 1–2 Myr).
- **Justification for the paper:** τ_rw is a calibrated *effective* rate, not a mineral precipitation rate.

**Result** (`output/calib_ls_000`–`019`, 2026-09-24). It converged in 2 rounds and 20 runs. None of the
finite-difference probes landed on the collapsed branch, so the upward-probe fix works.

| constant | before | fitted |
|---|---|---|
| `ALPHA_REF` | 356.4 (test) | **348.3** |
| `KD_MG_HT` | 2.864×10⁻² | **3.055×10⁻²** |
| `K_NA_CONT_REMOVAL` | 5.390×10⁻³ | **5.254×10⁻³** |
| `TAU_RW_REF` | 5 Myr | **39.2 Myr** |
| `K_CL_SUBDUCTION` | 1.962×10⁻⁴ | unchanged (analytic) |

The fitted Earth:

| Quantity | Value | Target |
|---|---|---|
| Na | 464 mM (−1.0 %) | 469 |
| Ca | 10.5 mM (+1.8 %) | 10.3 |
| Mg | 56.2 mM (+6.4 %) | 52.8 |
| Alk | 5.35 mM (×2.3) | 2.3 |
| C | 4.65 mM (×2.2) | 2.1 |
| net seafloor flux | 0.888 Teq/yr | 0.9 |
| reverse-weathering Mg sink | 0.0200 Tmol/yr | 0.02 |
| T | 293.8 K | |
| pCO₂ | 611 ppm | |

- **τ_rw = 39.2 Myr** is 7.8× slower than before. The Earth sink was 0.157 Tmol/yr at 5 Myr, so one rescaling
  landed on target; flux ∝ 1/τ_rw holds.
- **Earth's climate hardly moves** (T 293.9 → 293.8 K), as expected, because reverse weathering is under 1 % of its
  alkalinity sink.
- **The effect is on hot ocean worlds.** There reverse weathering is a CO₂ source, now 7.8× weaker (§36.2). The
  §37.18 pilot numbers above the RWR transition are therefore stale.
- **Pasted into `constants.py`** (all five).

### 37.22 The ocean pilot on the refitted constants, and the convergence cut-off

This reruns the §37.18 pilot, 24 runs: no Cl, blank ocean, 3 km, with `ALPHA_REF` = 348.3 and `TAU_RW_REF` =
39.2 Myr taken from `constants.py` (§37.21). Total wall time was 370 s. Seafloor T = 1.02 T_s − 16.7 K, with a
floor at 273.15 K.

| S | 1×: T_s (K) / pCO₂ (bar) / Da | vs §37.18 (α 450, τ_rw 5 Myr) |
|---|---|---|
| 0.4 | 247.2 / 0.83 / 0.023 (floor) | +3.1 K |
| 0.6 | 287.2 / 0.59 / 0.026 | +0.3 K |
| 0.8 | 294.5 / 0.128 / 0.090 | −2.7 K |
| 0.85 | 299.8 / 0.098 / 0.22 | −4.9 K |
| 0.9 | **310.1 / 0.106 / 1.45** | **−24.1 K** (was on the hot branch) |
| 0.95 | 341.2 / 1.10 / 45 | −4.8 K |
| 1.0 | 348.5 / 1.85 / 86 | −4.1 K |
| 1.1 | 358.5 / 3.62 / 194 | −4.2 K |

- **The RWR transition moves from S 0.85–0.9 to S 0.9–0.95.** Reverse weathering, which is a CO₂ source on hot
  worlds (§36.2), is 7.8× weaker. S = 0.9 now stays on the cool branch at Da ≈ 1.5, right at the
  kinetic/thermodynamic boundary.
- **The negative-feedback band is S ≈ 0.6–0.9** (T_s 287–310 K). pCO₂ flattens to its minimum (~0.1 bar) around
  S 0.85–0.9, just before the jump.
- **The hot branch is 4–5 K cooler,** at ~40 % lower pCO₂. It is still steady, not a runaway.
- **The cold end is 3 K warmer.** α fell from 450 to 348, so there is less weathering and more CO₂ at the floor.

Other arms:
- **0.3×:** still floor-bound for S ≤ 0.9, at pCO₂ 2.2 mbar (was 3.8). S 1.0 jumps to 334 K, 0.15 bar.
- **3×:** S 0.8 now settles hot (347 K, 7.5 bar) instead of hitting the CO₂ ceiling. S 0.6 settles at 319 K,
  3.2 bar. S 0.4 and 0.9–1.1 still hit the ceiling.
- **Crust 0.1×:** CO₂ ceiling at S 0.8 and 0.9. **Crust 10×:** floor-bound, ~0.85 mbar.

**Convergence cut-off: 1e-7 → 1e-6 mol/kgw.** The 10× crust S 0.9 run reached 2 Gyr without converging.
Trace Na at 0.26–0.6 µM flickered, and with the denominator floored at 1e-6 that read as ~5.7/Gyr of drift.
Species below 1 µM are now ignored in the drift test. Rerunning the whole pilot, 20 runs are bit-identical.
Four stop earlier: 2000 → 95, 433 → 333, 341 → 321 and 96 → 81 Myr. Their states differ by ≤ 0.12 % in major
species and ≤ 0.03 K. Wall time was unchanged (369 s).

Not yet explained: 1× S 0.4 and 0.6, and 0.3× S 1.1, converge only at 1.5–1.9 Gyr. They are cheap (≤ 20 s), so
this was left alone.

### 37.23 Shelf depth: 1000 → 140 m

`EARTH_SHELF_DEPTH = 140.0` m in `constants.py`, the mean depth of the shelf break (Harris et al. 2014, Mar.
Geol. 352, 4), replaces the uncited 1000 m in `planet.py`. It sets only the pressure at which shelf carbonate
precipitates.

**Earth A/B at the calibration setup** (land 0.3, 3.7 km, seeded ocean, current constants):

| shelf depth | T (K) | pCO₂ (ppm) | net seafloor flux (Teq/yr) | RW Mg sink (Tmol/yr) | Na / Ca / Mg (mM) |
|---|---|---|---|---|---|
| 1000 m | 293.79 | 610.7 | 0.8884 | 0.0200 | 464.2 / 10.490 / 56.19 |
| 140 m | 293.79 | 610.5 | 0.8901 | 0.0199 | 464.1 / 10.493 / 56.20 |

Every quantity moves by < 0.3 %, so the §37.21 calibration stands without a refit. Ocean worlds have no shelf, so
they are unaffected. Still open: the shelf uses the deep-water temperature, not the warmer shelf water.

### 37.24 The full rerun (24–25 Sep), and why its first figures zigzagged

The rerun is in `/data/pt426/sweep_output`: 5,491 ocean-world runs, all at pe −3, with the §37.21 constants, no
Cl, a blank ocean and the §37.19 convergence check. The old set was moved to `/data/pt426/sweep_output_seeded`.
The first figures had jagged lines, repeated Da = 1 circles and pH zigzags, all in the crust = 1 panels.

**Cause 1, plotting: the Cl sweep leaked into every figure.** `plot_results.load_data` didn't record the Cl
ratio or the run length, and `_sweep_mask` didn't filter on them. The 133 `basic_cl` runs (Cl 0.02, 4.5 Gyr)
share every other axis with `basic` at crust = 1, so each line alternated between the two runs at each S. This
affected the basic 1× panels, the 3 km depth line, the Mg/Si cross-section and the composition figures.
**Fix:**
- `load_data` records `cl_ratio` and `t_end_gyr`: `t_end_yr` from the output when present, else the
  `_tend` name tag, else 2 Gyr.
- `_ref_setup` pins both to their most-run values in `_sweep_mask` (`setup=False` keeps all of them), in the
  same way `_ref_chem` pins α.
- `planet.time_evolve` now writes `t_end_yr` to the output.

**Cause 2, plotting: three smaller selection bugs.**
- **Duplicate runs.** The α × outgassing plane wrote `out_1.0` and `out_10.0` beside `basic`'s `out_1` and
  `out_10`: 14 identical duplicate runs. `load_data` now drops duplicate configurations, and `OUTGASSING_PLANE`
  uses ints.
- **The α sensitivity figure plotted 0.01× outgassing.** `_best_operating_point` took the first of several
  tied (outgassing, crust) pairs. There every planet is CO₂-starved, so all α values sat on one line. Ties now
  go to the pair nearest (1×, 1×).
- **The depth figure dropped 3 km.** The continental baseline's 3.7 km ocean arm displaced it. The figure now
  uses `DEPTHS_SHOWN = (300, 1000, 3000, 20000, 50000)` when those depths are available.

**Cause 3, model: runs that aren't at steady state are drawn as if they were.**
- **39 `wall_timeout` runs.** One is the S = 0.55 dip on the 0.1× crust, 0.3× outgassing line, which stopped
  at 1186 Myr. `DA_TRUSTWORTHY` includes `wall_timeout`, so they're drawn as line points.
- **2 Gyr timeouts still moving.** 343 of the 583 2 Gyr timeouts have pCO₂ still changing by > 10 % over the
  last 500 Myr. By crust production: 292 of 400 at 0.01–0.03×, 9 of 84 at 0.1×, 12 of 36 at ≥ 1×. By depth:
  30 of 42 at ≥ 20 km. Relaxation time scales as depth / crust production, so these are 2 Gyr snapshots. The
  draft's "runs that reach 2 Gyr are treated as converged" doesn't hold for them.
- **Deep oceans oscillate.**
  - 20 km at S 0.55: pCO₂ 0.49–0.99 bar, period ~300–400 Myr.
  - 50 km at S 0.7: Ca 33–68 mM, period ~200 Myr, with a 2.3 bar pCO₂ spike at 1.3 Gyr.
  - Both reproduce at rtol 1e-5 with max_step 2 Myr (100× and 10× tighter), with no chemistry fallbacks, so
    this is model dynamics, not numerics.
  - A scan finds 31 oscillating runs above the 1 Pa climate floor, mostly at ≥ 20 km. Many more flicker at
    pCO₂ ≈ 0 in the CO₂-starved corner, which is climatically irrelevant.
  - The old seeded set shows both behaviours too: every deep run was unsettled at 2 Gyr, and some oscillated
    (20 km at S 0.55–0.6). So neither is a regression.

**Cause 4, model: 11 `solver_failure` runs.** All are at S 1.1, 1× outgassing, Mg/Si 1.75–2.0. PHREEQC's
weathering step fails ("Maximum iterations exceeded") from the blank, hot start at t = 0, so every derivative
is the outgassing-only fallback. LSODA takes a 20 Myr step on that constant derivative, then fails Newton
convergence at 40 Myr. It is a narrow corner whose neighbours are on the hot wall anyway.

**Also fixed: a data hazard in `rerun_wall_timeouts.py`.** It rebuilt names without the Cl ratio or run length,
and only checked that *some* file had the rebuilt name. So the one `basic_cl` wall-timeout would have rerun and
overwritten the converged `basic` `_cl0` run instead. Now:
- each rebuilt name must equal its own source file;
- Cl and `_tend` runs are skipped with a message.

A dry run lists 43 runs to redo, 5 of them continental.

**Termination counts (ocean worlds):**
- **New:** converged 2936, out_of_domain 1741, timeout 676, chemistry_void 49, wall_timeout 39,
  fallback_limit 39, solver_failure 11. The 27 fallback_limit runs in the α arm are at 10× and 30× `ALPHA_REF`.
- **Old:** timeout 6428, out_of_domain 4424, fallback_limit 28, chemistry_void 6, converged 8.

### 37.25 The continental baseline ran without Cl; the 'earth' sweep

**The problem.** The new `continental_baseline_ions` figure showed:
- Cl off the axis;
- Alk ≈ 590 mM and DIC ≈ 360 mM;
- Ca 0.7 mM;
- Na ≈ 490 mM.

`continental_baseline.py` runs through `parameter_sweep.run_simulation`, so it had inherited the sweep settings
(no Cl, blank ocean). The user had already moved it to 3.7 km. With no Cl⁻, the continental Na⁺ is balanced by
alkalinity. That makes a soda ocean, supersaturated with calcite, on the Ca-collapsed branch.

**Test at S = 1, land 0.3, 3.7 km** (mM):

| setup | T (K) | pCO₂ (ppm) | Alk | C | Ca | Mg | Na | Cl |
|---|---|---|---|---|---|---|---|---|
| blank, no Cl (as run) | 294.3 | 688 | 587 | 360 | 0.70 | 46 | 493 | 0 |
| seeded, Cl 0.02 (calibration setup) | 293.8 | 611 | 5.3 | 4.6 | 10.5 | 56.2 | 464 | 546 |
| blank, Cl 0.02, 4 Gyr | 294.4 | 698 | 288 | 183 | 0.71 | 47 | 501 | 308, still rising |

- **The calibration setup reproduces the calibration exactly.**
- **Cl alone isn't enough.** From a blank ocean, Cl is still rising at 4 Gyr and the ocean stays collapsed.
- **The climate barely moves,** because continental weathering sets Earth's carbon balance.

**Implemented:**
- **`parameter_sweep.run_simulation`** takes a per-run `seed`, and the resume guard now compares against it.
- **`continental_baseline.py`** has `SWEEP = 'earth'`: land 0.3, S 0.3–1.45 (24 runs), every other axis at
  Earth. It uses `EARTH_SEED` (seawater: Cl 546, SO₄ 28.2, K 10.2 mM), `EARTH_CL_RATIO` = 0.02 and
  `EARTH_T_END_GYR` = 4, as in `calibrate_earth.py`.
  - Names end `_land0.3_tend4`, distinct from the baseline's `_land0.3_cl0`.
  - `'all'` now includes it.
  - `run()` passes `cl`, `t_end_gyr` and `seed` through.
- **`plot_results.plot_continental_baseline`** draws `continental_baseline_{tp,chem,ions}` from these runs. It
  pins the Cl ratio and run length through `_arm(..., setup)`, and skips with a message if there are none.
  The other continental figures still use the sweep-setup baseline arms.

**Status.**
- The 24 runs were made in scratch and copied into `sweep_output` (no name clashes). The figures were redrawn.
- At S = 1: 293.79 K, 611 ppm, Na 464, Ca 10.5, Mg 56.2, Cl 546 mM.
- S 1.2–1.35 hit the 900 s cap on the hot branch (365–376 K). Rerun them with a longer
  `KAMINO_WALL_SHALLOW`, since `rerun_wall_timeouts.py` skips Cl runs.
- The earlier `SWEEP = 'all'` job stopped before its land-fraction series and full baseline line. Both 3.7 km
  arms have only the 9 coarse S points, and 9 half-written files remain.

### 37.26 The new sweep against the old one

Compared run by run: `/data/pt426/sweep_output` (new) against `/data/pt426/sweep_output_seeded` (old, 9–18 Sep).

| | old | new |
|---|---|---|
| α | 4.9 (anchored to a primary flux at pe +4, mostly Fe) | 348.3 (net 0.9 Teq/yr, Coogan & Dosso 2022) |
| η_HT / η_Na / τ_rw | 1.97e-2 / 6.1e-3 / 5 Myr | 3.05e-2 / 5.25e-3 / 39.2 Myr |
| Cl and initial ocean | Cl 0.02, seeded (Cl 546, SO₄ 23.45 mM) | no Cl, blank |
| crust rate in θ_r; sedimentation; floor | 1/50 Myr; old terrigenous + 0.3 m/Myr floor; 274 K | 1/130 Myr; dust + cosmic floor + porosity; 273.15 K |
| convergence | `r_avg` EMA | 50 Myr windowed drift |
| depth and composition sweeps | 0.1× outgassing, pe −3 and +4 | 1× outgassing, pe −3 only |

**Basic plane** (3 km, Earth crust, reverse weathering on, pe −3; 931 paired runs):

| | old | new |
|---|---|---|
| steady (converged or 2 Gyr timeout) | 438 (47 %) | 541 (58 %) |
| stopped at the CO₂ ceiling / hot wall | 392 / 98 | 294 / 88 |
| temperate, 273–320 K | 171 (18 %) | 126 (14 %) |
| steady with the seafloor at its floor | 193 | 370 |
| kinetic (Da < 1) and above the floor | 198 | 55 |
| median Da of steady runs, by crust rate | 0.01–0.13 (0.82 at 0.01×) | 3.8–16 |
| how steady runs ended | all 438 ran to 2 Gyr | 417 converged, median 716 Myr |

- **Climate.**
  - Runs steady in both are a median 6.4 K colder, with pCO₂ 0.8 dex lower (IQR −31 to 0 K).
  - By crust rate: 0.01× −44 K and −2.4 dex; 0.03× −27 K; 0.1× −7 K; ≥ 0.3× about 0.
  - By outgassing: 0.3–3× −29 to −35 K; 10× −5 K; ≤ 0.03× unchanged, because both sets are floored there.
  - 93 runs that hit the CO₂ ceiling in the old set now have a steady state.
- **The kinetic band** (Da < 1, above the floor) and the first thermodynamic S:
  - Old: the band existed at every outgassing down to 0.01×.
  - New: there is none at ≤ 0.03× outgassing, and none at 0.1× outgassing for crust ≥ 0.1×. It has moved to
    0.1–10× outgassing and narrowed.
  - At 1× outgassing and 1× crust: S 0.60–1.00 with the transition at 1.05, now S 0.55–0.85 with the transition
    at 0.90.
- **Chemistry.** Medians over runs steady in both, in mM:

  | | Cl | Ca | Mg | Alk | DIC | Si | SO₄ | Na | pH | salinity (g/kg) |
  |---|---|---|---|---|---|---|---|---|---|---|
  | old | 373 | 60 | 18 | 12 | 14 | 4.5 | 23.5 | 0.02 | 6.7 | 24 |
  | new | 0 | 0.14 | 12 | 23 | 21 | 2.9 | 0 | ≈ 0 | 9.2 | 1.7 |

**Earth arm** (land 0.3, seeded with Cl in both; the old one at 3 km, the new one at 3.7 km):
- **T is within 1.2 K for S ≥ 0.9.** At S = 1 it is 294.3 → 293.8 K and 686 → 611 ppm.
- **It is 2–10 K colder at S 0.35–0.85,** where seafloor weathering matters more.
- **S 0.35 now converges.** Before, it hit the CO₂ ceiling.
- **The S = 1 seafloor flux rose 27×,** from 0.033 to 0.89 Teq/yr, with Da 0.0035 → 0.084.
- **The ions moved towards their targets:** Mg 74 → 56 mM (target 52.8), Na 427 → 464 (469), and SO₄ is now
  an explicit 28.2 mM.

**Composition** (not directly comparable: the old set is at 0.1× outgassing, the new at 1×):
- **The Mg/Si pattern holds.** Mg/Si ≥ 1.75 is cold and on the floor in both, with identical T, because the
  floor decouples T from outgassing. Mg/Si 0.5 is hot in both.
- **At Mg/Si 0.8–1.5, S 0.6–0.8,** the new 1× runs are within 2–7 K of the old 0.1× runs. The larger α roughly
  offsets the 10× higher outgassing.
- **At S = 1** the new runs are on the hot branch (348 K), where the old ones were at 299 K.
- **The depth sweep** also moved from 0.1× to 1× outgassing.

**Why.**
1. **α is up about 70×, and it sets the climate differences.**
   - The old anchor left the model's net seafloor flux at ~0.04 Teq/yr, ~20× below Coogan & Dosso (§37.16).
   - At a given outgassing, weathering now draws CO₂ lower. Planets are colder, more of them sit on the seafloor
     floor (no feedback there), and some are CO₂-starved.
   - Da ∝ α, so planets become thermodynamic at lower S. Together these shift and shrink the negative-feedback
     band.
   - The stronger sink also balances more outgassing, so more runs have steady states.
   - The §37.17 sedimentation (−30–35 % ocean-world reactive area) and the 1/130 Myr crust rate in θ_r partly
     offset the α increase.
   - Isolated test (§37.16): at 1× outgassing, α 14.6 → 440 cooled S 0.4–0.8 by 40–60 K.
2. **No Cl and a blank ocean set the chemistry differences.**
   - The old seed put in 546 mM NaCl. Ocean worlds have no Na source, so the Na sink removed the Na. Cl has no
     fast sink (~4 Gyr), so it stayed. Charge balance then left a Ca–Mg chloride brine: Ca 60 mM, low
     alkalinity, pH 6.7.
   - Without Cl, weathering cations are balanced by alkalinity. The result is a dilute, alkaline Mg-bicarbonate
     ocean, with Ca held low by calcite saturation.
   - The climate effect of Cl is second-order (§37.9): ≤ 2–4 K in temperate states, 8–9 K near the inner edge.
3. **τ_rw 5 → 39 Myr** weakens reverse weathering, a CO₂ source on hot worlds. The hot branch is 4–5 K cooler,
   and the transition moves by about +0.05 in S (§37.22).
4. **Minor:** the floor 274 → 273.15 K, and the convergence check, which stops runs earlier without changing
   their states (≤ 0.8 %, §37.19). The old EMA almost never fired, so nearly every old run ran to 2 Gyr.
5. **Earth barely moves** because continental weathering sets its carbon balance. The seafloor increase shows
   only at low S.

**Caveat.** The new deep (≥ 20 km) and low-crust (0.01–0.03×) runs are 2 Gyr snapshots, and some oscillate
(§37.24). Numbers from those parts of the grid aren't steady states.
