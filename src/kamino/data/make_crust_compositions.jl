#=
make_crust_compositions.jl — primary-melt compositions on a (Mg/Si, dIW) grid, via MAGEMin.

Generates `crust_compositions.csv`, the lookup table that `crust_composition.oxide_composition`
interpolates. Two axes, both anchored on quantities that stellar spectroscopy constrains:

  mg_si     mantle molar Mg/Si, 0.5-2.0. Controls olivine vs orthopyroxene in the mantle
            (Guimond et al. 2024, section 3.1.1) and, through the melt, the normative
            feldspar/feldspathoid balance. Mantles go olivine-free below ~0.8 and
            orthopyroxene-free above ~1.6; the grid deliberately crosses both.

  delta_iw  CORE-FORMATION oxygen fugacity, log10 units relative to iron-wustite, -5 to -1.
            Sets how much iron the mantle kept as FeO instead of losing it to the core as metal,
            via the metal-silicate equilibrium dIW = 2 log10(a_FeO_sil / a_Fe_met). The axis is
            LOGARITHMIC in FeO: -5 -> 0.26 wt%, -2 -> 8.05 (Earth), -1 -> 24.1, which is already
            the ~25 wt% ceiling above which HGP18 is unreliable. Do not extend it upward.

            This is NOT the fO2 of a modern degassing melt. After core formation the mantle
            self-oxidises by Fe disproportionation (3FeO -> Fe0 + Fe2O3; Guimond et al. 2024
            section 2.4), reaching ~FMQ ~ IW+3.5 on Earth at unchanged total Fe. The derived
            modern value is written to the `delta_iw_melt` column for the outgassing model; the
            crust itself only ever sees the core-formation value, through mantle FeO.

Why MAGEMin and not pMELTS
--------------------------
pMELTS is calibrated for PERIDOTITE melting and cannot cover the stellar Mg/Si range: it fails
outright above molar Mg/Si ~1.6 (its solution models do not span the ferropericlase-bearing
assemblages stable there) and below ~0.8 it extrapolates badly, returning 69 wt% SiO2 rhyolites.
MAGEMin (Riel et al. 2022) uses the Holland, Green & Powell (2018) igneous dataset, which spans
ultramafic to granitic, and converges across the whole range. It stabilises ferropericlase and
nepheline at high Mg/Si on its own -- independent confirmation, from a different thermodynamic
dataset, that the CIPW desilication cascade in crust_composition.cipw_norm infers the right phase.

Method
------
  1. bulk mantle = McDonough & Sun (1995) pyrolite, modified on the two axes IN THIS ORDER so
     they stay orthogonal (the failure mode of section 23.4 was an axis that moved two things):
       a. set FeO from delta_iw, renormalising the non-Fe oxides to (100 - FeO) at their
          pyrolite proportions -- iron removed to the core leaves the mantle, so the rest
          rescales;
       b. within that budget, re-split MgO/SiO2 at fixed (MgO + SiO2) mass to hit mg_si,
          leaving Al/Ca/Na/Ti/Cr untouched;
  2. start at PSTART on the solid adiabat for a trial potential temperature;
  3. decompress ISENTROPICALLY to PEND. MAGEMin is a fixed-(P,T) Gibbs minimiser with no
     isentropic mode, so each step root-finds the temperature that holds entropy constant.
     This matters: prescribing the solid adiabat ignores the latent heat of melting and
     over-melts, by ~50 C at T_p = 1350;
  4. bisect the potential temperature until the melt fraction hits F_TARGET;
  5. batch melting -- the melt stays with the residue, so the final liquid is the pooled
     primary melt, and its composition is what the CIPW norm converts to mineralogy.

PEND is the melt SEGREGATION pressure, not the base of the crust. Under batch melting the liquid
re-equilibrates at every step, so carrying it to 0.2 GPa yields an over-equilibrated andesite.
1 GPa is a representative mean segregation depth beneath a ridge.

Environment
-----------
    julia -e 'using Pkg; Pkg.add("MAGEMin_C")'      # one-off
    julia make_crust_compositions.jl                 # production sweep (~5.5 h serial)
    julia make_crust_compositions.jl --slice 1/8     # one of 8 disjoint shards, for parallelism
    python merge_crust_slices.py                     # concatenate the shards, checking for holes
    julia make_crust_compositions.jl --probe         # 3-point timing/convergence check FIRST
    julia make_crust_compositions.jl --calibrate     # T_p scan at Earth against PRIMELT
    julia make_crust_compositions.jl --fixed-p 20    # isobaric batch melting, Guimond-comparable
    julia make_crust_compositions.jl --closure cpx-out --out sensitivity.csv

Citations for the paper
-----------------------
  * MAGEMin: Riel, Kaus, Green & Berlie (2022), G3 23, e2022GC010427.
  * Thermodynamic dataset: Holland, Green & Powell (2018), J. Petrol. 59, 881.
  * Pyrolite bulk mantle: McDonough & Sun (1995), Chem. Geol. 120, 223.
  * Melt fraction and the Mg/Si construction: Guimond et al. (2024), RiMG 90, 259.
=#

using MAGEMin_C, Printf

const OX = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","Cr2O3","H2O"]

# Melting path, kbar. 3 GPa is sub-solidus across this T_p range; see the PEND note above.
const PSTART, PEND, DP = 30.0, 10.0, 2.0

# Melt fraction, fixed across the whole grid.
#
# 0.20 is the value Guimond et al. (2024) adopt, and their justification is petrological rather
# than arbitrary: for typical Earth mantle it is the degree at which clinopyroxene is lost from
# the melting assemblage, beyond which melting is much less productive (Katz et al. 2003). Using
# their number makes this table directly comparable to their Figure 9.
#
# Holding F rather than T_p is deliberate. A mantle that cannot melt also cannot transport heat
# by magmatism, so it warms until melting carries the heat out -- mantle temperature is
# self-regulated, not free (Nature Comms Earth Environ 2022, 3, 261; the counterweight is
# Korenaga 2016, who argues the adjustment is too slow to fully self-regulate). Holding T_p fixed
# instead makes the Mg-rich end look like planets that "barely melt", which is an artifact of the
# closure, not a result.
#
# KNOWN LIMITATION, do not let it pass silently. F = 0.20 is at or PAST cpx-out at the Mg-rich
# end, which is exactly the condition that produces ultracalcic melts (CaO/Al2O3 up to 1.77
# against MORB's ~0.78 -- a clinopyroxenite where a picrite belongs; Medard et al. 2004 exclude
# such melts from a volatile-free fertile lherzolite source, which is what this is). The grid
# therefore reports CaO_Al2O3 per row and warns above ULTRACALCIC_WARN. Run `--closure cpx-out`
# for the alternative that solves F per composition instead.
const F_TARGET = 0.20
const ULTRACALCIC_WARN = 1.2

# REJECTED ALTERNATIVE -- constant homologous temperature T_p/T_solidus, which is the quantity the
# mantle-buffering literature actually identifies. It fails here for an instructive reason: the
# true multi-component solidus is set by the first infinitesimal melt, which minor alkalis
# control. Na2O stays at pyrolite's 0.36 wt% while SiO2 falls, so silica-poor bulks become
# nepheline-normative and their alkaline eutectic melts LOW -- the solidus collapses 1509 -> 1049 C
# from Mg/Si 1.25 to 2.0, driving T_p DOWN 400 C and extinguishing melting altogether. The solidus
# measures a trace-driven eutectic, not the bulk refractoriness that governs heat transport. (The
# literature definition uses a dry peridotite solidus parameterization, which coincides with the
# true solidus for Earth-like compositions but not across this range.) Do not retry it.

# Bracket and tolerance for the T_p root-find. F is monotonic in T_p, so bisection is safe. The
# bracket is wider than the constant-F = 0.117 version needed: F = 0.20 pushes every solution
# hotter, and the reduced (low-dIW, low-FeO) mantles are more refractory still -- less FeO-rich
# mantles melt up to ~100 C higher (Guimond et al. 2024, section 4.1).
const TP_LO, TP_HI, TP_TOL = 1150.0, 1900.0, 5.0

# Earth's anchor on each axis. dIW = -2 reproduces BSE FeO = 8.05 wt% by construction.
const MGSI_EARTH, DIW_EARTH = 1.25, -2.0

# Earth's mantle potential temperature is NOT set here -- it is solved, like every other grid
# point, as the T_p delivering F_TARGET. Recorded for reference: at F = 0.117 the calibration in
# section 24.2 gave T_p = 1325 C. F = 0.20 lands hotter; --calibrate reports the new value.
const PRIMELT_REF = (SiO2=48.76, MgO=11.27, Al2O3=17.05, CaO=11.91)  # PRIMELT primary melt, T_p 1350

# Molar Mg/Si grid. Non-uniform on purpose: dense where the assemblage changes -- olivine-out
# near 0.8, orthopyroxene-out and the nepheline/ferropericlase switch near 1.6.
const MGSI_GRID = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

# Core-formation dIW grid. Uniform 0.5 spacing = uniform in log FeO, which is what keeps the
# downstream linear interpolation honest (FeO spans 0.26 -> 24 wt% across it).
const DIW_GRID = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]

# Oxides written to the CSV, matching what cipw_norm consumes (plus Cr2O3, which it drops).
const CSV_OXIDES = ["SiO2","TiO2","Al2O3","Cr2O3","FeO","MgO","CaO","Na2O","K2O"]

# --- Mantle construction -----------------------------------------------------------------------

# McDonough & Sun (1995) pyrolite, NON-Fe oxides (wt%). Must stay identical to PYROLITE_NON_FE in
# crust_composition.py -- the two modules share this calibration.
const PYROLITE_NON_FE = Dict("SiO2"=>45.00, "Al2O3"=>4.45, "CaO"=>3.55, "MgO"=>37.80,
                             "K2O"=>0.029, "Na2O"=>0.36, "TiO2"=>0.201, "Cr2O3"=>0.384)
const M_OX = Dict("SiO2"=>60.084, "Al2O3"=>101.961, "CaO"=>56.077, "MgO"=>40.304, "FeO"=>71.844,
                  "K2O"=>94.196, "Na2O"=>61.979, "TiO2"=>79.866, "Cr2O3"=>151.990)
const NCAT = Dict("SiO2"=>1, "Al2O3"=>2, "CaO"=>1, "MgO"=>1,
                  "K2O"=>2, "Na2O"=>2, "TiO2"=>1, "Cr2O3"=>2)

const PYRO_SUM = sum(values(PYROLITE_NON_FE))
# Cation moles contributed by 1 wt% of the renormalised non-Fe budget.
const K_CATIONS = sum(PYROLITE_NON_FE[o] / PYRO_SUM * NCAT[o] / M_OX[o] for o in keys(PYROLITE_NON_FE))
const EARTH_MANTLE_FEO = 8.05

x_feo_from_wt(feo) = (feo / M_OX["FeO"]) / (feo / M_OX["FeO"] + (100.0 - feo) * K_CATIONS)

# Metal-silicate activity constant X_Fe(metal)/gamma_FeO, CALIBRATED so dIW = -2 gives 8.05 wt%.
# The first-principles value carries ~0.4 log units of slop (gamma_FeO 1.0 -> Earth dIW -2.35,
# 1.5 -> -1.99), so only the Earth anchor is meaningful. Mirrors _FEO_ACTIVITY_CONST in
# crust_composition.py; change both together.
const FEO_ACTIVITY_CONST = x_feo_from_wt(EARTH_MANTLE_FEO) / 10.0^(DIW_EARTH / 2)

"Mantle FeO (wt%) implied by a core-formation oxygen fugacity, from Fe + 1/2 O2 = FeO."
function feo_from_delta_iw(diw::Float64)
    x = FEO_ACTIVITY_CONST * 10.0^(diw / 2)
    0.0 < x < 1.0 || error("delta_iw=$diw implies X_FeO=$x, outside (0,1)")
    return 100.0 * x * K_CATIONS / ((1.0 - x) / M_OX["FeO"] + x * K_CATIONS)
end

"""
Pyrolite with FeO set from the redox axis and MgO/SiO2 re-split to the target molar Mg/Si.

Order matters and is the point: iron first (it changes the size of the silicate budget, because
Fe that went to the core is no longer in the mantle), then Mg/Si within what remains. Applied the
other way round the two axes would not be orthogonal.
"""
function mantle_composition(mgsi::Float64, feo_wt::Float64)
    mgsi > 0 || error("mg_si must be positive")
    0.0 <= feo_wt < 100.0 || error("feo_wt out of range")
    scale = (100.0 - feo_wt) / PYRO_SUM
    ox = Dict(o => v * scale for (o, v) in PYROLITE_NON_FE)
    tot = ox["MgO"] + ox["SiO2"]
    si  = tot / (1 + mgsi * M_OX["MgO"] / M_OX["SiO2"])
    ox["SiO2"], ox["MgO"] = si, tot - si
    # Order must match OX exactly.
    return [ox["SiO2"], ox["Al2O3"], ox["CaO"], ox["MgO"], feo_wt,
            ox["K2O"], ox["Na2O"], ox["TiO2"], 0.0, ox["Cr2O3"], 0.0]
end

# Ferric iron is deliberately OFF: the O component above is 0.0, i.e. all iron is ferrous. Fe3+ is
# second-order for melt major elements, and switching it on would perturb the Earth anchor that
# every downstream calibration rests on. Guimond et al. run their pMELTS grids the same way. The
# derived Fe3+/SigmaFe belongs in the outgassing handoff, not in this melting calculation.

"Core mass fraction implied by keeping `feo_wt` in the mantle, at solar bulk Fe/Mg. INFORMATIONAL."
function core_mass_fraction(mgsi::Float64, feo_wt::Float64)
    # Solar Fe/Mg molar ~0.83 (Lodders 2003). Kamino sets planet mass and radius independently,
    # so this column never feeds back into the model -- it is here for the paper's tables.
    X = mantle_composition(mgsi, feo_wt)
    n_mg = X[4] / M_OX["MgO"]
    n_fe_mantle = feo_wt / M_OX["FeO"]
    n_fe_core = 0.83 * n_mg - n_fe_mantle
    n_fe_core <= 0 && return 0.0            # mantle already holds all the iron: no metal core
    m_core = n_fe_core * 55.845
    return m_core / (m_core + 100.0)
end

# --- MAGEMin driving ---------------------------------------------------------------------------

minim(data, X, P, T) = single_point_minimization(P, T, data, X=X, Xoxides=OX, sys_in="wt")

"Solid mantle adiabat, used only to set the STARTING temperature of the isentrope."
adiabat_T(Tp, P_kbar) = (Tp + 273.15) * exp(3.0e-5 * P_kbar * 1e8 / (3300.0 * 1200.0)) - 273.15

"""
Temperature at pressure `P` whose entropy matches `S0`, by secant iteration.

MAGEMin has no isentropic mode, so the isentrope has to be tracked by hand. Entropy is monotonic
in T at fixed P, which makes this well behaved; the step is clamped so a flat secant near a phase
boundary cannot throw the iterate across the melting interval.
"""
function T_at_entropy(data, X, P, Tguess, S0)
    T1, T2 = Tguess, Tguess - 8.0
    s1 = minim(data, X, P, T1).entropy[1]
    for _ in 1:25
        s2 = minim(data, X, P, T2).entropy[1]
        abs(s2 - S0) < 1e-6 && return T2
        ds = s2 - s1
        abs(ds) < 1e-12 && return T2
        Tn = clamp(T2 - (s2 - S0) * (T2 - T1) / ds, T2 - 60.0, T2 + 60.0)
        T1, s1, T2 = T2, s2, Tn
    end
    return T2
end

"Isentropic decompression melting. Returns (MAGEMin output, final T)."
function isentropic_melt(data, X, Tp)
    T = adiabat_T(Tp, PSTART)
    out = minim(data, X, PSTART, T)
    S0 = out.entropy[1]
    P = PSTART
    while P - DP >= PEND - 1e-9
        P -= DP
        T = T_at_entropy(data, X, P, T, S0)
        out = minim(data, X, P, T)
    end
    return out, T
end

"Isobaric batch melting at fixed P, marching T until the melt fraction reaches `Ftarget`."
function isobaric_melt(data, X, P, Ftarget)
    lo, hi = 900.0, 2200.0
    while hi - lo > 1.0
        mid = (lo + hi) / 2
        _, F = melt_oxides(minim(data, X, P, mid))
        F < Ftarget ? (lo = mid) : (hi = mid)
    end
    T = (lo + hi) / 2
    return minim(data, X, P, T), T
end

"Melt composition in wt%, renormalised over the oxides the CIPW norm consumes, plus melt fraction."
function melt_oxides(out)
    i = findfirst(==("liq"), out.ph)
    i === nothing && return nothing, 0.0
    frac = out.frac_M_wt
    (frac === nothing || isnan(frac) || frac <= 0) && return nothing, 0.0
    comp = Dict(OX[j] => out.bulk_M_wt[j] * 100 for j in eachindex(OX))
    kept = Dict(o => get(comp, o, 0.0) for o in CSV_OXIDES)
    tot = sum(values(kept))
    tot <= 0 && return nothing, frac
    return Dict(o => v / tot * 100 for (o, v) in kept), frac
end

"Phases in the residue (everything but the liquid), as a ;-joined string."
residue(out) = join(filter(!=("liq"), out.ph), ";")

# --- Closures ----------------------------------------------------------------------------------

"""
Potential temperature giving melt fraction `Ftarget`, by bisection. NaN if unreachable.

F increases monotonically with T_p, so bisection is safe. A composition too refractory to reach
the target even at TP_HI is a result, not something to extrapolate through.
"""
function Tp_for_F(data, X, Ftarget)
    _, F_hi = melt_oxides(isentropic_melt(data, X, TP_HI)[1])
    F_hi < Ftarget && return NaN
    lo, hi = TP_LO, TP_HI
    while hi - lo > TP_TOL
        mid = (lo + hi) / 2
        _, F = melt_oxides(isentropic_melt(data, X, mid)[1])
        F < Ftarget ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
Potential temperature at which clinopyroxene is just exhausted from the residue.

The alternative closure to F_TARGET, and the one section 24.5 showed removes the ultracalcic
melts (CaO/Al2O3 falls to 0.84-1.01 across almost the whole grid). Melt fraction comes out
solved rather than assumed. NaN if cpx is absent at TP_LO (nothing to exhaust) or still present
at TP_HI (never exhausts in range) -- both are reported, never silently substituted.
"""
function Tp_for_cpx_out(data, X)
    has_cpx(Tp) = any(p -> p in ("cpx", "aug"), isentropic_melt(data, X, Tp)[1].ph)
    has_cpx(TP_LO) || return NaN
    has_cpx(TP_HI) && return NaN
    lo, hi = TP_LO, TP_HI
    while hi - lo > TP_TOL
        mid = (lo + hi) / 2
        has_cpx(mid) ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# --- Grid --------------------------------------------------------------------------------------

const CSV_HEADER = "mg_si,delta_iw,mantle_feo,T_p,melt_fraction,T_end,closure,residual_phases," *
                   join(CSV_OXIDES, ",") *
                   ",CaO_Al2O3,mg_number,delta_iw_melt,core_mass_fraction,warnings"

# Earth anchor for the self-oxidation bridge; mirrors DELTA_IW_SELF_OXIDATION in constants.py.
const SELF_OXIDATION = 5.5

"""
Validate against a published experiment: melt an EXPLICIT bulk composition at fixed P to a
fixed melt fraction, and report the liquid. Used to compare against piston-cylinder data
(Brugman et al. 2021), where both the bulk and the run conditions are prescribed.

    julia make_crust_compositions.jl --validate --bulk "42.0,4.85,3.76,40.04,8.23,0.0,0.21,0.0" \
                                     --fixed-p 15 --ftarget 0.05
"""
function validate(data, bulkspec, P, Ftarget, mgsi, diw)
    X = if bulkspec == ""
        mantle_composition(mgsi, feo_from_delta_iw(diw))
    else
        v = parse.(Float64, split(bulkspec, ","))
        length(v) == 8 || error("--bulk needs 8 values: SiO2,Al2O3,CaO,MgO,FeO,K2O,Na2O,TiO2")
        [v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], 0.0, 0.0, 0.0]
    end
    @printf("bulk (wt%%): SiO2 %.2f  Al2O3 %.2f  CaO %.2f  MgO %.2f  FeO %.2f  Na2O %.2f\n",
            X[1], X[2], X[3], X[4], X[5], X[7])
    @printf("isobaric at %.1f kbar, melting to F = %.3f\n", P, Ftarget)
    flush(stdout)
    out, T = isobaric_melt(data, X, P, Ftarget)
    m, F = melt_oxides(out)
    if m === nothing
        println("no melt"); return
    end
    @printf("T = %.0f C   F = %.4f   residue = %s\n", T, F, residue(out))
    @printf("melt: SiO2 %.2f  MgO %.2f  FeO %.2f  Al2O3 %.2f  CaO %.2f  Na2O %.2f  TiO2 %.2f\n",
            m["SiO2"], m["MgO"], m["FeO"], m["Al2O3"], m["CaO"], m["Na2O"], m["TiO2"])
    flush(stdout)
end

"One grid point. Returns (csv_row | nothing, status_string)."
function grid_point(data, mgsi, diw, closure, fixedP)
    feo = feo_from_delta_iw(diw)
    X = mantle_composition(mgsi, feo)
    warns = String[]
    feo > 25.0 && push!(warns, "mantle FeO $(round(feo,digits=1)) wt% above HGP18 reliability limit")

    local out, Tend, Tp
    if fixedP > 0
        out, Tend = isobaric_melt(data, X, fixedP, F_TARGET)
        Tp = NaN
    else
        Tp = closure == "cpx-out" ? Tp_for_cpx_out(data, X) : Tp_for_F(data, X, F_TARGET)
        isnan(Tp) && return nothing, (closure == "cpx-out" ?
            "no cpx-out bracket in [$TP_LO, $TP_HI] (absent at the cold end or never exhausted)" :
            "cannot reach F=$F_TARGET below $TP_HI C")
        out, Tend = isentropic_melt(data, X, Tp)
    end

    m, F = melt_oxides(out)
    m === nothing && return nothing, "no melt at the end of the path"

    ca_al = m["Al2O3"] > 0 ? m["CaO"] / m["Al2O3"] : NaN
    ca_al > ULTRACALCIC_WARN && push!(warns,
        "ultracalcic: CaO/Al2O3 $(round(ca_al,digits=2)) > $ULTRACALCIC_WARN (melting past cpx-out)")
    n_mg = m["MgO"] / M_OX["MgO"]
    mgnum = n_mg / (n_mg + m["FeO"] / M_OX["FeO"])

    row = @sprintf("%.4g,%.4g,%.4f,%s,%.6g,%.1f,%s,%s,", mgsi, diw, feo,
                   isnan(Tp) ? "NaN" : @sprintf("%.0f", Tp), F, Tend, closure, residue(out)) *
          join([@sprintf("%.6g", m[o]) for o in CSV_OXIDES], ",") *
          @sprintf(",%.4f,%.4f,%.2f,%.4f,%s", ca_al, mgnum, diw + SELF_OXIDATION,
                   core_mass_fraction(mgsi, feo), join(warns, ";"))
    return row, (isempty(warns) ? "ok" : join(warns, "; "))
end

"""
Shard `n` of `k`: every k-th Mg/Si value, so the shards are disjoint and together exhaustive.

Grid points cost ~130 s each, so the full 153-point grid is ~5.5 h in one process. Sharding by
PROCESS rather than by thread is deliberate: MAGEMin keeps mutable per-point workspaces inside the
`data` handle, and section 23.5 of the development history records what a poisoned solver looks
like -- a failure that does not raise, leaves the library unrecoverable, and silently degrades
every later calculation in the process. Separate processes cannot do that to each other.

Sharding on Mg/Si rather than on the flat index keeps each shard's cost similar: cost varies far
more along Mg/Si (89-156 s in the probe) than along dIW.
"""
function shard(grid, n, k)
    1 <= n <= k || error("--slice n/k needs 1 <= n <= k")
    return grid[n:k:end]
end

function sweep(data, outpath, closure, fixedP, mgsi_grid, diw_grid)
    @printf("MAGEMin melting: %d Mg/Si x %d dIW = %d points, closure=%s, %s\n",
            length(mgsi_grid), length(diw_grid), length(mgsi_grid)*length(diw_grid), closure,
            fixedP > 0 ? @sprintf("isobaric %.1f GPa", fixedP/10) :
                         @sprintf("isentropic %.1f -> %.1f GPa", PSTART/10, PEND/10))
    closure == "fixed-F" && @printf("  F held at %.3f; T_p solved per composition.\n", F_TARGET)
    flush(stdout)

    rows = String[CSV_HEADER]
    failures = String[]
    warned = 0
    t0 = time()
    for mgsi in mgsi_grid, diw in diw_grid
        row, status = try
            grid_point(data, mgsi, diw, closure, fixedP)
        catch e
            (nothing, "EXCEPTION $(typeof(e))")
        end
        if row === nothing
            @printf("  Mg/Si=%.2f dIW=%+.1f  SKIPPED: %s\n", mgsi, diw, status)
            push!(failures, @sprintf("Mg/Si=%.2f dIW=%+.1f: %s", mgsi, diw, status))
        else
            f = split(row, ",")
            @printf("  Mg/Si=%.2f dIW=%+.1f FeO=%5.2f T_p=%-4s F=%.3f SiO2=%5.2f MgO=%5.2f CaO=%5.2f Ca/Al=%.2f %s\n",
                    mgsi, diw, parse(Float64, f[3]), f[4], parse(Float64, f[5]),
                    parse(Float64, f[9]), parse(Float64, f[14]), parse(Float64, f[15]),
                    parse(Float64, f[18]), status == "ok" ? "" : "[$status]")
            status != "ok" && (warned += 1)
            push!(rows, row)
        end
        flush(stdout)
    end

    open(outpath, "w") do io
        println(io, "# Generated by make_crust_compositions.jl -- do not edit by hand.")
        println(io, "# closure=$closure  F_TARGET=$F_TARGET  path=" *
                    (fixedP > 0 ? "isobaric $(fixedP/10) GPa" : "isentropic $(PSTART/10)->$(PEND/10) GPa"))
        println(io, "# delta_iw is the CORE-FORMATION fO2 (sets mantle FeO); delta_iw_melt is the")
        println(io, "# derived modern-melt value for the outgassing model. core_mass_fraction")
        println(io, "# assumes solar Fe/Mg and is INFORMATIONAL -- kamino sets mass/radius directly.")
        for r in rows; println(io, r); end
    end

    @printf("\nWritten: %s  (%d rows, %d with warnings, %d failed) in %.1f min\n",
            outpath, length(rows)-1, warned, length(failures), (time()-t0)/60)
    if !isempty(failures)
        # Never let a partial grid pass as a complete one: crust_composition.load_crust_table
        # rejects a table that is not a full rectangle, and this is why.
        println("FAILURES (the table is NOT a complete grid and will be rejected downstream):")
        for f in failures; println("  ", f); end
    end
end

function calibrate(data)
    println("T_p scan at Earth's pyrolite (Mg/Si = $MGSI_EARTH, dIW = $DIW_EARTH), isentropic $(PSTART/10) -> $(PEND/10) GPa")
    @printf("%-7s %-8s %-7s %-8s %-8s %-8s %-8s %s\n",
            "T_p", "T_end", "F", "SiO2", "MgO", "Al2O3", "CaO", "misfit*")
    X = mantle_composition(MGSI_EARTH, feo_from_delta_iw(DIW_EARTH))
    for Tp in 1300.0:25.0:1500.0
        out, Tend = isentropic_melt(data, X, Tp)
        m, F = melt_oxides(out)
        if m === nothing
            @printf("%-7.0f %-8.1f no melt\n", Tp, Tend); flush(stdout); continue
        end
        misfit = abs(m["SiO2"]-PRIMELT_REF.SiO2) + abs(m["MgO"]-PRIMELT_REF.MgO) +
                 abs(m["Al2O3"]-PRIMELT_REF.Al2O3) + abs(m["CaO"]-PRIMELT_REF.CaO)
        @printf("%-7.0f %-8.1f %-7.3f %-8.2f %-8.2f %-8.2f %-8.2f %.2f\n",
                Tp, Tend, F, m["SiO2"], m["MgO"], m["Al2O3"], m["CaO"], misfit)
        flush(stdout)
    end
    println("\n*misfit = sum |melt - PRIMELT primary melt| over SiO2/MgO/Al2O3/CaO, wt%")
    println(" PRIMELT reference: SiO2 $(PRIMELT_REF.SiO2)  MgO $(PRIMELT_REF.MgO)  Al2O3 $(PRIMELT_REF.Al2O3)  CaO $(PRIMELT_REF.CaO)")
    println(" NOTE: at F = $F_TARGET the best-fit T_p is hotter, and the misfit LARGER, than the")
    println(" F = 0.117 calibration of section 24.2 (T_p 1325, misfit 2.22). That is the expected")
    println(" cost of adopting Guimond's melt fraction; quantify it, do not tune it away.")
end

"""
Run an explicit list of (Mg/Si, dIW) points under a chosen closure.

For spot-checking a closure against the literature without paying for a full grid:

    julia make_crust_compositions.jl --points "1.25,-2.0;1.6,-2.0" --closure cpx-out
"""
function points(data, spec, closure, fixedP)
    println("Closure = $closure")
    println(CSV_HEADER); flush(stdout)
    for pair in split(spec, ";")
        mgsi, diw = parse.(Float64, split(strip(pair), ","))
        t = time()
        row, status = try
            grid_point(data, mgsi, diw, closure, fixedP)
        catch e
            (nothing, "EXCEPTION $(typeof(e)): $e")
        end
        if row === nothing
            @printf("# Mg/Si=%.2f dIW=%+.1f  FAILED (%.0f s): %s\n", mgsi, diw, time()-t, status)
        else
            println(row)
            @printf("# ^ %.0f s  %s\n", time()-t, status)
        end
        flush(stdout)
    end
end

"Three-point convergence and timing check. Run this BEFORE committing to the full grid."
function probe(data)
    println("Probe: convergence + timing at the Earth anchor and both Mg/Si extremes.")
    flush(stdout)
    for (mgsi, diw) in [(MGSI_EARTH, DIW_EARTH), (0.5, DIW_EARTH), (2.0, DIW_EARTH)]
        t = time()
        row, status = try
            grid_point(data, mgsi, diw, "fixed-F", -1.0)
        catch e
            (nothing, "EXCEPTION $(typeof(e)): $e")
        end
        @printf("  Mg/Si=%.2f dIW=%+.1f  %.1f s  %s\n", mgsi, diw, time()-t,
                row === nothing ? "FAILED: $status" : "ok")
        row !== nothing && println("    ", row)
        flush(stdout)
    end
    @printf("\nFull grid is %d points. Coarsen DIW_GRID to 1.0 spacing if a point costs > ~2 min.\n",
            length(MGSI_GRID)*length(DIW_GRID))
end

function main()
    args = ARGS
    outpath = joinpath(@__DIR__, "crust_compositions.csv")
    i = findfirst(==("--out"), args); i !== nothing && (outpath = args[i+1])
    closure = "fixed-F"
    i = findfirst(==("--closure"), args); i !== nothing && (closure = args[i+1])
    fixedP = -1.0
    i = findfirst(==("--fixed-p"), args); i !== nothing && (fixedP = parse(Float64, args[i+1]))
    mgsi_grid = MGSI_GRID
    i = findfirst(==("--slice"), args)
    if i !== nothing
        n, k = parse.(Int, split(args[i+1], "/"))
        mgsi_grid = shard(MGSI_GRID, n, k)
        outpath = replace(outpath, ".csv" => @sprintf("_slice%d_of_%d.csv", n, k))
        @printf("Shard %d/%d: Mg/Si = %s -> %s\n", n, k, mgsi_grid, outpath)
    end

    data = Initialize_MAGEMin("ig", verbose=false)
    try
        i = findfirst(==("--points"), args)
        if "--validate" in args
            j = findfirst(==("--bulk"), args); bulk = j === nothing ? "" : args[j+1]
            j = findfirst(==("--ftarget"), args); ft = j === nothing ? F_TARGET : parse(Float64, args[j+1])
            j = findfirst(==("--mgsi"), args); mg = j === nothing ? MGSI_EARTH : parse(Float64, args[j+1])
            j = findfirst(==("--diw"), args); dw = j === nothing ? DIW_EARTH : parse(Float64, args[j+1])
            validate(data, bulk, fixedP > 0 ? fixedP : 15.0, ft, mg, dw)
        elseif "--calibrate" in args
            calibrate(data)
        elseif i !== nothing
            points(data, args[i+1], closure, fixedP)
        elseif "--probe" in args
            probe(data)
        else
            sweep(data, outpath, closure, fixedP, mgsi_grid, DIW_GRID)
        end
    finally
        Finalize_MAGEMin(data)
    end
end

main()
