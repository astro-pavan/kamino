#=
make_crust_compositions.jl — primary-melt compositions vs mantle Mg/Si, via MAGEMin.

Generates the lookup table that `crust_composition.oxide_composition` interpolates, replacing
the PRIMELT1 spreadsheet + SiO2-rescale proxy. The proxy scaled SiO2 only and renormalised,
which inflated Al2O3 (14 -> 20 wt%) and CaO (10 -> 14 wt%) as pure artifacts and pushed the
endpoints outside primary-melt space.

Why MAGEMin and not pMELTS
--------------------------
pMELTS (Ghiorso et al. 2002) is calibrated for PERIDOTITE melting and cannot reach the ends of
the stellar Mg/Si distribution: it fails outright above molar Mg/Si ~1.6 (the start state will
not converge at any temperature nudge) because its solution models do not span the
ferropericlase-bearing assemblages that appear there, and below ~0.8 it is extrapolating badly
outside its calibration, returning 69 wt% SiO2 rhyolites. MAGEMin (Riel et al. 2022) uses the
Holland, Green & Powell (2018) igneous dataset, which spans ultramafic to granitic, and
converges across the whole Hypatia range. It stabilises ferropericlase and nepheline at high
Mg/Si on its own — independent confirmation of the CIPW desilication cascade in
crust_composition.cipw_norm.

Method
------
  1. bulk mantle = McDonough & Sun (1995) pyrolite, with MgO/SiO2 re-split to the target molar
     Mg/Si at fixed (MgO + SiO2) mass, leaving Al/Ca/Fe/Na untouched — so Mg/Si is a genuinely
     orthogonal axis, which it is not in the proxy;
  2. start at P_START on the solid adiabat for the given potential temperature;
  3. decompress ISENTROPICALLY to P_END. MAGEMin is a fixed-(P,T) Gibbs minimiser with no
     isentropic mode, so each step root-finds the temperature whose entropy matches the start.
     This matters: prescribing the solid adiabat instead ignores the latent heat of melting and
     over-melts. At T_p = 1350 the isentrope arrives at 1309 C where the solid adiabat would
     give ~1362 C;
  4. batch melting — the melt stays with the residue, so the final liquid is the pooled primary
     melt, and its composition is what the CIPW norm converts to mineralogy.

P_END is the melt SEGREGATION pressure, not the base of the crust. Under batch melting the
liquid re-equilibrates at every step, so carrying it to 0.2 GPa yields an over-equilibrated
andesite. 1 GPa is a representative mean segregation depth beneath a ridge.

Environment
-----------
Julia's runtime must live on the data partition — the home partition is full, and juliaup
otherwise fails with "No space left on device":

    export PATH=/data/pt426/.juliaup/bin:\$PATH
    export JULIA_DEPOT_PATH=/home/pt426/data/julia_depot
    export JULIAUP_DEPOT_PATH=/home/pt426/data/julia_depot

Usage
-----
    julia make_crust_compositions.jl                      # production sweep
    julia make_crust_compositions.jl --calibrate          # T_p scan against Earth MORB
    julia make_crust_compositions.jl --tp 1350 --out other.csv

Citations for the paper
-----------------------
  * MAGEMin: Riel, Kaus, Green & Berlie (2022), G3 23, e2022GC010427.
  * Thermodynamic dataset: Holland, Green & Powell (2018), J. Petrol. 59, 881.
  * Pyrolite bulk mantle: McDonough & Sun (1995), Chem. Geol. 120, 223.
=#

using MAGEMin_C, Printf

const OX = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","Cr2O3","H2O"]

# Melting path, kbar. 3 GPa is sub-solidus across this T_p range; see the P_END note above.
const PSTART, PEND, DP = 30.0, 10.0, 2.0

# Melt fraction held constant across the Mg/Si sweep, instead of holding T_p constant.
#
# Physical basis: a mantle that cannot melt also cannot transport heat by magmatism, so it warms
# until melting carries the heat out -- mantle temperature is self-regulated, not free (Nature
# Comms Earth Environ 2022, 3, 261). Holding T_p fixed instead makes the Mg-rich end look like
# planets that "barely melt" (F ~ 0.003 at Mg/Si 2.0), which is an artifact of the closure: at
# constant F those planets simply run ~170 C hotter and produce a picritic melt (MgO ~16.5 wt%).
#
# The value is the melt fraction Earth's pyrolite produces at the calibrated T_p = 1325 C, so the
# Earth anchor is preserved exactly: Mg/Si = 1.25 still returns T_p = 1325 and its basalt.
#
# REJECTED ALTERNATIVE -- constant homologous temperature T_p/T_solidus, which is the quantity the
# buffering literature actually identifies. It fails here, badly, and for an instructive reason:
# the true multi-component solidus is set by the first infinitesimal melt, which minor alkalis
# control. Na2O stays at pyrolite's 0.36 wt% while SiO2 falls, so silica-poor bulks become
# nepheline-normative and their alkaline eutectic melts LOW -- the solidus collapses 1509 -> 1049 C
# from Mg/Si 1.25 to 2.0, driving T_p DOWN 400 C and extinguishing melting altogether. The solidus
# measures a trace-driven eutectic, not the bulk refractoriness that governs heat transport. (The
# literature definition uses a dry peridotite solidus parameterization, which coincides with the
# true solidus for Earth-like compositions but not across this range.) Do not retry it against the
# true solidus; if the homologous framing is wanted, define the reference against a fixed
# melt-fraction contour, which is near-identical to what is done here.
const F_TARGET = 0.117

# Bracket and tolerance for the T_p root-find. F is monotonic in T_p, so bisection is safe.
# The bracket is deliberately tight: every composition in MGSI_GRID solves to 1167-1523 C, and
# each bisection step costs a full isentropic path, so a 1100-2200 bracket spent about half the
# runtime resolving temperatures nothing needs. Widen it if MGSI_GRID is extended and a
# composition reports "cannot reach F".
const TP_LO, TP_HI, TP_TOL = 1150.0, 1650.0, 5.0

# Earth's mantle potential temperature (deg C), set by --calibrate. 1325 minimises the misfit to
# the PRIMELT primary melt (2.22 wt%, against 2.63 at 1350 and 2.96 at 1300) AND is the only
# value whose melt fraction, F = 0.117, lands in the 0.08-0.12 range of real MORB primary melts;
# 1350 gives F = 0.153, too melt-rich.
#
# This sits 25 C below the project's EARTH_MANTLE_POTENTIAL_TEMPERATURE of 1350. That is a model
# offset, not a claim about Earth: MAGEMin/HGP18 is slightly more melt-productive than PRIMELT at
# the same potential temperature, so matching the MELT -- which is what the weathering model
# actually consumes -- needs a marginally cooler label. Anchor on the melt, not the number.
const TP_EARTH = 1325.0

# Molar Mg/Si grid, spanning the Hypatia stellar range. Non-uniform on purpose: below ~1.3 the
# normative mineralogy barely moves, while Nepheline switches on and Albite collapses between
# 1.4 and 1.6, so the resolution goes there.
const MGSI_GRID = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.4, 1.45,
                   1.5, 1.55, 1.6, 1.7, 1.8, 1.9, 2.0]

# Oxides written to the CSV, matching what cipw_norm consumes.
const CSV_OXIDES = ["SiO2","TiO2","Al2O3","Cr2O3","FeO","MgO","CaO","Na2O","K2O"]

"McDonough & Sun pyrolite with MgO/SiO2 re-split to the target molar Mg/Si, total mass fixed."
function mantle_composition(mgsi::Float64)
    mgsi > 0 || error("mg_si must be positive")
    tot = 45.00 + 37.80
    si  = tot / (1 + mgsi * 40.304 / 60.084)
    return [si, 4.45, 3.55, tot - si, 8.05, 0.029, 0.36, 0.201, 0.0, 0.384, 0.0]
end

"Solid mantle adiabat, used only to set the STARTING temperature of the isentrope."
adiabat_T(Tp, P_kbar) = (Tp + 273.15) * exp(3.0e-5 * P_kbar * 1e8 / (3300.0 * 1200.0)) - 273.15

minim(data, X, P, T) = single_point_minimization(P, T, data, X=X, Xoxides=OX, sys_in="wt")

"""
Temperature at pressure `P` whose entropy matches `S0`, by secant iteration.

MAGEMin has no isentropic mode, so the isentrope has to be tracked by hand. Entropy is
monotonic in T at fixed P, which makes this well behaved; the step is clamped so a flat
secant near a phase boundary cannot throw the iterate across the melting interval.
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

"Isentropic decompression melting. Returns (MAGEMin output, final T) or (nothing, T) on failure."
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

"Melt composition in wt%, renormalised over the oxides the CIPW norm consumes, or nothing."
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

function calibrate(data)
    println("T_p scan at Earth's pyrolite (molar Mg/Si = 1.25), isentropic 3.0 -> 1.0 GPa")
    @printf("%-7s %-8s %-7s %-8s %-8s %-8s %-8s %s\n",
            "T_p", "T_end", "F", "SiO2", "MgO", "Al2O3", "CaO", "misfit*")
    X = mantle_composition(1.25)
    ref = (SiO2=48.76, MgO=11.27, Al2O3=17.05, CaO=11.91)   # PRIMELT at T_p = 1350
    for Tp in 1275.0:25.0:1425.0
        out, Tend = isentropic_melt(data, X, Tp)
        m, F = melt_oxides(out)
        if m === nothing
            @printf("%-7.0f %-8.1f no melt\n", Tp, Tend); continue
        end
        misfit = abs(m["SiO2"]-ref.SiO2) + abs(m["MgO"]-ref.MgO) +
                 abs(m["Al2O3"]-ref.Al2O3) + abs(m["CaO"]-ref.CaO)
        @printf("%-7.0f %-8.1f %-7.3f %-8.2f %-8.2f %-8.2f %-8.2f %.2f\n",
                Tp, Tend, F, m["SiO2"], m["MgO"], m["Al2O3"], m["CaO"], misfit)
    end
    println("\n*misfit = sum |melt - PRIMELT primary melt| over SiO2/MgO/Al2O3/CaO, wt%")
    println(" PRIMELT reference: SiO2 48.76  MgO 11.27  Al2O3 17.05  CaO 11.91")
end

"""
Potential temperature giving melt fraction `Ftarget`, by bisection.

F increases monotonically with T_p, so bisection is safe. Returns NaN if the target is
unreachable inside [TP_LO, TP_HI] -- a composition too refractory to reach Earth's melt
fraction even at 2200 C is a result, not something to extrapolate through.
"""
function Tp_for_F(data, X, Ftarget)
    lo, hi = TP_LO, TP_HI
    _, F_hi = melt_oxides(isentropic_melt(data, X, hi)[1])
    F_hi < Ftarget && return NaN
    while hi - lo > TP_TOL
        mid = (lo + hi) / 2
        _, F = melt_oxides(isentropic_melt(data, X, mid)[1])
        F < Ftarget ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

function sweep(data, _unused, outpath)
    @printf("MAGEMin isentropic melting at CONSTANT melt fraction F = %.3f, %d Mg/Si values, %.1f -> %.1f GPa\n",
            F_TARGET, length(MGSI_GRID), PSTART/10, PEND/10)
    println("  T_p is solved per composition, not held fixed -- see the F_TARGET note above.")
    flush(stdout)
    rows = String[]
    push!(rows, "T_p,mg_si,melt_fraction," * join(CSV_OXIDES, ","))
    nfail = 0
    for mgsi in MGSI_GRID
        X = mantle_composition(mgsi)
        local out, Tend, Tp
        try
            Tp = Tp_for_F(data, X, F_TARGET)
            if isnan(Tp)
                @printf("  Mg/Si=%.2f  cannot reach F=%.3f below %.0f C\n", mgsi, F_TARGET, TP_HI)
                nfail += 1; continue
            end
            out, Tend = isentropic_melt(data, X, Tp)
        catch e
            @printf("  Mg/Si=%.2f  FAILED (%s)\n", mgsi, typeof(e)); nfail += 1; continue
        end
        m, F = melt_oxides(out)
        if m === nothing
            @printf("  Mg/Si=%.2f  no melt at the end of the path\n", mgsi); nfail += 1; continue
        end
        push!(rows, @sprintf("%.0f,%.4g,%.6g,", Tp, mgsi, F) *
              join([@sprintf("%.6g", m[o]) for o in CSV_OXIDES], ","))
        @printf("  Mg/Si=%.2f T_p=%.0f (%+.0f) T_end=%.0f F=%.3f SiO2=%.2f MgO=%.2f Al2O3=%.2f CaO=%.2f  %s\n",
                mgsi, Tp, Tp - TP_EARTH, Tend, F, m["SiO2"], m["MgO"], m["Al2O3"], m["CaO"],
                join(out.ph, ","))
        flush(stdout)
    end
    open(outpath, "w") do io
        for r in rows; println(io, r); end
    end
    @printf("\nWritten: %s  (%d rows", outpath, length(rows)-1)
    nfail > 0 ? @printf(", %d failed)\n", nfail) : println(")")
end

function main()
    args = ARGS
    docalib = "--calibrate" in args
    Tp = TP_EARTH
    i = findfirst(==("--tp"), args);  i !== nothing && (Tp = parse(Float64, args[i+1]))
    outpath = joinpath(@__DIR__, "crust_compositions.csv")
    i = findfirst(==("--out"), args); i !== nothing && (outpath = args[i+1])

    data = Initialize_MAGEMin("ig", verbose=false)
    try
        docalib ? calibrate(data) : sweep(data, Tp, outpath)
    finally
        Finalize_MAGEMin(data)
    end
end

main()
