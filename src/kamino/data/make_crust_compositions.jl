#=
make_crust_compositions.jl -- primary-melt compositions on a (Mg/Si, dIW) grid, via MAGEMin.

Writes crust_compositions.csv, the table crust_composition.oxide_composition interpolates.
Isobaric batch melting at P_MELT, solving temperature for a fixed melt fraction; the liquid is
the crust composition.

    julia make_crust_compositions.jl                       # full grid
    julia make_crust_compositions.jl --points "1.25,-2.0"  # spot checks
    julia make_crust_compositions.jl --points "1.25,-2.0" --ftarget 0.12
    julia make_crust_compositions.jl --out other.csv

Batch melting is path-independent -- the liquid depends only on (bulk, P, T) -- so an isobaric
solve at the segregation pressure gives the same melt as tracking an isentrope down to it.

MAGEMin: Riel, Kaus, Green & Berlie (2022), G3 23, e2022GC010427.
Dataset: Holland, Green & Powell (2018), J. Petrol. 59, 881.
Pyrolite: McDonough & Sun (1995), Chem. Geol. 120, 223.
F = 0.20 and the Mg/Si construction: Guimond et al. (2024), RiMG 90, 259.
=#

using MAGEMin_C, Printf

const OX = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","Cr2O3","H2O"]
const CSV_OXIDES = ["SiO2","TiO2","Al2O3","Cr2O3","FeO","MgO","CaO","Na2O","K2O"]

# Melt segregation pressure (kbar), a representative mean depth beneath a ridge.
const P_MELT = 10.0
# Where clinopyroxene leaves the assemblage for Earth's mantle, beyond which melting is much
# less productive (Guimond et al. 2024, after Katz et al. 2003).
# Not const: --ftarget overrides it, for the F-sensitivity question the Earth-vs-MORB offset
# raises (development history 32.10). MAGEMin dominates the runtime, so the global costs nothing.
F_TARGET = 0.20
const T_LO, T_HI, T_TOL = 900.0, 2200.0, 1.0
const ULTRACALCIC_WARN = 1.2
const F_TOL = 0.05

# Non-uniform, dense where the assemblage changes: the source goes quartz-normative below the
# line (Mg/Si, dIW) = (0.5, -1.05) .. (0.8, -2.82), and opx-out and the nepheline/ferropericlase
# switch sit near Mg/Si 1.6. Earth's anchors (1.25, -2.0) must stay on the grid.
const MGSI_GRID = sort(unique([0.5:0.05:0.9; 1.0:0.1:1.2; 1.25; 1.3:0.05:1.8; 1.9; 2.0]))
const DIW_GRID = sort(unique([-5.0:0.5:-3.0; -2.9:0.1:-1.0]))

# --- Mantle construction -----------------------------------------------------------------------

# Must stay identical to PYROLITE_NON_FE in crust_composition.py; iron is set by the dIW axis.
const PYROLITE_NON_FE = Dict("SiO2"=>45.00, "Al2O3"=>4.45, "CaO"=>3.55, "MgO"=>37.80,
                             "K2O"=>0.029, "Na2O"=>0.36, "TiO2"=>0.201, "Cr2O3"=>0.384)
const M_OX = Dict("SiO2"=>60.084, "Al2O3"=>101.961, "CaO"=>56.077, "MgO"=>40.304, "FeO"=>71.844,
                  "K2O"=>94.196, "Na2O"=>61.979, "TiO2"=>79.866, "Cr2O3"=>151.990)
const NCAT = Dict("SiO2"=>1, "Al2O3"=>2, "CaO"=>1, "MgO"=>1,
                  "K2O"=>2, "Na2O"=>2, "TiO2"=>1, "Cr2O3"=>2)

const PYRO_SUM = sum(values(PYROLITE_NON_FE))
const K_CATIONS = sum(PYROLITE_NON_FE[o] / PYRO_SUM * NCAT[o] / M_OX[o] for o in keys(PYROLITE_NON_FE))
const EARTH_MANTLE_FEO, DIW_EARTH = 8.05, -2.0

x_feo_from_wt(feo) = (feo / M_OX["FeO"]) / (feo / M_OX["FeO"] + (100.0 - feo) * K_CATIONS)

# Calibrated so dIW = -2 reproduces BSE FeO exactly; mirrors _FEO_ACTIVITY_CONST in
# crust_composition.py. Only the Earth anchor is meaningful -- the absolute value is a label.
const FEO_ACTIVITY_CONST = x_feo_from_wt(EARTH_MANTLE_FEO) / 10.0^(DIW_EARTH / 2)

"Mantle FeO (wt%) from core-formation fO2, inverting Fe + 1/2 O2 = FeO. Logarithmic in dIW."
function feo_from_delta_iw(diw::Float64)
    x = FEO_ACTIVITY_CONST * 10.0^(diw / 2)
    0.0 < x < 1.0 || error("delta_iw=$diw implies X_FeO=$x, outside (0,1)")
    return 100.0 * x * K_CATIONS / ((1.0 - x) / M_OX["FeO"] + x * K_CATIONS)
end

"Pyrolite with FeO set from dIW, then MgO/SiO2 re-split to the target molar Mg/Si."
function mantle_composition(mgsi::Float64, feo_wt::Float64)
    mgsi > 0 || error("mg_si must be positive")
    0.0 <= feo_wt < 100.0 || error("feo_wt out of range")
    # Iron first, then Mg/Si within what remains: applied the other way the axes would not be
    # orthogonal, because iron removed to the core rescales the whole silicate budget.
    scale = (100.0 - feo_wt) / PYRO_SUM
    ox = Dict(o => v * scale for (o, v) in PYROLITE_NON_FE)
    tot = ox["MgO"] + ox["SiO2"]
    si = tot / (1 + mgsi * M_OX["MgO"] / M_OX["SiO2"])
    ox["SiO2"], ox["MgO"] = si, tot - si
    return [ox["SiO2"], ox["Al2O3"], ox["CaO"], ox["MgO"], feo_wt,   # order must match OX
            ox["K2O"], ox["Na2O"], ox["TiO2"], 0.0, ox["Cr2O3"], 0.0]
end

# Ferric iron is off (the O component is 0.0): Fe3+ is second-order for melt major elements and
# enabling it would perturb the Earth anchor every downstream calibration rests on.

# --- Melting -----------------------------------------------------------------------------------

minim(data, X, T) = single_point_minimization(P_MELT, T, data, X=X, Xoxides=OX, sys_in="wt")

"Liquid composition in wt% over CSV_OXIDES, plus melt fraction. (nothing, 0) if there is no melt."
function melt_oxides(out)
    findfirst(==("liq"), out.ph) === nothing && return nothing, 0.0
    frac = out.frac_M_wt
    (frac === nothing || isnan(frac) || frac <= 0) && return nothing, 0.0
    comp = Dict(OX[j] => out.bulk_M_wt[j] * 100 for j in eachindex(OX))
    kept = Dict(o => get(comp, o, 0.0) for o in CSV_OXIDES)
    tot = sum(values(kept))
    tot <= 0 && return nothing, frac
    return Dict(o => v / tot * 100 for (o, v) in kept), frac
end

residue(out) = join(filter(!=("liq"), out.ph), ";")

"Temperature giving F_TARGET at P_MELT, by bisection. NaN if unreachable below T_HI."
function T_for_F(data, X)
    melt_oxides(minim(data, X, T_HI))[2] < F_TARGET && return NaN
    lo, hi = T_LO, T_HI
    while hi - lo > T_TOL
        mid = (lo + hi) / 2
        melt_oxides(minim(data, X, mid))[2] < F_TARGET ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# --- Grid --------------------------------------------------------------------------------------

const CSV_HEADER = "mg_si,delta_iw,mantle_feo,T_melt,melt_fraction,residual_phases," *
                   join(CSV_OXIDES, ",") * ",CaO_Al2O3"

"One grid point. Returns (csv_row | nothing, status)."
function grid_point(data, mgsi, diw)
    feo = feo_from_delta_iw(diw)
    X = mantle_composition(mgsi, feo)
    T = T_for_F(data, X)
    isnan(T) && return nothing, "cannot reach F=$F_TARGET below $T_HI C"

    out = minim(data, X, T)
    m, F = melt_oxides(out)
    m === nothing && return nothing, "no melt at $T C"

    ca_al = m["Al2O3"] > 0 ? m["CaO"] / m["Al2O3"] : NaN
    warns = String[]
    # F = 0.20 melts past cpx-out at high Mg/Si, which is what makes those melts ultracalcic
    # (Medard et al. 2004 exclude them from the fertile lherzolite source this assumes).
    ca_al > ULTRACALCIC_WARN && push!(warns, @sprintf("ultracalcic CaO/Al2O3 %.2f", ca_al))
    # Catches a bisection that converged on a bracket end rather than on the target.
    abs(F - F_TARGET) > F_TOL && push!(warns, @sprintf("F=%.3f off target", F))

    row = @sprintf("%.4g,%.4g,%.4f,%.0f,%.6g,%s,", mgsi, diw, feo, T, F, residue(out)) *
          join([@sprintf("%.6g", m[o]) for o in CSV_OXIDES], ",") *
          @sprintf(",%.4f", ca_al)
    return row, (isempty(warns) ? "ok" : join(warns, "; "))
end

function sweep(data, outpath)
    @printf("MAGEMin: %d Mg/Si x %d dIW = %d points, isobaric %.1f GPa, F = %.2f\n",
            length(MGSI_GRID), length(DIW_GRID), length(MGSI_GRID)*length(DIW_GRID),
            P_MELT/10, F_TARGET)
    flush(stdout)

    rows, failures = String[CSV_HEADER], String[]
    t0 = time()
    for mgsi in MGSI_GRID, diw in DIW_GRID
        row, status = try
            grid_point(data, mgsi, diw)
        catch e
            (nothing, "EXCEPTION $(typeof(e))")
        end
        if row === nothing
            @printf("  Mg/Si=%.2f dIW=%+.1f  SKIPPED: %s\n", mgsi, diw, status)
            push!(failures, @sprintf("Mg/Si=%.2f dIW=%+.1f: %s", mgsi, diw, status))
        else
            f = split(row, ",")
            @printf("  Mg/Si=%.2f dIW=%+.1f FeO=%5.2f T=%-5s F=%.3f SiO2=%5.2f Ca/Al=%.2f %s\n",
                    mgsi, diw, parse(Float64, f[3]), f[4], parse(Float64, f[5]),
                    parse(Float64, f[7]), parse(Float64, f[16]), status == "ok" ? "" : "[$status]")
            push!(rows, row)
        end
        flush(stdout)
    end

    open(outpath, "w") do io
        println(io, "# Generated by make_crust_compositions.jl -- do not edit by hand.")
        println(io, "# isobaric $(P_MELT/10) GPa, F_TARGET=$F_TARGET, batch melting")
        println(io, "# delta_iw is the CORE-FORMATION fO2, which sets mantle FeO.")
        println(io, "# T_melt is the melting temperature at that pressure, NOT a potential temperature.")
        for r in rows; println(io, r); end
    end

    @printf("\nWritten: %s  (%d rows, %d failed) in %.1f min\n",
            outpath, length(rows)-1, length(failures), (time()-t0)/60)
    if !isempty(failures)
        # load_crust_table rejects a table that is not a full rectangle, and this is why.
        println("FAILURES (the table is NOT a complete grid and will be rejected downstream):")
        for f in failures; println("  ", f); end
    end
end

function points(data, spec)
    println(CSV_HEADER); flush(stdout)
    for pair in split(spec, ";")
        mgsi, diw = parse.(Float64, split(strip(pair), ","))
        t = time()
        row, status = try
            grid_point(data, mgsi, diw)
        catch e
            (nothing, "EXCEPTION $(typeof(e)): $e")
        end
        if row === nothing
            @printf("# Mg/Si=%.2f dIW=%+.1f FAILED (%.0f s): %s\n", mgsi, diw, time()-t, status)
        else
            println(row)
            @printf("# ^ %.0f s  %s\n", time()-t, status)
        end
        flush(stdout)
    end
end

function main()
    outpath = joinpath(@__DIR__, "crust_compositions.csv")
    i = findfirst(==("--out"), ARGS); i !== nothing && (outpath = ARGS[i+1])
    k = findfirst(==("--ftarget"), ARGS)
    k !== nothing && (global F_TARGET = parse(Float64, ARGS[k+1]))
    j = findfirst(==("--points"), ARGS)

    data = Initialize_MAGEMin("ig", verbose=false)
    try
        j === nothing ? sweep(data, outpath) : points(data, ARGS[j+1])
    finally
        Finalize_MAGEMin(data)
    end
end

main()
