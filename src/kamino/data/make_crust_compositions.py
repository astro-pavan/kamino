"""make_crust_compositions.py — primary-melt compositions on a (T_p, Mg/Si) grid via pMELTS.

Generates the lookup table that `crust_composition.oxide_composition` interpolates, replacing
the PRIMELT1 spreadsheet + SiO2-rescale proxy. The proxy scaled SiO2 only and renormalised,
which inflated Al2O3 (14 -> 20 wt%) and CaO (10 -> 14 wt%) as pure artifacts and pushed the
endpoints outside primary-melt space (38-58 wt% SiO2). Here both axes come from the same
thermodynamic calculation instead.

Method
------
Adiabatic decompression melting with pMELTS (Ghiorso et al. 2002), the MELTS calibration for
peridotite melting at 1-3 GPa:

  1. bulk mantle = McDonough & Sun (1995) pyrolite, with MgO/SiO2 re-split to the target molar
     Mg/Si at fixed (MgO + SiO2) mass, leaving Al/Ca/Fe/Na untouched -- so Mg/Si is a genuinely
     orthogonal axis, which it is not in the proxy;
  2. start sub-solidus at P_START on the solid adiabat for the given potential temperature;
  3. decompress isentropically (runMode 3) in P_STEP increments to P_END;
  4. batch melting -- melt stays with the residue, so the final liquid IS the pooled primary
     melt, and its composition is what the CIPW norm then converts to mineralogy.

Environment
-----------
Needs alphaMELTS for Python and libgsl.so.27. See ALPHAMELTS_DIR / the GSL note below; run
with `--check` to verify the environment before committing to a full grid.

  DO NOT symlink the system libgsl.so.25 (GSL 2.6) as libgsl.so.27. All 38 gsl_ symbols
  alphaMELTS needs are present, so it loads and runs -- then corrupts the heap
  (realloc(): invalid pointer / SIGABRT). Build GSL 2.7.1 from source; it takes ~1 minute:
      curl -O https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz && tar xzf gsl-2.7.1.tar.gz
      cd gsl-2.7.1 && ./configure --prefix=$PREFIX --disable-static && make -j8 && make install
      export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH

Usage
-----
    python make_crust_compositions.py --check
    python make_crust_compositions.py --workers 26
    python make_crust_compositions.py --tp 1300 1700 9 --mgsi 0.9 1.6 8

Citations for the paper
-----------------------
  * pMELTS: Ghiorso, Hirschmann, Reiners & Kress (2002), G3 3(5), 1030.
  * alphaMELTS for Python: Antoshechkina & Ghiorso (2018), AGU Fall Meeting.
  * Pyrolite bulk mantle: McDonough & Sun (1995), Chem. Geol. 120, 223.
"""

import argparse
import csv
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_CSV = os.path.join(DATA_DIR, "crust_compositions.csv")

ALPHAMELTS_DIR = os.environ.get(
    "ALPHAMELTS_DIR",
    "/data/pt426/alphamelts/alphamelts_py-2.3.1/alphamelts-py-2.3.1-linux",
)

# Where to find libgsl.so.27, tried in order. Set GSL_LIB_DIR to override.
# NB /data/pt426/melts-deps/lib holds a libgsl.so.27 that CANNOT work on this machine -- it is
# built against GLIBC 2.35 and RHEL 9 has 2.34 -- so it is deliberately not in this list.
GSL_CANDIDATES = [
    os.environ.get("GSL_LIB_DIR", ""),
    "/data/pt426/gsl-2.7.1/lib",
]

_GSL_SONAME = "libgsl.so.27"
# libgsl.so declares no DT_NEEDED on a CBLAS -- GSL expects the application to supply one at
# link time. When the loader pulls libgsl in as a dependency of libalphamelts.so it resolves
# CBLAS from the same directory, but preloading libgsl on its own leaves cblas_* undefined
# ("undefined symbol: cblas_ctrmv"), so it has to be loaded first and RTLD_GLOBAL.
_CBLAS_SONAME = "libgslcblas.so.0"


def preload_gsl():
    """Load libgsl.so.27 into the process before alphaMELTS is imported.

    libalphamelts.so links against it, and the dynamic loader resolves that dependency from
    LD_LIBRARY_PATH -- which cannot be set from inside a running process, because the loader
    caches it at startup. Loading the library ourselves with RTLD_GLOBAL puts it in the global
    symbol scope under its soname, so the later dlopen of libalphamelts.so finds it already
    resident and never searches the path. That makes this script work as plain
    `python make_crust_compositions.py`, with no environment setup by the caller.

    Returns the path loaded, or raises RuntimeError naming what to do about it.
    """
    import ctypes
    for directory in GSL_CANDIDATES:
        if not directory:
            continue
        path = os.path.join(directory, _GSL_SONAME)
        if not os.path.isfile(path):
            continue
        try:
            cblas = os.path.join(directory, _CBLAS_SONAME)
            if os.path.isfile(cblas):
                ctypes.CDLL(cblas, mode=ctypes.RTLD_GLOBAL)
            ctypes.CDLL(path, mode=ctypes.RTLD_GLOBAL)
            return path
        except OSError as exc:
            # Present but unloadable, e.g. built against a newer glibc than this machine has.
            last = f"{path}: {exc}"
            break
    else:
        last = f"{_GSL_SONAME} not found in {[d for d in GSL_CANDIDATES if d]}"
    try:
        import ctypes as _c
        _c.CDLL(_GSL_SONAME, mode=_c.RTLD_GLOBAL)   # already on the loader path?
        return _GSL_SONAME
    except OSError:
        pass
    raise RuntimeError(
        f"Could not load {_GSL_SONAME} ({last}).\n"
        f"Build GSL 2.7.1 (~1 minute, no dependencies) and point GSL_LIB_DIR at its lib/:\n"
        f"    curl -O https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz && tar xzf gsl-2.7.1.tar.gz\n"
        f"    cd gsl-2.7.1 && ./configure --prefix=$PREFIX --disable-static && make -j8 "
        f"&& make install\n"
        f"Do NOT symlink the system GSL 2.6 as {_GSL_SONAME}: every symbol resolves, it runs, "
        f"and then it corrupts the heap."
    )

PMELTS_MODE = 2  # 2 = pMELTS; the mode can only be set ONCE per process (MELTS limitation)

# Oxygen buffer. Without one, fO2 is unconstrained, pMELTS drives the Fe2O3 liquid component
# to a negative mole fraction and the melt comes back at ~36 wt% SiO2 -- not a silicate melt.
# FMQ is the usual reference for MORB-source mantle.
FO2_BUFFER = "FMQ"

# Melting path, in bars. 3 GPa is sub-solidus across this T_p range. P_END is the melt
# SEGREGATION pressure, not the base of the crust: this is batch melting, so the liquid keeps
# re-equilibrating with its residue at every step, and carrying that to 0.2 GPa yields
# 57 wt% SiO2 at F = 0.29 -- an over-equilibrated andesite, not a primary melt. Stopping at
# 1 GPa, a representative mean segregation depth beneath a ridge, gives 49 wt% SiO2 at
# F = 0.12 for Earth-like T_p, which is MORB-like.
P_START, P_END, P_STEP = 30000.0, 10000.0, 1000.0

# Solid adiabat T(P) = T_p * exp(alpha*P / (rho*Cp)), i.e. ~12 K/GPa over this range.
THERMAL_EXPANSIVITY = 3.0e-5   # 1/K
MANTLE_DENSITY = 3300.0        # kg/m3
HEAT_CAPACITY = 1200.0         # J/kg/K

ABSOLUTE_ZERO = 273.15

# McDonough & Sun (1995) primitive mantle, wt%. Molar Mg/Si = 1.25.
PYROLITE = {
    "SiO2": 45.00, "TiO2": 0.201, "Al2O3": 4.45, "Cr2O3": 0.384, "FeO": 8.05,
    "MnO": 0.135, "MgO": 37.80, "CaO": 3.55, "Na2O": 0.36, "K2O": 0.029, "P2O5": 0.021,
}

MOLAR_MASS = {"SiO2": 60.084, "MgO": 40.304}

TP_DEFAULT = [1300.0, 1350.0, 1400.0, 1450.0, 1500.0, 1550.0]

# Non-uniform on purpose. Below ~1.3 the normative mineralogy barely moves; above it Nepheline
# switches on and Albite collapses between 1.4 and 1.6, so that is where the resolution goes.
MGSI_DEFAULT = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60]

# pMELTS fails outright above Mg/Si ~1.6 (the start state will not converge at any nudge), which
# truncates the Mg-rich tail of the stellar Mg/Si distribution. Below ~0.8 it converges but the
# melts leave basalt space -- 69 wt% SiO2 (a rhyolite) at Mg/Si 0.5 -- and pMELTS is calibrated
# for PERIDOTITE melting, which a 62 wt% SiO2 bulk is not. Those rows are flagged, not dropped.
SILICIC_MELT_WARN = 57.0   # wt% SiO2; above this the melt is no longer basaltic

# A melt fraction this high is not a crust: the "melt" is essentially the bulk mantle, so the
# CIPW norm would return a peridotite. pMELTS gives F = 1.00 at T_p >= 1650 with this path.
MELT_FRACTION_WARN = 0.5

# Oxides carried through to the CSV, in the order cipw_norm expects to find them.
CSV_OXIDES = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO", "CaO",
              "Na2O", "K2O", "P2O5"]

# alphaMELTS' internal bulk-composition slot order, used to decode dispComposition.
MELTS_OXIDE_ORDER = ["sio2", "tio2", "al2o3", "fe2o3", "cr2o3", "feo", "mno", "mgo", "nio",
                     "coo", "cao", "na2o", "k2o", "p2o5", "h2o", "co2", "so3",
                     "cl2o-1", "f2o-1"]


def mantle_composition(mg_si):
    """Pyrolite with MgO/SiO2 re-split to the target molar Mg/Si, total mass conserved.

    Only the Mg-Si pair moves, so Al2O3/CaO/FeO/Na2O stay at their pyrolite values and the
    axis is orthogonal to everything else -- unlike the SiO2-rescale-plus-renormalise proxy,
    where changing Mg/Si silently enriched Al and Ca.
    """
    if mg_si <= 0:
        raise ValueError("mg_si must be positive")
    total = PYROLITE["SiO2"] + PYROLITE["MgO"]
    si_mass = total / (1.0 + mg_si * MOLAR_MASS["MgO"] / MOLAR_MASS["SiO2"])
    comp = dict(PYROLITE)
    comp["SiO2"] = si_mass
    comp["MgO"] = total - si_mass
    return comp


def molar_mg_si(comp):
    """Molar Mg/Si of a wt% oxide dict, for verifying the re-split hit its target."""
    return (comp["MgO"] / MOLAR_MASS["MgO"]) / (comp["SiO2"] / MOLAR_MASS["SiO2"])


def adiabat_temperature(T_p_celsius, pressure_bar):
    """Temperature (degC) on the solid mantle adiabat at `pressure_bar`, from potential temperature.

    An approximation: MELTS has no phase-suppression hook, so the isentrope cannot be run
    sub-solidus to define T_p self-consistently. Over 0-3 GPa this gives ~12 K/GPa, within the
    usual 10-16 K/GPa range for solid peridotite, and it only sets where the melting path
    STARTS -- the melt compositions themselves come from the isentropic calculation.
    """
    T_p_kelvin = T_p_celsius + ABSOLUTE_ZERO
    pressure_pa = pressure_bar * 1.0e5
    factor = THERMAL_EXPANSIVITY * pressure_pa / (MANTLE_DENSITY * HEAT_CAPACITY)
    import math
    return T_p_kelvin * math.exp(factor) - ABSOLUTE_ZERO


def _melt_one(args):
    """Run one (T_p, Mg/Si) decompression path. Executed in its own process -- see below.

    A failed calcEquilibriumState leaves the MELTS library unrecoverable ("Could not
    re-initialize MELTS library after failure"), poisoning every later calculation in the
    same process. One process per grid point makes each failure local.
    """
    T_p, mg_si = args[0], args[1]
    t_nudge = args[2] if len(args) > 2 else 0.0
    preload_gsl()
    sys.path.insert(0, ALPHAMELTS_DIR)

    # alphaMELTS writes *_tbl.txt output tables into the CWD, so every worker would otherwise
    # race on the same filenames (and litter wherever the script was launched from).
    import tempfile
    workdir = tempfile.mkdtemp(prefix="melts_")
    os.chdir(workdir)

    from meltsdynamic import MELTSdynamic

    comp = mantle_composition(mg_si)
    node = MELTSdynamic(PMELTS_MODE)
    node.engine.setSystemProperties("Log fO2 Path", FO2_BUFFER)
    for oxide, value in comp.items():
        node.engine.setBulkComposition(oxide, value)

    # Establish the starting state isothermally; this also fixes the entropy the isentropic
    # decompression then holds.
    node.engine.pressure = P_START
    node.engine.temperature = adiabat_temperature(T_p, P_START) + t_nudge
    node.engine.calcEquilibriumState(1, 0)
    if node.engine.status.failed:
        return T_p, mg_si, float("nan"), None, f"start state: {node.engine.status.message}"

    melt_fraction, melt = 0.0, None
    pressure = P_START
    while pressure - P_STEP >= P_END - 1e-9:
        pressure -= P_STEP
        node = node.addNodeAfter()
        node.engine.pressure = pressure
        try:
            node.engine.calcEquilibriumState(3, 0)   # 3 = isentropic
        except Exception as exc:
            return T_p, mg_si, float("nan"), None, f"{type(exc).__name__} at {pressure:.0f} bar"
        # calcEquilibriumState does NOT raise on failure -- it prints and leaves the library
        # unrecoverable, so the flag has to be read or the run silently degrades.
        if node.engine.status.failed:
            return (T_p, mg_si, float("nan"), None,
                    f"{node.engine.status.message} at {pressure:.0f} bar")
        if "liquid1" not in (node.engine.liquidNames or []):
            continue
        try:
            mass = node.engine.getProperty("mass", ["liquid1"])
            values = node.engine.getProperty("dispComposition", ["liquid1"])
        except Exception:
            continue
        if mass is None or mass != mass or mass <= 0.0:
            continue
        melt_fraction = mass / 100.0
        melt = {ox: float(values[MELTS_OXIDE_ORDER.index(ox.lower())]) for ox in CSV_OXIDES}

    if melt is None:
        return T_p, mg_si, 0.0, None, "no melt produced along the path"

    # Fe2O3 is not carried (the CIPW norm routes all iron to fayalite), so renormalise the
    # oxides that are, keeping the reported composition summing to 100 wt%.
    total = sum(melt.values())
    if total <= 0:
        return T_p, mg_si, melt_fraction, None, "melt composition summed to zero"
    melt = {ox: value / total * 100.0 for ox, value in melt.items()}
    return T_p, mg_si, melt_fraction, melt, None


def check_environment():
    """Verify alphaMELTS + GSL load and that pMELTS produces a sane pyrolite melt."""
    if not os.path.isdir(ALPHAMELTS_DIR):
        print(f"FAIL: ALPHAMELTS_DIR not found: {ALPHAMELTS_DIR}")
        return False
    print(f"alphaMELTS dir : {ALPHAMELTS_DIR}")
    try:
        print(f"GSL            : {preload_gsl()}")
    except RuntimeError as exc:
        print(f"FAIL: {exc}")
        return False
    try:
        sys.path.insert(0, ALPHAMELTS_DIR)
        from meltsdynamic import MELTSdynamic  # noqa: F401
    except OSError as exc:
        print(f"FAIL: cannot load the MELTS library: {exc}")
        return False

    T_p, mg_si = 1350.0, 1.25
    _, _, F, melt, error = _melt_one((T_p, mg_si))
    if error or melt is None:
        print(f"FAIL: reference calculation failed: {error}")
        return False
    print(f"reference run  : T_p={T_p:.0f}C Mg/Si={mg_si:.2f} -> F={F:.3f}, "
          f"SiO2={melt['SiO2']:.2f} MgO={melt['MgO']:.2f} Al2O3={melt['Al2O3']:.2f} "
          f"CaO={melt['CaO']:.2f} wt%")
    if not (40.0 <= melt["SiO2"] <= 56.0):
        print("WARN: melt SiO2 is outside the primary-melt range (40-56 wt%).")
    print("Environment OK.")
    return True


def _linspace(lo, hi, n):
    if n < 2:
        return [float(lo)]
    step = (float(hi) - float(lo)) / (n - 1)
    return [float(lo) + i * step for i in range(n)]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="Verify the environment with one reference run, then exit.")
    parser.add_argument("--tp", nargs=3, metavar=("LO", "HI", "N"), default=None,
                        help=f"Uniform T_p grid in degC. Default: {TP_DEFAULT}.")
    parser.add_argument("--tp-values", nargs="+", type=float, default=None,
                        help="Explicit T_p values, overriding --tp.")
    parser.add_argument("--mgsi", nargs=3, metavar=("LO", "HI", "N"), default=None,
                        help=f"Uniform Mg/Si grid. Default is non-uniform: {MGSI_DEFAULT}.")
    parser.add_argument("--mgsi-values", nargs="+", type=float, default=None,
                        help="Explicit Mg/Si values, overriding --mgsi.")
    parser.add_argument("--workers", type=int, default=8, help="Parallel processes (default 8).")
    parser.add_argument("--output", default=OUTPUT_CSV, help=f"CSV path (default {OUTPUT_CSV}).")
    args = parser.parse_args()

    if args.check:
        sys.exit(0 if check_environment() else 1)

    # Fail in the parent with one clear message rather than in every worker with an OSError
    # raised from inside ctypes, four frames deep in a concurrent.futures traceback.
    try:
        preload_gsl()
    except RuntimeError as exc:
        print(f"FAIL: {exc}")
        sys.exit(1)

    if args.tp_values:
        tp_values = sorted(float(v) for v in args.tp_values)
    elif args.tp:
        tp_values = _linspace(args.tp[0], args.tp[1], int(args.tp[2]))
    else:
        tp_values = list(TP_DEFAULT)
    if args.mgsi_values:
        mgsi_values = sorted(float(v) for v in args.mgsi_values)
    elif args.mgsi:
        mgsi_values = _linspace(args.mgsi[0], args.mgsi[1], int(args.mgsi[2]))
    else:
        mgsi_values = list(MGSI_DEFAULT)
    grid = [(tp, ms) for tp in tp_values for ms in mgsi_values]

    print(f"pMELTS decompression melting: {len(tp_values)} T_p x {len(mgsi_values)} Mg/Si "
          f"= {len(grid)} points, {args.workers} workers")
    print(f"  T_p   : {tp_values[0]:.0f} - {tp_values[-1]:.0f} degC")
    print(f"  Mg/Si : {mgsi_values[0]:.2f} - {mgsi_values[-1]:.2f}")
    print(f"  path  : {P_START/1e4:.1f} -> {P_END/1e4:.1f} GPa in {P_STEP:.0f} bar steps")
    print(f"  output: {args.output}")

    rows, failures = [], []
    # max_tasks_per_child=1 is REQUIRED, not a tuning knob. ProcessPoolExecutor reuses its
    # workers, but the MELTS calculation mode can only be set once per process and a failed
    # calcEquilibriumState leaves the library unrecoverable -- so every task after the first in
    # a given worker fails. Symptom: successes == max_workers exactly, regardless of grid size
    # (a 72-point grid on 8 workers returned 8 rows and 64 failures).
    with ProcessPoolExecutor(max_workers=args.workers,
                             mp_context=mp.get_context("spawn"),
                             max_tasks_per_child=1) as executor:
        futures = {executor.submit(_melt_one, point): point for point in grid}
        for done, future in enumerate(as_completed(futures), start=1):
            T_p, mg_si, F, melt, error = future.result()
            if melt is None:
                failures.append((T_p, mg_si, error))
                print(f"[{done}/{len(grid)}] FAILED T_p={T_p:.0f} Mg/Si={mg_si:.2f}: {error}")
                continue
            row = {"T_p": T_p, "mg_si": mg_si, "melt_fraction": F}
            row.update(melt)
            rows.append(row)
            flag = "  <-- F high, not a crust" if F > MELT_FRACTION_WARN else ""
            if melt["SiO2"] > SILICIC_MELT_WARN:
                flag += "  <-- silicic, outside pMELTS peridotite calibration"
            print(f"[{done}/{len(grid)}] T_p={T_p:.0f} Mg/Si={mg_si:.2f} F={F:.3f} "
                  f"SiO2={melt['SiO2']:.2f} MgO={melt['MgO']:.2f} "
                  f"Al2O3={melt['Al2O3']:.2f}{flag}")

    for nudge in (5.0, -5.0, 12.0):
        if not failures:
            break
        retry = [(T_p, mg_si, nudge) for T_p, mg_si, _ in failures]
        print(f"\nRetrying {len(retry)} failed point(s) with a {nudge:+.0f} K start nudge...")
        failures = []
        with ProcessPoolExecutor(max_workers=args.workers,
                                 mp_context=mp.get_context("spawn"),
                                 max_tasks_per_child=1) as executor:
            futures = {executor.submit(_melt_one, point): point for point in retry}
            for future in as_completed(futures):
                T_p, mg_si, F, melt, error = future.result()
                if melt is None:
                    failures.append((T_p, mg_si, error))
                    continue
                row = {"T_p": T_p, "mg_si": mg_si, "melt_fraction": F}
                row.update(melt)
                rows.append(row)
                print(f"  recovered T_p={T_p:.0f} Mg/Si={mg_si:.2f} F={F:.3f} "
                      f"SiO2={melt['SiO2']:.2f}")

    if not rows:
        print("No grid points succeeded; nothing written.")
        sys.exit(1)

    rows.sort(key=lambda r: (r["T_p"], r["mg_si"]))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["T_p", "mg_si", "melt_fraction"] + CSV_OXIDES)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                             for k, v in row.items()})

    print(f"\nWritten: {args.output}  ({len(rows)} rows)")
    silicic = [r for r in rows if r["SiO2"] > SILICIC_MELT_WARN]
    if silicic:
        print(f"  {len(silicic)} row(s) have melt SiO2 > {SILICIC_MELT_WARN:.0f} wt% -- not "
              f"basaltic, and pMELTS is extrapolating outside its peridotite calibration there. "
              f"Treat the low-Mg/Si end as indicative only.")
    hot = [r for r in rows if r["melt_fraction"] > MELT_FRACTION_WARN]
    if hot:
        print(f"  {len(hot)} row(s) have melt fraction > {MELT_FRACTION_WARN:.2f} -- at those T_p "
              f"the mantle is essentially all molten, so the 'melt' is the bulk composition and "
              f"the norm returns a peridotite, not a crust. Lower the T_p ceiling.")
    if failures:
        print(f"  {len(failures)} grid point(s) failed and were omitted:")
        for T_p, mg_si, error in failures[:10]:
            print(f"    T_p={T_p:.0f} Mg/Si={mg_si:.2f}: {error}")
        print("  Interpolation over a grid with holes is silently wrong -- fill or shrink "
              "the grid before using the table.")


if __name__ == "__main__":
    main()
