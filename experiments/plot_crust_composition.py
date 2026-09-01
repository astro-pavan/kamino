"""Normative crust mineralogy as a function of mantle Mg/Si.

Two figures:
  crust_mineralogy_mgsi.png  — stacked normative mineralogy vs mantle Mg/Si, comparing the old
                               SiO2-rescale proxy with the MAGEMin melts.
  crust_reactivity_mgsi.png  — the potential temperature the constant-F closure requires, and
                               what the resulting mineralogy delivers: cation charge per kg of
                               rock, and the charge-weighted and Na dissolution rate constants.

Reads the CSV written by src/kamino/data/make_crust_compositions.jl. That table is generated at
CONSTANT MELT FRACTION, not constant T_p — a mantle that cannot melt cannot transport heat by
magmatism, so it warms until it does — which means T_p VARIES along the Mg/Si axis and is itself
a result, not a control. Mg/Si is the only independent variable here, so nothing is faceted by
T_p (an earlier version did, and produced one single-point panel per composition).

Usage:
    python plot_crust_composition.py --csv <crust_compositions.csv> [--outdir output]
"""

import argparse
import csv
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

from kamino.crust_composition import cipw_norm, oxide_composition
from kamino.mineral_info import MINERAL_MOLAR_MASS
from kamino.chemistry import stoichiometry, elements, ION_CHARGE, get_k, al_idx, na_idx
from kamino.constants import EARTH_DELTA_IW, EARTH_MANTLE_MG_SI

_STYLE = os.path.join(os.path.dirname(__file__), "planetary-chem-paper.mplstyle")

# Fixed mineral order and colours -- petrological grouping (feldspars, feldspathoid, pyroxenes,
# olivines), so a mineral keeps its colour across every panel and figure. Hues are the project's
# own axes.prop_cycle (Paul Tol 'bright', published colourblind-safe), consumed unchanged so
# these figures match the rest of the paper.
# Every phase `cipw_norm` can emit, silica-rich -> silica-poor, each Fe/desilication endmember
# placed next to the parent it derives from. Quartz, Akermanite, Hedenbergite and Ferrosilite were
# added 2026-08-27: the norm has emitted Quartz since the silica-oversaturated crusts were made
# mass-conservative, Akermanite since larnite was rerouted through diopside, and the two
# Fe-pyroxenes unconditionally since the Fe-pyroxene adoption (development history 25 and 27).
# Before that they were missing here, so the stacked bands did not sum to 1 and the shortfall was
# invisible -- at reduced silica-rich compositions the Fe-pyroxenes alone reach ~39 wt%.
#
# Hues are the project's own axes.prop_cycle (Paul Tol 'bright', published colourblind-safe).
# There are more phases than Tol hues, so each derived endmember shares its PARENT's hue and is
# separated by hatching -- semantically exact, and it keeps the palette unextended.
MINERAL_ORDER = ["Quartz", "Anorthite", "Albite", "Nepheline",
                 "Diopside", "Hedenbergite", "Akermanite",
                 "Enstatite", "Ferrosilite", "Forsterite", "Fayalite"]
MINERAL_COLOUR = {
    "Quartz":       "#332288",
    "Anorthite":    "#4477AA",
    "Albite":       "#EE6677",
    "Nepheline":    "#228833",
    "Diopside":     "#CCBB44",
    "Hedenbergite": "#CCBB44",
    "Akermanite":   "#CCBB44",
    "Enstatite":    "#66CCEE",
    "Ferrosilite":  "#66CCEE",
    "Forsterite":   "#AA3377",
    "Fayalite":     "#BBBBBB",
}
MINERAL_HATCH = {"Hedenbergite": "...", "Akermanite": "///", "Ferrosilite": "..."}
LINE_COLOUR = "#4477AA"

# Pore conditions for the rate constants, matching the ocean-world state used throughout the
# Mg/Si analysis (T_pore 343 K, pH 6.6).
K_PRESSURE, K_TEMPERATURE, K_PH = 3.0e7, 343.0, 6.6

# Earth's calibrated potential temperature, drawn as a reference line on the T_p panel.
EARTH_TP = 1325.0

# dIW values for the redox panel. Grid nodes of the MAGEMin table, so nothing is interpolated.
DIW_CUT = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0]

LABEL_MIN_BAND = 0.07   # only direct-label a band thick enough to hold the text
LABEL_MAX_COUNT = 4     # legend carries identity; direct labels are selective (dataviz rule)
LABEL_EDGE_MARGIN = 2   # keep direct labels this many points clear of each end

MODEL_LABEL = "MAGEMin (HGP18)"


# The columns cipw_norm actually consumes. The CSV also carries provenance and derived columns
# (closure, residual_phases, warnings, mantle_feo, delta_iw_melt, ...), and the previous
# "everything except T_p/mg_si/melt_fraction" filter fed those to the norm as if they were
# oxides -- float("fixed-F") raises, so the script could not run at all.
OXIDE_COLUMNS = ("SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MgO", "CaO", "Na2O", "K2O")


def load_table(path, delta_iw=None):
    """Read one dIW slice of the CSV into (mg_si, T_p, melt_fraction, oxides), ordered by Mg/Si.

    The table is a 2-D (Mg/Si x dIW) grid. This function returns a single dIW CUT through it,
    because everything downstream plots against Mg/Si alone. Without that, all 153 rows collapsed
    onto the Mg/Si axis with nine values stacked at each point and the lines sawtoothed between
    redox states.
    """
    rows = [r for r in csv.DictReader(line for line in open(path)
                                      if not line.lstrip().startswith("#"))]
    if delta_iw is None:
        delta_iw = EARTH_DELTA_IW
    available = sorted({float(r["delta_iw"]) for r in rows})
    diw = min(available, key=lambda v: abs(v - delta_iw))
    rows = sorted((r for r in rows if np.isclose(float(r["delta_iw"]), diw)),
                  key=lambda r: float(r["mg_si"]))
    if not rows:
        raise SystemExit(f"no rows at delta_iw={delta_iw}; table has {available}")
    mg_si = np.array([float(r["mg_si"]) for r in rows])
    T_p = np.array([float(r["T_p"]) for r in rows])
    melt_fraction = np.array([float(r["melt_fraction"]) for r in rows])
    oxides = [{k: float(r[k]) for k in OXIDE_COLUMNS if k in r} for r in rows]
    return mg_si, T_p, melt_fraction, oxides, diw


def norm_series(oxide_list):
    """CIPW norm for each composition -> (n_minerals, n_points) weight-fraction array."""
    comps = []
    for oxides in oxide_list:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            comps.append(cipw_norm(oxides))
    return np.array([[c.get(m, 0.0) for c in comps] for m in MINERAL_ORDER])


def reactivity(oxides):
    """Cation charge per kg of rock, charge-weighted k, and k[Na] for one composition."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        comp = cipw_norm(oxides)
    totals = {e: 0.0 for e in elements}
    for mineral, fraction in comp.items():
        # No .get() fallback: a phase the norm emits but MINERAL_MOLAR_MASS or `stoichiometry`
        # does not know is a real inconsistency between this script and the model, and a 150 g/mol
        # guess would hide it behind a plausible number.
        if mineral not in MINERAL_MOLAR_MASS:
            raise KeyError(f"{mineral} is emitted by cipw_norm but missing from "
                           f"mineral_info.MINERAL_MOLAR_MASS")
        moles = fraction / MINERAL_MOLAR_MASS[mineral] * 1000.0
        for i, element in enumerate(elements):
            totals[element] += moles * stoichiometry[mineral][i]
    charge = sum(totals[e] * ION_CHARGE[i]
                 for i, e in enumerate(elements) if ION_CHARGE[i] > 0)
    weights = np.where(ION_CHARGE > 0, ION_CHARGE, 0.0)
    weights[al_idx] = 0.0        # Al precipitates as clay; it delivers no ocean alkalinity
    k = get_k(K_PRESSURE, K_TEMPERATURE, K_PH, comp)
    return charge, float(np.dot(weights, k)), float(k[na_idx])


def stacked_panel(ax, x, stack, title):
    """One stacked-area panel of normative mineralogy."""
    polys = ax.stackplot(x, stack, colors=[MINERAL_COLOUR[m] for m in MINERAL_ORDER],
                         edgecolor="white", linewidth=0.6)
    # Derived endmembers share their parent's hue, so they are separated by texture instead.
    for poly, mineral in zip(polys, MINERAL_ORDER):
        if mineral in MINERAL_HATCH:
            poly.set_hatch(MINERAL_HATCH[mineral])
            poly.set_edgecolor("white")
    # Direct-label at most the four thickest bands -- the legend carries the rest. Placement is
    # restricted to interior x so a label cannot overhang the panel edge or collide with the
    # y-axis ticks, which is what happens if a band peaks at the first or last point.
    cumulative = np.cumsum(stack, axis=0)
    m = LABEL_EDGE_MARGIN if len(x) > 2 * LABEL_EDGE_MARGIN + 1 else 1
    interior = slice(m, -m) if len(x) > 2 * m else slice(None)
    candidates = []
    for i, mineral in enumerate(MINERAL_ORDER):
        band = stack[i][interior]
        if band.size == 0:
            continue
        j = int(np.argmax(band)) + (m if len(x) > 2 * m else 0)
        if stack[i][j] >= LABEL_MIN_BAND:
            candidates.append((stack[i][j], i, j, mineral))
    for _, i, j, mineral in sorted(candidates, reverse=True)[:LABEL_MAX_COUNT]:
        y = cumulative[i, j] - stack[i][j] / 2.0
        ax.text(x[j], y, mineral, ha="center", va="center", fontsize=6,
                color="white", fontweight="bold", clip_on=True)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=8, pad=4)
    ax.set_xlabel("mantle Mg/Si (molar)")


def figure_mineralogy(mg_si, oxides, diw, outdir):
    """Mineralogy along Mg/Si at fixed dIW, and along dIW at fixed Mg/Si.

    The left panel used to show the pre-MAGEMin "SiO2-rescale proxy" for comparison. That proxy
    was removed with the old `oxide_composition(T_p, mg_si)` signature, so the call raised
    ValueError and this figure could not be produced at all. The redox axis is the meaningful
    second cut now, so the panel shows that instead.
    """
    fig, axes = plt.subplots(1, 2, figsize=(4.8, 2.8), sharey=True)

    stacked_panel(axes[0], mg_si, norm_series(oxides),
                  f"{MODEL_LABEL}\nvarying Mg/Si at $\\Delta$IW = {diw:+g}")
    axes[0].set_ylabel("normative weight fraction")

    diw_axis = np.array(DIW_CUT)
    diw_ox = [oxide_composition(EARTH_MANTLE_MG_SI, float(d)) for d in diw_axis]
    stacked_panel(axes[1], diw_axis, norm_series(diw_ox),
                  f"{MODEL_LABEL}\nvarying $\\Delta$IW at Mg/Si = {EARTH_MANTLE_MG_SI:g}")
    axes[1].set_xlabel("core-formation $\\Delta$IW")

    handles = [Patch(facecolor=MINERAL_COLOUR[m], edgecolor="white", linewidth=0.6,
                     hatch=MINERAL_HATCH.get(m), label=m) for m in MINERAL_ORDER]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, fontsize=6.5,
               handlelength=1.2, columnspacing=1.0, borderpad=0.2)
    path = os.path.join(outdir, "crust_mineralogy_mgsi.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


def figure_reactivity(mg_si, T_p, oxides, outdir):
    charge, k_charge, k_na = zip(*[reactivity(o) for o in oxides])
    fig, axes = plt.subplots(4, 1, figsize=(3.33, 6.4), sharex=True)
    series = [(T_p, "$T_p$ required\n($^\\circ$C)", False),
              (charge, "cation charge\n(eq kg$^{-1}$ rock)", False),
              (k_charge, "charge-weighted $k$\n(mol m$^{-2}$ s$^{-1}$)", False),
              (k_na, "$k$[Na]\n(mol m$^{-2}$ s$^{-1}$)", True)]
    for ax, (y, label, logscale) in zip(axes, series):
        ax.plot(mg_si, y, marker="o", markersize=3, linewidth=1.0, color=LINE_COLOUR)
        ax.set_ylabel(label)
        if logscale:
            ax.set_yscale("log")
        ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.5)
    axes[0].axhline(EARTH_TP, color="0.45", linestyle="--", linewidth=0.8)
    axes[0].text(mg_si.min(), EARTH_TP, " Earth", fontsize=6, color="0.35",
                 va="bottom", ha="left")
    axes[-1].set_xlabel("mantle Mg/Si (molar)")
    path = os.path.join(outdir, "crust_reactivity_mgsi.png")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


def main():
    global MODEL_LABEL
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, help="CSV from make_crust_compositions.jl")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "../output"))
    parser.add_argument("--model-label", default=MODEL_LABEL,
                        help="Name of the melting model that produced the CSV.")
    parser.add_argument("--delta-iw", type=float, default=EARTH_DELTA_IW,
                        help="dIW slice for the Mg/Si panels (default: Earth).")
    args = parser.parse_args()
    MODEL_LABEL = args.model_label

    plt.style.use(_STYLE)
    os.makedirs(args.outdir, exist_ok=True)
    mg_si, T_p, melt_fraction, oxides, diw = load_table(args.csv, args.delta_iw)
    print(f"Loaded {len(mg_si)} compositions at dIW {diw:+g}: "
          f"Mg/Si {mg_si.min():g}-{mg_si.max():g}, "
          f"T_p {T_p.min():.0f}-{T_p.max():.0f} degC, "
          f"F {melt_fraction.min():.3f}-{melt_fraction.max():.3f}")
    figure_mineralogy(mg_si, oxides, diw, args.outdir)
    figure_reactivity(mg_si, T_p, oxides, args.outdir)


if __name__ == "__main__":
    main()
