"""Diagnostic plots for the state of a single planet.

Answers three questions about one saved run, at one instant (default: the final state):
  what is precipitating, what is dissolving, and what are the sources and sinks.

    python experiments/plot_planet_state.py output/rg_final_earth.json
    python experiments/plot_planet_state.py output/rg_final_earth.json --element Alkalinity
    python experiments/plot_planet_state.py output/rg_final_earth.json --index 0   # initial state

Produces, next to the run JSON:
    <name>_state.png    2x3 snapshot dashboard  — the working diagnostic
    <name>_ledger.png   process x element flux ledger — the dense one-panel version
    <name>_sankey.png   source -> ocean -> sink flow for one element — the talk figure

All three read the same diagnostic dict (kamino.diagnostics.diagnose_run) and share one
process -> colour mapping, so a process is the same colour in every panel of every figure.
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from kamino.diagnostics import (
    diagnose_run, PLOT_ELEMENTS, EARTH_SEAWATER,
    SOURCE_TERMS, SINK_TERMS, TRANSFER_TERMS, ALL_TERMS,
)

_style = os.path.join(os.path.dirname(__file__), 'planetary-chem-paper.mplstyle')
if os.path.exists(_style):
    plt.style.use(_style)

# ---------------------------------------------------------------------------
# One colour system, shared by all three figures.
#
# Identity is carried by a fixed process order (never cycled), but the hue FAMILY also
# encodes polarity: sources warm, sinks cool, internal transfers grey. Every process is
# also directly labelled wherever it appears, so identity is never colour-alone.
# ---------------------------------------------------------------------------

PROCESS_COLORS = {
    # sources — warm ramp
    'outgassing':         '#8c2d04',
    'seafloor LT':        '#cc4c02',
    'seafloor HT':        '#ec7014',
    'continental':        '#fe9929',
    # sinks — cool ramp
    'precipitation':      '#084081',
    'shelf carbonate':    '#0868ac',
    'reverse weathering': '#2b8cbe',
    'biogenic':           '#4eb3d3',
    'Cl subduction':      '#7bccc4',
    'Na albitization':    '#a8ddb5',
    # internal transfer — neutral
    'HT Mg-Ca exchange':  '#969696',
}

# Diverging pair for the ledger: source (warm) / sink (cool) through a NEUTRAL grey
# midpoint. Not a rainbow, and no hue at zero.
LEDGER_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'source_sink', ['#08306b', '#4292c6', '#d9d9d9', '#f16913', '#7f2704'])

ELEMENT_LABELS = {
    'Alkalinity': 'Alk', 'C': 'C', 'Si': 'Si', 'Al': 'Al', 'Fe': 'Fe',
    'Ca': 'Ca', 'Mg': 'Mg', 'Na': 'Na', 'Cl': 'Cl',
}

MINERAL_LABELS = {
    'Calcite': 'Calcite', 'Siderite': 'Siderite', 'Kaolinite': 'Kaolinite',
    'Goethite': 'Goethite', 'SiO2(am)': 'SiO$_2$(am)', 'Halite': 'Halite',
    'Sepiolite(d)': 'Sepiolite', 'Saponite-Na': 'Saponite-Na', 'Greenalite': 'Greenalite',
}

GREY = '#969696'
INK = '#252525'


def _fmt(x, unit=''):
    if not np.isfinite(x):
        return 'n/a'
    if x == 0:
        return '0'
    if abs(x) >= 1e4 or abs(x) < 1e-2:
        return f'{x:.1e}{unit}'
    return f'{x:.2f}{unit}'


def _present_terms(diag, elements=None):
    """Process terms that actually move something, in canonical order."""
    els = elements or PLOT_ELEMENTS
    return [t for t in ALL_TERMS
            if t in diag['fluxes']
            and any(abs(diag['fluxes'][t][e]) > 0 for e in els)]


# ---------------------------------------------------------------------------
# Figure 1 — snapshot dashboard
# ---------------------------------------------------------------------------

def plot_dashboard(diag, out_path):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    s = diag['scalars']

    # (a) scalar state --------------------------------------------------------
    ax = axes[0, 0]
    ax.axis('off')
    Da = s['Da']
    regime = 'kinetic' if (np.isfinite(Da) and Da < 1) else 'thermodynamic'
    rows = [
        ('Surface temperature', f"{s['T_surface']:.1f} K"),
        ('Seafloor temperature', f"{s['T_seafloor']:.1f} K"),
        (r'$P_{\mathrm{CO_2}}$', f"{s['P_CO2_bar']:.2e} bar"),
        ('Surface pH', f"{s['pH_surface']:.2f}"),
        ('Seafloor pH', f"{s['pH_seafloor']:.2f}"),
        ('Salinity', f"{s['salinity']:.1f} g/kg"),
        ('Damkohler number', f"{Da:.2g}  ({regime})"),
        ('Supply efficiency', f"{100 * s['supply_efficiency']:.0f} %"),
        ('Land fraction', f"{s['land_fraction']:.2f}"),
        ('Termination', diag['termination']),
        ('Time', f"{diag['time_yr']:.2e} yr"),
    ]
    ax.set_title('(a)  Planet state', loc='left', fontweight='bold')
    for i, (k, v) in enumerate(rows):
        y = 0.94 - i * 0.088
        ax.text(0.0, y, k, transform=ax.transAxes, color=GREY, va='center')
        ax.text(1.0, y, v, transform=ax.transAxes, color=INK, va='center',
                ha='right', fontweight='bold')

    # (b) ocean composition ---------------------------------------------------
    ax = axes[0, 1]
    comp_els = [e for e in PLOT_ELEMENTS if e in EARTH_SEAWATER]
    vals = [max(diag['composition'][e] * 1e3, 1e-6) for e in comp_els]  # mmol/kg
    earth = [EARTH_SEAWATER[e] * 1e3 for e in comp_els]
    x = np.arange(len(comp_els))
    ax.bar(x, vals, width=0.62, color='#4292c6', zorder=3, label='Model')
    ax.scatter(x, earth, marker='_', s=260, linewidths=2.2, color=INK, zorder=5,
               label='Earth seawater')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([ELEMENT_LABELS[e] for e in comp_els])
    ax.set_ylabel('Concentration (mmol/kg)')
    ax.set_title('(b)  Ocean composition', loc='left', fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.35, zorder=0)
    ax.legend(frameon=False, fontsize=7, loc='upper left')

    # (c) saturation state — what is precipitating ----------------------------
    ax = axes[0, 2]
    si_ocean = diag['SI']['ocean']
    si_pore = diag['SI']['pore']
    mins = [m for m in MINERAL_LABELS if m in si_ocean or m in si_pore]
    y = np.arange(len(mins))
    ocean_vals = [si_ocean.get(m, np.nan) for m in mins]
    pore_vals = [si_pore.get(m, np.nan) for m in mins]
    colors = ['#0868ac' if (np.isfinite(v) and v > 0) else GREY for v in ocean_vals]
    ax.barh(y, [v if np.isfinite(v) else 0 for v in ocean_vals], height=0.6,
            color=colors, zorder=3, label='Ocean')
    ax.scatter([v if np.isfinite(v) else np.nan for v in pore_vals], y,
               marker='D', s=26, facecolors='none', edgecolors='#cc4c02',
               linewidths=1.3, zorder=5, label='Pore space')
    ax.axvline(0, color=INK, linewidth=1.0, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels([MINERAL_LABELS[m] for m in mins])
    ax.invert_yaxis()
    ax.set_xlabel('Saturation index')
    ax.set_title('(c)  Saturation — SI > 0 precipitates', loc='left', fontweight='bold')
    ax.grid(True, axis='x', linestyle='--', alpha=0.35, zorder=0)
    ax.legend(fontsize=7, loc='lower right', framealpha=0.92, facecolor='white',
              edgecolor='none')
    ax.margins(x=0.12)

    # (d) sources and sinks per element ---------------------------------------
    ax = axes[1, 0]
    terms = _present_terms(diag)
    els = [e for e in PLOT_ELEMENTS if any(abs(diag['fluxes'][t][e]) > 0 for t in terms)]
    x = np.arange(len(els))
    pos = np.zeros(len(els))
    neg = np.zeros(len(els))
    for t in terms:
        v = np.array([diag['fluxes'][t][e] for e in els])
        up = np.where(v > 0, v, 0.0)
        dn = np.where(v < 0, v, 0.0)
        ax.bar(x, up, bottom=pos, width=0.66, color=PROCESS_COLORS[t],
               label=t, zorder=3, linewidth=0.4, edgecolor='white')
        ax.bar(x, dn, bottom=neg, width=0.66, color=PROCESS_COLORS[t],
               zorder=3, linewidth=0.4, edgecolor='white')
        pos += up
        neg += dn
    ax.scatter(x, [diag['net'][e] for e in els], marker='D', s=26, color=INK,
               zorder=6, label='net')
    ax.axhline(0, color=INK, linewidth=1.0, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels([ELEMENT_LABELS[e] for e in els])
    ax.set_yscale('symlog', linthresh=1e-3)
    ax.set_ylabel('Flux (Tmol/yr)')
    ax.set_title('(d)  Sources (+) and sinks (-)', loc='left', fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.35, zorder=0)
    # Process legend is figure-level (see below) so it never sits on top of the bars.
    process_handles = [Line2D([0], [0], marker='s', color='w', markersize=8,
                              markerfacecolor=PROCESS_COLORS[t], label=t) for t in terms]
    process_handles.append(Line2D([0], [0], marker='D', color='w', markersize=6,
                                  markerfacecolor=INK, label='net'))

    # (e) dissolution / supply side -------------------------------------------
    ax = axes[1, 1]
    supply = [t for t in SOURCE_TERMS if t in diag['fluxes']]
    alk = [diag['fluxes'][t]['Alkalinity'] for t in supply]
    y = np.arange(len(supply))
    ax.barh(y, alk, height=0.6, color=[PROCESS_COLORS[t] for t in supply], zorder=3)
    for yi, v in zip(y, alk):
        ax.text(v, yi, f'  {_fmt(v)}', va='center', ha='left' if v >= 0 else 'right',
                fontsize=7, color=INK)
    ax.axvline(0, color=INK, linewidth=1.0, zorder=4)
    ax.set_yticks(y)
    ax.set_yticklabels(supply)
    ax.invert_yaxis()
    ax.set_xlabel('Alkalinity flux (Teq/yr)')
    ax.set_title('(e)  What is dissolving', loc='left', fontweight='bold')
    ax.grid(True, axis='x', linestyle='--', alpha=0.35, zorder=0)
    txt = (f"Da = {_fmt(Da)} ({regime}-limited)\n"
           f"supply efficiency = {100 * s['supply_efficiency']:.0f}%\n"
           f"reactive area = {_fmt(s['A_reactive'])}")
    ax.text(0.98, 0.96, txt, transform=ax.transAxes, ha='right', va='top',
            fontsize=7, color=GREY)
    ax.margins(x=0.18)  # room for the value labels at the bar ends

    # (f) imbalance / residence -----------------------------------------------
    ax = axes[1, 2]
    els_r = [e for e in PLOT_ELEMENTS if diag['composition'][e] > 0]
    rates = []
    for e in els_r:
        tau = diag['residence'][e]
        rates.append(1e9 / tau if np.isfinite(tau) and tau > 0 else 0.0)  # per Gyr
    x = np.arange(len(els_r))
    cols = ['#f16913' if diag['net'][e] > 0 else '#4292c6' for e in els_r]
    ax.bar(x, np.maximum(rates, 1e-6), width=0.62, color=cols, zorder=3)
    ax.axhline(0.1, color=INK, linestyle='--', linewidth=1.0, zorder=4)
    ax.text(0.02, 0.1, ' convergence threshold', transform=ax.get_yaxis_transform(),
            fontsize=6.5, color=GREY, va='bottom')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels([ELEMENT_LABELS[e] for e in els_r])
    ax.set_ylabel(r'$|F_{\mathrm{net}}|\,/\,$inventory (per Gyr)')
    ax.set_title('(f)  Imbalance — how far from steady state', loc='left', fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.35, zorder=0)
    ax.legend(handles=[
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#f16913', markersize=8,
               label='net source (accumulating)'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#4292c6', markersize=8,
               label='net sink (depleting)'),
    ], frameon=False, fontsize=6.5, loc='upper right')

    name = diag['config'].get('name', 'planet')
    fig.suptitle(f"{name}  —  {diag['termination']} at {diag['time_yr']:.2e} yr",
                 fontweight='bold')
    fig.legend(handles=process_handles, loc='outside lower center',
               ncol=min(len(process_handles), 6), frameon=False, fontsize=8)
    fig.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')


# ---------------------------------------------------------------------------
# Figure 2 — process x element ledger
# ---------------------------------------------------------------------------

def plot_ledger(diag, out_path):
    terms = _present_terms(diag)
    els = [e for e in PLOT_ELEMENTS if any(abs(diag['fluxes'][t][e]) > 0 for t in terms)]

    M = np.array([[diag['fluxes'][t][e] for e in els] for t in terms])

    scale = np.nanmax(np.abs(M)) if M.size else 1.0
    scale = scale if scale > 0 else 1.0
    linthresh = max(scale * 1e-4, 1e-9)
    norm = mcolors.SymLogNorm(linthresh=linthresh, vmin=-scale, vmax=scale, base=10)

    fig, ax = plt.subplots(figsize=(1.05 * len(els) + 4.6, 0.52 * len(terms) + 3.0))
    im = ax.imshow(M, cmap=LEDGER_CMAP, norm=norm, aspect='auto')

    ax.set_xticks(np.arange(len(els)))
    ax.set_xticklabels([ELEMENT_LABELS[e] for e in els])
    ax.set_yticks(np.arange(len(terms)))
    ax.set_yticklabels(terms)
    ax.set_xticks(np.arange(len(els) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(terms) + 1) - 0.5, minor=True)
    ax.grid(which='minor', color='white', linewidth=1.6)
    ax.tick_params(which='minor', length=0)

    # A cell's number is its own label — no colour-only reading required.
    for i in range(len(terms)):
        for j in range(len(els)):
            v = M[i, j]
            if v == 0:
                ax.text(j, i, '.', ha='center', va='center', color=GREY, fontsize=8)
                continue
            shade = abs(norm(v) - 0.5) * 2
            ax.text(j, i, f'{v:.2g}', ha='center', va='center', fontsize=6.5,
                    color='white' if shade > 0.62 else INK)

    # Separator between source block, sink block and transfers.
    n_src = sum(1 for t in terms if t in SOURCE_TERMS)
    n_snk = sum(1 for t in terms if t in SINK_TERMS)
    for boundary in (n_src, n_src + n_snk):
        if 0 < boundary < len(terms):
            ax.axhline(boundary - 0.5, color=INK, linewidth=1.6)

    cbar = fig.colorbar(im, ax=ax, pad=0.02, aspect=28)
    cbar.set_label('Flux (Tmol/yr;  Teq/yr for Alk)\n$-$ sink        $+$ source')

    net = [diag['net'][e] for e in els]
    ax.set_xlabel('net:  ' + '   '.join(f'{ELEMENT_LABELS[e]} {v:+.2g}'
                                        for e, v in zip(els, net)), fontsize=7, color=GREY)

    name = diag['config'].get('name', 'planet')
    ax.set_title(f'{name} — source/sink ledger at {diag["time_yr"]:.2e} yr',
                 fontweight='bold', loc='left')
    fig.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')


# ---------------------------------------------------------------------------
# Figure 3 — Sankey flow for one element
# ---------------------------------------------------------------------------

def _band(ax, x0, y0, x1, y1, w0, w1, color, alpha=0.78):
    """A smooth flow band from (x0, y0) with width w0 to (x1, y1) with width w1."""
    mx = (x0 + x1) / 2
    verts = [
        (x0, y0 - w0 / 2), (mx, y0 - w0 / 2), (mx, y1 - w1 / 2), (x1, y1 - w1 / 2),
        (x1, y1 + w1 / 2), (mx, y1 + w1 / 2), (mx, y0 + w0 / 2), (x0, y0 + w0 / 2),
        (x0, y0 - w0 / 2),
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor='white',
                           linewidth=0.8, alpha=alpha, zorder=3))


def plot_sankey(diag, element, out_path):
    f = {t: diag['fluxes'][t][element] for t in ALL_TERMS if t in diag['fluxes']}
    sources = {t: v for t, v in f.items() if v > 0}
    sinks = {t: -v for t, v in f.items() if v < 0}

    if not sources and not sinks:
        print(f'No {element} flux to draw — skipping Sankey.')
        return

    total = max(sum(sources.values()), sum(sinks.values()), 1e-30)
    gap = 0.06 * total  # vertical spacer between bands, in flux units

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.set_xlim(0, 10)
    ax.axis('off')

    box_x0, box_x1 = 4.3, 5.7
    box_h = total
    ax.add_patch(Rectangle((box_x0, -box_h / 2), box_x1 - box_x0, box_h,
                           facecolor='#deebf7', edgecolor=INK, linewidth=1.2, zorder=4))
    ax.text((box_x0 + box_x1) / 2, 0, 'OCEAN', ha='center', va='center',
            fontweight='bold', color=INK, zorder=5)

    def _stack(items, total_h):
        """Top-down y-centres for a set of bands filling total_h, with gaps."""
        n = len(items)
        span = sum(items.values()) + gap * max(n - 1, 0)
        y = span / 2
        out = {}
        for k, v in items.items():
            out[k] = (y - v / 2, v)
            y -= v + gap
        return out

    for k, (yc, w) in _stack(sources, box_h).items():
        _band(ax, 2.2, yc, box_x0, yc * (box_h / max(box_h, 1e-30)), w, w, PROCESS_COLORS[k])
        ax.text(2.1, yc, f'{k}\n{f[k]:.2g}', ha='right', va='center', fontsize=8, color=INK)

    for k, (yc, w) in _stack(sinks, box_h).items():
        _band(ax, box_x1, yc, 7.8, yc, w, w, PROCESS_COLORS[k])
        ax.text(7.9, yc, f'{k}\n{-f[k]:.2g}', ha='left', va='center', fontsize=8, color=INK)

    net = diag['net'][element]
    unit = 'Teq/yr' if element == 'Alkalinity' else 'Tmol/yr'
    ax.text(0.5, 0.02,
            f'sources {sum(sources.values()):.3g}   sinks {sum(sinks.values()):.3g}   '
            f'net {net:+.2g} {unit}',
            transform=ax.transAxes, ha='center', fontsize=8, color=GREY)

    name = diag['config'].get('name', 'planet')
    ax.set_title(f'{name} — {ELEMENT_LABELS.get(element, element)} budget '
                 f'({unit}), band width $\\propto$ flux',
                 fontweight='bold')
    ax.set_ylim(-box_h * 0.85, box_h * 0.85)
    fig.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run', help='Path to a run JSON written by Planet.time_evolve')
    ap.add_argument('--index', type=int, default=-1,
                    help='Time index to diagnose (default -1, the final state)')
    ap.add_argument('--element', default='C',
                    help='Element for the Sankey figure (default C; try Alkalinity, Ca, Mg)')
    ap.add_argument('--outdir', default=None, help='Where to write figures (default: alongside the run)')
    args = ap.parse_args()

    diag = diagnose_run(args.run, index=args.index)

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.run))
    os.makedirs(outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(args.run))[0]

    plot_dashboard(diag, os.path.join(outdir, f'{stem}_state.png'))
    plot_ledger(diag, os.path.join(outdir, f'{stem}_ledger.png'))
    plot_sankey(diag, args.element, os.path.join(outdir, f'{stem}_sankey.png'))


if __name__ == '__main__':
    main()
