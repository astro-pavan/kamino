"""Mantle -> melt -> minerals diagnostic: three pies per composition, with the melt temperature.

    python experiments/plot_crust_diagnostic.py                      # 3x3 grid + default individuals
    python experiments/plot_crust_diagnostic.py --points "0.5,-1.0"  # individuals for these cells
    python experiments/plot_crust_diagnostic.py --grid-only

Everything is read through kamino.crust_composition so the figure can never disagree with the
model: `mantle_composition` is the melting calculation's INPUT, `oxide_composition` the melt it
produced, and `mineral_composition` the assemblage the weathering law actually dissolves. Reading
left to right therefore follows one parcel of rock through the whole pipeline.

Design notes
------------
Pies rather than bars because the question is "where did each cation end up" -- three compositions
of the same 100 g, compared as shapes, with every slice above the label threshold carrying its own
name and percentage so identity is never colour alone. Slices below the threshold are named in a
note under the pie (individual figures) or in the companion CSV (the grid).

Oxides and minerals are TWO SEPARATE ENCODINGS and get two separate legends. They have to be:
the mineral hues encode mineral family (plot_crust_grid.COLORS), which no oxide palette can
mirror, so a colour shared between the two sides is a coincidence and must not be read as a
relationship. The first two pies are directly comparable to each other -- same oxides, same slot
order (PIPELINE_OXIDES), so mantle and melt can be differenced by eye -- and the third is read
against the mineral legend alone.

The oxide slot order's ring-adjacent pairs were validated the way plot_crust_grid's were --
Machado et al. (2009) at severity 1.0, OKLab dE x100: worst adjacent CVD dE 9.2, worst tritan
9.6, worst normal-vision 20.8. Do not reassign without re-checking.

Cr2O3 is dropped from both oxide pies and the remainder renormalised, because that is what the
norm does with it (no chromite in the phase set, and Cr is not a tracked ion).
"""

import argparse
import os
import sys
import warnings

import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

from kamino.constants import EARTH_DELTA_IW, EARTH_MANTLE_MG_SI
from kamino.crust_composition import (PIPELINE_OXIDES, load_crust_table, mantle_composition,
                                      mineral_composition, oxide_composition)
from plot_crust_grid import COLORS as MINERAL_COLORS
from plot_crust_grid import HATCHED as MINERAL_HATCH
from plot_crust_grid import (INK, INK_2, INK_3, MINERALS, SURFACE, draw_reference_pie,
                             hatch_ink)

# Oxide hues -- independent of the mineral palette, which encodes mineral family and cannot also
# carry oxide identity. Assigned by exhausting the permutations against the CVD checks.
OXIDE_COLORS = {'SiO2': '#2a78d6', 'TiO2': '#e34948', 'Al2O3': '#eda100', 'FeOt': '#008300',
                'MgO': '#e87ba4', 'CaO': '#4a3aa7', 'Na2O': '#eb6834', 'K2O': '#1baf7a'}
OXIDE_LABEL = {'SiO2': 'SiO$_2$', 'TiO2': 'TiO$_2$', 'Al2O3': 'Al$_2$O$_3$', 'FeOt': 'FeO$_t$',
               'MgO': 'MgO', 'CaO': 'CaO', 'Na2O': 'Na$_2$O', 'K2O': 'K$_2$O'}

# The 3x3 diagnostic grid. Rows run oxidised-at-the-top, matching plot_crust_grid.
GRID_MGSI = [0.5, 1.25, 2.0]
GRID_DIW = [-1.0, -2.0, -5.0]

# Individual figures: Earth, the quartz-out transition, and the four corners of the grid.
DEFAULT_POINTS = [(EARTH_MANTLE_MG_SI, EARTH_DELTA_IW), (0.5, -1.2), (0.5, -5.0), (2.0, -1.0)]

# Label thresholds in wt%. The grid needs a coarser one: 27 pies x 8 slices of labels at grid
# scale is unreadable, and the CSV carries the full numbers anyway.
THRESH_SINGLE, THRESH_GRID = 1.0, 5.0


def cell_data(mgsi, diw, t_interp, f_interp):
    """Mantle oxides, melt oxides, minerals (all wt%, fixed order) plus the melt diagnostics."""
    mantle = {ox: mantle_composition(mgsi, diw)[ox] for ox in PIPELINE_OXIDES}
    total = sum(mantle.values())                      # renormalise over the norm's oxide set
    mantle = {ox: v / total * 100.0 for ox, v in mantle.items()}
    melt = {ox: oxide_composition(mgsi, diw)[ox] for ox in PIPELINE_OXIDES}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        minerals = mineral_composition(mgsi, diw)
    minerals = {m: minerals.get(m, 0.0) * 100.0 for m in MINERALS}
    return mantle, melt, minerals, float(t_interp([[mgsi, diw]])[0]), float(f_interp([[mgsi, diw]])[0])


def _spread_labels(texts, min_dy):
    """Push overlapping radial labels apart vertically, per side of the pie.

    Two thin slices side by side put their labels at almost the same angle, and matplotlib will
    happily stack them on top of each other. Separating them per side preserves the left/right
    alignment that keeps the text clear of the wedges.
    """
    for side in (1, -1):
        column = sorted((t for t in texts if t.get_text()
                         and (1 if t.get_position()[0] >= 0 else -1) == side),
                        key=lambda t: t.get_position()[1])
        for above, below in zip(column[1:], column):
            x, y = above.get_position()
            floor = below.get_position()[1] + min_dy
            if y < floor:
                above.set_position((x, floor))


def draw_pie(ax, comp, colors, threshold, title, hatches=None, label_fmt=str, fontsize=7.4,
             min_dy=0.17):
    """One composition pie. Returns the names of the slices too small to label."""
    items = [(k, v) for k, v in comp.items() if v > 1e-6]
    labels = [f'{label_fmt(k)}\n{v:.1f}%' if v >= threshold else '' for k, v in items]
    # radius < 1 and labeldistance > 1 keep the labels clear of the wedges; matplotlib's own
    # per-angle horizontal alignment must NOT be overridden or long names run back over the pie.
    wedges, texts = ax.pie([v for _, v in items], labels=labels, startangle=90, counterclock=False,
                           radius=0.82, labeldistance=1.12, colors=[colors[k] for k, _ in items],
                           wedgeprops=dict(edgecolor=SURFACE, linewidth=1.1))
    for (name, _), w in zip(items, wedges):
        if hatches and name in hatches:
            w.set_hatch(hatches[name])
            w.set_edgecolor(hatch_ink(colors[name]))
    for t in texts:
        t.set(fontsize=fontsize, color=INK_2)
    _spread_labels(texts, min_dy)
    ax.set_xlim(-1.34, 1.34)          # headroom for the pushed-apart labels, so none is clipped
    ax.set_ylim(-1.30, 1.30)
    ax.set_title(title, fontsize=8.6, color=INK, pad=6)
    return [label_fmt(k) for k, v in items if v < threshold]


def draw_cell(axs, mgsi, diw, data, threshold, pie_titles=True, fontsize=7.4, min_dy=0.17):
    """The three pies for one composition, on three prepared axes."""
    mantle, melt, minerals, _, _ = data
    kw = dict(fontsize=fontsize, min_dy=min_dy)
    return {
        'mantle': draw_pie(axs[0], mantle, OXIDE_COLORS, threshold,
                           'Mantle oxides' if pie_titles else '',
                           label_fmt=OXIDE_LABEL.get, **kw),
        'melt': draw_pie(axs[1], melt, OXIDE_COLORS, threshold,
                         'Primary melt oxides' if pie_titles else '',
                         label_fmt=OXIDE_LABEL.get, **kw),
        'minerals': draw_pie(axs[2], minerals, MINERAL_COLORS, threshold,
                             'Normative minerals' if pie_titles else '',
                             hatches=MINERAL_HATCH, **kw),
    }


def header(mgsi, diw, data):
    """The one-line state summary that carries the melt temperature."""
    _, _, _, t_melt, f = data
    return (f'Mg/Si = {mgsi:g}    $\\Delta$IW = {diw:+g}    '
            f'$T_\\mathrm{{melt}}$ = {t_melt:.0f} $^\\circ$C    F = {f:.3f}')


def oxide_legend():
    """Handles for the two oxide pies."""
    return ([plt.Rectangle((0, 0), 1, 1, facecolor=OXIDE_COLORS[o], edgecolor=SURFACE, lw=0.8)
             for o in PIPELINE_OXIDES], [OXIDE_LABEL[o] for o in PIPELINE_OXIDES])


def mineral_legend():
    """Handles for the mineral pie, kept separate: the two encodings are unrelated."""
    return ([plt.Rectangle((0, 0), 1, 1, facecolor=MINERAL_COLORS[m],
                           edgecolor=hatch_ink(MINERAL_COLORS[m]) if m in MINERAL_HATCH else SURFACE,
                           lw=0.8, hatch=MINERAL_HATCH.get(m)) for m in MINERALS], list(MINERALS))


def figure_single(mgsi, diw, data, out):
    fig, axs = plt.subplots(1, 3, figsize=(11.2, 4.5))
    fig.patch.set_facecolor(SURFACE)
    small = draw_cell(axs, mgsi, diw, data, THRESH_SINGLE)
    fig.text(0.5, 0.945, header(mgsi, diw, data), ha='center', fontsize=11, color=INK)
    fig.text(0.5, 0.895, 'mantle $\\rightarrow$ 20% batch melt at 1 GPa $\\rightarrow$ CIPW norm',
             ha='center', fontsize=8.4, color=INK_3)
    note = '    '.join(f'{pie}: ' + ', '.join(names) for pie, names in small.items() if names)
    if note:
        fig.text(0.5, 0.035, f'below {THRESH_SINGLE:g}% and unlabelled --- {note}',
                 ha='center', fontsize=7.4, color=INK_3)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.80, bottom=0.10, wspace=0.28)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=220, facecolor=SURFACE)
    plt.close(fig)
    return f'{out}.png'


def figure_grid(cells, out):
    # Pie axes are square, so a cell's height must track its 3-pie width or the rows fill with
    # dead space; the title and legend get their own bands rather than overlapping the top row.
    fig = plt.figure(figsize=(20.5, 10.4))
    fig.patch.set_facecolor(SURFACE)
    title_row, body, legend_row = fig.subfigures(3, 1, height_ratios=[0.135, 1.0, 0.075])
    for sub in (title_row, body, legend_row):
        sub.patch.set_facecolor(SURFACE)
    # Row 0 carries the per-pie titles, so it needs a little more height than the other two.
    panels = body.subfigures(len(GRID_DIW), len(GRID_MGSI), wspace=0.01, hspace=0.02,
                             height_ratios=[1.12, 1.0, 1.0])

    for iy, diw in enumerate(GRID_DIW):
        for ix, mgsi in enumerate(GRID_MGSI):
            sf = panels[iy][ix]
            sf.patch.set_facecolor(SURFACE)
            axs = sf.subplots(1, 3)
            draw_cell(axs, mgsi, diw, cells[(mgsi, diw)], THRESH_GRID,
                      pie_titles=(iy == 0), fontsize=6.6, min_dy=0.30)
            sf.suptitle(header(mgsi, diw, cells[(mgsi, diw)]), fontsize=9.4, color=INK,
                        y=0.97 if iy == 0 else 0.99)
            sf.subplots_adjust(left=0.035, right=0.965, top=0.74 if iy == 0 else 0.84,
                               bottom=0.02, wspace=0.32)

    # MORB in the title band, left of the heading: the same reference the mineral grid carries,
    # so the two figures are read against the same measured crust. It goes on the title SUBFIGURE
    # -- an axes added to the parent is painted over by the subfigures' own patches -- and its
    # width fraction is derived from the band's aspect so the pie comes out round.
    band_h_in = fig.get_size_inches()[1] * 0.135 / 1.21
    pie_in = 0.72
    draw_reference_pie(title_row,
                       [0.012, 0.20, pie_in / fig.get_size_inches()[0], pie_in / band_h_in],
                       label_size=6.4, title_size=7.4)
    title_row.text(0.5, 0.60, 'Mantle, primary melt and normative crust across the composition grid',
                   ha='center', fontsize=14, color=INK)
    title_row.text(0.5, 0.14, f'Pies left to right per cell: bulk mantle, 20% batch melt at 1 GPa, '
                   f'CIPW-normative crust.  Slices below {THRESH_GRID:g}% are unlabelled '
                   f'(the companion CSV carries every value).',
                   ha='center', fontsize=9, color=INK_2)
    # Two legends, not one: the oxide and mineral palettes are independent encodings, and a
    # single strip invites reading a shared colour as a relationship between them.
    ox_h, ox_lab = oxide_legend()
    mn_h, mn_lab = mineral_legend()
    left, right = legend_row.subfigures(1, 2, width_ratios=[1.0, 1.45])
    for sub in (left, right):
        sub.patch.set_facecolor(SURFACE)
    left.legend(ox_h, ox_lab, loc='center', ncol=4, frameon=False, fontsize=8.2,
                handlelength=1.1, columnspacing=1.3, title='Oxides (wt%)',
                title_fontproperties=dict(size=8.6, weight='bold'))
    right.legend(mn_h, mn_lab, loc='center', ncol=6, frameon=False, fontsize=8.2,
                 handlelength=1.1, columnspacing=1.3, title='Minerals (wt%)',
                 title_fontproperties=dict(size=8.6, weight='bold'))

    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return f'{out}.png'


def main(points, outdir, grid_only):
    _, axes = load_crust_table()
    grid = (axes['mg_si'], axes['delta_iw'])
    t_interp = RegularGridInterpolator(grid, axes['T_melt'], bounds_error=True)
    f_interp = RegularGridInterpolator(grid, axes['melt_fraction'], bounds_error=True)
    os.makedirs(outdir, exist_ok=True)
    rows = []

    def record(mgsi, diw, data):
        mantle, melt, minerals, t_melt, f = data
        rows.append(dict(mg_si=mgsi, delta_iw=diw, T_melt=t_melt, melt_fraction=f,
                         **{f'mantle_{k}': v for k, v in mantle.items()},
                         **{f'melt_{k}': v for k, v in melt.items()},
                         **{f'mineral_{k}': v for k, v in minerals.items()}))

    if not grid_only:
        for mgsi, diw in points:
            data = cell_data(mgsi, diw, t_interp, f_interp)
            record(mgsi, diw, data)
            tag = f'{mgsi:g}'.replace('.', 'p') + '_' + f'{diw:+g}'.replace('.', 'p')
            print('wrote', figure_single(mgsi, diw, data, os.path.join(outdir, f'crust_diag_{tag}')))

    cells = {(mgsi, diw): cell_data(mgsi, diw, t_interp, f_interp)
             for diw in GRID_DIW for mgsi in GRID_MGSI}
    for (mgsi, diw), data in cells.items():
        if grid_only or (mgsi, diw) not in points:
            record(mgsi, diw, data)
    print('wrote', figure_grid(cells, os.path.join(outdir, 'crust_diag_grid')))

    csv = os.path.join(outdir, 'crust_diag.csv')
    pd.DataFrame(rows).sort_values(['mg_si', 'delta_iw']).to_csv(csv, index=False)
    print(f'wrote {csv} ({len(rows)} cells)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default=None,
                    help='individual figures, "mg_si,delta_iw;..." (default: Earth + three others)')
    ap.add_argument('--outdir', default='output')
    ap.add_argument('--grid-only', action='store_true')
    a = ap.parse_args()
    pts = DEFAULT_POINTS if a.points is None else [
        tuple(float(v) for v in p.split(',')) for p in a.points.split(';') if p.strip()]
    main(pts, a.outdir, a.grid_only)
