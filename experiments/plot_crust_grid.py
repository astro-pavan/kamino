"""Small-multiples map of crust mineralogy over the two composition axes.

One pie per (mantle Mg/Si, core-formation dIW) cell of crust_compositions.csv, showing the
CIPW-normative assemblage the weathering model actually consumes. Reads the table through
kamino.crust_composition so the figure can never disagree with the model.

    python experiments/plot_crust_grid.py [--out output/crust_grid]

Design notes
------------
Pie charts are a weak form for comparing values, and they are used here only because the
question is "how does the assemblage change across the grid", where each panel is read as a
shape rather than as eight numbers. Three things carry identity so it is never colour alone:

  * slice order is IDENTICAL in every panel, so angular position encodes mineral;
  * a legend is always present, and each panel is direct-labelled with its melt SiO2;
  * mineral weight fractions are written to a companion CSV as the table view.

Colours are the eight documented categorical hues, assigned so that RING-ADJACENT pairs (the
ones that touch, which is the pairlist a fixed-order pie actually exposes) clear the CVD floor
under protanopia, deuteranopia AND tritanopia: worst pair dE 9.7 (target 8), worst
normal-vision dE 20.8 (floor 15). The natural slot order fails tritan at 4.9 on
Anorthite/Diopside; swapping the Anorthite and Forsterite hues fixes it using the same eight
documented values, with no re-stepping. Do not reorder the minerals without re-checking.
"""

import argparse
import os
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

from kamino.constants import EARTH_DELTA_IW, EARTH_MANTLE_MG_SI
from kamino.crust_composition import CRUST_TABLE, cipw_norm

# Fixed slice order, silica-rich -> silica-poor, with the albite->nepheline cascade pair
# adjacent. Changing this changes which colour pairs touch; see the module docstring.
MINERALS = ['Quartz', 'Albite', 'Nepheline', 'Anorthite',
            'Diopside', 'Akermanite', 'Enstatite', 'Forsterite', 'Fayalite']
COLORS = {'Quartz': '#2a78d6', 'Albite': '#eb6834', 'Nepheline': '#1baf7a',
          'Anorthite': '#4a3aa7', 'Diopside': '#e87ba4', 'Akermanite': '#e87ba4',
          'Enstatite': '#008300', 'Forsterite': '#eda100', 'Fayalite': '#e34948'}
# Nine minerals, eight documented hues. Rather than invent a ninth, Akermanite takes
# Diopside's hue with a hatch -- composite encoding, and semantically exact: akermanite IS
# the diopside desilication product (norm step 5c). The hue ring is therefore unchanged and
# its CVD validation still holds; the two are adjacent, so the boundary is a texture edge.
HATCHED = {'Akermanite': '///'}

INK, INK_2, INK_3 = '#1a1a19', '#55554e', '#8a8a80'
SURFACE, MUTED = '#ffffff', '#e8e8e2'

# Subsampled so the panels stay legible; every value is an exact grid node (no interpolation).
MGSI_SHOW = [0.5, 0.7, 0.9, 1.0, 1.1, 1.25, 1.4, 1.6, 1.8, 2.0]
DIW_SHOW = [-5.0, -4.0, -3.0, -2.0, -1.0]


def assemblage(row):
    """(weight-fraction dict, is_mass_violating) for one table row."""
    ox = {'SiO2': row.SiO2, 'TiO2': row.TiO2, 'Al2O3': row.Al2O3, 'FeOt': row.FeO,
          'MgO': row.MgO, 'CaO': row.CaO, 'Na2O': row.Na2O, 'K2O': row.K2O}
    # cipw_norm RAISES when the phase set cannot express the melt (a residual silica deficit
    # after both desilication cascades). That is a real exclusion, not a warning, so catch it
    # and mark the cell rather than letting it kill the figure.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return cipw_norm(ox, emit_quartz=True), False
    except ValueError:
        return {}, True


def main(table, out):
    df = pd.read_csv(table, comment='#')
    nx, ny = len(MGSI_SHOW), len(DIW_SHOW)

    fig, axes = plt.subplots(ny, nx, figsize=(1.15 * nx + 2.6, 1.16 * ny + 1.9))
    fig.patch.set_facecolor(SURFACE)

    records = []
    for iy, diw in enumerate(reversed(DIW_SHOW)):        # oxidised at the top
        for ix, mgsi in enumerate(MGSI_SHOW):
            ax = axes[iy][ix]
            ax.set_facecolor(SURFACE)
            for s in ax.spines.values():
                s.set_visible(False)
            ax.set_xticks([]); ax.set_yticks([])

            sel = df[np.isclose(df.mg_si, mgsi) & np.isclose(df.delta_iw, diw)]
            if sel.empty:
                continue
            row = sel.iloc[0]
            c, bad = assemblage(row)
            records.append(dict(mg_si=mgsi, delta_iw=diw, mantle_feo=row.mantle_feo,
                                melt_SiO2=row.SiO2, mass_violating=bad,
                                **{m: c.get(m, 0.0) for m in MINERALS}))

            if bad:
                # Mass-violating: the norm assigned more SiO2 than the melt contains, so the
                # assemblage is not a real composition. Show the cell as excluded rather than
                # drawing a plausible-looking pie of wrong numbers.
                ax.add_patch(Circle((0, 0), 1.0, facecolor=MUTED, edgecolor='none'))
                ax.plot([-.45, .45], [-.45, .45], color=INK_3, lw=1.4, zorder=3)
                ax.plot([-.45, .45], [.45, -.45], color=INK_3, lw=1.4, zorder=3)
                ax.set_xlim(-1.30, 1.30); ax.set_ylim(-1.42, 1.14)
                ax.set_aspect('equal')
                continue

            keep = [(m, c.get(m, 0.0)) for m in MINERALS if c.get(m, 0.0) > 1e-4]
            wedges, _ = ax.pie([v for _, v in keep], colors=[COLORS[m] for m, _ in keep],
                               startangle=90, counterclock=False, radius=1.0,
                               wedgeprops=dict(edgecolor=SURFACE, linewidth=1.1))
            for wedge, (m, _) in zip(wedges, keep):
                if m in HATCHED:
                    wedge.set_hatch(HATCHED[m])
                    wedge.set_edgecolor(INK)
            # Selective direct label: melt silica is the scalar that separates the regimes.
            # The Earth cell folds its name into the same label rather than adding a second one
            # above the pie, which collided with the row above once the rows were tightened.
            is_earth = np.isclose(mgsi, EARTH_MANTLE_MG_SI) and np.isclose(diw, EARTH_DELTA_IW)
            ax.text(0, -1.28, f'Earth  {row.SiO2:.0f}' if is_earth else f'{row.SiO2:.0f}',
                    ha='center', va='center', fontsize=6.6,
                    color=INK if is_earth else INK_2,
                    fontweight='bold' if is_earth else 'normal')
            ax.set_xlim(-1.30, 1.30); ax.set_ylim(-1.42, 1.14)

            if is_earth:
                ax.add_patch(Circle((0, 0), 1.16, fill=False, edgecolor=INK, lw=1.6, zorder=4))

    for ix, mgsi in enumerate(MGSI_SHOW):
        axes[0][ix].set_title(f'{mgsi:g}', fontsize=8.5, color=INK, pad=7)
    for iy, diw in enumerate(reversed(DIW_SHOW)):
        feo = df[np.isclose(df.delta_iw, diw)].mantle_feo.iloc[0]
        axes[iy][0].set_ylabel(f'{diw:+g}', fontsize=8.5, color=INK,
                               rotation=0, ha='right', va='center', labelpad=12)
        axes[iy][-1].text(1.75, 0, f'{feo:.2f}', fontsize=7.5, color=INK_2,
                          ha='left', va='center', transform=axes[iy][-1].transData)

    fig.text(0.5, 0.962, 'Normative crust mineralogy across the composition grid',
             ha='center', fontsize=13, color=INK)
    fig.text(0.5, 0.929,
             'MAGEMin primary melts at 20% melting, through the CIPW norm.  '
             'Number under each pie is melt SiO$_2$ (wt%).  '
             'Hatched = Akermanite (proxied kinetics).',
             ha='center', fontsize=8.5, color=INK_2)
    fig.text(0.5, 0.086, 'Mantle molar Mg/Si', ha='center', fontsize=10, color=INK)
    fig.text(0.018, 0.5, 'Core-formation $\\Delta$IW', va='center', rotation=90,
             fontsize=10, color=INK)
    fig.text(0.982, 0.5, 'Mantle FeO (wt%)', va='center', rotation=-90,
             fontsize=9, color=INK_2)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=COLORS[m],
                             edgecolor=INK if m in HATCHED else SURFACE, lw=0.8,
                             hatch=HATCHED.get(m)) for m in MINERALS]
    labels = list(MINERALS)
    if any(r['mass_violating'] for r in records):
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=MUTED, edgecolor='none'))
        labels.append('beyond norm')
    fig.legend(handles, labels, loc='lower center', ncol=len(labels), frameon=False,
               fontsize=8, bbox_to_anchor=(0.5, 0.002), handlelength=1.1, columnspacing=1.3)

    fig.subplots_adjust(left=0.072, right=0.945, top=0.885, bottom=0.150,
                        wspace=0.10, hspace=0.02)

    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=220, facecolor=SURFACE)
    pd.DataFrame(records).to_csv(f'{out}.csv', index=False)
    print(f'wrote {out}.png / .pdf and {out}.csv ({len(records)} cells, '
          f'{sum(r["mass_violating"] for r in records)} excluded)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=CRUST_TABLE)
    ap.add_argument('--out', default='output/crust_grid')
    a = ap.parse_args()
    main(a.table, a.out)
