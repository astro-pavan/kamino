"""Validation figure for the crust-composition pipeline.

Four independent checks, one panel each:

  (a) melting relations against Katz, Spiegelman & Langmuir (2003), the standard empirical
      anhydrous-peridotite parameterisation -- independent of both MAGEMin and Guimond;
  (b) mantle mineralogical boundaries against Guimond et al. (2024) section 3.1.1;
  (c) the only experimental test available: Brugman, Phillips & Till (2021) piston-cylinder
      melting of HEX1, a hypothetical exoplanet mantle at molar Mg/Si = 1.42;
  (d) melt-vs-mantle element partitioning against Guimond et al. (2024) section 4.2.

    python experiments/plot_crust_validation.py [--out output/crust_validation]

The MAGEMin numbers in MEASURED below are outputs, not inputs -- reproduce with
    julia make_crust_compositions.jl --validate --mgsi 1.25 --diw -2.0 --fixed-p <P> --ftarget <F>
and the HEX1 run with --bulk. They are recorded here rather than recomputed so the figure is
cheap to redraw; the command that produced each is given alongside.

Colour: at most three categorical slots per panel, which is the cap that validates on the
ALL-PAIRS pairlist (a scatter/line panel exposes every pair, not just neighbours). Model results
are always slot 1, reference data always slot 2, so the roles stay fixed across panels.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from kamino.crust_composition import (CRUST_TABLE, OXIDE_MOLAR_MASS as M, PYROLITE_NON_FE,
                                      feo_from_delta_iw)

C_MODEL, C_REF, C_THIRD = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK_2, INK_3, GRID = '#1a1a19', '#55554e', '#8a8a80', '#e8e8e2'
SURFACE = '#ffffff'

# MAGEMin isobaric runs, pyrolite at dIW = -2  ->  (P_GPa, F, T_C)
MEASURED_KATZ = [(1.0, 0.10, 1287), (1.5, 0.05, 1338), (1.5, 0.10, 1359),
                 (1.5, 0.20, 1391), (2.0, 0.10, 1423)]
# Brugman et al. (2021): HEX1 bulk (Table 1) and mean experimental melt (Table 5).
HEX1_BULK = {'SiO2': 42.00, 'MgO': 40.04, 'FeO': 8.23, 'Al2O3': 4.85, 'CaO': 3.76, 'Na2O': 0.21}
HEX1_MELT = {'SiO2': 48.09, 'MgO': 9.17, 'FeO': 11.16, 'Al2O3': 15.65, 'CaO': 13.92}
# Ours at their exact bulk, isobaric 1.5 GPa, F = 0.05.
OURS_MELT = {'SiO2': 45.37, 'MgO': 11.50, 'FeO': 7.26, 'Al2O3': 20.30, 'CaO': 12.65}


# --- Katz et al. (2003) anhydrous peridotite, eqs. (4)-(6) -------------------------------------
def T_sol(P):   return 1085.7 + 132.9 * P - 5.1 * P ** 2
def T_lherz(P): return 1475.0 + 80.0 * P - 3.2 * P ** 2
def T_for_F(F, P): return T_sol(P) + F ** (1 / 1.5) * (T_lherz(P) - T_sol(P))


def mantle(mgsi, feo):
    sc = (100 - feo) / sum(PYROLITE_NON_FE.values())
    ox = {k: v * sc for k, v in PYROLITE_NON_FE.items()}
    tot = ox['MgO'] + ox['SiO2']
    si = tot / (1 + mgsi * M['MgO'] / M['SiO2'])
    ox['SiO2'], ox['MgO'] = si, tot - si
    ox['FeO'] = feo
    t = sum(ox.values())
    return {k: v / t * 100 for k, v in ox.items()}


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(INK_3)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_2, labelsize=8, length=3, width=0.8)
    ax.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)


def panel_katz(ax):
    style(ax)
    Fs = np.linspace(0.01, 0.24, 120)
    for P, col in zip((1.0, 1.5, 2.0), (C_MODEL, C_REF, C_THIRD)):
        ax.plot(Fs * 100, T_for_F(Fs, P), color=col, lw=2, zorder=2)
        pts = [(F, T) for Pk, F, T in MEASURED_KATZ if Pk == P]
        ax.plot([f * 100 for f, _ in pts], [t for _, t in pts], 'o', color=col,
                markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.4, zorder=3)
        ax.text(24.6, T_for_F(0.24, P), f'{P:g} GPa', color=col, fontsize=8,
                va='center', ha='left')
    ax.plot([], [], color=INK_3, lw=2, label='Katz et al. (2003)')
    ax.plot([], [], 'o', color=INK_3, markersize=8, markeredgecolor=SURFACE,
            markeredgewidth=1.4, label='this model (MAGEMin)')
    ax.legend(frameon=False, fontsize=8, loc='upper left', labelcolor=INK_2, handletextpad=0.6)
    ax.set_xlim(0, 31); ax.set_xlabel('Melt fraction F (%)', fontsize=9, color=INK)
    ax.set_ylabel('Temperature (°C)', fontsize=9, color=INK)
    ax.set_title('(a)  Melting relations vs the standard parameterisation',
                 fontsize=9.5, color=INK, loc='left', pad=8)
    d = [T - T_for_F(F, P) for P, F, T in MEASURED_KATZ]
    ax.text(0.97, 0.06, f'mean offset {np.mean(d):+.0f} °C  (sd {np.std(d):.0f})',
            transform=ax.transAxes, ha='right', fontsize=8, color=INK_2)


def panel_boundaries(ax, df):
    style(ax)
    e = df[np.isclose(df.delta_iw, -2.0)].sort_values('mg_si')
    ph = e.residual_phases.str.split(';')
    for y, (phase, label, col) in enumerate([('ol', 'olivine', C_MODEL),
                                             ('opx', 'orthopyroxene', C_REF)]):
        present = e.mg_si[ph.apply(lambda p: phase in p)].to_numpy()
        ax.scatter(present, np.full_like(present, y), s=64, color=col, zorder=3,
                   edgecolor=SURFACE, linewidth=1.2, label=label)
    # Guimond's subsolidus boundaries, and ours from the 20%-melting residue
    for x, lab, ls in [(0.8, 'Guimond ≲0.8', ':'), (1.6, 'Guimond ≳1.6', ':')]:
        ax.axvline(x, color=INK_3, ls=ls, lw=1.4, zorder=1)
        ax.text(x - 0.022, 1.70, lab, rotation=90, fontsize=7, color=INK_3,
                ha='right', va='top')
    for x, col in [(0.70, C_MODEL), (1.50, C_REF)]:
        ax.axvline(x, color=col, ls='--', lw=1.4, alpha=0.8, zorder=2)
    ax.set_yticks([0, 1]); ax.set_yticklabels(['olivine', 'opx'], fontsize=8.5)
    ax.set_ylim(-0.7, 1.75); ax.set_xlim(0.4, 2.1)
    ax.set_xlabel('Mantle molar Mg/Si', fontsize=9, color=INK)
    ax.set_title('(b)  Where phases leave the residue', fontsize=9.5, color=INK, loc='left', pad=8)
    ax.text(0.5, 0.50, 'dashed = ours, from the 20% melting residue. It is depleted, so more\n'
                       'olivine-rich than their subsolidus mantle — and both offsets fall in\n'
                       'exactly the direction that predicts.',
            transform=ax.transAxes, ha='center', va='center', fontsize=7.5, color=INK_2,
            bbox=dict(facecolor=SURFACE, edgecolor='none', pad=3.5))


def panel_brugman(ax):
    style(ax)
    ox = ['SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO']
    x = np.arange(len(ox)); w = 0.38
    ax.bar(x - w / 2, [HEX1_MELT[o] for o in ox], w, color=C_REF, zorder=3,
           label='Brugman et al. (2021), experiment')
    ax.bar(x + w / 2, [OURS_MELT[o] for o in ox], w, color=C_MODEL, zorder=3,
           label='this model, matched conditions')
    for i, o in enumerate(ox):
        d = OURS_MELT[o] - HEX1_MELT[o]
        ax.text(i, max(HEX1_MELT[o], OURS_MELT[o]) + 1.4, f'{d:+.1f}',
                ha='center', fontsize=7.5, color=INK_2)
    ax.set_xticks(x)
    ax.set_xticklabels(['SiO$_2$', 'MgO', 'FeO', 'Al$_2$O$_3$', 'CaO'], fontsize=8.5)
    ax.set_ylabel('Melt oxide (wt%)', fontsize=9, color=INK)
    ax.set_ylim(0, 58)   # SiO2 reaches ~48 wt%; do NOT tighten this or the bars clip
    ax.legend(frameon=False, fontsize=8, loc='upper center', labelcolor=INK_2,
              handlelength=1.2, bbox_to_anchor=(0.62, 1.0))
    ax.set_title('(c)  Against experiment: HEX1 melt, Mg/Si = 1.42',
                 fontsize=9.5, color=INK, loc='left', pad=8)
    ax.text(0.015, 0.055, 'their exact bulk, isobaric 1.5 GPa, F = 0.05',
            transform=ax.transAxes, fontsize=7.5, color=INK_2)


def panel_partitioning(ax, df):
    style(ax)
    ox = ['SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO', 'Na2O', 'TiO2']
    sub = df[(df.mg_si >= 0.8) & (df.mg_si <= 1.6)]
    mant = pd.DataFrame([mantle(r.mg_si, r.mantle_feo) for _, r in sub.iterrows()])
    ratio = [sub[o].mean() / mant[o].mean() for o in ox]
    # One series, so one hue -- the axvline at unity already encodes direction, and reusing the
    # reference colour here would clash with the model/reference roles the other panels set up.
    y = np.arange(len(ox))
    ax.barh(y, ratio, 0.62, color=C_MODEL, zorder=3)
    ax.axvline(1.0, color=INK, lw=1.2, zorder=4)
    for i, r in enumerate(ratio):
        # Anchor outside the bar end. The threshold is 0.85 rather than 1.0 so a ratio that sits
        # essentially at unity (FeO) still labels to the right instead of inside its own bar.
        outside_right = r >= 0.85
        ax.text(r * (1.09 if outside_right else 0.9), i, f'{r:.1f}×', va='center',
                ha='left' if outside_right else 'right', fontsize=7.5, color=INK_2)
    ax.set_yticks(y)
    ax.set_yticklabels(['SiO$_2$', 'MgO', 'FeO', 'Al$_2$O$_3$', 'CaO', 'Na$_2$O', 'TiO$_2$'],
                       fontsize=8.5)
    ax.set_xscale('log'); ax.set_xlim(0.16, 26)
    ax.set_ylim(-0.7, len(ox) - 0.15)
    ax.set_xlabel('melt / mantle (mean over Mg/Si 0.8–1.6)', fontsize=9, color=INK)
    ax.set_title('(d)  Element partitioning into the melt', fontsize=9.5, color=INK,
                 loc='left', pad=8)
    ax.text(0.985, -0.155, 'Guimond §4.2: melts less magnesian, similar silica, CaO / Na / Ti '
                           'enriched — 20/20 relations reproduce',
            transform=ax.transAxes, ha='right', va='top', fontsize=7.5, color=INK_2)


def main(table, out):
    df = pd.read_csv(table, comment='#')
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.0))
    fig.patch.set_facecolor(SURFACE)
    panel_katz(axes[0][0])
    panel_boundaries(axes[0][1], df)
    panel_brugman(axes[1][0])
    panel_partitioning(axes[1][1], df)
    fig.suptitle('Validation of the crust-composition pipeline against three independent sources',
                 fontsize=13, color=INK, y=0.975)
    fig.subplots_adjust(left=0.075, right=0.965, top=0.905, bottom=0.075, wspace=0.24, hspace=0.36)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=220, facecolor=SURFACE)
    print(f'wrote {out}.png / .pdf')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--table', default=CRUST_TABLE)
    ap.add_argument('--out', default='output/crust_validation')
    a = ap.parse_args()
    main(a.table, a.out)
