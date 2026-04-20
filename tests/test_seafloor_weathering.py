import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt

from kamino.kamino_chem.ocean_chemistry import *
from kamino.constants import YR, EARTH_ATM


def plot_weathering_vs_pressure():
    """
    Plot seafloor weathering fluxes as a function of pore pressure
    at fixed T = 280 K and P_CO2 = 280 ppm * 1 bar.
    """
    T = 280.0                          # K
    P_CO2 = 280e-6 * EARTH_ATM        # Pa  (~28 Pa)

    # Pressure range: 1 bar (surface) to ~500 bar (deep ocean/seafloor)
    P_range = np.logspace(np.log10(EARTH_ATM), np.log10(500 * EARTH_ATM), 60)

    # Trace initial ocean composition (same as Planet initialisation)
    b_ocean = np.zeros(len(elements))
    for name, val in [('Alkalinity', 1e-3), ('C', 1e-3), ('Ca', 5e-4),
                      ('Mg', 5e-5), ('Si', 1e-4), ('Al', 1e-9), ('Fe', 1e-9)]:
        idx = int(np.where(elements == name)[0][0])
        b_ocean[idx] = val

    fluxes = []
    Da_vals = []
    supply_vals = []

    for P in P_range:
        flux, diag = get_weathering_flux(P, T, P_CO2, b_ocean)
        fluxes.append(flux.copy())
        Da_vals.append(diag['Da'])
        supply_vals.append(diag['supply_efficiency'])

    fluxes = np.array(fluxes)       # shape (n_P, n_elements)
    P_bar = P_range / EARTH_ATM

    fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    plt.subplots_adjust(hspace=0.05)

    # --- Panel 1: per-element fluxes ---
    ax1 = axes[0]
    for j, elem in enumerate(elements):
        f = fluxes[:, j] * YR   # mol / kgw / yr  (normalised by ocean mass elsewhere)
        if np.any(np.abs(f) > 0):
            ax1.plot(P_bar, f, label=elem, linewidth=1.5)
    ax1.set_yscale('symlog', linthresh=1e-10)
    ax1.set_ylabel('Weathering flux (mol kgw$^{-1}$ yr$^{-1}$)', fontsize=11)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.set_title(f'Seafloor weathering vs. pressure  (T = {T} K, $P_{{CO_2}}$ = 280 ppm)', fontsize=12)

    # --- Panel 2: Damköhler number ---
    ax2 = axes[1]
    ax2.plot(P_bar, Da_vals, color='tab:purple', linewidth=1.5)
    ax2.axhline(1, color='k', linestyle='--', linewidth=0.8, alpha=0.6, label='Da = 1')
    ax2.set_ylabel('Damköhler number (Da)', fontsize=11)
    ax2.set_yscale('log')
    ax2.legend(fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.5)

    # --- Panel 3: supply efficiency ---
    ax3 = axes[2]
    ax3.plot(P_bar, np.array(supply_vals) * 100, color='tab:green', linewidth=1.5)
    ax3.set_ylabel('Supply efficiency (%)', fontsize=11)
    ax3.set_xlabel('Pore pressure (bar)', fontsize=11)
    ax3.set_xscale('log')
    ax3.grid(True, linestyle='--', alpha=0.5)

    plt.savefig('weathering_vs_pressure.pdf', bbox_inches='tight')
    print('Saved weathering_vs_pressure.pdf')
    plt.show()


if __name__ == '__main__':
    plot_weathering_vs_pressure()