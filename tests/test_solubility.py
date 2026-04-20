import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

import numpy as np
import matplotlib.pyplot as plt
from phreeqc import Phreeqc

T_C = 288.0 - 273.15  # 14.85 °C

_here = os.path.dirname(__file__)
DB_PATH = os.path.abspath(os.path.join(_here, '..', 'ocean_chem.dat'))

# For each gas: which element to dissolve at 1 mM, and what pe to set.
#
# pe matters because PHREEQC distributes elements across oxidation states at
# thermodynamic equilibrium.  Without a suitable pe, oxidising conditions
# (default pe ≈ 4) convert H2S/SO2 → SO4²⁻, CH4/CO → CO2, and N2 → NH4⁺,
# making their SI values meaningless.  The pe values here keep each element in
# the oxidation state relevant to the gas being tested.
#
# Rough stability boundaries at pH 7 (from redox half-reactions):
#   CO2 stable vs CH4 : pe > -2     → use pe = 4 (default) for CO2
#   CH4 stable vs CO2 : pe < -2     → use pe = -12 for CH4
#   CO  stable        : intermediate → use pe = -8  for CO
#   H2S/SO2 stable vs SO4 : pe < -4 → use pe = -8  for H2S, SO2
#   N2  stable vs NH4+: pe < -5     → use pe = -5  for N2
#   NH3/NH4+ acid-base only (pe-insensitive in this range): default
#   NO2 → NO3⁻ (oxidising product): default
#   O2  stable vs H2O : pe > 14     → use pe = 14  for O2
GAS_CONFIG: dict[str, dict] = {
    'CO2(g)': {'element': 'C',      'conc': 1e-3, 'pe': None},
    'SO2(g)': {'element': 'S',      'conc': 1e-3, 'pe': -8},
    'H2S(g)': {'element': 'S',      'conc': 1e-3, 'pe': -8},
    'CH4(g)': {'element': 'C(-4)',  'conc': 1e-3, 'pe': -12},
    'CO(g)':  {'element': 'C(2)',   'conc': 1e-3, 'pe': -8},
    'N2(g)':  {'element': 'N(0)',   'conc': 1e-3, 'pe': -5},
    'NH3(g)': {'element': 'N(-3)',  'conc': 1e-3, 'pe': None},
    'NO2(g)': {'element': 'N(5)',   'conc': 1e-3, 'pe': None},
    'O2(g)':  {'element': 'O(0)',   'conc': 1e-3, 'pe': 14},
}

REACTIVE_GASES = ['CO2(g)', 'SO2(g)', 'H2S(g)', 'NH3(g)', 'NO2(g)']
INERT_GASES    = ['N2(g)', 'CH4(g)', 'CO(g)', 'O2(g)']

PH_RANGE = np.linspace(4, 10, 40)


def make_input(gas: str, pH: float) -> str:
    """
    SOLUTION block with pH fixed and one element dissolved at 1 mM.
    No EQUILIBRIUM_PHASES — we just do speciation and read SI.
    pGas_eq = 10^SI(gas) atm is the gas partial pressure that would be in
    equilibrium with this solution.
    """
    cfg = GAS_CONFIG[gas]
    lines = [
        'SOLUTION 1',
        f'    temp    {T_C:.4f}',
        f'    pH      {pH:.6f}',
    ]
    if cfg['pe'] is not None:
        lines.append(f'    pe      {cfg["pe"]:.2f}')
    lines += [
        '    units   mol/kgw',
        f'    {cfg["element"]}   {cfg["conc"]:.6e}',
        '',
        'SELECTED_OUTPUT',
        '    -pH',
        f'    -saturation_indices  {gas}',
        '',
    ]
    return '\n'.join(lines)


def run_gas(gas: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (pH_array, SI_array).  SI = log10(p_gas / atm).
    """
    p = Phreeqc()
    if p.load_database(DB_PATH) == 1:
        raise RuntimeError(f"Failed to load database:\n{p.get_error_string()}")

    si_key = f'si_{gas}'
    phs, sis = [], []

    for pH in PH_RANGE:
        try:
            if p.run_string(make_input(gas, pH)) == 1:
                raise RuntimeError(p.get_error_string())
            out = p.get_selected_output()
            si = float(out[si_key][-1])
        except Exception as e:
            print(f'  Warning: {gas} pH={pH:.2f}: {e}')
            si = np.nan
        phs.append(pH)
        sis.append(si)

    return np.array(phs), np.array(sis)


def run_all() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    results = {}
    for gas in REACTIVE_GASES + INERT_GASES:
        print(f'Running {gas}...')
        results[gas] = run_gas(gas)
    return results


def plot_results(results: dict):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        'Equilibrium gas pressure vs ocean pH  (T = 288 K, [element] = 1 mmol/kgw)',
        fontsize=12,
    )

    for i, gas in enumerate(REACTIVE_GASES):
        ax = axes[i // 3, i % 3]
        phs, sis = results[gas]
        mask = np.isfinite(sis)
        ax.plot(phs[mask], sis[mask], 'o-', markersize=3, linewidth=1)
        ax.set_xlabel('pH')
        ax.set_ylabel('log₁₀(p / atm)')
        cfg = GAS_CONFIG[gas]
        pe_note = f'  pe={cfg["pe"]}' if cfg['pe'] is not None else ''
        ax.set_title(gas.replace('(g)', '') + pe_note)
        ax.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.4)
        ax.grid(True, alpha=0.3)

    ax_inert = axes[1, 2]
    for gas in INERT_GASES:
        phs, sis = results[gas]
        mask = np.isfinite(sis)
        cfg = GAS_CONFIG[gas]
        pe_note = f' (pe={cfg["pe"]})' if cfg['pe'] is not None else ''
        label = gas.replace('(g)', '') + pe_note
        ax_inert.plot(phs[mask], sis[mask], 'o-', markersize=3, linewidth=1, label=label)
    ax_inert.set_xlabel('pH')
    ax_inert.set_ylabel('log₁₀(p / atm)')
    ax_inert.set_title('Inert gases')
    ax_inert.axhline(0, color='k', linewidth=0.5, linestyle='--', alpha=0.4)
    ax_inert.legend(fontsize=8)
    ax_inert.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(_here, '..', 'solubility_plot.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    results = run_all()
    plot_results(results)
