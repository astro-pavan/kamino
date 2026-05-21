import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

plt.style.use('experiments/planetary-chem-paper.mplstyle')

DEFAULT_OUTPUT_PATH = '/data/pt426/kamino_experiments/'

# (name, label, colour, pK1, pK2, is_base)
#
# is_base=True  → base gas (e.g. NH3): pK1 is pKa of conjugate acid (NH4+).
#                 Enhancement is large at LOW pH (gas converts to ionic form).
# is_base=False → acidic gas: ionises with increasing pH, so enhancement
#                 grows with pH.
# pK1=None      → no ionisation; H_eff/H0 = 1 at all pH (dashed reference).
#
# Gases with no simple acid-base speciation (N2O, O3, NO2):
#   N2O — inert; H_eff/H0 ≈ 1
#   NO2 — disproportionates in water (2 NO2 + H2O → HNO3 + HNO2); not
#          described by a single Henry + pKa, omitted.
#   O3  — reacts irreversibly with water; omitted.
#   HCl — pKa ≈ −7 (strong acid); effectively fully ionised at any ocean pH
#          so H_eff/H0 ≈ 10^(pH+7) → handled separately.
_GASES = [
    ('CO2', r'CO$_2$',  'tab:blue',   6.35,  10.33, False),
    ('SO2', r'SO$_2$',  'tab:orange', 1.85,   7.20, False),
    ('H2S', r'H$_2$S',  'tab:olive',  7.00,  13.90, False),
    ('HCN', r'HCN',     'tab:red',    9.21,   None, False),
    ('NH3', r'NH$_3$',  'tab:purple', 9.25,   None, True),
    ('HCl', r'HCl',     'tab:cyan',  -7.00,   None, False),  # strong acid
    ('CH4', r'CH$_4$',  'tab:green',  None,   None, False),  # inert reference
    ('N2O', r'N$_2$O',  'silver',     None,   None, False),  # inert reference
    ('CO',  r'CO',      'tab:brown',  None,   None, False),  # inert reference
]


def _h_eff_ratio(pH: np.ndarray, pK1, pK2, is_base: bool) -> np.ndarray:
    """Return H_eff / H0 (solubility enhancement from acid-base speciation).

    For an acidic gas AH with dissociation constants K1, K2:
        [total dissolved] = [AH] (1 + K1/[H+] + K1 K2/[H+]^2)
        H_eff / H0        = 1 + K1/[H+] + K1 K2/[H+]^2

    For a basic gas B (e.g. NH3) with conjugate acid pKa:
        [total dissolved] = [B] (1 + [H+]/Ka)
        H_eff / H0        = 1 + [H+] / Ka   (large at low pH)
    """
    if pK1 is None:
        return np.ones_like(pH)

    h = 10.0 ** (-pH)

    if is_base:
        Ka = 10.0 ** (-pK1)
        return 1.0 + h / Ka

    K1 = 10.0 ** (-pK1)
    ratio = 1.0 + K1 / h
    if pK2 is not None:
        K2 = 10.0 ** (-pK2)
        ratio += K1 * K2 / h ** 2
    return ratio


def plot_solubility_enhancement(output_path):
    """Plot H_eff/H0 vs ocean pH for a range of climatically relevant gases.

    H_eff/H0 is the factor by which acid-base speciation increases a gas's
    effective solubility above its bare Henry's law value.  A value of 10^N
    means the ocean holds 10^N × more of that gas than a non-ionising liquid
    would at the same partial pressure.  The curves are purely analytical
    (pKa values at 25 °C, 1 atm); absolute solubilities depend on H0.
    """
    pH = np.linspace(3.0, 12.0, 500)

    fig, ax = plt.subplots(figsize=(7, 5))

    for name, label, color, pK1, pK2, is_base in _GASES:
        ratio = _h_eff_ratio(pH, pK1, pK2, is_base)
        log_r = np.log10(np.maximum(ratio, 1e-30))

        lw = 1.8 if pK1 is not None else 1.0
        ls = '-'  if pK1 is not None else '--'
        ax.plot(pH, log_r, label=label, color=color, linewidth=lw, linestyle=ls)

    ax.axhline(0, color='k', linewidth=0.6, linestyle=':', alpha=0.4)
    ax.axvline(8.1, color='k', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.text(8.15, 0.3, 'modern\nocean', fontsize=7, va='bottom', alpha=0.7,
            transform=ax.get_xaxis_transform())

    # Mark pKa values as vertical ticks on the curves (optional annotation)
    for name, label, color, pK1, pK2, is_base in _GASES:
        if pK1 is not None and 3 < pK1 < 12:
            ratio_at_pK1 = _h_eff_ratio(np.array([pK1]), pK1, pK2, is_base)[0]
            ax.plot(pK1, np.log10(ratio_at_pK1), 'o', color=color,
                    markersize=4, zorder=5)
        if pK2 is not None and 3 < pK2 < 12:
            ratio_at_pK2 = _h_eff_ratio(np.array([pK2]), pK1, pK2, is_base)[0]
            ax.plot(pK2, np.log10(ratio_at_pK2), 's', color=color,
                    markersize=4, zorder=5)

    ax.set_xlabel('Ocean pH')
    ax.set_ylabel(r'$\log_{10}(H_\mathrm{eff}\,/\,H_0)$')
    ax.set_title('Solubility enhancement from ocean acid-base speciation  (25 °C)',
                 fontsize=10)
    ax.set_xlim(3, 12)
    ax.legend(ncol=3, fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.3)

    out_path = os.path.join(output_path, 'gas_solubility_pH.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Solubility enhancement of atmospheric gases vs ocean pH.')
    parser.add_argument('--path', default=DEFAULT_OUTPUT_PATH)
    args = parser.parse_args()

    os.makedirs(args.path, exist_ok=True)
    plot_solubility_enhancement(args.path)
