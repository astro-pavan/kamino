from kamino.constants import *
from kamino.kamino_chem.mineral_rates import K_FUNCTIONS
from kamino.kamino_chem.mineral_info import *

from phreeqc import Phreeqc
import numpy as np
import numpy.typing as npt
import re
import pandas as pd

database_path = 'ocean_chem.dat'

elements = np.array([
    'Alkalinity',
    'C',
    'Si',
    'Al',
    'Fe',
    'Ca',
    'Mg',
    'Na',
    #'K',
    #'S',
    #'N',
    #'F',
    #'Cl'
])

element_string = ' '.join(elements[1:])

minerals = []

try:
    with open(database_path, 'r') as f:
        in_phases = False
        for line in f:
            if line.startswith('PHASES'):
                in_phases = True
                continue
            if line.startswith('RATES'):
                in_phases = False
                continue
            if in_phases:
                if line.startswith('\t') or line.startswith('  '):
                    continue
                if line.strip() == '':
                    continue
                if line[0].isalpha():
                    mineral_name = line.split()[0]
                    minerals.append(mineral_name)
except FileNotFoundError:
    print(f"Database file {database_path} not found")

available_mineral_string = ' '.join(minerals)

# Parses a PHREEQC database file to extract the stoichiometry and alkalinity generation of each mineral upon dissolution.

# Mapping PHREEQC aqueous species back to elemental buckets
species_to_element = {
    'Ca+2': 'Ca',
    'Mg+2': 'Mg',
    'Na+': 'Na',
    # 'K+': 'K',
    'Fe+2': 'Fe',
    'Fe+3': 'Fe',
    'Al+3': 'Al',
    'SiO2': 'Si',
    'HCO3-': 'C',
    'CO3-2': 'C',
    #'SO4-2': 'S',
    #'HS-': 'S',
    #'Cl-': 'Cl',
    #'F-': 'F'
}

# Contribution of 1 mole of aqueous species to total alkalinity
alk_contributions = {
    'HCO3-': 1.0,
    'CO3-2': 2.0,
    'H+': -1.0,
    'OH-': 1.0,
    'HS-': 1.0
}

def parse_stoichiometry(filepath: str) -> dict[str, npt.NDArray[np.float64]]:
    """Parse mineral reaction stoichiometry from a PHREEQC database PHASES section.

    Works with both reaction formats found in PHREEQC databases:
      - Indented reactions (ocean_chem.dat / SUPCRT92 style)
      - Unindented reactions (hydrothermal.dat / SupPHREEQC bl-1kb.dat style)
    """

    def parse_side(side_str: str, stoich: dict, multiplier: float, is_lhs: bool = False) -> None:
        terms = [t.strip() for t in re.split(r'\s\+\s', side_str.strip()) if t.strip()]
        if is_lhs and terms:
            terms = terms[1:]  # first LHS term is the mineral formula itself
        for term in terms:
            parts = term.split()
            try:
                coeff, species = (float(parts[0]), parts[1]) if len(parts) > 1 else (1.0, parts[0])
            except (ValueError, IndexError):
                coeff, species = 1.0, parts[0]
            if species in species_to_element:
                stoich[species_to_element[species]] += coeff * multiplier
            if species in alk_contributions:
                stoich['Alkalinity'] += coeff * alk_contributions[species] * multiplier

    result: dict[str, npt.NDArray[np.float64]] = {}
    in_phases = False
    current_mineral = None

    with open(filepath, 'r') as f:
        for raw_line in f:
            line = raw_line.split('#')[0].strip()
            if not line:
                continue
            if line == 'PHASES':
                in_phases = True
                continue
            elif line in ('RATES', 'SOLUTION_SPECIES', 'END'):
                in_phases = False
                continue
            if not in_phases:
                continue

            if '=' not in line and not line.startswith('-') and not raw_line[0].isspace():
                current_mineral = line.split()[0]
            elif current_mineral and '=' in line and 'log_k' not in line and 'delta_H' not in line:
                stoich: dict[str, float] = {str(el): 0.0 for el in elements}
                lhs_str, rhs_str = line.split('=', 1)
                parse_side(lhs_str, stoich, multiplier=-1.0, is_lhs=True)
                parse_side(rhs_str, stoich, multiplier=1.0)
                result[current_mineral] = np.array([round(v, 6) for v in stoich.values()])
                current_mineral = None

    return result

stoichiometry: dict[str, npt.NDArray[np.float64]] = {}
stoichiometry.update(parse_stoichiometry(database_path))
stoichiometry.update(parse_stoichiometry('hydrothermal.dat'))

p_LT = Phreeqc()
if p_LT.load_database("ocean_chem.dat") == 1:
    print(p_LT.get_error_string())

p_HT = Phreeqc()
p_HT.load_database("hydrothermal.dat")

def solve_solution(P: float, T: float, b: npt.NDArray[np.float64], pH: float | None=None, P_CO2: float | None=None, precipitating_minerals: list[str]=[], equilibriating_minerals: list[str]=[], high_temperature: bool=False, fO2: float=0, trace_approximation: bool=False, verbose: bool=False):

    solution_lines: list[str] = [
        # KNOBS: increase iteration limit for robustness with extreme-K minerals at low T
        'KNOBS',
        '    -iterations 2000',
        '',
        'SOLUTION 1',
        f'    pressure  {P / EARTH_ATM:.4f}',
        f'    temp      {T + ABSOLUTE_ZERO:.4f}',
        '    units     mol/kgw'
    ]

    if pH is not None:
        solution_lines.append(f'    pH     {pH:5f}')

    for element, x in zip(elements, b):
        molality = 0 if x < 1e-9 and trace_approximation else x
        solution_lines.append(f'    {element}    {molality:.15e}')

    solution_lines.append('')

    equilibrium_lines: list[str] = []

    if equilibriating_minerals or precipitating_minerals or P_CO2 is not None or fO2 > 0:

        equilibrium_lines.append('EQUILIBRIUM_PHASES 1')

        if P_CO2 is not None:
            equilibrium_lines.append(f'    CO2(g)  {np.log10(P_CO2 / EARTH_ATM):.4f}  {1e6}')

        if fO2 > 0:
            equilibrium_lines.append(f'    O2(g)   {np.log10(fO2 / EARTH_ATM):.4f}  {1e6}')

        for phase in equilibriating_minerals:
            equilibrium_lines.append(f'    {phase}  {0.0}  {100.0}')

        equilibrium_lines.append('')

        for phase in precipitating_minerals:
            equilibrium_lines.append(f'    {phase}  {0.0}  {0.0}')

        equilibrium_lines.append('')

    output_lines: list[str] = [
        'SELECTED_OUTPUT',
        '    -pH',
        '    -pe',
        '    -alkalinity',
        f'    -totals {element_string}',
        f'    -saturation_indices {hydrothermal_mineral_string if high_temperature else available_mineral_string}'
        ]
    
    X_combined = equilibriating_minerals + precipitating_minerals
    if X_combined:
        output_lines.append('    -equilibrium_phases ' + ' '.join(X_combined))
    
    input_lines = solution_lines + equilibrium_lines + output_lines + ['\n']
    input_string = '\n'.join(input_lines)

    if verbose:
        print(input_string)

    p = p_HT if high_temperature else p_LT 

    if p.run_string(input_string) == 1:
        print(p.get_error_string())

    output_dict = p.get_selected_output()
    if output_dict and verbose:
        df = pd.DataFrame(output_dict)
        df = df.loc[:, (df != -999.999).any()]
        pd.set_option('display.float_format', '{:.2e}'.format)
        print(df.T)

    if output_dict:
        return output_dict
    else:
        raise ValueError

def get_k(P: float, T: float, pH: float, composition: dict[str, float]) -> npt.NDArray[np.float64]:
    
    k = np.zeros(elements.shape)

    for mineral, mineral_fraction in composition.items():
        k_eff = K_FUNCTIONS[mineral](T, pH)
        k += k_eff * stoichiometry[mineral] * mineral_fraction

    return k

def get_b_eq(P: float, T: float, P_CO2: float, composition: dict[str, float], b_input: npt.NDArray[np.float64] | None=None, precipitating_minerals: list[str]=[], high_temperature: bool=False, fO2: float=0) -> tuple[npt.NDArray[np.float64], float]:

    b_eq = np.zeros(elements.shape)

    # Anorthite and Forsterite have extreme PHREEQC log_k values at seafloor T (~31 at 1.85°C)
    # due to SUPCRT92 temperature extrapolation. Equilibrating them from dilute water causes
    # numerical divergence. Exclude them for the low-T database only; the high-T database
    # (bl-1kb.dat) has proper analytic fits and does not have this problem.
    # Anorthite drives Al to molar levels without an Al-sink (Kaolinite is not in
    # hydrothermal.dat). Exclude it for both databases. Forsterite is only problematic
    # with SUPCRT92 low-T extrapolation; keep it for the HT database.

    _problematic = {'Anorthite'} if high_temperature else {'Anorthite', 'Forsterite'}
    # Under oxic conditions Fayalite (Fe2SiO4) is thermodynamically unstable —
    # it oxidises to Goethite/Magnetite during early alteration and is absent as
    # a reactive phase. Equilibrating it with O2(g) at low T drives PHREEQC to a
    # non-physical state (pH > 20) because log_k ~21 at 2°C gives enormous Fe2+
    # that the O2/H2O redox couple cannot stably handle.
    if fO2 > 0:
        _problematic.add('Fayalite')

    equilibrium_minerals = [m for m in composition.keys() if m not in _problematic]

    b = b_input if b_input is not None else np.array([])

    output = solve_solution(P, T, b, P_CO2=P_CO2, equilibriating_minerals=equilibrium_minerals, precipitating_minerals=precipitating_minerals, high_temperature=high_temperature, fO2=fO2)

    for i, element in enumerate(elements):
        output_key = 'Alk(eq/kgw)' if element == 'Alkalinity' else f'{element}(mol/kgw)'
        b_eq[i] = output[output_key][-1]

    pH = float(output['pH'][-1])

    return b_eq, pH

def aqueous_flux(P: float, T: float, P_CO2: float, J: float, A_reactive: float, composition: dict[str, float]) -> npt.NDArray[np.float64]:

    b_eq, pH = get_b_eq(P, T, P_CO2, composition)
    k = get_k(P, T, pH, composition)
    
    return A_reactive / ((1 / k) + (A_reactive / (J * b_eq)))

def get_precipitation_flux(P: float, T: float, b: npt.NDArray[np.float64], precipitating_minerals: list[str], equilibrium_minerals: list[str]=[], fO2: float=0) -> npt.NDArray[np.float64]:

    output = solve_solution(P, T, b, precipitating_minerals=precipitating_minerals, equilibriating_minerals=equilibrium_minerals, fO2=fO2)

    aqueous_fluxes = np.zeros(elements.shape)

    for i, element in enumerate(elements):
        output_key = 'Alk(eq/kgw)' if element == 'Alkalinity' else f'{element}(mol/kgw)'
        aqueous_fluxes[i] = (float(output[output_key][-1]) - float(output[output_key][0]))

    return aqueous_fluxes

def reactive_area(T: float, pH: float, rate: float, alpha: float, clog: bool=True, cover: bool=True, sedimentation_rate: float | None = None) -> float:

    T_ref = 280 # Pore space temperature
    T_c = 7 # Activation energy 92 kJ
    pH_ref = 8.1
    t_clog_ref = 20e6 * YR # Coogan & Gillis (2018), Fig 6
    beta = 0 # no pH dependence
    h_cover = 100 # m
    S_ref = 5 / (1e6 * YR) # 5 m / Myr
    S_min = 0.3 / (1e6 * YR) # 0.3 m / Myr background (dust)

    t_clog = t_clog_ref * np.exp(- (T - T_ref) / T_c) * (pH / pH_ref) ** beta
    S = max(sedimentation_rate, S_min) if sedimentation_rate is not None else S_ref
    t_cover = h_cover / S

    if clog and cover:
        t_reduce = (1 / t_clog + 1 / t_cover) ** -1
    elif clog:
        t_reduce = t_clog
    elif cover:
        t_reduce = t_cover
    else:
        return alpha

    return alpha * rate * t_reduce * (1 - np.exp(- 1 / (rate * t_reduce)))

rate_ref = 1 / (50e6 * YR)
J_ref = 1.4e15 / YR # kg / yr
A_seafloor = 0.7 * 4 * np.pi * R_EARTH ** 2
J_ref_normalised = J_ref / A_seafloor

def get_weathering_flux(P: float, T: float, P_CO2: float, b_input: npt.NDArray[np.float64], alpha: float | None=None, rate: float | None=None, J: float | None=None, cover: bool=True, clog: bool=False, high_temperature: bool=False, precipitating_minerals: list[str] | None=None, fO2: float=0, sedimentation_rate: float | None=None) -> tuple[npt.NDArray[np.float64], dict[str, float]]:

    molar_mass = 0.216 # kg / mol

    if alpha is None:
        alpha = ALPHA_REF

    if rate is None:
        rate = rate_ref

    if J is None:
        J = J_ref_normalised * (rate / rate_ref) # hydrothermal flux proportional to crust production rate

    if len(b_input) == 0:
        b_input = np.zeros(elements.shape)

    if high_temperature:
        composition = hydrothermal_composition
        default_precipitating = []          # Kaolinite not in hydrothermal.dat
    else:
        composition = basalt_composition
        default_precipitating = ['Kaolinite', 'Goethite']

    if precipitating_minerals is None:
        precipitating_minerals = default_precipitating

    # Pass b_input to get_b_eq so PHREEQC equilibrates starting from the input
    # ocean composition, not dilute water.  b_eq - b_input then gives the net
    # compositional change due to rock-water interaction.
    b_eq_primary, pH = get_b_eq(P, T, P_CO2, composition, b_input=b_input,
                                 precipitating_minerals=precipitating_minerals,
                                 high_temperature=high_temperature, fO2=fO2)
    k_primary = get_k(P, T, pH, composition)

    A_reactive = reactive_area(T, pH, rate, alpha, clog, cover, sedimentation_rate=sedimentation_rate)

    # Flux formula: F = A*(b_eq - b_in) / (b_eq/k + A/J)
    # Correct first-order kinetics box-model generalisation for non-zero b_input.
    # Equivalent to A/(1/k + A/(J*b_eq)) when b_input=0 (existing behaviour unchanged),
    # but the denominator uses b_eq (not delta_b), so it stays positive and gives the
    # correct sign when b_eq[i] < b_input[i] (e.g. Mg removal by Chrysotile/Talc).
    safe_k = np.where(k_primary != 0, k_primary, np.inf)
    F_primary = np.where(k_primary != 0,
                         A_reactive * (b_eq_primary - b_input) / (b_eq_primary / safe_k + A_reactive / J),
                         0.0)

    Da_primary = (k_primary * A_reactive * molar_mass) / J
    b_pore = b_input + (b_eq_primary - b_input) * (1 - np.exp(-Da_primary))

    # C has zero stoichiometry in all basalt minerals, so Da[C]=0 and b_pore[C]=0.
    # Use the equilibrium C concentration so PHREEQC has a consistent carbonate system.
    C_idx = int(np.where(elements == 'C')[0][0])
    b_pore[C_idx] = b_eq_primary[C_idx]

    flux = F_primary

    d_b_carb = get_precipitation_flux(P, T, b_pore, carbonate_minerals, [])
    flux += J * d_b_carb

    supply_efficiency = 1 - np.exp(-Da_primary[0])

    weathering_diagnostics = {
        'Da': Da_primary[0],
        'A_reactive': A_reactive,
        'supply_efficiency': supply_efficiency
    }

    return flux, weathering_diagnostics

# ---------------------------------------------------------------------------
# Continental silicate weathering (Walker-Hays-Kasting parameterization)
# ---------------------------------------------------------------------------

# Reference alkalinity flux per unit land area at modern Earth conditions.
# Calibrated so that 0.3 land fraction gives ~8 Tmol eq/yr total, which
# balances modern volcanic outgassing after accounting for seafloor weathering.
_CONT_LAND_AREA_EARTH = 0.3 * 4 * np.pi * R_EARTH ** 2   # m²
EARTH_CONTINENTAL_WEATHERING_REF = (8e12 / YR) / _CONT_LAND_AREA_EARTH  # mol_eq / m² / s

def get_continental_weathering_flux(
    T: float,
    P_CO2: float,
    F_alk_ref: float = EARTH_CONTINENTAL_WEATHERING_REF,
    T_ref: float = 288.0,
    P_CO2_ref: float = EARTH_ATM * 280e-6,
    beta: float = 0.3,
    T_e: float = 17.0,
) -> npt.NDArray[np.float64]:
    """Walker-Hays-Kasting continental silicate weathering parameterization.

    Returns flux per unit land area [mol/m²/s] using CaSiO3 stoichiometry:
        CaSiO3 + 2CO2 + H2O -> Ca2+ + 2HCO3-

    Parameters
    ----------
    T        : surface temperature [K]
    P_CO2    : atmospheric CO2 partial pressure [Pa]
    F_alk_ref: reference alkalinity flux at (T_ref, P_CO2_ref) [mol_eq/m²/s]
    beta     : CO2 sensitivity exponent (default 0.3)
    T_e      : temperature e-folding scale [K] (default 17 K ~ 70 kJ/mol)
    """
    f = (P_CO2 / P_CO2_ref) ** beta * np.exp((T - T_ref) / T_e)
    F_alk = F_alk_ref * f   # mol_eq / m² / s

    flux = np.zeros(len(elements))
    _alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
    _C_idx   = int(np.where(elements == 'C')[0][0])
    _Ca_idx  = int(np.where(elements == 'Ca')[0][0])

    # Per mol CaSiO3 dissolved: 1 mol Ca2+, 2 mol HCO3-
    # => 2 eq Alk, 2 mol C, 1 mol Ca per mol CaSiO3
    flux[_alk_idx] = F_alk        # eq / m² / s
    flux[_C_idx]   = F_alk        # mol C / m² / s  (1 HCO3- = 1 C = 1 eq Alk)
    flux[_Ca_idx]  = F_alk / 2    # mol Ca / m² / s (2 eq per mol Ca2+)

    return flux


def get_P_CO2(P: float, T: float, b: npt.NDArray[np.float64]):

    output = solve_solution(P, T, b, precipitating_minerals=['CO2(g)'])

    si_CO2 = float(output['si_CO2(g)'][-1])
    P_CO2 = EARTH_ATM * 10 ** si_CO2

    return P_CO2

from scipy.optimize import least_squares

# print("Calibrating weathering model to modern Earth baseline (1 Tmol/yr Alkalinity)...")

_T_ref_calib = 280
_P_ref_calib = 1000 * 10 * 3000
_P_CO2_ref_calib = EARTH_ATM * 280e-6
_seafloor_flux = 1e12 / YR
_seafloor_flux_normalised = _seafloor_flux / A_seafloor

# Residual function to find the alpha that produces the correct Alkalinity flux (Index 0)
def _calibration_residual(a_array):
    a_val = a_array[0]
    flux, _ = get_weathering_flux(
        _P_ref_calib, _T_ref_calib, _P_CO2_ref_calib, 
        np.array([]), alpha=a_val, rate=rate_ref, J=J_ref_normalised
    )
    alkalinity_flux = flux[0] 
    return (alkalinity_flux - _seafloor_flux_normalised) / _seafloor_flux_normalised

# Run the optimization
_root = least_squares(_calibration_residual, [100.0])
ALPHA_REF = float(_root.x[0])
# print(f"Calibration complete: ALPHA_REF = {ALPHA_REF:.2e}")

