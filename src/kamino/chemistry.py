import os

from kamino.constants import *
from kamino.mineral_rates import K_FUNCTIONS
from kamino.mineral_info import *

from phreeqc import Phreeqc
import numpy as np
import numpy.typing as npt
import re
import pandas as pd

_data_dir = os.path.join(os.path.dirname(__file__), 'data')
database_path = os.path.join(_data_dir, 'ocean_chem.dat')
_hydrothermal_path = os.path.join(_data_dir, 'hydrothermal.dat')

elements = np.array([
    'Alkalinity',
    'C',
    'Si',
    'Al',
    'Fe',
    'Ca',
    'Mg',
    #'Na',
    #'Cl',
    #'K',
    #'S',
    #'N',
    #'F',
])

alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
ca_idx  = int(np.where(elements == 'Ca')[0][0])
mg_idx  = int(np.where(elements == 'Mg')[0][0])
c_idx   = int(np.where(elements == 'C')[0][0])
si_idx = int(np.where(elements == 'Si')[0][0])

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
    # 'Na+': 'Na',
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
stoichiometry.update(parse_stoichiometry(_hydrothermal_path))

p_LT = Phreeqc()
if p_LT.LoadDatabase(database_path) == 1:
    print(p_LT.GetErrorString())

p_HT = Phreeqc()
p_HT.LoadDatabase(_hydrothermal_path)

def solve_solution(P: float, T: float, b: npt.NDArray[np.float64], pH: float | None=None, P_CO2: float | None=None, precipitating_minerals: list[str]=[], equilibriating_minerals: list[str]=[], high_temperature: bool=False, fO2: float=0, trace_approximation: bool=True, precipitation_SI: float=0, verbose: bool=False):

    solution_lines: list[str] = [
        #'KNOBS',
        #'    -iterations 2000',
        #'    -convergence_tolerance 1e-9',
        #'    -step_size 0.5',   # half-step Newton: more stable near charge-balance limit
        #'',
        'SOLUTION 1',
        f'    pressure  {P / EARTH_ATM:.4f}',
        f'    temp      {max(T + ABSOLUTE_ZERO, 0.01):.4f}',  # LLNL database valid from 0.01°C
        '    units     mol/kgw'
    ]

    if pH is not None:
        solution_lines.append(f'    pH     {pH:5f}')

    for element, x in zip(elements, b):
        molality = 1e-9 if x < 1e-9 and trace_approximation else x
        solution_lines.append(f'    {element}    {molality:.15e}')

    solution_lines.append('')

    equilibrium_lines: list[str] = []

    if equilibriating_minerals or precipitating_minerals or P_CO2 is not None or fO2 > 0:

        equilibrium_lines.append('EQUILIBRIUM_PHASES 1')

        if P_CO2 is not None:
            P_CO2 = np.maximum(P_CO2, 1e-2)
            equilibrium_lines.append(f'    CO2(g)  {np.log10(P_CO2 / EARTH_ATM):.4f}  {1e6}') # type: ignore

        if fO2 > 0:
            equilibrium_lines.append(f'    O2(g)   {np.log10(fO2 / EARTH_ATM):.4f}  {1e6}')

        for phase in equilibriating_minerals:
            equilibrium_lines.append(f'    {phase}  {0.0}  {100.0}')

        equilibrium_lines.append('')

        for phase in precipitating_minerals:
            equilibrium_lines.append(f'    {phase}  {precipitation_SI:.4f}  {0.0}')

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

    if p.RunString(input_string) == 1:
        # print()
        # print(input_string)
        # print(p.GetErrorString())
        raise ChemistryError

    output_dict = p.GetSelectedOutput()
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

    # Carbonate minerals are excluded: CO2(g) already constrains carbonate chemistry
    _excluded = set(carbonate_minerals)

    # These silicate minerals cause problems for the equilibrium calculations at low T
    if not high_temperature:
        _excluded |= {'Anorthite', 'Forsterite'}
    else:
        _excluded.add('Anorthite')

    # These minerals cause problems with oxic conditions
    if fO2 > 0:
        _excluded.add('Fayalite')

    equilibrium_minerals = [m for m in composition.keys() if m not in _excluded]

    b = b_input if b_input is not None else np.array([])

    output = solve_solution(P, T, b, P_CO2=P_CO2, equilibriating_minerals=equilibrium_minerals, precipitating_minerals=precipitating_minerals, high_temperature=high_temperature, fO2=fO2)

    for i, element in enumerate(elements):
        output_key = 'Alk(eq/kgw)' if element == 'Alkalinity' else f'{element}(mol/kgw)'
        b_eq[i] = output[output_key][-1]

    pH = float(output['pH'][-1])

    return b_eq, pH

def get_gas_partial_pressure(P: float, T: float, b: npt.NDArray[np.float64], gases: list[str]) -> list[float]:

    output = solve_solution(P, T, b)

    P_gases = []

    for gas in gases:
        si = float(output[f'si_{gas}(g)'][-1])
        P_gases.append(EARTH_ATM * 10 ** si)

    return P_gases

def get_P_CO2(P: float, T: float, b: npt.NDArray[np.float64]) -> float:
    return get_gas_partial_pressure(P, T, b, ['CO2'])[0]

def get_pH(P: float, T: float, b: npt.NDArray[np.float64]) -> float:
    output = solve_solution(P, T, b)
    pH = float(output['pH'][-1])
    return pH

class ChemistryError(Exception):
    pass