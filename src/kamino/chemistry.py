import os

from kamino.constants import ABSOLUTE_ZERO, EARTH_ATM
from kamino.mineral_rates import K_FUNCTIONS
from kamino.mineral_info import (
    MINERAL_MOLAR_MASS,
    carbonate_minerals,
    clay_minerals,
    evaporite_minerals,
    hydrothermal_mineral_string,
    ht_secondary_minerals,
    primary_minerals,
    reverse_weathering_minerals,
    silica_minerals,
)

# The SIT/ThermoChimie base names gases conventionally (O2, H2, N2, CH4); the former Pitzer
# runtime database used ThermoChimie's alternative spellings (Oxg, Hdg, Ntg, Mtg). This is the
# same class of trap as the silica naming in section 11, and it is silent: an absent gas phase
# makes PHREEQC reject the whole SELECTED_OUTPUT block rather than warn.
OXYGEN_GAS_PHASE = 'O2(g)'

_NON_CONVERGENCE_WARNINGS = (
    'Maximum iterations exceeded',
    'Numerical method failed',
    'has not converged',
)

DISSOLVE_ONLY_PRIMARY = False

from phreeqc import Phreeqc
import numpy as np
import numpy.typing as npt
import re
import pandas as pd

class ChemistryError(Exception):
    """PHREEQC failed to load a database or to converge on a solution."""

_data_dir = os.path.join(os.path.dirname(__file__), 'data')
# Runtime LT database. Switched from hybrid_ocean.dat (Pitzer) to lt_weathering_sit.dat (SIT)
# -- see section 25.12/25.14 of the development history. Two reasons: hybrid_ocean.dat is not
# reproducible (its generator, make_hybrid.py, is not in the repository), and pitzer.dat carries
# NO aluminium parameterisation at all, which matters because the feldspar/nepheline/kaolinite
# equilibria this model turns on are written against Al+3. Every ocean state the model reaches
# sits at I = 0.7-2.8 mol/kgw, inside SIT's validity.
# Runtime LT database, built by data/make_database.py and therefore reproducible from the
# databases bundled with the `phreeqc` package -- unlike the previous hybrid_ocean.dat, whose
# generator (make_hybrid.py) was never committed. See development history 25.12/25.14.
#
# The Pitzer build is the default: it is ~8x faster per solve than the SIT build, and the two
# agree on the alkalinity weathering flux to within 10-20% and on the feedback strength
# d ln(F)/dT to within ~15%, because the model sits in the kinetically-limited regime where
# F -> A_r*k and the equilibrium concentrations the activity model sets divide out.
#
# Override to compare, without editing code:
#   KAMINO_LT_DATABASE=src/kamino/data/lt_weathering_sit.dat python ...
database_path = os.environ.get(
    'KAMINO_LT_DATABASE', os.path.join(_data_dir, 'lt_weathering_pitzer.dat'))
_hydrothermal_path = os.path.join(_data_dir, 'hydrothermal.dat')

elements = np.array([
    'Alkalinity',
    'C',
    'Si',
    'Al',
    'Fe',
    'Ca',
    'Mg',
    'Na',
    'Cl',
    'S',    # sulfate — fixed background for charge balance; not a dynamic ODE variable
    #'K',
    #'N',
    #'F',
])

alk_idx = int(np.where(elements == 'Alkalinity')[0][0])
ca_idx  = int(np.where(elements == 'Ca')[0][0])
mg_idx  = int(np.where(elements == 'Mg')[0][0])
c_idx   = int(np.where(elements == 'C')[0][0])
si_idx  = int(np.where(elements == 'Si')[0][0])
al_idx  = int(np.where(elements == 'Al')[0][0])
fe_idx  = int(np.where(elements == 'Fe')[0][0])
na_idx  = int(np.where(elements == 'Na')[0][0])
cl_idx  = int(np.where(elements == 'Cl')[0][0])
so4_idx = int(np.where(elements == 'S')[0][0])

# Conservative-ion charge per element (signed), for deriving alkalinity as the charge
# balance of the tracked ions: Alk = sum(ION_CHARGE * concentration). Fe is +2 (Fe2+,
# soluble/conservative under the anoxic conditions where it is mobile; oxic Fe
# precipitates as Goethite and is removed with its charge). Al nominal +3 (trace).
# C/Si are 0 (the carbonate system C is solved via alkalinity+DIC, not counted here).
ION_CHARGE = np.zeros(elements.shape)
ION_CHARGE[ca_idx] = 2.0
ION_CHARGE[mg_idx] = 2.0
ION_CHARGE[na_idx] = 1.0
ION_CHARGE[fe_idx] = 2.0
ION_CHARGE[al_idx] = 3.0
ION_CHARGE[cl_idx] = -1.0
ION_CHARGE[so4_idx] = -2.0

# PHREEQC reports Alkalinity under its own key, not as a -totals element.
element_string = ' '.join(elements[elements != 'Alkalinity'])

# Parses a PHREEQC database file to extract the stoichiometry and alkalinity generation of each mineral upon dissolution.

# Mapping PHREEQC aqueous species back to elemental buckets
species_to_element = {
    'Ca+2': 'Ca',
    'Mg+2': 'Mg',
    'Na+': 'Na',
    'Cl-': 'Cl',
    # 'K+': 'K',
    'Fe+2': 'Fe',
    'Fe+3': 'Fe',
    'Al+3': 'Al',
    'SiO2': 'Si',
    'H4SiO4': 'Si',
    'H4(SiO4)': 'Si',   # ThermoChimie/SIT spelling of the same species
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

# Species that appear in PHASES reactions but are deliberately NOT tracked. Anything
# appearing in a reaction that is in neither this set nor the two maps above raises,
# rather than silently contributing zero to the element and alkalinity budgets.
IGNORED_SPECIES = {
    'H2O', 'e-',                                  # solvent / electron
    'SO4-2',                                      # sulfur chemistry intentionally decoupled:
                                                  # S is pinned as a charge-balance background
                                                  # (Planet.dY_dt sets F_net[so4_idx] = 0)
    'K+',                                         # K not in `elements`
    'B(OH)3', 'Ba+2', 'Sr+2',                     # B/Ba/Sr not in `elements`
    'CO2', 'Oxg', 'Hdg', 'Ntg', 'Mtg', 'HSg-',    # gas phases, not aqueous buckets
    'Al(OH)4-',                                   # aluminate. Only in the LT Anorthite reaction,
                                                  # which the HT (Al+3-basis) entry overrides.
                                                  # That basis is the one consistent with Kaolinite,
                                                  # so Anorthite + Kaolinite nets to the textbook
                                                  # Alk = 2. Do not "fix" without checking that.
}

_COEFF_RE = re.compile(r'\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$')

# A coefficient fused to the species it multiplies, with no separating space: Thermoddem writes
# '2.000H+' and '1.000H4SiO4' where llnl writes '2 H+'. PHREEQC accepts both. Anchored at the
# start and requiring a following letter or '(' so it cannot bite a bare number or a species
# whose name legitimately begins with a digit.
_FUSED_RE = re.compile(r'^(\d+(?:\.\d+)?)(?=[A-Za-z(])(.+)$')

def _iter_terms(side: str):
    """Yield (coeff, species) for one side of a reaction.

    Handles '2 X', 'X', negative coefficients written as '- X' (which appear in llnl.dat-derived
    reactions, e.g. 'CaSiO3 + 2 H+ = Ca+2 - H2O + H4SiO4'), and the fused Thermoddem form
    '2.000H+' (see _FUSED_RE). Species names may themselves end in '+'/'-' (HCO3-, Ca+2); only a
    bare '+'/'-' token is an operator.
    """
    sign, coeff = 1.0, None
    for tok in side.split():
        if tok == '+':
            sign, coeff = 1.0, None
        elif tok == '-':
            sign, coeff = -1.0, None
        elif _COEFF_RE.match(tok):
            coeff = float(tok)
        else:
            fused = _FUSED_RE.match(tok)
            if fused:
                coeff, tok = float(fused.group(1)), fused.group(2)
            yield sign * (1.0 if coeff is None else coeff), tok
            sign, coeff = 1.0, None

def _db_phase_names(filepath: str) -> set[str]:
    """Every phase name in a database's PHASES section, gases included."""
    names, in_phases = set(), False
    with open(filepath) as f:
        for raw in f:
            line = raw.split('#')[0].strip()
            if not line:
                continue
            if line == 'PHASES':
                in_phases = True
            elif line in ('RATES', 'SOLUTION_SPECIES', 'END', 'PITZER', 'SIT'):
                in_phases = False
            elif in_phases and '=' not in line and not line.startswith('-') and not raw[0].isspace():
                names.add(line.split()[0])
    return names


def parse_stoichiometry(filepath: str, keep: set[str] | None = None) -> dict[str, npt.NDArray[np.float64]]:
    """Parse mineral reaction stoichiometry from a PHREEQC database PHASES section.

    Works with both reaction formats found in PHREEQC databases:
      - Indented reactions (ocean_chem.dat / SUPCRT92 style)
      - Unindented reactions (hydrothermal.dat / SupPHREEQC bl-1kb.dat style)

    `keep`, if given, restricts parsing to those phase names. The SIT database is the full
    ThermoChimie set (~1772 phases spanning actinides, ore metals and organics); without a filter
    the unmapped-species check below would demand a species map for chemistry this model has no
    business knowing about, and PHREEQC would be asked for 1772 saturation indices per solve.

    Raises
    ------
    ChemistryError
        If a reaction contains an aqueous species that is neither mapped to an
        element (species_to_element / alk_contributions) nor explicitly listed in
        IGNORED_SPECIES. Silently scoring an unknown species as zero is how
        element and alkalinity budgets go quietly wrong.
    """
    unmapped: set[str] = set()

    def parse_side(side_str: str, stoich: dict, multiplier: float, is_lhs: bool = False) -> None:
        terms = list(_iter_terms(side_str))
        if is_lhs and terms:
            terms = terms[1:]  # first LHS term is the mineral formula itself
        for coeff, species in terms:
            known = False
            if species in species_to_element:
                stoich[species_to_element[species]] += coeff * multiplier
                known = True
            if species in alk_contributions:
                stoich['Alkalinity'] += coeff * alk_contributions[species] * multiplier
                known = True
            if not known and species not in IGNORED_SPECIES:
                unmapped.add(species)

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
            elif line in ('RATES', 'SOLUTION_SPECIES', 'END', 'PITZER'):
                in_phases = False
                continue
            if not in_phases:
                continue

            if '=' not in line and not line.startswith('-') and not raw_line[0].isspace():
                current_mineral = line.split()[0]
            elif current_mineral and '=' in line and 'log_k' not in line and 'delta_H' not in line:
                if keep is not None and current_mineral not in keep:
                    current_mineral = None      # not a phase this model uses; skip its reaction
                    continue
                stoich: dict[str, float] = {str(el): 0.0 for el in elements}
                lhs_str, rhs_str = line.split('=', 1)
                parse_side(lhs_str, stoich, multiplier=-1.0, is_lhs=True)
                parse_side(rhs_str, stoich, multiplier=1.0)
                result[current_mineral] = np.array([round(v, 6) for v in stoich.values()])
                current_mineral = None

    if unmapped:
        raise ChemistryError(
            f"{filepath}: PHASES reactions contain species that are neither mapped to "
            f"an element nor listed in IGNORED_SPECIES: {sorted(unmapped)}. "
            f"Add them to species_to_element/alk_contributions to track them, or to "
            f"IGNORED_SPECIES to deliberately exclude them."
        )

    return result

# Phases the model can use as a crust component or precipitating sink. Everything else in the
# SIT database is irrelevant here and is not parsed -- see parse_stoichiometry's `keep`.
_MODEL_MINERALS = (set(MINERAL_MOLAR_MASS) | set(carbonate_minerals) | set(clay_minerals)
                   | set(silica_minerals) | set(reverse_weathering_minerals)
                   | set(evaporite_minerals) | set(ht_secondary_minerals))

_stoich_LT = parse_stoichiometry(database_path, keep=_MODEL_MINERALS)
_stoich_HT = parse_stoichiometry(_hydrothermal_path)

# Some minerals appear in BOTH databases with reactions written on different bases
# (different formula units, or different aqueous Al species). Merging the two dicts
# would let load order silently pick a winner, so every mineral whose two entries
# disagree must declare which database it takes its stoichiometry from.
#
# The stoichiometry vector is combined with a rate constant from mineral_rates.py
# (see get_k), so it must be on the same formula-unit basis as that rate constant.
STOICHIOMETRY_SOURCE: dict[str, str] = {
    # HT writes Mg2Si2O6 (double formula unit); LT writes MgSiO3. K_FUNCTIONS['Enstatite']
    # is a rate per mole of MgSiO3, so LT is the consistent basis. Enstatite is only used
    # in the LT path anyway — get_weathering_flux swaps it for Diopside at HT.
    'Enstatite': 'LT',
    # LT writes Anorthite against aluminate (Al(OH)4-), HT against Al+3. The Al+3 basis is
    # the one Kaolinite is written on, so Anorthite dissolution followed by Kaolinite
    # precipitation nets to the textbook Alk = 2. Keep HT.
    'Anorthite': 'HT',
    # Gas phase, never used as a crust or precipitating mineral, so the vector is unused.
    # Declared only to satisfy the collision check below.
    'CO2(g)': 'HT',
}

stoichiometry: dict[str, npt.NDArray[np.float64]] = {**_stoich_LT, **_stoich_HT}

# Only phases the model actually uses need a declared source; the SIT database shares hundreds
# of irrelevant names with hydrothermal.dat and flagging those would be noise.
_undeclared = sorted(
    m for m in set(_stoich_LT) & set(_stoich_HT) & _MODEL_MINERALS
    if not np.array_equal(_stoich_LT[m], _stoich_HT[m]) and m not in STOICHIOMETRY_SOURCE
)
if _undeclared:
    raise ChemistryError(
        f"Minerals defined in both databases with DIFFERENT stoichiometry, but not "
        f"declared in STOICHIOMETRY_SOURCE: {_undeclared}. Merge order would silently "
        f"pick a winner. Declare which database each should use."
    )

for _mineral, _src in STOICHIOMETRY_SOURCE.items():
    stoichiometry[_mineral] = (_stoich_LT if _src == 'LT' else _stoich_HT)[_mineral]

# Minerals PHREEQC can report saturation indices for.
#
# The SIT database carries ~1772 phases, almost all irrelevant here (actinides, cements,
# zeolites, ore minerals). Requesting saturation indices for all of them would make PHREEQC
# compute 1772 SIs on every single solve. Restrict the request to the phases the model can
# actually use as a crust component or a precipitating sink -- section 19 flagged this list as
# unfiltered back when it was 84 phases; at 1772 it stops being cosmetic.
# Gas phases whose saturation indices the code reads back (get_ocean_state, and
# get_gas_partial_pressure for any gas the caller asks for). These are NOT minerals and are not
# in _MODEL_MINERALS, but they must still be requested or the SI lookup KeyErrors.
# Gas saturation indices to request. CO2 is the only one read back in the default path
# (get_ocean_state, get_P_CO2); the others are requested when present so the fO2-buffered path
# and caller-supplied gas lists resolve. Which names exist is database-dependent -- SIT uses the
# conventional O2/H2/N2/CH4, the Pitzer database uses ThermoChimie's Oxg/Hdg/Ntg/Mtg -- so the
# list is intersected with what the loaded database actually defines. Requesting an absent phase
# makes PHREEQC reject the whole SELECTED_OUTPUT block.
_SI_GAS_CANDIDATES = ['CO2(g)', 'H2O(g)', 'O2(g)', 'Oxg(g)']
_REQUIRED_GASES = ['CO2(g)']

minerals = [m for m in _stoich_LT if m in _MODEL_MINERALS]
# Hedenbergite and Ferrosilite were exempted here while they sat behind an off-by-default flag.
# They are now always emitted by the CIPW norm, so a database lacking them would silently drop the
# crust's iron out of the reactive assemblage instead of failing -- they must be required.
_missing = sorted(_MODEL_MINERALS - set(_stoich_LT) - set(_stoich_HT))
if _missing:
    raise ChemistryError(
        f"Minerals the model references are absent from the LT database {database_path}: "
        f"{_missing}. Add them to make_database.ADDED_PHASES and regenerate.")
_db_phases = _db_phase_names(database_path)
_SI_GASES = [g for g in _SI_GAS_CANDIDATES if g in _db_phases]
_missing_gases = [g for g in _REQUIRED_GASES if g not in _db_phases]
if _missing_gases:
    raise ChemistryError(
        f"Gas phases whose saturation indices the model reads back are absent from "
        f"{database_path}: {_missing_gases}.")
# OXYGEN_GAS_PHASE is only used when fO2 > 0; point it at whichever name this database uses.
OXYGEN_GAS_PHASE = next((g for g in ('O2(g)', 'Oxg(g)') if g in _db_phases), OXYGEN_GAS_PHASE)
available_mineral_string = ' '.join(minerals + _SI_GASES)

def _load_database(path: str) -> Phreeqc:
    """Load a PHREEQC database, raising if it fails rather than continuing blind."""
    p = Phreeqc()
    if p.LoadDatabase(path) == 1:
        raise ChemistryError(f"PHREEQC failed to load database {path}:\n{p.GetErrorString()}")
    return p

p_LT = _load_database(database_path)
p_HT = _load_database(_hydrothermal_path)

def _knobs_block() -> list[str]:
    """PHREEQC KNOBS block: Newton solver settings.

    On stock defaults, carbon-rich ocean compositions (C > 0.1 mol/kgw, i.e. pCO2 above
    a few bar) routinely exhaust PHREEQC's iteration limit and hand back a non-converged
    iterate, which solve_solution then has to reject as a ChemistryError. That is a
    solver-tuning failure, not a thermodynamic one: on 200 inputs captured from a failing
    high-pCO2 trajectory, the defaults converged 144/200 while the settings below
    converged 200/200. Over a full run at S=0.7 they cut ChemistryErrors 2785 -> 393.

    -step_size / -pe_step_size are the load-bearing pair; raising -iterations alone only
    reaches ~90%. Smaller Newton steps are what the stiff, high-ionic-strength carbonate
    systems need.
    """
    return [
        'KNOBS',
        '    -iterations 1000',
        '    -diagonal_scale true',
        '    -step_size 5',
        '    -pe_step_size 5',
        '',
    ]

def _solution_block(P: float, T: float, b: npt.NDArray[np.float64],
                    pH: float | None, trace_approximation: bool) -> list[str]:
    """PHREEQC SOLUTION block: bulk composition at (P, T)."""
    lines = [
        'SOLUTION 1',
        f'    pressure  {P / EARTH_ATM:.4f}',
        f'    temp      {max(T + ABSOLUTE_ZERO, 0.01):.4f}',  # LLNL database valid from 0.01°C
        '    units     mol/kgw',
    ]
    if pH is not None:
        lines.append(f'    pH     {pH:5f}')
    for element, x in zip(elements, b):
        molality = 1e-9 if x < 1e-9 and trace_approximation else x
        lines.append(f'    {element}    {molality:.15e}')
    lines.append('')
    return lines

def _equilibrium_block(P_CO2: float | None, fO2: float,
                       equilibriating_minerals: list[str],
                       equilibriating_amounts: dict[str, float] | None,
                       precipitating_minerals: list[str],
                       precipitation_SI: float,
                       dissolve_only_primary: bool = False) -> list[str]:
    """PHREEQC EQUILIBRIUM_PHASES block.

    Equilibriating phases are given a finite amount (they may dissolve); precipitating
    phases are given amount 0 (they may only precipitate out of solution).

    When ``dissolve_only_primary`` is enabled, primary igneous minerals
    (mineral_info.primary_minerals) also get PHREEQC's `dissolve_only` modifier, making
    them a source only.
    """
    if not (equilibriating_minerals or precipitating_minerals or P_CO2 is not None or fO2 > 0):
        return []

    lines = ['EQUILIBRIUM_PHASES 1']

    if P_CO2 is not None:
        P_CO2 = np.maximum(P_CO2, 1e-2)
        lines.append(f'    CO2(g)  {np.log10(P_CO2 / EARTH_ATM):.4f}  {1e6}') # type: ignore

    if fO2 > 0:
        lines.append(f'    {OXYGEN_GAS_PHASE}  {np.log10(fO2 / EARTH_ATM):.4f}  {1e6}')

    for phase in equilibriating_minerals:
        amount = equilibriating_amounts.get(phase, 100.0) if equilibriating_amounts else 100.0
        modifier = '  dissolve_only' if (dissolve_only_primary and phase in primary_minerals) else ''
        lines.append(f'    {phase}  {0.0}  {amount:.6f}{modifier}')

    lines.append('')

    for phase in precipitating_minerals:
        lines.append(f'    {phase}  {precipitation_SI:.4f}  {0.0}')

    lines.append('')
    return lines

def _output_block(high_temperature: bool, reported_phases: list[str]) -> list[str]:
    """PHREEQC SELECTED_OUTPUT block: what to read back out."""
    lines = [
        'SELECTED_OUTPUT',
        '    -pH',
        '    -pe',
        '    -alkalinity',
        f'    -totals {element_string}',
        f'    -saturation_indices {hydrothermal_mineral_string if high_temperature else available_mineral_string}',
    ]
    if reported_phases:
        lines.append('    -equilibrium_phases ' + ' '.join(reported_phases))
    return lines

def solve_solution(P: float, T: float, b: npt.NDArray[np.float64], pH: float | None=None, P_CO2: float | None=None, precipitating_minerals: list[str]=[], equilibriating_minerals: list[str]=[], equilibriating_amounts: dict[str, float] | None=None, high_temperature: bool=False, fO2: float=0, trace_approximation: bool=True, precipitation_SI: float=0, verbose: bool=False, dissolve_only_primary: bool | None=None):

    # Per-call override of the module-level DISSOLVE_ONLY_PRIMARY toggle (None -> use global).
    do_primary = DISSOLVE_ONLY_PRIMARY if dissolve_only_primary is None else dissolve_only_primary

    input_lines = (
        _knobs_block()
        + _solution_block(P, T, b, pH, trace_approximation)
        + _equilibrium_block(P_CO2, fO2, equilibriating_minerals, equilibriating_amounts,
                             precipitating_minerals, precipitation_SI, do_primary)
        + _output_block(high_temperature, equilibriating_minerals + precipitating_minerals)
        + ['\n']
    )
    input_string = '\n'.join(input_lines)

    if verbose:
        print(input_string)

    p = p_HT if high_temperature else p_LT

    def _attempt() -> tuple[bool, str]:
        """Submit input_string once. Returns (converged, non_convergence_warnings)."""
        if p.RunString(input_string) == 1:
            raise ChemistryError(
                f"PHREEQC returned an error:\n{p.GetErrorString()}\n"
                f"--- input ---\n{input_string}"
            )
        warnings = p.GetWarningString() or ''
        return (not any(s in warnings for s in _NON_CONVERGENCE_WARNINGS)), warnings

    converged, warnings = _attempt()
    if not converged:
        converged, warnings = _attempt()   # one retry, from PHREEQC's post-failure state
    if not converged:
        raise ChemistryError(
            f"PHREEQC failed to converge (returned a non-converged iterate, twice):\n"
            f"{warnings.strip()}\n--- input ---\n{input_string}"
        )

    output_dict = p.GetSelectedOutput()
    if not output_dict:
        # Deliberately not a ChemistryError: Planet.dY_dt catches those and degrades
        # to outgassing-only. Empty output means the request itself was malformed,
        # which should fail loudly rather than silently alter the physics.
        raise ValueError(
            f"PHREEQC returned no selected output.\n--- input ---\n{input_string}"
        )

    if verbose:
        df = pd.DataFrame(output_dict)
        df = df.loc[:, (df != -999.999).any()]
        pd.set_option('display.float_format', '{:.2e}'.format)
        print(df.T)

    return output_dict

def get_k(P: float, T: float, pH: float, composition: dict[str, float]) -> npt.NDArray[np.float64]:
    
    k = np.zeros(elements.shape)

    for mineral, mineral_fraction in composition.items():
        k_eff = K_FUNCTIONS[mineral](T, pH)
        k += k_eff * stoichiometry[mineral] * mineral_fraction

    return k

def get_b_eq(P: float, T: float, P_CO2: float, composition: dict[str, float], b_input: npt.NDArray[np.float64] | None=None, precipitating_minerals: list[str]=[], high_temperature: bool=False, fO2: float=0, water_rock_ratio: float | None=None, dissolve_only: bool | None=None) -> tuple[npt.NDArray[np.float64], float]:

    b_eq = np.zeros(elements.shape)

    # Carbonate minerals are excluded: CO2(g) already constrains carbonate chemistry
    _excluded = set(carbonate_minerals)

    # These minerals cause problems with oxic conditions
    if fO2 > 0:
        _excluded.add('Fayalite')

    equilibrium_minerals = [m for m in composition.keys() if m not in _excluded]

    b = b_input if b_input is not None else np.array([])

    # If a water-to-rock mass ratio is given, compute physically constrained per-mineral mole amounts: amount_i = (1000 g/kgw / w_r) * fraction_i / M_i
    # This prevents PHREEQC from dissolving unrealistically large amounts of crust.
    eq_amounts: dict[str, float] | None = None
    if water_rock_ratio is not None:
        rock_g_per_kgw = 1000.0 / water_rock_ratio
        eq_amounts = {
            m: rock_g_per_kgw * composition[m] / MINERAL_MOLAR_MASS.get(m, 150.0)
            for m in equilibrium_minerals
        }

    # At HT conditions, the greenschist secondary assemblage precipitates: Quartz buffers Si, Clinochlore (Mg-chlorite) is the Mg sink.
    if high_temperature:
        precipitating_minerals = list(precipitating_minerals) + ht_secondary_minerals

    output = solve_solution(P, T, b, P_CO2=P_CO2, equilibriating_minerals=equilibrium_minerals, equilibriating_amounts=eq_amounts, precipitating_minerals=precipitating_minerals, high_temperature=high_temperature, fO2=fO2, dissolve_only_primary=dissolve_only)

    for i, element in enumerate(elements):
        output_key = 'Alk(eq/kgw)' if element == 'Alkalinity' else f'{element}(mol/kgw)'
        b_eq[i] = output[output_key][-1]

    pH = float(output['pH'][-1])

    return b_eq, pH

def get_gas_partial_pressure(P: float, T: float, b: npt.NDArray[np.float64], gases: list[str], pH: float | None=None) -> list[float]:

    output = solve_solution(P, T, b, pH=pH)

    P_gases = []

    for gas in gases:
        si = float(output[f'si_{gas}(g)'][-1])
        P_gases.append(EARTH_ATM * 10 ** si)

    return P_gases

def get_ocean_state(P: float, T: float, b: npt.NDArray[np.float64]) -> tuple[float, float]:
    output = solve_solution(P, T, b)
    P_CO2 = EARTH_ATM * 10 ** float(output['si_CO2(g)'][-1])
    pH = float(output['pH'][-1])
    return P_CO2, pH

def get_P_CO2(P: float, T: float, b: npt.NDArray[np.float64]) -> float:
    return get_gas_partial_pressure(P, T, b, ['CO2'])[0]

def get_pH(P: float, T: float, b: npt.NDArray[np.float64]) -> float:
    output = solve_solution(P, T, b)
    pH = float(output['pH'][-1])
    return pH