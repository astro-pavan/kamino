import subprocess
import numpy as np
import pandas as pd
from typing import Union
import tempfile
import sys
from pathlib import Path
import importlib.resources

from kamino.constants import *

try:
    from . import build_phreeqc
except ImportError:
    import build_phreeqc

debug_mode = False

available_elements = ['Si', 'Al', 'Fe', 'Ca', 'Mg', 'Na', 'K', 'C', 'S', 'N', 'F', 'Cl']
available_element_string = 'Si Al Fe Ca Mg Na K C S N F Cl'

available_minerals = ['Albite', 'Anorthite', 'K-Feldspar', 'Enstatite', 'Diopside', 'Forsterite', 'Fayalite', 'Quartz', 'Magnetite', 'Calcite', 'Aragonite', 'Kaolinite', 'Goethite']

insoluble_minerals = ['Magnetite', 'Goethite']

available_mineral_string = ''
for mineral in available_minerals:
    available_mineral_string += mineral + ' '

def get_phreeqc_resources():
    """
    Locates the PHREEQC executable and database.
    Returns: (path_to_executable, path_to_database)
    """
    # 1. Locate the package directory using modern importlib
    # This points to src/kamino/ocean_chemistry
    try:
        pkg_path = importlib.resources.files("kamino.ocean_chemistry")
    except (ImportError, AttributeError):
        # Fallback for older python or non-package execution
        pkg_path = Path(__file__).parent

    # 2. Define expected paths based on your build script structure
    # Expected: ocean_chemistry/phreeqc_bin/bin/phreeqc
    bin_dir = pkg_path / "phreeqc_bin"
    
    # Executable logic (handling Windows/Linux)
    exe_name = "phreeqc.exe" if sys.platform == "win32" else "phreeqc"
    
    # Check common locations (bin/phreeqc or just phreeqc)
    executable = bin_dir / "bin" / exe_name

    if not executable.exists(): # type: ignore
        print(f"PHREEQC executable not found at {executable}")
        print("Attempting to compile PHREEQC automatically... (This may take a minute)")
        
        try:
            # Call the build function from your build_phreeqc.py script
            build_phreeqc.build()
            print("Compilation complete.")
        except Exception as e:
            raise RuntimeError(
                f"Automatic build failed. Please check you have a C++ compiler (g++ or clang) installed.\nError: {e}"
            )

        # Re-check after build
        if not executable.exists(): # type: ignore
             raise FileNotFoundError(f"Build appeared to succeed, but {executable} is still missing.")
        
    database = bin_dir / "basic_v2.dat"
    if not database.exists(): # type: ignore
        # Search recursively if not immediately found
        found = list(bin_dir.rglob("basic_v2.dat")) # type: ignore
        if found:
            database = found[0]
            
    return executable, database

executable_path, database_path = get_phreeqc_resources()

def solution_block(P: float, T: float, composition: dict[str, float], pH: Union[float, None], trace_approximation: bool=False) -> list[str]:
    """
    Generates a SOLUTION block for a PHREEQC input file.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    composition : dict[str, float]
        Molality of each element in the solution in mol/kgw.
    pH : Union[float, None]
        pH of the solution. If None, then the solution will be charge balanced.
    trace_approximation : bool, optional
        Sets molality of any element below 1e-9 to 0 if True, by default False.

    Returns
    -------
    list[str]
        List of text lines for input file.
    """

    if 'Alkalinity' in composition.keys():
        pH_line = ''
    else:
        pH_line = f'    pH        {pH:.8f}' if pH is not None else '    pH        7.0 charge'

    lines: list[str] = [
        'SOLUTION 1',
        f'    pressure  {P / EARTH_ATM:.4f}',
        f'    temp      {T + ABSOLUTE_ZERO:.4f}',
        pH_line,
        f'    units     mol/kgw'
    ]

    for k in composition.keys():
        molality = 0 if composition[k] < 1e-9 and trace_approximation else composition[k]# any concentration too low is brought up to a trace amount of 1e-9
        lines.append(f'    {k}    {molality:.15e}')

    lines.append('')

    return lines

def equilbrium_phases_block(phases: list[str], amounts: Union[list[float], None], saturation_indexes: Union[list[float], None]) -> list[str]:
    """
    Generates a EQUILBRIUM_PHASES block for a PHREEQC input file.

    Parameters
    ----------
    phases : list[str]
        A list of phases to include.
    amounts : Union[list[float], None]
        The amount of each phase in mol. If None, then all values are set to zero.
    saturation_indexes : Union[list[float], None]
        The saturation index of each phase in mol. If None, then all values are set to zero.

    Returns
    -------
    list[str]
        List of text lines for input file.
    """


    lines: list[str] = ['EQUILIBRIUM_PHASES 1']

    amounts = list(np.zeros(len(phases))) if amounts is None else amounts
    saturation_indexes = list(np.zeros(len(phases))) if saturation_indexes is None else saturation_indexes

    for i, phase in enumerate(phases):
        lines.append(f'    {phase}  {saturation_indexes[i]}  {amounts[i]}')

    lines.append('')

    return lines

def kinetics_block(phases: list[str], amounts: list[float], dt: float, specific_surface_area: float, bad_steps_max: int=100, dissolution_only: bool=False) -> list[str]:
    """
    Generates a KINETICS block for a PHREEQC input file.

    Parameters
    ----------
    phases : list[str]
        A list of phases to include.
    amounts : list[str]
        The amount of each phase in mol.
    dt : float
        Length of time to run kinetics calculation in s.
    specific_surface_area : float
        Specific surface area of the mineral phases in m^2/kg.
    dissolution_only : bool, optional
        If True, will return a dissolution rate of 0 if precipitation is occuring. 

    Returns
    -------
    list[str]
        List of text lines for input file.
    """

    lines: list[str] = ['KINETICS 1', '']

    for i, phase in enumerate(phases):
        if phase not in insoluble_minerals:
            lines.append(f'    {phase}')
            lines.append(f'    -m0 {amounts[i]}')
            lines.append(f'    -parms  0 {specific_surface_area * 1000} 0 {2 if dissolution_only else 0}')
            lines.append('')

    lines.append('    -cvode true')
    lines.append(f'    -steps {dt}')
    lines.append(f'    -bad_step_max {bad_steps_max}')
    lines.append('')

    return lines

def output_block(saturation_indexes: list[str]=[], equilibrium_phases: list[str]=[], kinetic_reactants: list[str]=[], activities: list[str]=[]) -> list[str]:
    """
    Generates a SELECTED_OUTPUT block for a PHREEQC input file.

    Parameters
    ----------
    saturation_indexes : list[str]
        A list of the phase saturation indices to output, default [].
    equilibrium_phases : list[str]
        A list of the phase amount changes to output for equilbrium phase calculation, default [].
    kinetic_reactants : list[str]
        A list of the phase amount changes to output for kinetics calculation, default [].
    activities : list[str]
        A list of activities indices to output, default [].

    Returns
    -------
    list[str]
        List of text lines for input file.
    """

    lines: list[str] = [
        'SELECTED_OUTPUT',
        '    -file output.txt',
        f'    -totals {available_element_string}',
        '    -alkalinity'
        ]
    
    if equilibrium_phases:
        lines.append('    -equilibrium_phases ' + ' '.join(equilibrium_phases))
    if saturation_indexes:
        lines.append('    -saturation_indices ' + ' '.join(saturation_indexes))
    if kinetic_reactants:
        lines.append('    -kinetic_reactants ' + ' '.join(kinetic_reactants))
    if activities:
        lines.append('    -activities ' + ' '.join(activities))

    lines.append('')

    return lines

def knobs_block() -> list[str]:
    """
    Generates a KNOBS block for a PHREEQC input file.

    Returns
    -------
    list[str]
        List of text lines for input file.
    """
    
    lines: list[str] = [
        'KNOBS',
        f'    -step_size 5.0',
        f'    -pe_step_size 2.0',
        f'    -diagonal_scale true',
        ''
        ]
    
    return lines

def run_PHREEQC(lines: list[str], single_output: bool=False) -> pd.DataFrame:
    """
    Runs PHREEQC with given input blocks.

    Parameters
    ----------
    lines : list[str]
        List of text lines conatining the input blocks.
    single_output : bool
        If the calculation is only going to solve the solution and not run any further calculations, default False. 

    Returns
    -------
    pd.DataFrame
        Dataframe containing the selected output.

    Raises
    ------
    PHREEQCError
        Raise if the reaction calculation fails or if there is no calculation run.
    """

    if not executable_path.exists(): # type: ignore
        raise FileNotFoundError(f"PHREEQC executable not found at {executable_path}. Did you run the build script?")
    
    if not database_path.exists(): # type: ignore
        raise FileNotFoundError(f"PHREEQC database (basic_v2.dat) not found in {database_path.parent}.") # type: ignore

    # Create a temporary directory that is deleted automatically when the 'with' block ends
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path("tests") if debug_mode else Path(temp_dir)
        
        input_file = temp_path / "input"
        output_file = temp_path / "output.txt"

        full_input = [f'DATABASE {database_path.as_posix()}', ''] + lines + ['END'] # type: ignore
        with open(input_file, 'w') as f:
            f.write('\n'.join(full_input))

        # We set cwd=temp_path so PHREEQC runs entirely inside the temp folder
        try:
            subprocess.run(
                [str(executable_path), "input"],
                cwd=temp_path,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True # Raises CalledProcessError if PHREEQC crashes
            )
        except subprocess.CalledProcessError:
            raise PHREEQCError("PHREEQC binary crashed or returned a non-zero exit code.")

        if not output_file.exists():
            raise PHREEQCError("PHREEQC did not generate an output file.")

        try:
            output = pd.read_csv(output_file, sep=r'\s+') # type: ignore
        except pd.errors.EmptyDataError:
             raise PHREEQCError("PHREEQC output file was empty.")

        if len(output) <= 1 and not single_output:
            raise PHREEQCError("PHREEQC failed calulation.")

    return output

def get_output_pH(output_df: pd.DataFrame) -> float:
    """
    Returns pH in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.

    Returns
    -------
    float
        pH of the resulting solution.
    """
    return float(output_df.at[1, 'pH']) # type: ignore

def get_output_composition(output_df: pd.DataFrame) -> dict[str, float]:
    """
    Returns solution composition (molalities) in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.

    Returns
    -------
    dict[str, float]
        Element molalities in resulting solution in mol/kgw.
    """
    new_comp: dict[str, float] = {}
    for element in available_elements:
        new_comp[element] = float(output_df.at[1, element]) # type: ignore
    return new_comp

def get_output_saturation_indexes(output_df: pd.DataFrame, get_initial: bool=False) -> dict[str, float]:
    """
    Returns phase saturation indices in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.
    get_initial : bool
        Whether to the the initial value from the dataframe, default False

    Returns
    -------
    dict[str, float]
        Phase saturation indices.
    """
    saturation_indexes: dict[str, float] = {}
    phases = [key for key in output_df.columns if key.startswith('si_')]
    for phase in phases:
        i = 0 if get_initial else 1
        saturation_indexes[phase[3:]] = float(output_df.at[i, phase]) # type: ignore
    return saturation_indexes

def get_output_equilbrium_phases(output_df: pd.DataFrame) -> dict[str, float]:
    """
    Returns change in amounts for each phase in the equilibrium phase calculation in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.

    Returns
    -------
    dict[str, float]
        Change in amounts for each phase in mol.
    """
    delta_moles: dict[str, float] = {}
    phases = [key for key in output_df.columns if key.startswith('d_')]
    for phase in phases:
        delta_moles[phase[2:]] = float(output_df.at[1, phase]) # type: ignore
    return delta_moles

def get_output_kinetics_phases(output_df: pd.DataFrame) -> dict[str, float]:
    """
    Returns change in amounts for each phase in the kinetics calculation in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.

    Returns
    -------
    dict[str, float]
        Change in amounts for each phase in mol.
    """
    delta_moles: dict[str, float] = {}
    non_zero_minerals = [key for key in output_df.columns if key.startswith('dk_')]
    for mineral in non_zero_minerals:
        if mineral not in insoluble_minerals:
            delta_moles[mineral[3:]] = float(output_df.at[1, mineral]) # type: ignore
    return delta_moles

def get_output_activities(output_df: pd.DataFrame, get_initial: bool=False) -> dict[str, float]:
    """
    Returns solution activities in selected output file.

    Parameters
    ----------
    output_df : pd.DataFrame
        DataFrame from selected output file.

    Returns
    -------
    dict[str, float]
        Activities for each phase in mol/kgw.
    """
    activities: dict[str, float] = {}
    phases = [key for key in output_df.columns if key.startswith('la_')]
    for phase in phases:
        i = 0 if get_initial else 1
        activities[phase[3:]] = 10 ** float(output_df.at[i, phase]) # type: ignore
    return activities   

def solve_chemistry(P: float, T: float, composition: dict[str, float], pH: Union[float, None]) -> pd.DataFrame:
    """
    Solves for the chemical state of a solution.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    composition : dict[str, float]
        Molality of each element in the solution in mol/kgw.
    pH : Union[float, None]
        pH of the solution. If None, then the solution will be charge balanced.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the selected output of the chemical calculation.
    """
    
    input_lines = solution_block(P, T, composition, pH) + output_block(saturation_indexes=available_minerals, activities=['H+', 'HCO3-', 'CO3-2'])

    return run_PHREEQC(input_lines, single_output=True)

def equilbriate_phases(P: float, T: float, composition: dict[str, float], pH: Union[float, None], phase_amounts: dict[str, float]) -> pd.DataFrame:
    """
    Dissolves moles of mineral or gas phases in a solution.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    composition : dict[str, float]
        Molality of each element in the solution in mol/kgw.
    pH : Union[float, None]
        pH of the solution. If None, then the solution will be charge balanced.
    phase_amounts : dict[str, float]
        The amount of each phase in mol.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the selected output of the chemical calculation.
    """

    phases = list(phase_amounts.keys())
    amounts = list(phase_amounts.values())
    
    input_lines = solution_block(P, T, composition, pH) + equilbrium_phases_block(phases, amounts, None) + output_block(equilibrium_phases=phases, saturation_indexes=available_minerals)

    return run_PHREEQC(input_lines)

def kinetics(P: float, T: float, composition: dict[str, float], pH: Union[float, None], phase_amounts: dict[str, float], dt: float, specific_surface_area: float=0.01, bad_steps_max: int=100) -> pd.DataFrame:
    """
    Dissolves minearl phases with chemical kinetic calculations.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    composition : dict[str, float]
        Molality of each element in the solution in mol/kgw.
    pH : Union[float, None]
        pH of the solution. If None, then the solution will be charge balanced.
    phase_amounts : dict[str, float]
        The amount of each phase in mol.
    dt : float
        Kinetics calculatiion time step in s.
    specific_surface_area : float
        Specific surface area of dissolving/precipitating phase in m^2/kg, by default 0.01.
    bad_steps_max : int
        Maximum number of bad steps before the PHREEQC solver stops, by deafult 100.

    Returns
    -------
    pd.DataFrame
        Dataframe containing the selected output of the chemical calculation.
    """

    phases = list(phase_amounts.keys())
    amounts = list(phase_amounts.values())
    
    input_lines = knobs_block() + \
        solution_block(P, T, composition, pH) + \
        kinetics_block(phases, amounts, dt, specific_surface_area, bad_steps_max) + \
        output_block(kinetic_reactants=phases, saturation_indexes=available_minerals, activities=['H+', 'HCO3-', 'CO3-2'])

    return run_PHREEQC(input_lines)

class PHREEQCError(Exception):
    pass
