import numpy as np
import numpy.typing as npt

from kamino.chemistry import solve_solution, elements

def get_precipitation(P: float, T: float, b: npt.NDArray[np.float64], precipitating_minerals: list[str], equilibrium_minerals: list[str]=[], fO2: float=0, precipitation_timescale: float=0) -> tuple[npt.NDArray[np.float64], float, dict[str, float]]:

    output = solve_solution(P, T, b, precipitating_minerals=precipitating_minerals, equilibriating_minerals=equilibrium_minerals, fO2=fO2)

    aqueous_fluxes = np.zeros(elements.shape)

    for i, element in enumerate(elements):
        output_key = 'Alk(eq/kgw)' if element == 'Alkalinity' else f'{element}(mol/kgw)'
        aqueous_fluxes[i] = (float(output[output_key][-1]) - float(output[output_key][0]))

    pH = float(output['pH'][-1])

    si_dict = {}
    for min_name in precipitating_minerals:
        si_key = f'si_{min_name}'
        if si_key in output:
            si_dict[min_name] = float(output[si_key][0]) # Index 0 is the initial state of the solution BEFORE precipitation occurs

    # enforce always negative flux
    aqueous_fluxes = np.minimum(aqueous_fluxes, 0)

    if precipitation_timescale > 0:

        aqueous_fluxes = aqueous_fluxes / precipitation_timescale

    return aqueous_fluxes, pH, si_dict