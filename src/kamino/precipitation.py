import numpy as np
import numpy.typing as npt

from kamino.chemistry import solve_solution, elements, stoichiometry

def get_precipitation(P: float, T: float, b: npt.NDArray[np.float64], precipitating_minerals: list[str], equilibrium_minerals: list[str]=[], fO2: float=0, precipitation_timescale: float=0) -> tuple[npt.NDArray[np.float64], float, dict[str, float]]:

    output = solve_solution(P, T, b, precipitating_minerals=precipitating_minerals, equilibriating_minerals=equilibrium_minerals, fO2=fO2)

    pH = float(output['pH'][-1])

    # Pre-precipitation SI for each mineral (index 0 = state before equilibration)
    si_dict = {}
    for min_name in precipitating_minerals:
        si_key = f'si_{min_name}'
        if si_key in output:
            si_dict[min_name] = float(output[si_key][0])

    # Reconstruct per-mineral aqueous fluxes with sigmoid smoothing at SI=0.
    # Each mineral's flux is scaled by smooth(SI) = 0.5*(1+tanh(SI*5)), which is C-inf
    # everywhere. This removes the kink in dY/dt at SI=0 that causes LSODA to take
    # tiny steps when a mineral sits near saturation at quasi-steady state.
    aqueous_fluxes = np.zeros(elements.shape)

    for min_name in precipitating_minerals:
        d_key = f'd_{min_name}'
        if d_key not in output:
            continue
        moles_prec = float(output[d_key][-1])  # moles precipitated per kgw (>= 0)
        if moles_prec == 0.0:
            continue

        stoich_vec = stoichiometry.get(min_name)
        if stoich_vec is None:
            continue

        si = si_dict.get(min_name, 0.0)
        smooth = 0.5 * (1.0 + np.tanh(si * 5.0))

        # stoich_vec gives elements released per mol dissolved; negate for precipitation
        aqueous_fluxes -= moles_prec * stoich_vec * smooth

    # Safety clamp: precipitation can only remove elements, never add them
    aqueous_fluxes = np.minimum(aqueous_fluxes, 0)

    if precipitation_timescale > 0:
        aqueous_fluxes = aqueous_fluxes / precipitation_timescale

    return aqueous_fluxes, pH, si_dict
