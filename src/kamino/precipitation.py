import numpy as np
import numpy.typing as npt

from kamino.chemistry import solve_solution, elements, stoichiometry
from kamino.mineral_info import PRECIPITATE_MOLAR_MASS, PRECIPITATE_DENSITY

def get_precipitation_by_mineral(P: float, T: float, b: npt.NDArray[np.float64], precipitating_minerals: list[str], equilibrium_minerals: list[str]=[], fO2: float=0, precipitation_timescale: float=0, pe: float | None=None) -> tuple[dict[str, npt.NDArray[np.float64]], float, dict[str, float], dict[str, float]]:
    """As get_precipitation, but keeps each mineral's aqueous flux vector separate.

    Returns {mineral: flux_vector} rather than their sum, so a diagnostic can attribute
    an element sink to the phase responsible for it, plus {mineral: molar rate} -- the
    moles of the MINERAL itself, which the aqueous vectors cannot express (they carry
    only tracked elements, so a phase's H and O are invisible in them and Halite's mass
    is split across two entries). sediment_volume_rate consumes it.

    Note the fluxes here are NOT clamped at zero: the "precipitation only removes
    elements" clamp is applied by get_precipitation to the SUM, which allows one
    mineral's positive contribution to an element to cancel another's negative one
    before clipping. Clamping per mineral would not be equivalent. Individual vectors
    may therefore carry a small positive entry; only their sum is guaranteed <= 0.
    """

    output = solve_solution(P, T, b, precipitating_minerals=precipitating_minerals, equilibriating_minerals=equilibrium_minerals, fO2=fO2, pe=pe)

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
    fluxes_by_mineral: dict[str, npt.NDArray[np.float64]] = {}
    moles_by_mineral: dict[str, float] = {}

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
        flux = -moles_prec * stoich_vec * smooth

        rate = moles_prec * smooth
        if precipitation_timescale > 0:
            flux = flux / precipitation_timescale
            rate = rate / precipitation_timescale

        fluxes_by_mineral[min_name] = flux
        moles_by_mineral[min_name] = rate

    return fluxes_by_mineral, pH, si_dict, moles_by_mineral


def sum_precipitation(fluxes_by_mineral: dict[str, npt.NDArray[np.float64]]) -> npt.NDArray[np.float64]:
    """Clamped element-flux sum over minerals: precipitation can only remove elements.

    The clamp is applied to the SUM, so one mineral's positive contribution to an element
    can cancel another's negative one before clipping -- clamping per mineral would not be
    equivalent. Kept as one definition because Planet applies it to the fast and reverse-
    weathering assemblages separately.
    """
    aqueous_fluxes = np.zeros(elements.shape)
    for flux in fluxes_by_mineral.values():
        aqueous_fluxes += flux
    return np.minimum(aqueous_fluxes, 0)


def sediment_volume_rate(moles_by_mineral: dict[str, float]) -> float:
    """Volume of solid sediment produced per kg of seawater, in m^3 / kgw (per second if the
    molar rates carry a precipitation timescale).

    Each phase is converted by its OWN molar mass and density. The previous formulation mapped
    the whole carbon sink onto calcite and the whole silicon sink onto quartz-density silica,
    which both under-counted (it ignored the clays, evaporites and reverse-weathering phases
    entirely) and mis-counted (a mole of Sepiolite(d) carries 6 mol Si but occupies 287 cm^3,
    not 6 x 22.7 cm^3 of SiO2(am)).

    A phase with no entry in the tables raises rather than being silently skipped: a missing
    density would quietly remove that phase from the burial rate, which is the failure this
    function exists to fix.
    """
    volume = 0.0
    for mineral, moles in moles_by_mineral.items():
        try:
            M = PRECIPITATE_MOLAR_MASS[mineral]      # g/mol
            rho = PRECIPITATE_DENSITY[mineral]       # kg/m^3
        except KeyError:
            raise KeyError(
                f"{mineral!r} precipitates but has no entry in mineral_info."
                f"PRECIPITATE_MOLAR_MASS / PRECIPITATE_DENSITY; add one so it is "
                f"counted in the sedimentation rate."
            ) from None
        volume += moles * (M * 1e-3) / rho
    return volume


def get_precipitation(P: float, T: float, b: npt.NDArray[np.float64], precipitating_minerals: list[str], equilibrium_minerals: list[str]=[], fO2: float=0, precipitation_timescale: float=0, pe: float | None=None) -> tuple[npt.NDArray[np.float64], float, dict[str, float]]:

    fluxes_by_mineral, pH, si_dict, _ = get_precipitation_by_mineral(
        P, T, b, precipitating_minerals, equilibrium_minerals, fO2, precipitation_timescale, pe=pe,
    )

    return sum_precipitation(fluxes_by_mineral), pH, si_dict
