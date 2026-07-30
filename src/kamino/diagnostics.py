"""Replay a saved Planet run and extract a full source/sink budget for one instant.

time_evolve only persists (t, Y). Everything else the model computes — the individual
flux terms, the mineral saturation indices, the Damkohler number — is derived state that
lives for one call of Planet.dY_dt and is then overwritten. This module rebuilds the
Planet from its saved config, evaluates dY_dt once on a chosen state, and returns the
derived quantities in plotting units.

Fluxes are reported in Tmol/yr (Teq/yr for Alkalinity), integrated over the whole ocean,
with the model's sign convention: positive = source to the ocean, negative = sink.
"""

import os
import json
import numpy as np

from kamino.constants import YR, M_EARTH, R_EARTH
from kamino.chemistry import elements, ChemistryError
from kamino.planet import Planet
from kamino.precipitation import get_precipitation_by_mineral
from kamino.mineral_info import (
    carbonate_minerals, clay_minerals, silica_minerals,
    reverse_weathering_minerals, evaporite_minerals,
)

# Element order is fixed by kamino.chemistry.elements; S is pinned background, not evolved.
PLOT_ELEMENTS = [e for e in elements if e != 'S']

# Modern seawater, mol/kgw — for the "vs Earth" reference markers.
EARTH_SEAWATER = {
    'Alkalinity': 2.3e-3,
    'C':          2.0e-3,
    'Si':         0.1e-3,
    'Ca':        10.3e-3,
    'Mg':        52.8e-3,
    'Na':       480.0e-3,
    'Cl':       550.0e-3,
}

# Which processes are inherently sources vs sinks, for ordering the ledger rows.
SOURCE_TERMS = ['outgassing', 'seafloor LT', 'seafloor HT', 'continental']
SINK_TERMS   = ['precipitation', 'shelf carbonate', 'reverse weathering',
                'biogenic', 'Cl subduction', 'Na albitization']
TRANSFER_TERMS = ['HT Mg-Ca exchange']
ALL_TERMS = SOURCE_TERMS + SINK_TERMS + TRANSFER_TERMS


def planet_from_config(config: dict) -> Planet:
    """Rebuild a Planet from the config block of a run JSON.

    Planet.__init__ writes a config JSON as a side effect. We are only borrowing the
    object to evaluate dY_dt, so that file is deleted again: leaving it behind would
    litter output/ with '*_diagnostic.json' files that plot_results.load_data globs up.
    """
    planet = Planet(
        mass=float(config.get('mass', M_EARTH)),
        radius=float(config.get('radius', R_EARTH)),
        background_pressure=float(config['background_pressure']),
        instellation=float(config['instellation']),
        crust_production_rate=float(config['crust_production_rate']),
        outgassing=float(config['outgassing']),
        ocean_depth=float(config['ocean_depth']),
        land_fraction=float(config.get('land_fraction', 0.0)),
        crust_composition=config['crust_composition'],
        reverse_weathering=bool(config.get('reverse_weathering', True)),
        alpha=float(config.get('alpha', 1.43)),
        f_HT=float(config.get('f_HT', 0.0)),
        cl_outgassing_ratio=float(config.get('cl_outgassing_ratio', 0.02)),
        tau_prec=float(config.get('tau_prec', 100e3 * YR)),
        tau_rw=float(config.get('tau_rw', 5e6 * YR)),
        f_bio=float(config.get('f_bio', 0.0)),
        name=str(config.get('name', 'planet')) + '_diagnostic',
        climate_model=str(config.get('climate_model', 'analytic')),
    )

    try:
        os.remove(planet._output_filename)
    except OSError:
        pass

    return planet


def diagnose(planet: Planet, Y: np.ndarray, t: float = 0.0) -> dict:
    """Evaluate the model at state Y and return its full source/sink budget.

    Returns a dict with:
      fluxes        {process: {element: Tmol/yr}}   signed, + source / - sink
      net           {element: Tmol/yr}              sum over processes
      minerals      {reservoir: {mineral: {element: Tmol/yr}}}  per-phase precipitation
      SI            {reservoir: {mineral: SI}}      saturation index, > 0 = precipitating
      composition   {element: mol/kgw}              ocean composition
      scalars       T, P_CO2, pH, Da, salinity, ...
      residence     {element: yr}                   b / |net|, the imbalance timescale
    """
    dYdt = planet.dY_dt(t, np.asarray(Y, dtype=float))
    if not hasattr(planet, '_flux_terms'):
        raise ChemistryError('dY_dt did not reach the chemistry block for this state')

    # mol/kgw/s -> Tmol/yr over the whole ocean
    to_Tmol_yr = planet.ocean_water_mass * YR / 1e12

    b_ocean = np.maximum(np.asarray(Y, dtype=float)[2:-1], 0.0)
    state = planet._state

    fluxes = {
        term: {el: float(vec[i]) * to_Tmol_yr
               for i, el in enumerate(elements) if el in PLOT_ELEMENTS}
        for term, vec in planet._flux_terms.items()
    }
    net = {el: sum(fluxes[term][el] for term in fluxes) for el in PLOT_ELEMENTS}

    # Per-mineral precipitation, recomputed with the same arguments dY_dt used.
    minerals: dict[str, dict[str, dict[str, float]]] = {'ocean': {}, 'pore': {}}
    fast = planet.fast_ocean_precipitating_minerals
    slow = planet.rw_ocean_precipitating_minerals

    def _by_mineral(P, T, mins, tau):
        if not mins:
            return {}
        try:
            per_min, _, _ = get_precipitation_by_mineral(
                P, T, b_ocean, precipitating_minerals=mins, precipitation_timescale=tau)
        except ChemistryError:
            return {}
        return {m: {el: float(v[i]) * to_Tmol_yr
                    for i, el in enumerate(elements) if el in PLOT_ELEMENTS}
                for m, v in per_min.items()}

    minerals['ocean'].update(_by_mineral(state['P_pore'], state['T_seafloor'], fast, planet.tau_prec))
    minerals['ocean'].update(_by_mineral(state['P_pore'], state['T_seafloor'], slow, planet.tau_rw))

    composition = {el: float(b_ocean[i]) for i, el in enumerate(elements) if el in PLOT_ELEMENTS}

    # Salinity: mass of dissolved ions per kg water (Alkalinity is a charge balance, not mass).
    # composition is mol/kgw and molar_mass is g/mol, so the product is already g/kgw.
    molar_mass = {'C': 61.0, 'Si': 60.1, 'Al': 27.0, 'Fe': 55.8,
                  'Ca': 40.1, 'Mg': 24.3, 'Na': 23.0, 'Cl': 35.45}
    salinity = sum(composition.get(el, 0.0) * m for el, m in molar_mass.items())  # g/kg

    residence = {}
    for el in PLOT_ELEMENTS:
        n = abs(net[el])
        b_tot = composition[el] * planet.ocean_water_mass / 1e12  # Tmol
        residence[el] = (b_tot / n) if n > 0 else np.inf

    scalars = {
        'T_surface':   state['T_surface'],
        'T_seafloor':  state['T_seafloor'],
        'P_CO2_bar':   state['P_CO2'] / 1e5,
        'pH_surface':  state['pH_surface'],
        'pH_seafloor': state['pH_seafloor'],
        'Da':          state['Da'],
        'supply_efficiency': state['supply_efficiency'],
        'A_reactive':  state['A_reactive'],
        'salinity':    salinity,
        'land_fraction': planet.land_fraction,
        'f_HT':        planet.f_HT,
        'f_bio':       planet.f_bio,
    }

    return {
        'fluxes':      fluxes,
        'net':         net,
        'minerals':    minerals,
        'SI':          {'ocean': state['ocean_SI'], 'pore': state['pore_SI']},
        'composition': composition,
        'scalars':     scalars,
        'residence':   residence,
        'to_Tmol_yr':  to_Tmol_yr,
    }


def diagnose_run(path: str, index: int = -1) -> dict:
    """Load a run JSON, rebuild the planet, and diagnose the state at `index` (default: final)."""
    with open(path) as f:
        d = json.load(f)

    y = np.array(d['data']['y'], dtype=float)
    t = np.array(d['data']['time'], dtype=float)
    Y = y[:, index]

    planet = planet_from_config(d)
    result = diagnose(planet, Y, t=float(t[index]))
    result['config'] = d
    result['time_yr'] = float(t[index]) / YR
    result['termination'] = d.get('termination', 'unknown')
    return result
