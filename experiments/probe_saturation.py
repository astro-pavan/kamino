"""Recompute per-mineral saturation indices and precipitation fluxes at a run's FINAL ocean state.

Usage:  PROBE_GLOB='/path/to/*.json' python experiments/probe_saturation.py

Diagnostic behind §26. The SI tells you whether a precipitation timescale is physics or numerics:
a phase sitting at SI ~ 0 has reached equilibrium, so its flux is set by solute supply and its
timescale is a free numerical parameter; a phase held at SI >> 0 is kinetically limited, so its
timescale IS the flux and changing it changes the answer. Check this before tuning any relaxation
timescale in the model.

The two precipitation calls in planet.dY_dt are independent PHREEQC equilibrations on the
same b_ocean, so tau_prec and tau_rw never compete inside one solve -- the separation that
§21 describes acts THROUGH the shared ocean state over time. So the question this probes is
whether shrinking the separation from 50x to 7x left the reverse-weathering phases (chiefly
Sepiolite(d)) supersaturated at the attractor.
"""
import os, sys, json, glob, warnings
os.environ.setdefault('JAX_PLATFORMS', 'cpu')
os.environ.setdefault('MPLCONFIGDIR', os.path.expanduser('~/.cache/kamino-mpl'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
warnings.simplefilter('ignore')
import numpy as np
from kamino.precipitation import get_precipitation_by_mineral
from kamino.chemistry import elements
from kamino.mineral_info import (carbonate_minerals, clay_minerals, silica_minerals,
                                 reverse_weathering_minerals)
from kamino.constants import YR

FAST = carbonate_minerals + clay_minerals + silica_minerals
RW = reverse_weathering_minerals
G = 9.81

for f in sorted(glob.glob(os.environ['PROBE_GLOB'])):

    d = json.load(open(f))
    y = np.array(d['data']['y'])
    P_CO2, P_H2O = y[0][-1], y[1][-1]
    b = y[2:2 + len(elements), -1]
    # planet.py: T_seafloor = max(1.02*T_surface - 16.7, 274); P_pore adds the water column
    T_sf = max(1.02 * d['T'] - 16.7, 274.0)
    P_pore = (d['background_pressure'] + P_CO2 + P_H2O) + 1000 * G * d['ocean_depth']

    fast, _, si_f = get_precipitation_by_mineral(
        P_pore, T_sf, b, precipitating_minerals=FAST, precipitation_timescale=d['tau_prec'])
    rw, _, si_r = get_precipitation_by_mineral(
        P_pore, T_sf, b, precipitating_minerals=RW, precipitation_timescale=d['tau_rw'])

    print(f"\n{os.path.basename(f)}  mg/si={d['mantle_mg_si']}  "
          f"tau_prec={d['tau_prec']/YR:.3g} yr  tau_rw={d['tau_rw']/YR:.3g} yr  "
          f"ratio={d['tau_rw']/d['tau_prec']:.1f}x  T={d['T']:.2f}")
    si = {**si_f, **si_r}
    for name, flux in list(fast.items()) + list(rw.items()):
        tag = 'RW ' if name in RW else 'fast'
        # flux is per-element; report the largest-magnitude component in mmol/kgw/yr
        mag = np.abs(flux).max() * 1e3 * YR
        if mag < 1e-9 and abs(si.get(name, -99)) > 5:
            continue
        print(f"   {tag} {name:14s} SI={si.get(name, float('nan')):+7.3f}  |flux|max={mag:.4g} mmol/kgw/yr")
