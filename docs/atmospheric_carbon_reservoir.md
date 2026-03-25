# Separate Atmospheric Carbon Reservoir — Implementation Plan

## Motivation

The current model lumps the atmosphere and ocean into a single carbon reservoir.
`P_CO2` is derived from ocean DIC via PHREEQC and tracked as a state variable
(`Y[1]`). This causes a bookkeeping error for Earth-like planets with land:

- Continental silicate weathering (WHK) consumes CO2 **from the atmosphere** and
  delivers HCO3⁻ to the ocean via rivers.
- In the lumped model, `F_cont_C` adds C directly to ocean DIC, which **raises**
  P_CO2 — the opposite of the correct CO2 drawdown signal.
- This creates a structural carbon imbalance: outgassing (7.5 Tmol C/yr) plus
  continental weathering (8 Tmol C/yr) = 15.5 Tmol C/yr source, requiring
  ~15.5 Tmol C/yr of CaCO3 burial to balance. The Alk budget cannot support
  this (only ~9 Tmol Alk/yr supplied vs 31 Tmol/yr needed).

On a **fully ocean planet** (`land_fraction = 0`) this issue does not arise
because `F_cont = 0`. The fix is only needed when `land_fraction > 0`.

In the fast-exchange limit (ocean-atmosphere equilibration timescale << model
timestep), adding a separate `C_atm` reservoir is mathematically equivalent to
setting `F_cont_C = 0` in the lumped model — but the architecture is physically
correct and extensible.

## Diagnostic results (from `earth_debug.py`)

At reference conditions (T=288 K, P_CO2=280 ppm, Earth-like ocean):

| Flux | Alk (Tmol/yr) | C (Tmol/yr) | Ca (Tmol/yr) |
|---|---|---|---|
| F_out (volcanic) | 0.00 | 7.50 | 0.00 |
| F_cont (WHK ref) | 8.00 | 8.00 | 4.00 |
| F_diss (seafloor) | 1.80 | 0.00 | 0.43 |
| F_prec (current) | -1.96 | -1.96 | ~-1.4 |
| **F_NET** | **+5.87** | **+13.54** | **+3.45** |

Needed burial for C balance: **-15.5 Tmol/yr** (actual: -2.0 Tmol/yr, 8× deficit).

With `F_cont_C = 0` (Fix A) and halved outgassing (Fix B):
- Needed burial: -3.75 Tmol/yr (achievable with SI(Calcite) > 0.9, Ca ~1 mmol/kgw)
- Alk balance closes to within ~2.3 Tmol/yr (coverable by seafloor weathering)

## Required code changes

### `planet2.py` — major changes

#### 1. New constant
```python
M_CO2_ATM = 0.044  # kg/mol, molar mass of CO2
```

#### 2. State vector: replace Y[1]
`Y[1]` changes from `P_CO2` (Pa) to `C_atm` (mol CO2 in atmosphere,
normalised by `ocean_water_mass` → mol/kgw for dimensional consistency).

Conversion in both directions:
```python
# mol/kgw → Pa
P_CO2_atm = Y[1] * self.ocean_water_mass * self.gravity * M_CO2_ATM / self.surface_area

# Pa → mol/kgw  (for initialisation)
C_atm_init = P_CO2_init * self.surface_area / (self.gravity * M_CO2_ATM * self.ocean_water_mass)
```

#### 3. Initialisation
```python
# Currently:
Y[1] = 280e-6 * self.P_background

# Replace with:
P_CO2_init = 280e-6 * EARTH_ATM  # 280 ppm in Pa
Y[1] = P_CO2_init * self.surface_area / (self.gravity * M_CO2_ATM * self.ocean_water_mass)
```

#### 4. `solve_climate_from_chemistry`
- Signature unchanged externally; internally derive `P_CO2_atm` from `C_atm`
  before the Newton solve.
- The climate model (`get_T_surface`) and PHREEQC both receive `P_CO2_atm`.
- The function no longer needs the `P_CO2_est` argument (P_CO2 comes from
  `C_atm`, not from the previous ocean DIC state).

```python
def solve_climate_from_chemistry(self, b_ocean, C_atm):
    P_CO2_atm = C_atm * self.ocean_water_mass * self.gravity * M_CO2_ATM / self.surface_area
    # ... Newton solve as before, using P_CO2_atm in place of P_CO2_est ...
```

#### 5. `_compute_fluxes_and_derivatives` — flux routing changes

**Outgassing — split subaerial/submarine by land fraction:**
```python
# Subaerial fraction goes to atmosphere (continental/island-arc volcanism)
# Submarine fraction goes directly to ocean
f_subaerial = self.land_fraction   # proxy; 0 on ocean planet → all submarine
f_submarine = 1.0 - f_subaerial

F_out_to_atm   = self.outgassing_flux * f_subaerial / self.ocean_water_mass
F_out_to_ocean = self.outgassing_flux * f_submarine / self.ocean_water_mass
```

**Continental weathering — route C through atmosphere:**
```python
# get_continental_weathering_flux returns Alk, C, Ca as before.
# C is redirected: atmosphere loses it, ocean gains it from rivers.
F_cont_per_area  = get_continental_weathering_flux(T_new, P_CO2_atm)
F_cont_ocean     = (F_cont_per_area * self.land_area) / self.ocean_water_mass
F_cont_C_to_atm  = -F_cont_ocean[_C_idx]   # C removed from atmosphere (negative)
# F_cont_ocean already has the C going into the ocean — routing is unchanged
# for the ocean; only the atmosphere budget changes.
```

**Gas exchange — couples C_atm and ocean DIC:**
```python
# PHREEQC gives equilibrium P_CO2 for current ocean DIC
P_CO2_ocean = get_P_CO2(P_atm, T_new, b_ocean)

# Atmosphere P_CO2 from state variable
C_atm = Y[1]
P_CO2_atm = C_atm * self.ocean_water_mass * self.gravity * M_CO2_ATM / self.surface_area

# Exchange flux (positive = ocean degasses to atmosphere)
# tau_exchange ~ 1000 yr (same order as current tau_atm)
F_exchange_C = ((P_CO2_ocean - P_CO2_atm) / (self.gravity * M_CO2_ATM / self.surface_area)
                / (tau_exchange * self.ocean_water_mass))
# Adds to C_atm, removes from ocean C (or vice versa)
```

**`dYdt` construction:**
```python
dYdt[1] = (                     # dC_atm/dt
    F_out_to_atm[_C_idx]        # subaerial volcanic source
    + F_cont_C_to_atm           # weathering consumption (negative)
    + F_exchange_C              # gas exchange (positive = ocean→atm)
)

# Ocean C: submarine outgassing + river HCO3- + gas exchange + burial/dissolution
F_net_ocean = F_out_to_ocean + F_cont_ocean + F_diss + F_prec
F_net_ocean[_C_idx] -= F_exchange_C   # exchange removes C from ocean when +ve
dYdt[2:] = F_net_ocean
```

#### 6. CSV / JSON output
- Add `C_atm_mol` column to timeseries (= `Y[1] * ocean_water_mass`)
- Add `P_CO2_ppm` as a derived diagnostic (= `P_CO2_atm / EARTH_ATM * 1e6`)
- Summary JSON `P_CO2_Pa` field: compute from `C_atm` at final state

### `ocean_chemistry.py` — no required changes

`get_continental_weathering_flux` returns the correct stoichiometry (Alk, C, Ca).
The routing of where that C goes (withdrawn from atmosphere, delivered to ocean)
is handled entirely in `planet2.py`.

Optional cosmetic improvement: split into `get_continental_ocean_flux` (Alk, C,
Ca → ocean) and `get_continental_atm_flux` (C → atmosphere, negative sign) for
explicitness — but this is not required.

### `earth_validation.py` — trivial changes

`b_init` and planet parameters are unchanged. Update the output section to derive
`P_CO2_ppm` from `C_atm` via the conversion formula rather than reading
`P_CO2_Pa` directly, or keep `P_CO2_Pa` as a derived field in the summary JSON.

## Scope summary

| File | Nature of change | Size |
|---|---|---|
| `planet2.py` | Y vector, init, solve_climate, _compute_fluxes, dYdt routing, output | Major |
| `ocean_chemistry.py` | None required (optional cosmetic split) | None / Minor |
| `earth_validation.py` | Output reading only | Trivial |

## Fast-exchange equivalence

In the limit `tau_exchange → 0`, gas exchange is instantaneous and C_atm is
always determined by ocean DIC via PHREEQC equilibrium. In this limit:

- `F_cont` moves C from atmosphere to ocean with zero net change to total C.
- Subaerial outgassing → atmosphere → immediately partitions ~98% to ocean
  (ocean holds ~50× more C). Effectively the same as adding to ocean DIC.

This is **mathematically equivalent** to setting `F_cont_C = 0` in the current
lumped model. The separate reservoir adds physical clarity and extensibility but
does not change the geological-timescale dynamics.

## Future extensions enabled by this architecture

- **Organic carbon burial**: add `F_organic_burial` as a `C_atm` sink
  proportional to `land_fraction` (or productivity proxy). This is the missing
  ~3.5 Tmol C/yr that currently prevents Earth budget closure.
- **Slow exchange scenarios**: snowball Earth, stagnant-lid planets where
  ocean-atmosphere equilibration is slow relative to geological processes.
- **Explicit atmospheric mass tracking**: relevant for CO2-dominated atmospheres
  (Venus-like) where the atmosphere holds more C than the ocean.

## Related fixes still needed regardless of this change

1. **Initial Ca**: raise from 7×10⁻⁵ to ~10⁻³ mol/kgw so SI(Calcite) > 0.9
   and precipitation acts as an effective thermostat.
2. **Outgassing calibration**: `EARTH_OUTGASSING` represents total volcanic CO2;
   ~half is balanced by organic burial (absent here). Either halve `F_out_C` for
   inorganic-only scenarios, or implement organic burial (see above).
3. **PHREEQC charge balance**: initial conditions must satisfy
   2×Ca + 2×Mg + Na < Alk to avoid "non-carbonate alkalinity > total alkalinity"
   errors as Ca accumulates.
