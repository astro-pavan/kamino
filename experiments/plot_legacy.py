"""Compatibility shim for sweep output that predates the current model.

`plot_results.py` assumes the CURRENT run schema and nothing else. Everything needed to read
older output lives here, as a one-way SCHEMA UPGRADE rather than a second set of plotting code:
`upgrade()` rewrites an old dataframe into the current vocabulary, after which every figure in
`plot_results` renders it unchanged.

    python experiments/plot_results.py --path <old_sweep_dir> --legacy

Keeping this separate matters because the legacy branches were scattered through the plotting
code as extra members of set literals and extra arms of if/elif chains, where they were easy to
mistake for current behaviour and impossible to test independently. The current sweep output
(2045 runs, 2026-08) exercises NONE of them.

What changed, and how it is mapped
----------------------------------
1. TERMINATIONS. Old runs encoded the OUTCOME in the termination; the model now records only why
   integration stopped, and the outcome is read from the final state (see Planet.time_evolve).

       snowball    -> out_of_domain, domain_wall 'cold'
       hothouse    -> out_of_domain, domain_wall 'hot'
       co2_ceiling -> out_of_domain, domain_wall 'co2_high'
       acid_ocean  -> out_of_domain, domain_wall 'co2_high'   (see note)
       co2_floor   -> out_of_domain, domain_wall 'co2_low'    (see note)

   `acid_ocean` was not a CO2 wall, but like the ceiling it stopped the run mid-evolution, so its
   fate is unknown; 'co2_high' is the current label that classifies it that way. `co2_floor` was
   previously counted as a snowball outcome; under the current classifier it becomes "unknown
   fate", which is the more honest call for a run that stopped while still evolving. This is a
   DELIBERATE reclassification, not a faithful translation.

2. CRUST COMPOSITION. Old sweeps named a fixed mineralogy in the run name (`_comp_basalt_49`)
   instead of deriving it from mantle Mg/Si and dIW. `basalt_49` was the reference crust, so it
   maps onto the current reference (Mg/Si 1.25, dIW -2). The others have no Mg/Si or dIW at all
   -- they were hand-written mineral dictionaries -- so they are given NaN on both axes, which
   excludes them from the reference-crust mask and from the composition figures. Use
   `plot_named_compositions` to see them.

3. CRUST CARBONATE. `crust_carbonate_content` was a discontinued experiment (carbonate already in
   the subducting crust). Runs with a non-zero value are dropped: the current model has no such
   term, so they are not comparable to anything else in the frame.

4. SURFACE TEMPERATURE. `T` in the JSON is planet.py's `self._T`, a side effect set on EVERY
   call to dY_dt -- including Jacobian probes and event root-finding trials -- so in old output
   it is not guaranteed to correspond to the accepted final state that `P_CO2` comes from. Near
   a domain wall that mattered: a 1% Jacobian perturbation could flip between a real ~389 K
   state and the analytic model's literal 400.0 "no equilibrium found" sentinel, an 11 K error.
   `Planet.time_evolve` now re-evaluates T on the accepted final state before writing, so
   current output needs no correction (measured over the 2026-08 sweep: mean |dT| 0.0013 K, max
   0.30 K, and every run above 0.01 K is a `wall_timeout`, where the expired deadline prevents
   the re-evaluation from running at all).

   `recompute_T` restores the old behaviour for files written before that fix. It is applied
   only under --legacy, and only to runs whose `climate_model` is 'analytic' or unrecorded: it
   evaluates the ANALYTIC climate model, so applying it to a run made with the clima
   interpolator would silently substitute a different model. That is exactly why it no longer
   runs unconditionally in plot_results.

5. T_p / mg_si_ratio. The oldest runs stored `mantle_potential_temperature` and a `mg_si_ratio`
   that was a multiplier on 1.23 rather than a Mg/Si. Their mineralogy cannot be reconstructed by
   the current pipeline; `plot_results._crust_composition_of` raises on them rather than
   substituting Earth, and that guard stays there because it protects a per-run lookup.
"""

import numpy as np
import pandas as pd

# Old termination -> (current termination, domain_wall). See the module docstring.
TERMINATION_MAP = {
    'snowball':    ('out_of_domain', 'cold'),
    'hothouse':    ('out_of_domain', 'hot'),
    'co2_ceiling': ('out_of_domain', 'co2_high'),
    'acid_ocean':  ('out_of_domain', 'co2_high'),
    'co2_floor':   ('out_of_domain', 'co2_low'),
}

# Named crust compositions from pre-MAGEMin sweeps, with the SiO2 content in the name.
COMP_LABELS = {
    'komatiite_42': 'Komatiite (42% SiO₂)',
    'komatiite_44': 'Komatiite (44% SiO₂)',
    'basalt_47':    'Basalt (47% SiO₂)',
    'basalt_49':    'Basalt (49% SiO₂)',
    'basalt_51':    'Basalt (51% SiO₂)',
}
COMP_ORDER = ['komatiite_42', 'komatiite_44', 'basalt_47', 'basalt_49', 'basalt_51']

# The named composition that the current (Mg/Si, dIW) reference reproduces.
REFERENCE_COMP = 'basalt_49'


def comp_name_of(run_name):
    """Named crust composition encoded in a run name, or '' for current-schema runs."""
    import re
    m = re.search(r'_comp_(.+?)(?:_rw|_fht_|$)', str(run_name))
    return m.group(1) if m else ''


def recompute_T(instellation, P_CO2_bar, land_fraction=0.0):
    """Surface temperature from (instellation, final P_CO2) using the ANALYTIC climate model.

    Imports from kamino at call time so this module stays importable without the model on the
    path. See item 4 of the module docstring for when this is the right thing to do.
    """
    from kamino.climate.analytic import get_T_surface_analytic
    from kamino.constants import SOLAR_CONSTANT
    from kamino.planet import OCEAN_ALBEDO, LAND_ALBEDO

    albedo = LAND_ALBEDO * land_fraction + OCEAN_ALBEDO * (1 - land_fraction)
    P_CO2_Pa = float(np.clip(P_CO2_bar * 1e5, 0.0, 1e6))
    return float(get_T_surface_analytic(instellation * SOLAR_CONSTANT, P_CO2_Pa, albedo, False))


def upgrade(df, verbose=True):
    """Rewrite a legacy dataframe into the current schema, in place on a copy.

    Returns the upgraded frame. Safe to call on current-schema data: it is then a no-op apart
    from adding the `comp_name` column, which stays empty.
    """
    df = df.copy()
    notes = []

    # -- 1. terminations ---------------------------------------------------------------------
    if 'domain_wall' not in df.columns:
        df['domain_wall'] = None
    for old, (new, wall) in TERMINATION_MAP.items():
        hit = df['termination'] == old
        n = int(hit.sum())
        if not n:
            continue
        df.loc[hit, 'termination'] = new
        # Only fill a wall that is not already recorded; never overwrite real data.
        missing = hit & df['domain_wall'].isna()
        df.loc[missing, 'domain_wall'] = wall
        notes.append(f"{n} run(s) {old!r} -> {new!r} / domain_wall {wall!r}")

    # -- 2. named crust compositions --------------------------------------------------------
    if 'comp_name' not in df.columns:
        df['comp_name'] = df['name'].map(comp_name_of)
    named = df['comp_name'].astype(bool)
    if named.any():
        is_ref = named & (df['comp_name'] == REFERENCE_COMP)
        other = named & ~is_ref
        # basalt_49 IS the current reference crust, so it keeps the reference axes it was
        # already given by load_data's defaults. The rest have no position on those axes.
        df.loc[other, ['mg_si', 'delta_iw']] = np.nan
        notes.append(f"{int(is_ref.sum())} run(s) {REFERENCE_COMP!r} treated as the reference "
                     f"crust; {int(other.sum())} other named-composition run(s) given NaN "
                     f"Mg/Si and dIW")

    # -- 3. surface temperature ---------------------------------------------------------------
    # Only where the analytic model was used; see item 4 of the module docstring.
    if {'instellation', 'P_CO2'} <= set(df.columns):
        model = df['climate_model'] if 'climate_model' in df.columns else 'analytic'
        ok = pd.Series(model, index=df.index).fillna('analytic').eq('analytic') \
            & df['P_CO2'].notna()
        if ok.any():
            before = df.loc[ok, 'T'].to_numpy(dtype=float)
            df.loc[ok, 'T'] = [
                recompute_T(float(r.instellation), float(r.P_CO2),
                            float(getattr(r, 'land_fraction', 0.0) or 0.0))
                for r in df.loc[ok].itertuples()]
            moved = np.nanmax(np.abs(df.loc[ok, 'T'].to_numpy(dtype=float) - before))
            notes.append(f"T recomputed from (instellation, P_CO2) for {int(ok.sum())} run(s); "
                         f"largest change {moved:.3g} K")
        skipped = int((~ok).sum())
        if skipped:
            notes.append(f"{skipped} run(s) left with their stored T (non-analytic climate "
                         f"model or no P_CO2)")

    # -- 4. crust carbonate ------------------------------------------------------------------
    if 'crust_carbonate' in df.columns:
        drop = df['crust_carbonate'].fillna(0.0) != 0.0
        if drop.any():
            notes.append(f"{int(drop.sum())} run(s) with crust_carbonate_content != 0 dropped "
                         f"(no equivalent term in the current model)")
            df = df[~drop]

    if verbose:
        if notes:
            print("Legacy upgrade:")
            for n in notes:
                print(f"  {n}")
        else:
            print("Legacy upgrade: nothing to change (data already uses the current schema).")
    return df


def plot_named_compositions(df, output_path, split_panels=False, show_markers=False):
    """The pre-MAGEMin crust figure: one line per NAMED composition, coloured by melt SiO2.

    Kept here rather than as a branch inside `plot_crust_composition`, which now handles only the
    (Mg/Si, dIW) axes. Imports from plot_results at call time so the two modules do not import
    each other at module scope.
    """
    import matplotlib.colors as mcolors
    import cmasher as cmr
    import plot_results as pr

    if 'comp_name' not in df.columns:
        df = df.assign(comp_name=df['name'].map(comp_name_of))
    pool = df[df['comp_name'].astype(bool) & df['reverse_weathering']
              & (df['ocean_depth'] == 3000) & (df['land_fraction'] == 0.0)].copy()
    if pool['comp_name'].nunique() < 2:
        print("No named-composition sweep in this data -- skipping.")
        return

    picked = pr._best_operating_point(pool, 'comp_name', 'Named composition')
    if picked is None:
        return
    subset = pr._add_diag_columns(picked[0], output_path)

    values = [c for c in COMP_ORDER if c in set(subset['comp_name'])]
    numeric = [int(c.split('_')[-1]) for c in values]      # SiO2 wt% from the name
    norm = mcolors.Normalize(vmin=42, vmax=53)
    cmap = cmr.gem
    pr._faceted_lines(subset, 'comp_name', values, [cmap(norm(n)) for n in numeric],
                      cmap, norm, 'SiO₂ content (%)', 'lines_named_composition', output_path,
                      split_panels=split_panels, show_markers=show_markers,
                      ticks=numeric, ticklabels=[f'{v}%' for v in numeric])
