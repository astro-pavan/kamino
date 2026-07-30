"""Generate a seafloor-crust mineral composition from planetary parameters.

Pipeline (a single physical knob -> a `crust_composition` dict):

    mantle potential temperature T_p  [+ optional mantle Mg/Si]
      --import_primelt_spreadsheet()--> T_p -> oxide interpolator
      --oxide_composition()-----------> major-element oxide vector (wt%)
      --mineral_composition()---------> CIPW norm -> DB-valid mineral fractions

The oxide anchors are the batch primary magmas from the PRIMELT1 spreadsheet
(Herzberg & O'Hara 2002; Herzberg et al. 2007, G-Cubed) for four natural suites
spanning ambient MORB (Siqueiros, T_p~1365 C) to komatiite (Gorgona, T_p~1606 C).
Interpolating between them gives a continuous, internally consistent primary-melt
composition as a function of T_p; the CIPW norm then converts that melt to the
normative mineral assemblage the weathering model consumes.

The Mg/Si knob is a first-order proxy: a planet whose mantle is more Mg-rich (higher
Mg/Si) melts to a more silica-undersaturated crust. We inject this by scaling the
melt SiO2 relative to the reference (Earth-like) anchors before running the norm.
This captures the silica-saturation *direction* (olivine vs quartz normative), not
the buffered nonlinearity a full melting model would give.
"""

import os
import warnings

import numpy as np
from scipy.interpolate import interp1d

from kamino.mineral_info import MINERAL_MOLAR_MASS

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
PRIMELT_SPREADSHEET = os.path.join(_DATA_DIR, 'ggge967-sup-0002-primelt1.xls')

# Sheets in the PRIMELT1 workbook (note the workbook's 'Icleand' misspelling).
PRIMELT_SHEETS = ['Siqueiros', 'OJP', 'Icleand', 'Gorgona']

# Major-element oxides carried through the pipeline (the ones the CIPW norm uses).
PIPELINE_OXIDES = ['SiO2', 'TiO2', 'Al2O3', 'FeOt', 'MgO', 'CaO', 'Na2O', 'K2O']

# Reference mantle molar Mg/Si the PRIMELT anchors implicitly correspond to. The
# anchors are terrestrial primary magmas, so their source is ~Earth/solar (Mg/Si~1).
MG_SI_REF = 1.0

# Oxide molar masses (g/mol) — convert a weight-percent oxide analysis to cation
# moles for the CIPW norm and the Mg/Si adjustment.
OXIDE_MOLAR_MASS = {
    'SiO2':   60.084,
    'TiO2':   79.866,
    'Al2O3': 101.961,
    'Fe2O3': 159.688,
    'FeO':    71.844,
    'FeOt':   71.844,   # total iron reported as FeO
    'MnO':    70.937,
    'MgO':    40.304,
    'CaO':    56.077,
    'Na2O':   61.979,
    'K2O':    94.196,
    'P2O5':  141.945,
}

# Fe2O3 -> FeO mass-conversion factor (2 * M_FeO / M_Fe2O3): combine split iron into
# total FeO, since the CIPW norm routes all iron to fayalitic olivine.
_FE2O3_TO_FEO = 2 * OXIDE_MOLAR_MASS['FeO'] / OXIDE_MOLAR_MASS['Fe2O3']

# Module-level cache built by import_primelt_spreadsheet().
_INTERPOLATORS: dict[str, interp1d] | None = None
_ANCHORS: dict | None = None


def import_primelt_spreadsheet(path: str = PRIMELT_SPREADSHEET) -> tuple[dict[str, interp1d], dict]:
    """Read the PRIMELT1 batch primary magmas and build the T_p -> oxide interpolator.

    For each suite sheet, the batch-primary output row (the one whose potential-
    temperature column is labelled 'TP') gives a primary-magma oxide analysis and its
    T_p. Iron is combined to total FeO. A per-oxide linear interpolator over T_p is
    built (linear extrapolation outside the anchor range) and cached module-level so
    oxide_composition()/mineral_composition() can be called without re-reading the file.

    Returns (interpolators, anchors); also stores them in module globals.
    """
    import pandas as pd  # local import: pandas/xlrd only needed to (re)build the interpolator

    # Column indices of the PRIMELT output block (0-based); 16 holds T_p.
    COL = {'SiO2': 1, 'TiO2': 2, 'Al2O3': 3, 'Fe2O3': 5, 'FeO': 6,
           'MgO': 8, 'CaO': 9, 'Na2O': 10, 'K2O': 11, 'TP': 16}

    def _num(x: object) -> float:
        try:
            return float(x)
        except (TypeError, ValueError):
            return np.nan

    T_p_list: list[float] = []
    oxide_rows: list[dict[str, float]] = []
    for sheet in PRIMELT_SHEETS:
        df = pd.read_excel(path, sheet, header=None)
        # The batch-primary output row is directly below the header whose T_p column
        # reads exactly 'TP' (as opposed to 'TP_AFM' for the accumulated-fractional fit).
        hdr = next(r for r in range(df.shape[0]) if str(df.iat[r, COL['TP']]).strip() == 'TP')
        v = df.iloc[hdr + 1]
        feot = _num(v[COL['FeO']]) + _FE2O3_TO_FEO * _num(v[COL['Fe2O3']])
        oxide_rows.append({
            'SiO2': _num(v[COL['SiO2']]), 'TiO2': _num(v[COL['TiO2']]),
            'Al2O3': _num(v[COL['Al2O3']]), 'FeOt': feot,
            'MgO': _num(v[COL['MgO']]), 'CaO': _num(v[COL['CaO']]),
            'Na2O': _num(v[COL['Na2O']]), 'K2O': _num(v[COL['K2O']]),
        })
        T_p_list.append(_num(v[COL['TP']]))

    order = np.argsort(T_p_list)  # interp1d needs ascending x
    T_p = np.asarray(T_p_list)[order]
    oxides = {ox: np.asarray([oxide_rows[i][ox] for i in order]) for ox in PIPELINE_OXIDES}

    interpolators = {
        ox: interp1d(T_p, oxides[ox], kind='linear', bounds_error=False, fill_value='extrapolate')
        for ox in PIPELINE_OXIDES
    }

    global _INTERPOLATORS, _ANCHORS
    _INTERPOLATORS = interpolators
    _ANCHORS = {'T_p': T_p, 'oxides': oxides, 'suites': [PRIMELT_SHEETS[i] for i in order]}
    return interpolators, _ANCHORS


def oxide_composition(T_p: float, mg_si_ratio: float = MG_SI_REF) -> dict[str, float]:
    """Interpolated primary-magma oxides (wt%) at potential temperature `T_p` (deg C).

    `mg_si_ratio` is the planet's mantle molar Mg/Si (Earth reference = MG_SI_REF = 1.0).
    A higher Mg/Si source melts to a more silica-poor crust, so the melt SiO2 is scaled
    by MG_SI_REF / mg_si_ratio (first-order silica-saturation proxy) before renormalising
    to 100 wt%. mg_si_ratio = MG_SI_REF leaves the composition unchanged. The overall
    renormalisation does not affect the CIPW norm (which works on ratios) — only the
    SiO2 rescaling changes the resulting mineralogy.
    """
    if _INTERPOLATORS is None:
        import_primelt_spreadsheet()

    T_p_min, T_p_max = float(_ANCHORS['T_p'].min()), float(_ANCHORS['T_p'].max())
    if not (T_p_min <= T_p <= T_p_max):
        warnings.warn(
            f'T_p={T_p} deg C is outside the PRIMELT anchor range '
            f'[{T_p_min:.0f}, {T_p_max:.0f}]; extrapolating.', stacklevel=2)

    oxides = {ox: max(float(_INTERPOLATORS[ox](T_p)), 0.0) for ox in PIPELINE_OXIDES}

    if mg_si_ratio <= 0:
        raise ValueError('mg_si_ratio must be positive')
    oxides['SiO2'] *= MG_SI_REF / mg_si_ratio

    total = sum(oxides.values())
    return {ox: v / total * 100.0 for ox, v in oxides.items()}


def mineral_composition(T_p: float, mg_si_ratio: float = MG_SI_REF, **cipw_kwargs) -> dict[str, float]:
    """Full pipeline: T_p (+ mantle Mg/Si) -> weight-fraction crust mineral dict.

    Interpolates the primary-magma oxides for `T_p`, applies the `mg_si_ratio`
    adjustment, and runs the CIPW norm. Extra keyword arguments (emit_quartz,
    kfeldspar, verbose) pass through to cipw_norm. The returned dict is a valid
    `crust_composition` for get_weathering_flux / get_b_eq.
    """
    oxides = oxide_composition(T_p, mg_si_ratio=mg_si_ratio)
    return cipw_norm(oxides, **cipw_kwargs)


def cipw_norm(oxides: dict[str, float], emit_quartz: bool = False, kfeldspar: bool = False,
              verbose: bool = False) -> dict[str, float]:
    """Convert a bulk-rock major-element analysis into a crust mineral composition.

    Takes major-element oxides in weight-percent and returns a weight-fraction
    dictionary of igneous ("normative") minerals, restricted to the phases that
    exist in the model's PHREEQC database (see hydrothermal_minerals /
    MINERAL_MOLAR_MASS):

        Forsterite, Fayalite, Enstatite, Diopside, Anorthite, Albite
        (+ Wollastonite only as a Ca-overflow fallback, + Quartz if emit_quartz,
         + K-Feldspar only if kfeldspar=True)

    The result is normalised to sum to 1 and can be passed straight in as a
    `crust_composition` (like basalt_49 etc.). Values are WEIGHT fractions, which
    is how get_b_eq consumes a composition (grams of mineral per kgw / molar mass
    -> moles; see chemistry.py `water_rock_ratio` handling).

    This is a compact CIPW-style norm (Cross-Iddings-Pirsson-Washington) that
    allocates cations to phases in the standard order feldspar -> diopside ->
    pyroxene/olivine, using the silica balance to split the ferromagnesian
    remainder between orthopyroxene (Enstatite) and olivine (Forsterite).

    Because the database has no magnetite/ilmenite/apatite/nepheline and no
    Fe-bearing pyroxene, it makes three deliberate simplifications:
      * TiO2 and P2O5 are dropped (P removes a little Ca as apatite first);
      * ALL iron is routed to fayalitic olivine (the only Fe silicate available),
        so pass total iron as FeO/FeOt (FeO and Fe2O3 are summed);
      * Diopside and Enstatite are pure Mg endmembers.
    These match how the existing hardcoded compositions are built. K is likewise not
    a tracked ion (no sink for it), so K-Feldspar is off by default (see kfeldspar).

    Parameters
    ----------
    oxides : dict[str, float]
        Major-element oxides in wt% (missing oxides default to 0). Recognised
        keys: SiO2, TiO2, Al2O3, Fe2O3, FeO, FeOt, MnO, MgO, CaO, Na2O, K2O, P2O5.
    emit_quartz : bool
        If the rock is silica-oversaturated, keep the excess SiO2 as 'Quartz'.
        Default False -> excess silica is dropped with a warning (Quartz doubles
        as an HT precipitation buffer, so it is normally kept out of the crust).
    kfeldspar : bool
        Allocate K to K-Feldspar. Default False -> K is dropped (no K sink in the
        model, and the LT database lacks a K-Feldspar phase). Only set True if a
        K-Feldspar phase has been added to hybrid_ocean.dat.
    verbose : bool
        Print the cation budget and any saturation warnings.

    Returns
    -------
    dict[str, float]  weight-fraction crust composition (sums to 1).
    """
    g = lambda k: float(oxides.get(k, 0.0))

    # oxide wt% -> cation moles (per 100 g; only ratios matter, we normalise later)
    m_Al2O3 = g('Al2O3') / OXIDE_MOLAR_MASS['Al2O3']
    m_Fe2O3 = g('Fe2O3') / OXIDE_MOLAR_MASS['Fe2O3']
    m_FeO   = (g('FeO') + g('FeOt')) / OXIDE_MOLAR_MASS['FeO']
    m_MnO   = g('MnO')   / OXIDE_MOLAR_MASS['MnO']
    m_P2O5  = g('P2O5')  / OXIDE_MOLAR_MASS['P2O5']

    Si = g('SiO2') / OXIDE_MOLAR_MASS['SiO2']
    Al = 2 * m_Al2O3
    Fe = m_FeO + m_MnO + 2 * m_Fe2O3        # all iron (+Mn) as divalent -> fayalite
    Mg = g('MgO') / OXIDE_MOLAR_MASS['MgO']
    Ca = g('CaO') / OXIDE_MOLAR_MASS['CaO'] - (10.0 / 3.0) * m_P2O5   # apatite removes Ca
    Na = 2 * (g('Na2O') / OXIDE_MOLAR_MASS['Na2O'])
    K  = 2 * (g('K2O')  / OXIDE_MOLAR_MASS['K2O'])
    Ca = max(Ca, 0.0)

    warn_msgs: list[str] = []
    moles: dict[str, float] = {}

    def take(mineral: str, amount: float) -> None:
        if amount > 1e-12:
            moles[mineral] = moles.get(mineral, 0.0) + amount

    # 1. K-Feldspar  KAlSi3O8 -- OFF by default. The LT runtime database
    #    (hybrid_ocean.dat) has no K-Feldspar phase, and the model deliberately does
    #    not carry K as a tracked ion (no sink for it), so K is dropped and its Al/Si
    #    pass to plagioclase. MORB is K-poor (~0.1% K2O -> ~1% orthoclase) so this is
    #    a small approximation; warn if handed a genuinely K-rich rock.
    if kfeldspar:
        kf = min(K, Al)
        take('K-Feldspar', kf); Al -= kf; Si -= 3 * kf
    elif g('K2O') > 0.5:
        warn_msgs.append(f'K2O={g("K2O"):.2g} wt% dropped (K not modelled); '
                         f'pass kfeldspar=True if the LT database has a K-Feldspar phase')
    # 2. Albite  NaAlSi3O8
    ab = min(Na, Al)
    take('Albite', ab); Al -= ab; Si -= 3 * ab
    # 3. Anorthite  CaAl2Si2O8  (uses remaining Al)
    an = min(Al / 2.0, Ca)
    take('Anorthite', an); Ca -= an; Al -= 2 * an; Si -= 2 * an
    if Al > 1e-9:
        warn_msgs.append(f'peraluminous: {Al:.3g} mol excess Al (corundum-normative) dropped')
    # 4. Diopside  CaMgSi2O6  (remaining Ca paired with Mg)
    di = min(Ca, Mg)
    take('Diopside', di); Ca -= di; Mg -= di; Si -= 2 * di
    # 4b. Ca overflow (Mg exhausted before Ca) -> Wollastonite fallback
    if Ca > 1e-9:
        take('Wollastonite', Ca); Si -= Ca
        warn_msgs.append(f'Ca-oversaturated: {Ca:.3g} mol Ca to Wollastonite (Mg too low for diopside)')
        Ca = 0.0
    # 5. Ferromagnesian remainder: all Fe -> Fayalite; Mg -> Enstatite/Forsterite by silica
    take('Fayalite', Fe / 2.0); Si -= 0.5 * Fe          # Fe2SiO4 : 0.5 Si per Fe
    if Mg > 1e-12:
        if Si >= Mg:                    # silica-saturated: all orthopyroxene
            take('Enstatite', Mg); Si -= Mg
        elif Si <= 0.5 * Mg:            # silica-deficient: all olivine
            take('Forsterite', Mg / 2.0); Si -= 0.5 * Mg
            warn_msgs.append('silica-deficient (feldspathoid-normative); all Mg to olivine')
        else:                          # mixed opx + olivine set by silica balance
            en = 2 * Si - Mg
            take('Enstatite', en)
            take('Forsterite', (Mg - en) / 2.0)
            Si -= en + 0.5 * (Mg - en)
    # 6. Excess silica -> Quartz (optional)
    if Si > 1e-9:
        if emit_quartz:
            take('Quartz', Si)
        else:
            warn_msgs.append(f'silica-oversaturated: {Si:.3g} mol SiO2 excess dropped '
                             f'(set emit_quartz=True to keep as Quartz)')

    if verbose:
        print(f'cations  Si={g("SiO2")/OXIDE_MOLAR_MASS["SiO2"]:.4f}  Al={2*m_Al2O3:.4f}  '
              f'Fe={Fe:.4f}  Mg={g("MgO")/OXIDE_MOLAR_MASS["MgO"]:.4f}  '
              f'Ca={g("CaO")/OXIDE_MOLAR_MASS["CaO"]:.4f}  Na={Na:.4f}  K={K:.4f}')
        for msg in warn_msgs:
            print('  warning:', msg)
    for msg in warn_msgs:
        warnings.warn(msg, stacklevel=2)

    # mineral moles -> weight fractions
    weights = {m: n * MINERAL_MOLAR_MASS[m] for m, n in moles.items() if n > 0}
    total = sum(weights.values())
    if total <= 0:
        raise ValueError('CIPW norm produced no minerals; check the oxide input')
    return {m: w / total for m, w in weights.items()}
