# Tectonic Activity

XXXX

# Crust Composition

The oceanic crust composition is calculated from two parameters that stellar spectroscopy can constrain: the mantle molar Mg/Si ratio, and the oxygen fugacity of core formation expressed relative to the iron-wustite buffer (dIW). Mg/Si controls the olivine-to-pyroxene balance, while dIW sets how much iron the mantle retained as FeO rather than losing to the core as metal, and hence the iron content of the crust. Earth corresponds to Mg/Si = 1.25 and dIW = -2, which reproduces a bulk silicate Earth FeO content of 8.05 wt%.

A pyrolite mantle (McDonough & Sun 1995) is modified on these two axes and melted by isentropic decompression from 3.0 to 1.0 GPa using MAGEMin (Riel et al. 2022) with the Holland, Green & Powell (2018) igneous dataset. The melt fraction is held at 20%, the degree at which clinopyroxene leaves the melting assemblage for Earth's mantle (Guimond et al. 2024), and the mantle potential temperature is solved for each composition rather than prescribed. The resulting primary melt compositions are tabulated over a grid and interpolated at runtime, then converted to a normative mineral assemblage using the CIPW norm.

**See `docs/crust_composition.md` for the full methods reference**, including the dIW-to-FeO mapping, the desilication cascade, validation against Guimond et al. (2024), and the known limitations.

# Outgassing

Carbon and chlorine enter into the ocean via outgassing. This is assumed to be at a constant rate. The cholrine is assumed to enter the ocean as HCl and thus adds a unit of alkalinity per unit of cholrine. 