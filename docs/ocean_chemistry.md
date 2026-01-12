# Ocean Chemistry

## Preciptation

When calcite (CaCO3) precipitates out of the ocean, it removes dissolved carbon and twice as much alkalinity (as a carbonate ion is removed which is worth twice much alkalinity as bicarbonate). This calcite settles at the ocean floor or in the seafloor pore space and acts as a carbon sink.

The precipitation rate of calcite in the ocean is calculated using [``PHREEQC``](https://www.usgs.gov/software/phreeqc-version-3). ``PHREEQC`` uses a modofed version of the Kinec_v3 database (containing a reduced number of phases) to calculate the ocean chemistry. ``PHREEQC`` can calculate the concentrations of all phases in the seawater with just temperature, pressure, dissolved inorganic carbon (DIC), alkalinity and the, molarity of the divalent cations (in the simplest case, only Calcium). These concentrations are used with the rate laws in the Kinec_v3 database to calculate the precipitation rate.

To get the preciptation rate, a very small concentration starts in the solution ($10^{-10}\, \mathrm{mol}/\mathrm{kgw}$). The specific surface area is set at 0.01 $\mathrm{m} \mathrm{kg}^{-1}$ (NEED CITATION FOR THESE). The kitectics calulation is run for 1 second to get an instantaneous reaction rate. The reaction rate per unit mass of water from ``PHREEQC`` output is multiplied by the mass of the ocean to get the overal reaction rate (NOTE: this is a simplification and needs to be examined). 

The pressure and temperature of the precipiation reaction in the ocean is set at the seafloor pressure and temperature. 

By setting the ocean precipitaion rate with the seafloor conditions, if the bottom of the ocean is dissolving calcite instead of precipitating it, the precipiation rate is set to zero. This happens if the calcite compenstation depth is above the seafloor and hence no carbon can leave the ocean, stopping the carbon cycle until the ocean is saturated enough for calcite to precipitate again.

## Dissolution of CO2

``PHREEQC`` is also used to calculate the amount of CO2 in the atmosphere. ``PHREEQC`` calculates the saturation index of CO2 from the surface pressure, temperature, DIC, alkalinity and divalent cation molalities. The saturation index can be used to calculate CO2 partial pressure:

$$ P_{CO2} = P_{\mathrm{surface}} 10^{\mathrm{SI}_{CO2}} $$