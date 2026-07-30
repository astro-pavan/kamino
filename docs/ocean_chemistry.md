# Chemistry functions

Geochemical calculations are primarily carried out by PHREEQC. 

$b_{eq} = f_{\mathrm{eq,PHREEQC}}(P, T, P_{\mathrm{CO2}}, b_{in})$

$\mathrm{pH} = f_{\mathrm{pH,PHREEQC}}(P, T, b_{in})$

$P_\mathrm{CO2} = f_{\mathrm{pH,PHREEQC}}(P, T, b_{in})$

$k = f_{\mathrm{kinetics}}(T, \mathrm{pH})$

$k$ is always calculated with $b_{eq}$ and uses the pH calculated from the equilbrium calculation.

The precipitation rate is calculated as follows:

$F_{\mathrm{prec}} = \frac{b_{eq} - b_{in}}{\tau_{prec}},$

where $F_{\mathrm{prec}}$ is the precipiation rate per unit mass of water, $b_{eq}$ is the equilbrium concentration when all supersaturated minerals have precipiated, $b_{in}$ is the input concentrations, and $\tau_{prec}$ is the precipiation timescale and controls the strength of the precipitation.

# Low Temperature Water rock interactions

The low temperature weathering law is calculated as follows. The water enters the low temperature pore space at seawater concentrations. The primary baslat minerals dissolve, modifyoing the pore space concentrations. In this pore space water, secondary minerals precipitate out, such as clays and carbonate minerals. The rest of the water then re enters the ocean.

## Primary Dissolution

The water rock interactions are based on the weathering formula from Maher & Chamberlain (2014), which models a weathering flux in both the kinetic and thermodynamic limits. 

The water rock reaction starts by modelling the dissolution of the primary minerals, given by the following formula,

$F_p=\frac{A_r(b_{eq} - b_{in})}{\frac{b_{eq}}{k_p} + \frac{A_r}{J}},$

where $F_p$ is the flux per unit area of seafloor from the dissolution, $A_r$ is the reactive area per unit area of seafloor, $b_{in}$ is the concentration from the incoming seawater, $b_{eq}$ is the equilbrium conecntration of the reaction, $k_p$ is the reaction rate per unit area of seafloor and $J$ is the hydrothermal flux per unit area seafloor. $b_{eq}$ and $k_p$ are both calculated by PHREEQC.

## Secondary Mineral Formation

After the primary dissolution, the pore space water concentration ($b_p$) is calculated as follows:

$b_p = b_in + (b_{eq} - b_{in})\left(1 - e^{-\mathrm{Da}}\right),$

where $\mathrm{Da}$ is the Damkohler coefficient of the weathering reaction and is given by,

$\mathrm{Da} = \frac{k_p A_r}{J b_{eq}}.$

The secondary precipitation is calculated from the concentration of the pore water. 

## Seafloor Reactive Area

The ion fluxes from 

# Precipitation

# Other Processes

Cl subduction

# Summary by Ion

| Element | Source | Sink |
| :--- | :--- | :--- |
| **Alk** | weathering | carbonate precipitation |
| **C** | outgassing | carbonate precipitation, (+ biogenic burial) |
| **Ca** | weathering, hydrothermal systems | carbonate precipitation |
| **Mg** | weathering | hydrothermal system, reverse weathering |
| **Fe** | weathering | clay formation |
| **Na** | weathering | reverse weathering, hydrothermal system, carbonate precipitation |
| **Si** | weathering | silica precipitation, reverse weathering, (+ biogenic burial) |
| **Al** | weathering | clay formation |
| **Cl** | outgassing | subduction |