# Model Overview

## The Carbon Cycle

The planet model solves the carbon cycle equations based on those in [Krissansen-Totton & Catling (2017)](https://www.nature.com/articles/ncomms15423) without continental weathering. They are:

$$ \frac{dC_o}{dt} = \frac{1}{M_o} \left(-J(C_o - C_p) - F_{\mathrm{precip,o}} + F_{\mathrm{outgas}} \right)$$

$$ \frac{dC_p}{dt} = \frac{1}{M_p} \left(+J(C_o - C_p) - F_{\mathrm{precip,p}} \right)$$

$$ \frac{dA_o}{dt} = \frac{1}{M_o} \left(-J(A_o - A_p) - 2F_{\mathrm{precip,o}} \right)$$

$$ \frac{dA_p}{dt} = \frac{1}{M_p} \left(+J(A_o - A_p) - 2F_{\mathrm{precip,p}} + F_{\mathrm{diss}}\right)$$

$$ \frac{db_{\mathrm{Ca},o}}{dt} = \frac{1}{M_o} \left(-J(b_{\mathrm{Ca},o} - b_{\mathrm{Ca},p}) - F_{\mathrm{precip,o}} \right)$$

$$ \frac{db_{\mathrm{Ca},p}}{dt} = \frac{1}{M_p} \left(+J(b_{\mathrm{Ca},o} - b_{\mathrm{Ca},p}) - F_{\mathrm{precip,p}} + \frac{1}{2}F_{\mathrm{diss}}\right)$$

where $C_o$ and $C_p$ is the dissolved inorganic carbon molality in the ocean and seafloor pore space respectively; $A_o$ and $A_p$ is the alkalinity in the ocean and pore space; $b_{\mathrm{Ca},o}$ and $b_{\mathrm{Ca},p}$ is the calcium molality in the ocean and pore space; $M_o$ and $M_p$ is the mass of the water in the ocean and pore space; and J is the flux of water through the pore space.

The fluxes in the carbon cycle equations are given in the table below (with $P$ as pressure, $T$ as temperature, $x_{CO2}$ as CO2 volume mixing ratio).

| Flux | Description | Input parameters | Calculated with |
| :--: | :---------: | :--------------: | :-------------: |
| $F_{\mathrm{precip},o}$ | Ocean precipitation rate of calcite | $P_{\mathrm{seafloor}}$, $T_{\mathrm{seafloor}}$, $A_o$, $C_o$, $b_{\mathrm{Ca},o}$ | ocean_chemistry |
| $F_{\mathrm{precip},p}$ | Pore space precipitation rate of calcite | $P_{\mathrm{seafloor}}$, $T_{\mathrm{pore}}$, $A_p$, $C_p$, $b_{\mathrm{Ca},p}$ | ocean_chemistry |
| $F_{\mathrm{diss}}$ | Alkalinity production rate from basalt dissolution | $P_{\mathrm{seafloor}}$, $T_{\mathrm{pore}}$, $x_{CO2}$ | seafloor_weathering |
| $F_{\mathrm{outgas}}$ | Volcanic outgassing | Constant | ... | 

NOTE: for every two units of alkalinity produced by basalt dissoltion one divalent cation is produced.

To improve solver stability, the preciptation fluxes have a maximum value of 

$$ F_{\mathrm{precip,max},o} = \frac{C_o M_o}{\tau}, $$ $$ F_{\mathrm{precip,max},p} = \frac{C_p M_p}{\tau}, $$

where $\tau$ is a precipation timescale (set at 0.001 $\mathrm{yr}$). 

The temperature structure of the ocean is set as a linear gradient in the absence of a better parameterisation. The relationship between seafloor and pore space temperature is from Krissansen-Totton & Catling (2017). The seafloor and pore space pressure and temperature are thus given by 

$$ P_{\mathrm{seafloor}} = P_{\mathrm{surface}} + \rho g d, $$ $$ T_{\mathrm{seafloor}} = T_{\mathrm{surface}} - ad, $$ $$ T_{\mathrm{pore}} = T_{\mathrm{seafloor}} + 9\,\mathrm{K}, $$

where $\rho$ is the density of water, $g$ is the planet's surface gravity, $d$ is the ocean depth and, $a$ is an ocean temperature gradient (NOTE: currently set at 0.033 $\mathrm{K}/\mathrm{m}$ for a 3 km deep ocean).

## Solving Climate

To close the carbon cycle equations, the climate model from speedy_climate is needed to relate $x_{CO2}$ and $T_{\mathrm{surface}}$

The climate emulator gives the surface temperature as a function instellation, surface pressure, CO2 mixing ratio, and H2O mixing ratio. Instellation and surface pressure is held constant for a planet model run. The H2O mixing ratio can be set with, 

$$x_{H2O}=RH \cdot \frac{e_s(T_{\mathrm{surface}})}{P_\mathrm{surface}},$$

where $RH$ is the relative humidity (set at 0.77) and $e_s(T)$ is the water saturation vapour pressure. This models increased evaporation and increased water vapour in the atmosphere with increasing temperatures. 

The CO2 mixing ratio can be calculated with the partial pressure of CO2 ($x_{CO2}=P_{CO2}/P_{\mathrm{surface}}$) which is calculated in ocean_chemistry as a function of $T_{\mathrm{surface}}$, $P_{\mathrm{surface}}$, $A_o$, $C_o$, $b_{\mathrm{Ca},o}$. This closes the system of equations and allows it to be solved.