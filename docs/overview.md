# Model Overview #

The planet model solves the carbon cycle equations based on those in [Krissansen-Totton & Catling (2017)](https://www.nature.com/articles/ncomms15423) without continental weathering. They are:

$$ \frac{dC_o}{dt} = \frac{1}{M_o} \left(-J(C_o - C_p) - F_{\mathrm{precip,o}} + F_{\mathrm{outgas}} \right)$$

$$ \frac{dC_p}{dt} = \frac{1}{M_p} \left(+J(C_o - C_p) - F_{\mathrm{precip,p}} \right)$$

$$ \frac{dA_o}{dt} = \frac{1}{M_o} \left(-J(A_o - A_p) - 2F_{\mathrm{precip,o}} \right)$$

$$ \frac{dA_p}{dt} = \frac{1}{M_p} \left(+J(A_o - A_p) - 2F_{\mathrm{precip,p}} + F_{\mathrm{diss}}\right)$$

$$ \frac{db_{\mathrm{Ca},o}}{dt} = \frac{1}{M_o} \left(-J(b_{\mathrm{Ca},o} - b_{\mathrm{Ca},p}) - F_{\mathrm{precip,o}} \right)$$

$$ \frac{db_{\mathrm{Ca},p}}{dt} = \frac{1}{M_p} \left(+J(b_{\mathrm{Ca},o} - b_{\mathrm{Ca},p}) - F_{\mathrm{precip,p}} + \frac{1}{2}F_{\mathrm{diss}}\right)$$

where $C_o$ and $C_p$ is the dissolved inorganic carbon molality in the ocean and seafloor pore space respectively; $A_o$ and $A_p$ is the alkalinity in the ocean and pore space; $b_{\mathrm{Ca},o}$ and $b_{\mathrm{Ca},p}$ is the calcium molality in the ocean and pore space; $M_o$ and $M_p$ is the mass of the water in the ocean and pore space; J is the flux of water through the pore space; $F_{\mathrm{precip,o}}$ and $F_{\mathrm{precip,p}}$ is the precipitation rate of calcite in the ocean and pore space; $F_{\mathrm{diss}}$ is the alkalinity production rate from basalt dissolution (NOTE: for every two units of alkalinity produced by basalt dissoltion one divalent cation is produced); and $F_{outgas}$ is the volcanic outgassing rate. 