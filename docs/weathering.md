# Seafloor Weathering

The planet model uses the silicate weathering model from [Hakim et al. (2021)](https://iopscience.iop.org/article/10.3847/PSJ/abe1b8/meta) which is based off the weathering model from Maher & Chamberlain (2014). Unlike the seafloor weathering model from [Krissansen-Totton & Catling (2017)](https://www.nature.com/articles/ncomms15423), the H21 model takes into account the thermodynamic limit, where the ions in the water reach equilbrium with the weathering rocks in a short enouygh time that the water becomes saturated and the supply limit, where the rate is limited by the amount of reactants. The weathering model is given by

$$ w = q \frac{C_{\mathrm{eq}}}{1 + \frac{q}{D_w}} $$

Where $w$ is the weathering rate per unit area, q is the fluid flow rate (or runoff), $C_{\mathrm{eq}}$ is the equilbrium concentration and $D_w$ is the Damköhler coefficient, which gives a net reaction rate. $D_w$ is calculated by

$$ D_w = \frac{L (1-\phi) \rho X_r A_{\mathrm{sp}}}{C_{\mathrm{eq}} \left(k_{\mathrm{eff}}^{-1} + mA_{\mathrm{sp}}t_s \right)} ,$$

where $L$ is the flow path length, $\phi$ is the rock porosity, $\rho$ is the weathered rock density, $X_r$ is the fraction of fresh minerals in the rock, $A_{\mathrm{sp}}$ is the specific surface area of the weathered rock, $k_{\mathrm{eff}}$ is the reaction rate of the rock dissolution, $m$ is the molar mass of the rock, and $t_s$ is the age of the rock.

$ C_{\mathrm{eq}} $ is the alkalinity concentration in the maximum weathering case, where the seawater becomes saturated with the weathering reaction products. Alkalinity is the sum of the bicarbonate and double the carbonate concentrations (carbonate is doubled as it has a -2 charge).  $ C_{\mathrm{eq}} $ is calculated by solving for the activities of all of the species in the dissoultion reactions of the basalt minerals (Wollastonite, Enstatite, Ferrosilite, Anorthite and Albite) and the water-carbonate system at different values of pressure, temperature and atmospheric CO2 mixing ratio.

$ k_{\mathrm{eff}} $ is the reaction rate of the rock dissolution reaction and depends on temperature and pH. The pH value is set at the pH of the saturated state (calculated at the same time as $ C_{\mathrm{eq}} $) As basalt is made up of multiple minerals, the lowest reaction rate is set as the overal reaction rate of the rock dissolution rate.

$ C_{\mathrm{eq}} $ and $ k_{\mathrm{eff}} $ are calculated with ``CHILI`` (CHemical weatherIng model based on LIthology) (Hakim et al. 2021).

### Weathering Parameters 

| Parameter | Units | Definition | Default Value |
| :-------: | :---: | :--------: | :-----------: |
| $q$ | $ \mathrm{m} \, \mathrm{yr}^{-1} $ | Runoff | 0.05 |
| $L$ | meters (m) | Flow path length | 100 |
| $\phi$ | ... | Porosity | 0.1 |
| $ \rho $ | $ \mathrm{kg} \mathrm{m}^{-3} $ | Rock density | 2700 | 
| $ X_r $ | ... |Fresh mineral fraction | 1 |
| $ A_{\mathrm{sp}} $ | $\mathrm{m} \, \mathrm{kg}^{-1} $ | Specific surface area | 100 |
| $m$ | $ \mathrm{kg} \, \mathrm{mol}^{-1} $ | Mean molar mass of rock | 0.216 |
| $ t_s $ | $ \mathrm{yr} $ | Rock age | $ 50 \times 10^{6} $ |

NOTE: These values were taken from Graham and Pierrehumbert (2020) and Hakim et al. (2021)
