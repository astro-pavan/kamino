# Seafloor Weathering #

The planet model uses the silicate weathering model from [Hakim et al. (2021)](https://iopscience.iop.org/article/10.3847/PSJ/abe1b8/meta) which is based off the weathering model from Maher & Chamberlain (2014). Unlike the seafloor weathering model from [Krissansen-Totton & Catling (2017)](https://www.nature.com/articles/ncomms15423), the H21 model takes into account the thermodynamic limit, where the ions in the water reach equilbrium with the weathering rocks in a short enouygh time that the water becomes saturated and the supply limit, where the rate is limited by the amount of reactants. The weathering model is given by

$$ w = q \frac{C_{\mathrm{eq}}}{1 + \frac{q}{D_w}} $$

Where $w$ is the weathering rate per unit area, q is the fluid flow rate (or runoff), $C_{\mathrm{eq}}$ is the equilbrium concentration and $D_w$ is the Damköhler coefficient, which gives a net reaction rate. 

$$ D_w = \frac{L(1-\phi)\rho X_r A_{\mathrm{sp}}}{C_{\mathrm{eq}} \left(k_{\mathrm{eff}}^{-1} + mA_{\mathrm{sp}}t_s \right)} ,$$