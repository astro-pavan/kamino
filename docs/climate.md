# Climate

The climate calculations are done by a modified version of [``HELIOS``](https://github.com/astro-pavan/HELIOS), a GPU accelerated 1D radiative convective climate solver. ``HELIOS`` runs very fast and is a good choice for parameter sweeps but was designed for general exoplanets, so needs some Earth-like atmospheric processes to be added to the code to accurately model Earth-like planets with water vapour in their atmospheres. To speed up the overall model further, an emulator was trained on a selection of 1000 HELIOS climate runs.

## HELIOS Modifications

### Rainout

As water vapour rises in the atmosphere, it eventually reaches a point where the temperature is low enough that it condenses and rains back to the surface. Hence, rather than setting a constant H2O volume mixing ratio through out the atmosphere, which leads to a much higher greenhouse effect than that on Earth, a rainout scheme (described below) is applied. 

If the set partial pressure from the H2O is higher than the H2O saturation vapour pressure multiplied by a constant value of relative humidity (set at 0.77 NEED CITATION), then the H2O partial pressure (and hence the H20 mixing ratio) is set at that maximum. This can be expressed as:

$$ x_{H2O} = \mathrm{min}\left(x_{H2O,\mathrm{surf}}, \, RH \cdot e_{s}(T)\right), $$

where $x_{H2O,\mathrm{surf}}$ is the surface H2O volume mixing ratio, $RH$ is the relative humidity, and $e_{s}(T)$ is the saturation vapour pressure as a function of temperature, which is calculated with the August-Roche-Magnus formula. This models the H2O condensation.

If the H2O volume mixing ratio in a layer is greater than the layer below it, its H2O volume mixing ratio is lowered to the value from the layer below. This models the fact water vapour cannot rise above the height it begins to rainout.

### Moist Convection

``HELIOS`` models convection by setting any convectively unstable layer (where temperature profile is above the adiabat), to the dry adiabat. 

As air containing water vapour rises adiabatically, the water vapour eventually condenses, releasing latent heat. This condensation leads to the moist adiabat being different to the dry adiabat. The moist adiabat can be calculated with the Saturated Adiabatic Lapse Rate:

$$ \Gamma_{\mathrm{w}} = g \frac{\left( 1 + \frac{H_{\mathrm{v}} r}{R_{\mathrm{sd}} T} \right)}{\left( c_{\mathrm{pd}} + \frac{H_{\mathrm{v}}^2 r}{R_{\mathrm{sw}} T^2} \right)}, $$

where $g$ is surface gravity, $H_{\mathrm{v}}$ is the heat of vaporization of water (2501000 J/kg), $R_{\mathrm{sd}}$ is the specific gas constant of dry air (287 J/kg·K), $T$ is the temperature, $R_{\mathrm{sw}}$ is the specific gas constant of water vapour (461.5 J/kg·K), $c_{\mathrm{pd}}$ is the specific heat of dry air at constant pressure (1003.5 J/kg·K) and $r$ is the the mixing ratio of the mass of water vapour to the mass of dry air. $r$ is given by:

$$ r = \frac{\frac{R_{\mathrm{sd}}}{R_{\mathrm{sw}}}e_s(T)}{p-e_s(T)}, $$ 

where $p$ is the pressure of the saturated air. 

This moist adiabat replaces the dry adiabat used by ``HELIOS``.

### Partial Clouds

Clouds play an important role in keeping Earth-like planets cool. However complete cloud cover is unlikely (NOTE: Is this generally the case or just applicable to Earth?), around 60% of Earth is covered in cloud. Complete cloud cover also cools the the planet signifcantly to signifiacntly colder temperatures than Earth's for Earth-like conditions. Hence partial clouds are an important part of modelling a an Earth-like planet.

Partial clouds are modelled by having two raditaive convective columns. The radiation calculation is performed for each column, one cloudy and one cloud free. The fluxes are then combined together:

$$ F_{\mathrm{avg},i} = f F_{\mathrm{cloudy},i} + (1 - f) F_{\mathrm{clear},i} $$

where $F_{\mathrm{avg},i}$ is the average flux in the layer, $F_{\mathrm{cloudy},i}$ is flux in the layer from the cloudy column, $F_{\mathrm{clear},i}$ is flux in the layer from the clear column and $f$ is the fraction of the planet covered by cloud.

### Cloud Destruction

High sea surface temperatures lead to the destruction of stratocumulus clouds (Schneider et al. 2019). This begins at 298 K and Stratoculmulus clouds are completely destroyed at 305 K. A simple cloud destruction scheme is implemented, where a cloud destruction factor $d_{\mathrm{cloud}}$ is calculated as:

$$
d_{\mathrm{cloud}}(T) = 
\begin{cases}
0, & T < 298\,\mathrm{K} \\
\dfrac{T - 298}{305 - 298}, & 298\,\mathrm{K} \leq T < 305\,\mathrm{K} \\
1, & T \geq 305\,\mathrm{K}
\end{cases}
$$

where $T$ is the surface temperature. The cloud fraction is multiplied by $(1 - d_{\mathrm{cloud}})$.

## Emulating HELIOS

Each climate run with ``HELIOS`` takes approximately 10 seconds with a desktop computer with a GPU. In solving the carbon cycle equations, the climate state needs to be solved on each iteration, making each run potentially very time consuming. To speed up the climate calculation, an emulator is used to calculate the surface temperature, the value needed for the carbon cycle calulations. The emulator is trained with Gaussian processes. Gaussian processes are chosen over a simple interpolator as there are large amount of input parameters and the climate state has some non-linear behaviours. Gaussian processes are chosen over neural networks as they can be trained on approximately 1000 data points, which takes a few hours to generate with a desktop computer.

To reduce the dimensionality of the problem, only some parameters are used as input parameters. Other parameters are fixed with Earth-like values and others are set for the type of planet that is being modelled. 

Two sizes of planets are tested Earth-like planets (with masses and radii similar to Earth) and super-Earths with masses around XX and Radii around XX. Two modes of planet rotation are modelled, normal Earth-like rapid rotation (where we can consider the heat from the star averaged out over the planet's surface) and tidally locked planets (with one face constantly facing the host star). Tidally locked planets typically orbit red dwarf stars with much lower stellar effective temperatures.

The input parameters for the climate emulator are given in the table below.

| Parameter | Units | Definition | Value |
| :-------: | :---: | :--------: | :-----------: |
| $S$ | $\mathrm{W}\mathrm{m}^{2}$ | Instellation | Varied |
| $P_{\mathrm{surface}}$ | $Pa$ | Surface pressure | Varied | 
| $x_{CO2}$ | ... | CO2 volume mixing ratio | Varied |
| $x_{H2O}$ | ... | H2O volume mixing ratio | Varied |
| $M$ | $\mathrm{kg}$ | Mass of the planet | $5.97 \times 10^{24}$ (Earth-like), $? \times 10^{24}$ (Super-Earth) |
| $R$ | $\mathrm{m}$ | Radius of planet | $6.37 \times 10^{6}$ (Earth-like), $? \times 10^{6}$, (Super-Earth)|
| $T^*_{\mathrm{eff}}$ | $\mathrm{K}$ |Stellar effective temperture | $5780\,\mathrm{K}$ (rapid-rotator), $3120\,\mathrm{K}$ (tidally-locked) |
| $f_{\mathrm{circ}}$ | ... | Circulation factor | 0.25 (rapid-rotator), 0.666 (tidally-locked) |
| $a$ | ... | Surface albedo | 0.05 |
| $f$ | ... | Cloud cover fraction | 0.55 |
| $RH$ | ... | Relative humidity | 0.77 |
