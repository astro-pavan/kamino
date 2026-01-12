# KAMINO

#### Coupled ocean-atmosphere model for exoplanets (new version)

![Picture of ocean exoplanet](Oceanworld.jpeg)

### About

KAMINO is a model for simulating the climate stability of ocean exoplanets, defined as planets that are entirely covered in oceans. It has three subpackages:

- **speedy_climate**: Uses [``HELIOS``](https://github.com/exoclime/HELIOS), a GPU accelerated 1D radiative convective climate model to produce a sample of climate runs across the relevant parameter space. It then trains an emulator using Gaussian processes on this sample to greatly speed up climate calculations. The repository contains two sets of climate runs so you don't have to perform the parameter sweep yourself.

- **seafloor_weathering**: Uses the weathering model from [Hakim et al. (2021)](https://iopscience.iop.org/article/10.3847/PSJ/abe1b8/meta) to calculate the production of alkalinity from seafloor weathering.

- **ocean_chemistry**: Uses [``PHREEQC``](https://www.usgs.gov/software/phreeqc-version-3), an aqueous geochemical solver, to calculate the ocean chemistry including the precipitation rate of the dissolved cations.

### Installation instructions

```bash
git clone --recursive https://github.com/astro-pavan/kamino.git
cd kamino
pip install -e .
```
