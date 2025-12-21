M_EARTH = 5.972e24 # kg
R_EARTH = 6371000 # m

R_JUPITER = 69911000 # m
R_SUN = 695700000 # m
AU = 1.495978707e11 # m

YR = 3.16e7 # s

STEFAN_BOLTZMANN = 5.670374419e-8 # W m^-2 K~^-4

G = 6.674e-11 # m^3 kg^-1 s^-2

IDEAL_GAS_CONSTANT = 8.314 # J K^-1 mol^-1
EARTH_ATM = 101325 # Pa
SOLAR_CONSTANT = 1360 # W / m^2
ABSOLUTE_ZERO = -273.15 # K

from typing import Dict

SPECTRAL_TYPE_DATA: Dict[str, Dict[str, float]] = {
    "O2":  {"Radius": 12.0, "Mass": 100.0, "Luminosity": 800000, "Temperature": 50000},
    "O6":  {"Radius": 9.8,  "Mass": 35.0,  "Luminosity": 180000, "Temperature": 38000},
    "B0":  {"Radius": 7.4,  "Mass": 18.0,  "Luminosity": 20000,  "Temperature": 30000},
    "B5":  {"Radius": 3.8,  "Mass": 6.5,   "Luminosity": 800,    "Temperature": 16400},
    "A0":  {"Radius": 2.5,  "Mass": 3.2,   "Luminosity": 80,     "Temperature": 10800},
    "A5":  {"Radius": 1.7,  "Mass": 2.1,   "Luminosity": 20,     "Temperature": 8620},
    "F0":  {"Radius": 1.3,  "Mass": 1.7,   "Luminosity": 6,      "Temperature": 7240},
    "F5":  {"Radius": 1.2,  "Mass": 1.3,   "Luminosity": 2.5,    "Temperature": 6540},
    "G0":  {"Radius": 1.05, "Mass": 1.1,   "Luminosity": 1.26,   "Temperature": 5920},
    "G2":  {"Radius": 1.0,  "Mass": 1.0,   "Luminosity": 1.0,    "Temperature": 5780},
    "G5":  {"Radius": 0.93, "Mass": 0.93,  "Luminosity": 0.79,   "Temperature": 5610},
    "K0":  {"Radius": 0.85, "Mass": 0.78,  "Luminosity": 0.4,    "Temperature": 5240},
    "K5":  {"Radius": 0.74, "Mass": 0.69,  "Luminosity": 0.16,   "Temperature": 4410},
    "M0":  {"Radius": 0.51, "Mass": 0.6,   "Luminosity": 0.072,  "Temperature": 3800},
    "M5":  {"Radius": 0.18, "Mass": 0.15,  "Luminosity": 0.0027, "Temperature": 3120},
    "M8":  {"Radius": 0.11, "Mass": 0.08,  "Luminosity": 0.0004, "Temperature": 2650},
    "L1":  {"Radius": 0.09, "Mass": 0.07,  "Luminosity": 0.00017,"Temperature": 2200}
} 
# Zombeck, Martin V. (1990). Handbook of Space Astronomy and Astrophysics (2nd ed.). Cambridge University Press. ISBN 978-0-521-34787-7. 
# In Solar units (Mass, Radius and Luminosity) and K (Temperature)

GAS_SPECIES_MMW = {
    'H2O' : 0.018,
    'CO2' : 0.044,
    'N2' : 0.028,
    'H2' : 0.002,
    'O2' : 0.032
} # kg / mol

MINERAL_MMW = {
    'Albite': 0.263,
    'Anorthite': 0.278,
    'Diopside':  0.217,
    'Enstatite': 0.201,
    'K-Feldspar': 0.278,
    'Forsterite': 0.141,
    'Fayalite': 0.204,
    'Magnetite': 0.232,
    'Quartz': 0.060,
    'Kaolinite': 0.258,
    'Goethite': 0.089,
    'Calcite': 0.100,
    'Aragonite': 0.100
} # kg / mol

MINERAL_SPECIFIC_GRAVITY = {
    'Albite': 2.62,
    'Anorthite': 2.73,
    'Diopside':  3.40,
    'Enstatite': 3.20,
    'K-Feldspar': 2.56,
    'Forsterite': 3.27,
    'Fayalite': 4.39,
    'Magnetite': 5.15,
    'Quartz': 2.62,
    'Kaolinite': 2.60,
    'Goethite': 3.80,
    'Calcite': 2.71,
    'Aragonite': 2.93
} # DIMENSIONLESS

WATER_MOLALITY = 55.510 # mol / kg