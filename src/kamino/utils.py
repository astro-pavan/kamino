import numpy as np

from kamino.constants import *

# Raise an exception on invalid floating point operations
np.seterr(invalid='raise')

def smooth_max(a, b, epsilon=1e-24):
    """
    Returns a smooth approximation of max(a, b).
    epsilon controls the 'roundness' of the corner. 
    Smaller epsilon = sharper corner.
    """
    return 0.5 * ((a + b) + np.sqrt((a - b)**2 + epsilon))

def smooth_min(a, b, epsilon=1e-24):
    """
    Returns a smooth approximation of min(a, b).
    """
    return 0.5 * ((a + b) - np.sqrt((a - b)**2 + epsilon))

def smooth_ReLU(x, epsilon=1e-20):
    return 0.5 * (x + np.sqrt(x**2 + epsilon))

def smooth_heaviside(x, k=1e10):
    return 0.5 *((k*x) / (np.sqrt(1 + (k * x) ** 2)) + 1)

def august_roche_magnus_formula(T: float) -> float:
    """
    Implementation of August-Roche-Magnus formula to calculate H2O partial pressure.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        H2O partial pressure in Pa
    """

    T_celsius = T + ABSOLUTE_ZERO
    return 610.94 * np.exp((17.625 * T_celsius)/(T_celsius + 243.04))