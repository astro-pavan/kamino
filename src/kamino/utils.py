import numpy as np

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