def get_T_ocean(T_surface: float, depth: float, a: float=0.0033, b: float=0) -> float:
    """
    Calculates temperature as a function of ocean depth from a linear function.

    Parameters
    ----------
    T_surface : float
        Ocean surface temperature in K.
    depth : float
        Ocean depth in m.
    a : float
        Temperature gradient in K/m, default 0.033.
    b : float
        Temperature offset in K, default 0.

    Returns
    -------
    float
        Ocean temperature in K.
    """

    return T_surface - a * depth + b

def get_T_ocean_KT18(T_surface: float) -> float:
    """
    Calculates temperature as a function of ocean depth from a the parameterization from KT18.

    Parameters
    ----------
    T_surface : float
        Ocean surface temperature in K.

    Returns
    -------
    float
        Ocean temperature in K.
    """

    return 1.02 * T_surface - 16.7