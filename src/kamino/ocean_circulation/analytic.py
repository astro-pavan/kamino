def get_T(T_surface: float, depth: float, a: float=0.033, b: float=0) -> float:
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