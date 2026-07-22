from math import sqrt

def wind_speed(u: float, v: float) -> float:
    """
    Returns the scalar wind speed, the square root of u squared plus v squared.
    """
    return sqrt(u**2 + v**2)
