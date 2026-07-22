from math import sqrt

from math import sqrt, atan2, degrees

def wind_speed(u: float, v: float) -> float:
    """
    Returns the scalar wind speed, the square root of u squared plus v squared.
    """
    return sqrt(u**2 + v**2)

def wind_direction(u: float, v: float) -> float:
    """
    Returns the meteorological direction in degrees the wind blows FROM,
    computed as (270 - degrees(atan2(v, u))) modulo 360.
    """
    return (270 - degrees(atan2(v, u))) % 360
