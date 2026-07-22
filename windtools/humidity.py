import math


def dewpoint(temperature_c, relative_humidity_pct):
    # Magnus-Tetens approximation. Returns dewpoint in degrees Celsius.
    b, c = 17.67, 243.5
    gamma = math.log(relative_humidity_pct / 100.0) + (b * temperature_c) / (c + temperature_c)
    return (c * gamma) / (b - gamma)
