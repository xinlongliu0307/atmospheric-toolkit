import numpy as np


def zw3_index(lons: np.ndarray, z: np.ndarray) -> dict:
    """
    Zonal wave three index by Fourier decomposition along a latitude circle.

    Extracts the wavenumber-3 harmonic of z(lon): amplitude of the
    component A * cos(3 * (lon - phase)), with phase reported in degrees
    in [0, 120), since a wave-3 pattern is periodic every 120 degrees.
    Follows the Fourier formulation of ZW3 (cf. Goyal et al. 2022),
    generalising the fixed-longitude index of Raphael (2004).

    Parameters:
        lons: 1D array of evenly spaced longitudes (degrees, full circle).
        z:    1D array of the field along the circle, same length as lons.

    Returns:
        dict with keys amplitude (same units as z) and phase (degrees).
    """
    if len(lons) != len(z):
        raise ValueError("lons and z must have the same length")
    n = len(z)
    z3 = np.fft.fft(z)[3]
    amplitude = 2.0 * np.abs(z3) / n
    phase = (-np.degrees(np.angle(z3)) / 3.0) % 120.0
    return {"amplitude": float(amplitude), "phase": float(phase)}
