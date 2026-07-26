import numpy as np


def find_asl(lons: np.ndarray, lats: np.ndarray, slp: np.ndarray) -> dict:
    """
    Locate the Amundsen Sea Low in a 2D sea-level-pressure field.

    Searches only the conventional ASL sector, 170-298 degrees east and
    60-80 degrees south, and returns the minimum-pressure grid point.

    Parameters:
        lons: 1D array of longitudes (degrees east).
        lats: 1D array of latitudes (degrees north, negative south).
        slp:  2D array of sea-level pressure, shape (len(lats), len(lons)).

    Returns:
        dict with keys lon, lat, central_pressure, and
        relative_central_pressure (central pressure minus the sector mean).
    """
    lon_mask = (lons >= 170.0) & (lons <= 298.0)
    lat_mask = (lats >= -80.0) & (lats <= -60.0)
    sector = slp[np.ix_(lat_mask, lon_mask)]
    sector_lons = lons[lon_mask]
    sector_lats = lats[lat_mask]

    j, i = np.unravel_index(np.argmin(sector), sector.shape)
    central_pressure = float(sector[j, i])
    return {
        "lon": float(sector_lons[i]),
        "lat": float(sector_lats[j]),
        "central_pressure": central_pressure,
        "relative_central_pressure": central_pressure - float(np.mean(sector)),
    }


def asl_timeseries(lons: np.ndarray, lats: np.ndarray, slp_stack: np.ndarray) -> dict:
    """
    Apply find_asl to each time step in a 3D stack of sea-level-pressure fields.

    Parameters:
        lons: 1D array of longitudes (degrees east).
        lats: 1D array of latitudes (degrees north, negative south).
        slp_stack: 3D array of pressure, shape (time, len(lats), len(lons)).

    Returns:
        dict with keys lon, lat, central_pressure, and
        relative_central_pressure, each a 1D array of length time.

    Raises:
        ValueError: if slp_stack is not 3D.
    """
    if slp_stack.ndim != 3:
        raise ValueError("slp_stack must be a 3D array")

    nt = slp_stack.shape[0]
    out = {k: np.empty(nt) for k in
           ("lon", "lat", "central_pressure", "relative_central_pressure")}
    for t in range(nt):
        step = find_asl(lons, lats, slp_stack[t])
        for k in out:
            out[k][t] = step[k]
    return out
