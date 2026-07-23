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
