import warnings
import numpy as np


def orca_band_profile(nav_lat: np.ndarray, nav_lon: np.ndarray,
                      field: np.ndarray, lat_south: float, lat_north: float,
                      target_lons: np.ndarray) -> np.ndarray:
    """
    Convert a curvilinear ORCA-grid field into a regular longitude profile.

    Tropical ORCA rows have constant latitude, so rows within the band are
    selected from the first column of nav_lat and averaged directly. Source
    longitudes are converted to 0-360 and rolled to ascending order before
    interpolation. Land gaps are propagated rather than interpolated
    across, so a target longitude adjacent to a missing source value is
    returned as nan.

    Parameters:
        nav_lat, nav_lon: 2D coordinate arrays.
        field: 2D data array of the same shape.
        lat_south, lat_north: latitude band limits (degrees north).
        target_lons: 1D array of output longitudes (degrees east, 0-360).

    Returns:
        1D array of the same shape as target_lons.

    Raises:
        ValueError: if the three 2D arrays do not share a shape.
    """
    if not (nav_lat.shape == nav_lon.shape == field.shape):
        raise ValueError("nav_lat, nav_lon and field must have the same shape")

    rows = (nav_lat[:, 0] >= lat_south) & (nav_lat[:, 0] <= lat_north)
    if not np.any(rows):
        return np.full(np.shape(target_lons), np.nan)

    with warnings.catch_warnings():
        # An all-land column averages to nan, which is the intended result.
        warnings.simplefilter("ignore", RuntimeWarning)
        profile = np.nanmean(np.where(np.isfinite(field[rows]),
                                      field[rows], np.nan), axis=0)

    lons = np.asarray(nav_lon[rows][0], dtype=float) % 360.0
    order = np.argsort(lons, kind="stable")
    lons, profile = lons[order], profile[order]

    out = np.interp(np.asarray(target_lons, dtype=float), lons,
                    np.nan_to_num(profile, nan=0.0))

    # Propagate gaps: a target is nan if its nearest source cell is nan.
    nearest = np.searchsorted(lons, np.asarray(target_lons, dtype=float))
    nearest = np.clip(nearest, 1, len(lons) - 1)
    left, right = nearest - 1, nearest
    gap = np.isnan(profile[left]) | np.isnan(profile[right])
    out[gap] = np.nan
    return out
