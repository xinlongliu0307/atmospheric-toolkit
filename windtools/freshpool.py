import numpy as np


def fresh_pool_edge(lons: np.ndarray, sss: np.ndarray,
                    threshold: float) -> dict:
    """
    Locate the eastern edge of the western Pacific freshwater pool.

    Scans a longitude profile of sea surface salinity from west to east and
    finds every upward crossing of the threshold isohaline, ignoring pairs
    bounded by nan so that satellite retrieval gaps are not read as edges.
    The easternmost crossing is reported, linearly interpolated between the
    two bracketing grid points.

    Parameters:
        lons: 1D array of longitudes (degrees east), ascending.
        sss: 1D array of sea surface salinity, same length as lons.
        threshold: isohaline defining the pool edge (e.g. 34.8).

    Returns:
        dict with keys edge (degrees east, nan if none), valid, reason
        ('ok' or 'no_crossing'), n_crossings, gradient (salinity change per
        degree at the crossing), and at_boundary (True when the edge lies
        within 2 degrees of either end, indicating truncation rather than
        measurement).

    Raises:
        ValueError: if lons and sss have different lengths.
    """
    lons = np.asarray(lons, dtype=float)
    sss = np.asarray(sss, dtype=float)
    if lons.shape != sss.shape:
        raise ValueError("lons and sss must have the same length")

    edge = float("nan")
    gradient = float("nan")
    n_crossings = 0

    for i in range(1, len(lons)):
        a, b = sss[i - 1], sss[i]
        if np.isnan(a) or np.isnan(b):
            continue
        if a < threshold <= b:
            n_crossings += 1
            span = lons[i] - lons[i - 1]
            edge = float(lons[i - 1] + (threshold - a) * span / (b - a))
            gradient = float((b - a) / span)

    valid = not np.isnan(edge)
    at_boundary = bool(valid and (edge - lons[0] <= 2.0
                                  or lons[-1] - edge <= 2.0))
    return {
        "edge": edge,
        "valid": bool(valid),
        "reason": "ok" if valid else "no_crossing",
        "n_crossings": int(n_crossings),
        "gradient": gradient,
        "at_boundary": at_boundary,
    }
