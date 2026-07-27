import numpy as np


def dipole_node(lons: np.ndarray, anomaly: np.ndarray) -> dict:
    """
    Locate the node of a west-positive, east-negative dipole in a profile.

    Fits a single sinusoid across the window and reports the descending
    zero crossing of the fit, together with the fraction of profile
    variance the fit explains. A dipole is reported as present only when
    that fraction reaches 0.5, so noisy or structureless profiles return
    no node rather than an arbitrary one.

    Known behaviour: for step-like dipoles the recovered node is
    compressed toward the window centre by a factor of two (an analytic
    property of the wavenumber-1 fit). The estimate is monotonic in the
    true node, so it attenuates a regression slope rather than biasing
    its sign. For pure sinusoidal dipoles the recovery is exact.

    Known limitation: detection weakens toward the window edges. For an
    ideal step dipole at fractional position p across the window,
    r_squared = 2 sin^2(pi p) / (pi^2 p (1-p)), peaking at 8/pi^2 = 0.811
    at the centre and crossing the 0.5 threshold near p = 0.23 and 0.77.
    Smoothed steps sit above this curve, since low-passing suppresses the
    higher harmonics more than the fundamental; a tanh of width 8 across a
    130 degree window reaches about 0.86 at the centre.
    Step dipoles with nodes outside the central ~54 percent of the window
    are therefore reported as invalid. This is intended: a step near the
    edge is a plateau rather than a dipole.

    IMPORTANT LIMITATION. The validity flag distinguishes large-scale from
    small-scale structure, not dipoles from monopoles. A centred Gaussian
    of width 30-40 degrees on a 130 degree window returns r_squared of
    0.96-1.00 and valid=True despite containing no dipole. More generally,
    on a field whose decorrelation length is a substantial fraction of the
    window, red noise passes at the same rate as real data. Use this
    estimator only where the window is long compared with the field's
    decorrelation length, and validate against a matched red-noise null
    before interpreting any node position.

    Parameters:
        lons: 1D array of longitudes (degrees east), monotonic.
        anomaly: 1D array of anomaly values, same length as lons.

    Returns:
        dict with keys node (degrees east, nan if invalid), contrast
        (peak-to-peak amplitude of the fit), r_squared, and valid.

    Raises:
        ValueError: if lons and anomaly have different lengths.
    """
    lons = np.asarray(lons, dtype=float)
    anomaly = np.asarray(anomaly, dtype=float)
    if lons.shape != anomaly.shape:
        raise ValueError("lons and anomaly must have the same length")

    span = float(lons[-1] - lons[0])
    x = 2.0 * np.pi * (lons - lons[0]) / span
    y = anomaly - anomaly.mean()

    design = np.column_stack([np.cos(x), np.sin(x)])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    a, b = float(coef[0]), float(coef[1])

    fit = design @ coef
    ss_tot = float(np.sum(y ** 2))
    ss_res = float(np.sum((y - fit) ** 2))
    r_squared = 0.0 if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot

    contrast = 2.0 * float(np.hypot(a, b))
    valid = bool(r_squared >= 0.5)
    if valid:
        x_node = (np.arctan2(b, a) + 0.5 * np.pi) % (2.0 * np.pi)
        node = float(lons[0] + span * x_node / (2.0 * np.pi))
    else:
        node = float("nan")

    return {"node": node, "contrast": contrast,
            "r_squared": float(r_squared), "valid": valid}
