# windtools.stats
import numpy as np
from scipy import stats
def deseason(x: np.ndarray, months: np.ndarray) -> np.ndarray:
    """
    Subtract the mean of x for each calendar month to remove seasonal cycles.

    Parameters
    ----------
    x : np.ndarray
        The input array containing values over time.
    months : np.ndarray
        An array of the same length as x, indicating the calendar month (1-12)
        corresponding to each value in x.

    Returns
    -------
    np.ndarray
        A flat array with the same shape as x, with seasonal means subtracted.

    Raises
    ------
    ValueError
        If the lengths of x and months do not match.
    """
    if len(x) != len(months):
        raise ValueError("Lengths of x and months must be equal.")
    x = np.asarray(x, dtype=float)
    deseasoned = np.zeros_like(x)
    for month in range(1, 13):
        mask = months == month
        if not np.any(mask):
            continue
        deseasoned[mask] = x[mask] - x[mask].mean()
    return deseasoned


def circ_mean(values: np.ndarray, period: float) -> float:
    """
    Calculate the circular mean of an array.

    Parameters
    ----------
    values : np.ndarray
        The input array containing angular values to average over a circle.
    period : float
        The full angle mapped onto a circle (e.g., 2 * pi or 360).

    Returns
    -------
    float
        The circular mean value, in the range [0, period).
    """
    angles = 2 * np.pi * np.asarray(values, dtype=float) / period
    complex_mean = np.exp(1j * angles).mean()
    angle_mean = np.angle(complex_mean)
    result = (angle_mean * period / (2 * np.pi)) % period
    # A tiny negative angle rounds up to exactly `period` under the modulo,
    # which would violate the documented [0, period) range.
    if result >= period:
        result = 0.0
    return float(result)


def phase_anomaly(phases: np.ndarray, months: np.ndarray, period: float) -> np.ndarray:
    """
    Calculate the wrapped difference between each phase and its calendar month's circular mean.

    Parameters
    ----------
    phases : np.ndarray
        The input array containing angular values to analyze over a circle.
    months : np.ndarray
        An array of the same length as phases, indicating the calendar month (1-12)
        corresponding to each value in phases.
    period : float
        The full angle mapped onto a circle (e.g., 2 * pi or 360).

    Returns
    -------
    np.ndarray
        A flat array with the same shape as phases, containing wrapped phase anomalies.
    """
    phases = np.asarray(phases, dtype=float)
    anomalies = np.zeros_like(phases)
    for month in range(1, 13):
        mask = months == month
        if not np.any(mask):
            continue
        ref = circ_mean(phases[mask], period)
        anomalies[mask] = (phases[mask] - ref + period / 2) % period - period / 2
    return anomalies


def effective_n(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the autocorrelation-adjusted sample size.

    Parameters
    ----------
    x : np.ndarray
        The first input array for correlation analysis.
    y : np.ndarray
        The second input array for correlation analysis, of the same length as x.

    Returns
    -------
    float
        The effective sample size adjusted for autocorrelation.
    """
    n = len(x)
    r1x = _lag_1_autocorr(x)
    r1y = _lag_1_autocorr(y)
    neff = n * (1 - r1x * r1y) / (1 + r1x * r1y)
    return max(3, min(n, neff))
def corr_neff(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Calculate Pearson correlation and effective sample size with a t-distribution p-value.

    Parameters
    ----------
    x : np.ndarray
        The first input array for correlation analysis.
    y : np.ndarray
        The second input array for correlation analysis, of the same length as x.

    Returns
    -------
    tuple
        A tuple containing the Pearson correlation coefficient, effective sample size,
        and a two-sided p-value from a t-distribution.
    """
    r, _ = stats.pearsonr(x, y)
    neff = effective_n(x, y)
    df = max(1.0, float(neff) - 2.0)
    denominator = max(1e-12, 1 - r ** 2)
    p_value = stats.t.sf(np.abs(r * np.sqrt(df / denominator)), df) * 2
    return r, neff, p_value
def _lag_1_autocorr(values: np.ndarray) -> float:
    """
    Compute the lag-1 autocorrelation of an array.

    Parameters
    ----------
    values : np.ndarray
        The input array for autocorrelation analysis.

    Returns
    -------
    float
        The lag-1 autocorrelation coefficient, computed non-circularly.
    """
    n = len(values)
    return ((values[:-1] - values.mean()) * (values[1:] - values.mean())).sum() / (n - 1) / values.var(ddof=0)