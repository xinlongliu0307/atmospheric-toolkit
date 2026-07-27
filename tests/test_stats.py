import numpy as np
import pytest
from windtools.stats import (deseason, circ_mean, phase_anomaly,
                             effective_n, corr_neff)

MONTHS = np.tile(np.arange(1, 13), 20)


def test_deseason_removes_a_pure_seasonal_cycle():
    cycle = np.array([0.0, 5.0, 9.0, 12.0, 11.0, 7.0,
                      2.0, -3.0, -8.0, -10.0, -7.0, -2.0])
    x = cycle[MONTHS - 1]
    assert deseason(x, MONTHS) == pytest.approx(np.zeros_like(x), abs=1e-12)


def test_deseason_recovers_known_anomalies():
    rng = np.random.default_rng(0)
    anom = rng.normal(0, 1, len(MONTHS))
    for m in range(1, 13):
        s = MONTHS == m
        anom[s] -= anom[s].mean()
    x = 100.0 + 5.0 * MONTHS + anom
    assert deseason(x, MONTHS) == pytest.approx(anom, abs=1e-10)


def test_deseason_output_has_zero_mean_in_every_month():
    rng = np.random.default_rng(1)
    out = deseason(rng.normal(0, 1, len(MONTHS)), MONTHS)
    for m in range(1, 13):
        assert abs(np.asarray(out)[MONTHS == m].mean()) < 1e-12


def test_deseason_returns_a_flat_array_aligned_with_the_input():
    rng = np.random.default_rng(2)
    out = np.asarray(deseason(rng.normal(0, 1, len(MONTHS)), MONTHS))
    assert out.shape == MONTHS.shape


def test_deseason_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        deseason(np.zeros(10), MONTHS)


def test_circ_mean_of_identical_values():
    assert circ_mean(np.array([30.0, 30.0, 30.0]), 120.0) == pytest.approx(30.0, abs=1e-9)


def test_circ_mean_handles_the_wrap():
    # Naive arithmetic mean of 119 and 1 is 60, wrong by half a cycle.
    assert circ_mean(np.array([119.0, 1.0]), 120.0) == pytest.approx(0.0, abs=1e-9)


def test_circ_mean_is_in_range():
    rng = np.random.default_rng(3)
    for _ in range(50):
        m = circ_mean(rng.uniform(0, 120, 20), 120.0)
        assert 0.0 <= m < 120.0


def test_phase_anomaly_returns_a_flat_array_aligned_with_the_input():
    rng = np.random.default_rng(4)
    out = np.asarray(phase_anomaly(rng.uniform(0, 120, len(MONTHS)), MONTHS, 120.0))
    assert out.shape == MONTHS.shape


def test_phase_anomaly_is_zero_for_constant_phase():
    p = np.full(len(MONTHS), 42.0)
    assert phase_anomaly(p, MONTHS, 120.0) == pytest.approx(np.zeros_like(p), abs=1e-9)


def test_phase_anomaly_wraps_correctly():
    # Two January values straddling the wrap: the circular mean is 0, so
    # the anomalies are -1 and +1, not +59 and -59.
    out = phase_anomaly(np.array([119.0, 1.0]), np.array([1, 1]), 120.0)
    assert out == pytest.approx(np.array([-1.0, 1.0]), abs=1e-9)


def test_phase_anomaly_is_bounded_by_half_the_period():
    rng = np.random.default_rng(5)
    out = phase_anomaly(rng.uniform(0, 120, len(MONTHS)), MONTHS, 120.0)
    assert np.all(np.abs(np.asarray(out)) <= 60.0 + 1e-9)


def _ar1(rng, rho, n):
    e = rng.normal(0, 1, n)
    out = np.empty(n)
    out[0] = e[0]
    for i in range(1, n):
        out[i] = rho * out[i - 1] + np.sqrt(1 - rho ** 2) * e[i]
    return out


def test_effective_n_approaches_n_for_white_noise():
    rng = np.random.default_rng(6)
    x, y = rng.normal(0, 1, 2000), rng.normal(0, 1, 2000)
    assert effective_n(x, y) > 0.8 * len(x)


def test_effective_n_shrinks_for_autocorrelated_series():
    rng = np.random.default_rng(7)
    x, y = _ar1(rng, 0.9, 2000), _ar1(rng, 0.9, 2000)
    assert effective_n(x, y) < 0.3 * len(x)


def test_effective_n_uses_a_non_circular_lag_one():
    # np.roll wraps the series, which is wrong. A strong linear ramp has
    # lag-1 autocorrelation near +1 non-circularly, but the wrap-around
    # pair destroys that, so a circular implementation gives a much larger
    # effective n.
    x = np.arange(500.0)
    assert effective_n(x, x) < 20.0


def test_effective_n_is_bounded():
    rng = np.random.default_rng(8)
    x = np.cumsum(rng.normal(0, 1, 500))
    assert 3.0 <= effective_n(x, x) <= 500.0


def test_corr_neff_handles_perfect_correlation_without_dividing_by_zero():
    rng = np.random.default_rng(9)
    x = rng.normal(0, 1, 300)
    r, neff, p = corr_neff(x, x)
    assert r == pytest.approx(1.0, abs=1e-12)
    assert np.isfinite(p) and 0.0 <= p <= 1.0


def test_corr_neff_p_exceeds_the_naive_p_when_autocorrelated():
    from scipy import stats as sps
    rng = np.random.default_rng(10)
    x, y = _ar1(rng, 0.9, 500), _ar1(rng, 0.9, 500)
    _, _, p = corr_neff(x, y)
    assert p > sps.pearsonr(x, y)[1]


def test_circ_mean_never_returns_the_period_itself():
    # Values whose circular mean sits on the wrap must return 0, not the
    # period, since the documented range is [0, period).
    for pair in ([119.0, 1.0], [110.0, 10.0], [90.0, 30.0]):
        m = circ_mean(np.array(pair), 120.0)
        assert 0.0 <= m < 120.0


def test_phase_anomaly_handles_a_partial_year_without_warnings():
    import warnings
    months = np.array([1, 1, 7, 7])
    phases = np.array([10.0, 20.0, 100.0, 110.0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = phase_anomaly(phases, months, 120.0)
    assert np.all(np.isfinite(out))


def test_deseason_promotes_integer_input_to_float():
    months = np.tile(np.arange(1, 13), 3)
    x = np.arange(len(months))
    out = deseason(x, months)
    assert np.asarray(out).dtype.kind == "f"
