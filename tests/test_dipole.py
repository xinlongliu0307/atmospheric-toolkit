import numpy as np
import pytest
from windtools.dipole import dipole_node

LONS = np.arange(160.0, 291.0, 1.0)
SPAN = LONS[-1] - LONS[0]
CENTRE = 0.5 * (LONS[0] + LONS[-1])


def sinusoidal_dipole(node_lon, amplitude=1.0):
    # cos(x - phi) whose descending zero crossing lies at node_lon.
    x = 2.0 * np.pi * (LONS - LONS[0]) / SPAN
    phi = 2.0 * np.pi * (node_lon - LONS[0]) / SPAN - 0.5 * np.pi
    return amplitude * np.cos(x - phi)


def step_dipole(node_lon, amplitude=1.0):
    # Positive west of the node, negative east of it.
    return amplitude * np.tanh((node_lon - LONS) / 8.0)


def test_recovers_sinusoidal_node_exactly():
    # Exact for the shape the estimator fits. Catches phase-convention
    # errors such as reporting the maximum instead of the zero crossing.
    for target in (190.0, 225.0, 265.0):
        r = dipole_node(LONS, sinusoidal_dipole(target))
        assert r["valid"] is True, f"failed at {target}"
        assert r["node"] == pytest.approx(target, abs=1.0)


def test_step_dipole_node_is_compressed_toward_centre():
    # Analytic property: a wavenumber-1 fit to a step-like dipole recovers
    # the node at the midpoint of the true node and the window centre.
    # Restricted to interior nodes; see the edge test below.
    for target in (200.0, 215.0, 240.0, 250.0):
        r = dipole_node(LONS, step_dipole(target))
        assert r["valid"] is True, f"failed at {target}"
        assert r["node"] == pytest.approx(0.5 * (target + CENTRE), abs=3.0)


def test_step_dipole_near_window_edge_is_invalid():
    # r2 = 2 sin^2(pi p) / (pi^2 p (1-p)) crosses 0.5 near p = 0.23 and
    # 0.77, which is 190E and 260E for this window. A step near the edge
    # is a plateau, not a dipole, and must not be reported as one.
    for target in (175.0, 185.0, 265.0, 278.0):
        assert dipole_node(LONS, step_dipole(target))["valid"] is False, \
            f"edge case at {target} was wrongly accepted"


def test_ideal_step_r_squared_matches_the_analytic_value():
    # An ideal step is a square wave: power falls across the odd harmonics
    # as 1/n^2, so the wavenumber-1 share is 8/pi^2. This is the analytic
    # claim in the docstring, and it holds only for an unsmoothed step.
    ideal = np.sign(CENTRE - LONS)
    assert dipole_node(LONS, ideal)["r_squared"] == pytest.approx(
        8.0 / np.pi ** 2, abs=0.03)


def test_smoothed_step_exceeds_the_ideal_step_value():
    # Smoothing low-passes the square wave, suppressing harmonic 3 more
    # than harmonic 1 and so raising the wavenumber-1 share above 8/pi^2.
    assert dipole_node(LONS, step_dipole(CENTRE))["r_squared"] > 8.0 / np.pi ** 2


def test_step_dipole_r_squared_peaks_at_the_centre():
    # The property that matters for the validity criterion: detection is
    # strongest at the window centre and weakens toward the edges.
    centre = dipole_node(LONS, step_dipole(CENTRE))["r_squared"]
    for offset in (25.0, 45.0):
        for target in (CENTRE - offset, CENTRE + offset):
            assert dipole_node(LONS, step_dipole(target))["r_squared"] < centre


def test_step_dipole_node_is_monotonic():
    # Compression is acceptable for regression only if it is monotonic.
    nodes = [dipole_node(LONS, step_dipole(t))["node"]
             for t in (200.0, 215.0, 225.0, 235.0, 250.0)]
    assert all(np.diff(nodes) > 0)


def test_flat_profile_is_invalid():
    r = dipole_node(LONS, np.zeros_like(LONS))
    assert r["valid"] is False
    assert np.isnan(r["node"])


def test_narrow_monopole_is_invalid():
    # A narrow centred bump spreads its power across many wavenumbers,
    # so the wavenumber-1 fit explains little of it.
    r = dipole_node(LONS, np.exp(-((LONS - 225.0) / 8.0) ** 2))
    assert r["valid"] is False


def test_noise_rarely_yields_a_valid_node():
    rng = np.random.default_rng(0)
    valid = sum(dipole_node(LONS, rng.normal(0, 1, len(LONS)))["valid"]
                for _ in range(300))
    assert valid / 300 < 0.10


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        dipole_node(LONS, np.zeros(len(LONS) - 1))
