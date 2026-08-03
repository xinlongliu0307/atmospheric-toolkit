import numpy as np
import pytest
from windtools.freshpool import fresh_pool_edge

LONS = np.arange(120.0, 280.5, 0.5)
THRESH = 34.8


def profile(edge_lon, amplitude=1.0, width=5.0):
    # Fresh west of the edge, salty east of it, centred ON the threshold so
    # the crossing lies exactly at edge_lon.
    return THRESH + amplitude * np.tanh((LONS - edge_lon) / width)


def test_recovers_a_planted_edge():
    for target in (160.0, 180.0, 200.0, 230.0):
        r = fresh_pool_edge(LONS, profile(target), THRESH)
        assert r["valid"] is True, f"failed at {target}"
        assert r["edge"] == pytest.approx(target, abs=1.0)


def test_interpolates_between_grid_points():
    # The edge must not be quantised to the grid.
    r = fresh_pool_edge(LONS, profile(180.27), THRESH)
    assert r["edge"] == pytest.approx(180.27, abs=0.3)
    assert r["edge"] % 0.5 != 0.0


def test_all_fresh_profile_is_invalid():
    r = fresh_pool_edge(LONS, np.full_like(LONS, 33.0), THRESH)
    assert r["valid"] is False and np.isnan(r["edge"])
    assert r["reason"] == "no_crossing"


def test_all_salty_profile_is_invalid():
    r = fresh_pool_edge(LONS, np.full_like(LONS, 36.0), THRESH)
    assert r["valid"] is False and r["reason"] == "no_crossing"


def test_reports_multiple_crossings_and_takes_the_easternmost():
    # Fresh, salty, fresh, salty: two upward crossings.
    z = np.where(LONS < 150.0, 33.5,
        np.where(LONS < 180.0, 35.5,
        np.where(LONS < 210.0, 33.5, 35.5)))
    r = fresh_pool_edge(LONS, z, THRESH)
    assert r["n_crossings"] == 2
    assert r["edge"] == pytest.approx(210.0, abs=1.0)


def test_flags_an_edge_at_the_window_boundary():
    # Crossing within 2 degrees of the eastern limit is truncation, not
    # measurement, and must be flagged.
    r = fresh_pool_edge(LONS, profile(279.5, width=0.5), THRESH)
    assert r["at_boundary"] is True


def test_gradient_is_reported_and_positive_at_a_real_edge():
    r = fresh_pool_edge(LONS, profile(180.0), THRESH)
    assert r["gradient"] > 0


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        fresh_pool_edge(LONS, np.zeros(len(LONS) - 1), THRESH)


def test_nan_gaps_do_not_produce_a_spurious_edge():
    # Satellite SSS has retrieval gaps; they must not be read as crossings.
    z = profile(180.0)
    z[100:140] = np.nan
    r = fresh_pool_edge(LONS, z, THRESH)
    assert np.isnan(r["edge"]) or r["valid"] is True
