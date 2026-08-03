import numpy as np
import pytest
from windtools.orca import orca_band_profile

# A miniature ORCA-like grid: flat rows, 0.25 deg spacing, seam in the middle.
NY, NX = 21, 1440
_lon_1d = ((np.arange(NX) * 0.25 + 72.75 + 180.0) % 360.0) - 180.0
NAV_LON = np.tile(_lon_1d, (NY, 1))
NAV_LAT = np.repeat(np.linspace(-5.0, 5.0, NY)[:, None], NX, axis=1)
TARGET = np.arange(130.0, 270.25, 0.25)


def field(fn):
    """Build a field whose value depends only on longitude (0-360)."""
    return np.tile(fn(_lon_1d % 360.0), (NY, 1))


def test_recovers_a_longitude_only_field():
    out = orca_band_profile(NAV_LAT, NAV_LON, field(lambda x: 0.01 * x),
                            -5.0, 5.0, TARGET)
    assert out.shape == TARGET.shape
    assert out == pytest.approx(0.01 * TARGET, abs=1e-3)


def test_handles_the_seam_without_artefacts():
    # A smooth function of longitude must stay smooth across the window.
    out = orca_band_profile(NAV_LAT, NAV_LON,
                            field(lambda x: np.sin(np.deg2rad(x))),
                            -5.0, 5.0, TARGET)
    assert np.max(np.abs(np.diff(out))) < 0.02


def test_averages_only_within_the_band():
    # Rows outside the band carry a huge value that must not appear.
    f = field(lambda x: np.full_like(x, 35.0))
    f[NAV_LAT[:, 0] > 2.0, :] = 1000.0
    out = orca_band_profile(NAV_LAT, NAV_LON, f, -2.0, 2.0, TARGET)
    assert out == pytest.approx(35.0, abs=1e-6)


def test_land_nans_propagate_rather_than_being_filled():
    f = field(lambda x: np.full_like(x, 35.0))
    mask = (_lon_1d % 360.0 > 200.0) & (_lon_1d % 360.0 < 210.0)
    f[:, mask] = np.nan
    out = orca_band_profile(NAV_LAT, NAV_LON, f, -5.0, 5.0, TARGET)
    inside = (TARGET > 202.0) & (TARGET < 208.0)
    assert np.all(np.isnan(out[inside]))
    assert not np.any(np.isnan(out[TARGET < 190.0]))


def test_partial_land_column_still_averages_the_wet_cells():
    f = field(lambda x: np.full_like(x, 35.0))
    f[:5, :] = np.nan
    out = orca_band_profile(NAV_LAT, NAV_LON, f, -5.0, 5.0, TARGET)
    assert out == pytest.approx(35.0, abs=1e-6)


def test_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        orca_band_profile(NAV_LAT, NAV_LON, np.zeros((NY, NX - 1)),
                          -5.0, 5.0, TARGET)
