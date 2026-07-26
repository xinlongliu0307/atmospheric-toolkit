import numpy as np
import pytest
from windtools.asl import asl_timeseries

LONS = np.arange(150.0, 310.0, 1.0)
LATS = np.arange(-85.0, -49.0, 1.0)
CENTRES = [(200.0, -70.0), (245.0, -72.0), (280.0, -68.0)]


def stack():
    lon2d, lat2d = np.meshgrid(LONS, LATS)
    fields = []
    for clon, clat in CENTRES:
        low = 25.0 * np.exp(-(((lon2d - clon) / 12.0) ** 2 + ((lat2d - clat) / 6.0) ** 2))
        fields.append(995.0 - low)
    return np.stack(fields)


def test_tracks_each_timestep():
    result = asl_timeseries(LONS, LATS, stack())
    assert result["lon"].shape == (3,)
    assert result["lon"] == pytest.approx([c[0] for c in CENTRES], abs=1.0)
    assert result["lat"] == pytest.approx([c[1] for c in CENTRES], abs=1.0)


def test_all_relative_pressures_negative():
    result = asl_timeseries(LONS, LATS, stack())
    assert np.all(result["relative_central_pressure"] < 0)


def test_rejects_wrong_dimensionality():
    with pytest.raises(ValueError):
        asl_timeseries(LONS, LATS, stack()[0])
