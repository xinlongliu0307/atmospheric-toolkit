import numpy as np
import pytest
from windtools.asl import find_asl


def synthetic_field(centre_lon, centre_lat, depth=25.0):
    # 1-degree grid over the whole ASL search neighbourhood.
    lons = np.arange(150.0, 310.0, 1.0)
    lats = np.arange(-85.0, -49.0, 1.0)
    lon2d, lat2d = np.meshgrid(lons, lats)
    background = 995.0
    low = depth * np.exp(-(((lon2d - centre_lon) / 12.0) ** 2 + ((lat2d - centre_lat) / 6.0) ** 2))
    return lons, lats, background - low


def test_finds_planted_minimum_inside_sector():
    lons, lats, slp = synthetic_field(centre_lon=245.0, centre_lat=-72.0)
    result = find_asl(lons, lats, slp)
    assert result["lon"] == pytest.approx(245.0, abs=1.0)
    assert result["lat"] == pytest.approx(-72.0, abs=1.0)
    assert result["central_pressure"] == pytest.approx(970.0, abs=0.5)


def test_relative_central_pressure_is_negative():
    # Relative central pressure = centre minus sector mean; a real low
    # must sit below its sector mean.
    lons, lats, slp = synthetic_field(centre_lon=245.0, centre_lat=-72.0)
    result = find_asl(lons, lats, slp)
    assert result["relative_central_pressure"] < 0


def test_minimum_outside_sector_is_ignored():
    # A low planted west of the sector boundary must not be reported;
    # the reported centre must lie inside the sector bounds regardless.
    lons, lats, slp = synthetic_field(centre_lon=155.0, centre_lat=-72.0)
    result = find_asl(lons, lats, slp)
    assert 170.0 <= result["lon"] <= 298.0
    assert -80.0 <= result["lat"] <= -60.0
