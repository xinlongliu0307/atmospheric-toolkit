import pytest
from windtools.speed import wind_speed


def test_wind_speed_nonzero():
    assert wind_speed(3, 4) == 5.0


def test_wind_speed_zero_wind():
    assert wind_speed(0, 0) == 0.0


def test_wind_speed_zero_v():
    assert wind_speed(3, 0) == 3.0
