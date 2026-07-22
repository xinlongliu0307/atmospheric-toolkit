from windtools.diagnostics import DIAGNOSTICS
from windtools.speed import wind_direction


def test_registry_has_base_diagnostics():
    assert {"wind_speed", "theta"} <= set(DIAGNOSTICS)


def test_wind_direction_westerly():
    # u=10, v=0 blows toward the east, so it comes FROM the west: 270 degrees.
    assert wind_direction(10, 0) == 270.0


def test_wind_direction_in_diagnostics():
    assert "wind_direction" in DIAGNOSTICS
