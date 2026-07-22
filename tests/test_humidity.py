import pytest
from windtools.humidity import dewpoint


def test_dewpoint_reference_value():
    # T=25 C, RH=50 % gives a Magnus dewpoint near 13.87 C.
    assert dewpoint(25.0, 50.0) == pytest.approx(13.87, abs=0.1)


def test_dewpoint_saturation_equals_temperature():
    # At 100 % relative humidity the dewpoint equals the air temperature.
    assert dewpoint(20.0, 100.0) == pytest.approx(20.0, abs=0.1)
