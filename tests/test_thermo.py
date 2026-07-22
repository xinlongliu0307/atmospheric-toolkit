import pytest
from metpy.calc import potential_temperature
from metpy.units import units
from windtools.thermo import theta


def _metpy_theta(pressure_hpa, temperature_k):
    result = potential_temperature(
        pressure_hpa * units.hPa, temperature_k * units.kelvin
    )
    return result.to("kelvin").magnitude


@pytest.mark.parametrize(
    "pressure_hpa, temperature_k",
    [(1000.0, 288.15), (850.0, 281.65), (500.0, 265.15), (300.0, 240.0)],
)
def test_theta_matches_metpy(pressure_hpa, temperature_k):
    # Hand-rolled theta must match MetPy's reference within a small tolerance.
    expected = _metpy_theta(pressure_hpa, temperature_k)
    assert theta(pressure_hpa, temperature_k) == pytest.approx(expected, rel=1e-3)


def test_theta_reference_pressure_is_identity():
    # At 1000 hPa, potential temperature equals the actual temperature.
    assert theta(1000.0, 288.15) == pytest.approx(288.15, rel=1e-6)
