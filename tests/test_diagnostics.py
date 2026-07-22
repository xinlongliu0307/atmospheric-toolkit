from windtools.diagnostics import DIAGNOSTICS


def test_registry_has_base_diagnostics():
    assert {"wind_speed", "theta"} <= set(DIAGNOSTICS)
