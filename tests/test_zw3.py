import numpy as np
import pytest
from windtools.zw3 import zw3_index


LONS = np.arange(0.0, 360.0, 2.5)  # 144 points, evenly spaced, full circle


def wave(k, amplitude, phase_deg, mean=0.0):
    # amplitude * cos(k * (lon - phase)) + mean
    return mean + amplitude * np.cos(np.deg2rad(k * (LONS - phase_deg)))


def test_recovers_planted_wave3_amplitude():
    z = wave(k=3, amplitude=50.0, phase_deg=0.0)
    result = zw3_index(LONS, z)
    assert result["amplitude"] == pytest.approx(50.0, rel=1e-6)


def test_recovers_planted_wave3_phase():
    z = wave(k=3, amplitude=50.0, phase_deg=30.0)
    result = zw3_index(LONS, z)
    assert result["phase"] == pytest.approx(30.0, abs=1e-6)


def test_phase_is_reported_modulo_120():
    # A wave-3 pattern repeats every 120 degrees, so phase 150 == phase 30.
    z = wave(k=3, amplitude=50.0, phase_deg=150.0)
    result = zw3_index(LONS, z)
    assert result["phase"] == pytest.approx(30.0, abs=1e-6)
    assert 0.0 <= result["phase"] < 120.0


def test_orthogonality_wave2_gives_zero_wave3():
    # A pure wave-2 field must contain no wave-3 signal.
    z = wave(k=2, amplitude=80.0, phase_deg=10.0)
    result = zw3_index(LONS, z)
    assert result["amplitude"] == pytest.approx(0.0, abs=1e-9)


def test_zonal_mean_does_not_leak_into_wave3():
    # Adding a large constant (the zonal mean) must not change the index.
    z = wave(k=3, amplitude=50.0, phase_deg=30.0, mean=5400.0)
    result = zw3_index(LONS, z)
    assert result["amplitude"] == pytest.approx(50.0, rel=1e-6)
    assert result["phase"] == pytest.approx(30.0, abs=1e-6)


def test_mixed_field_isolates_wave3():
    # Wave 1 + wave 3 + mean: only the wave-3 part must be reported.
    z = wave(k=1, amplitude=100.0, phase_deg=45.0, mean=5400.0) \
        + wave(k=3, amplitude=25.0, phase_deg=60.0)
    result = zw3_index(LONS, z)
    assert result["amplitude"] == pytest.approx(25.0, rel=1e-6)
    assert result["phase"] == pytest.approx(60.0, abs=1e-6)


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        zw3_index(LONS, np.zeros(len(LONS) - 1))
