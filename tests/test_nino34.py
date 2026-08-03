import numpy as np
import pytest
from windtools.nino34 import parse_nino34

SAMPLE = """ 1948        2026
 1948 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99 -99.99
 1950  24.55  25.06  25.87  26.28  26.18  26.46  26.29  25.88  25.74  25.69  25.47  25.29
 1951  25.24  25.71  26.90  27.58  27.92  27.73  27.60  27.02  27.23  27.20  27.25  26.91
 1952  25.00  25.00  25.00  25.00  25.00  25.00  25.00  25.00  25.00  25.00 -99.99 -99.99
   -99.99
  Nino 3.4 Mean using NOAA ERSST v6 from NCEI
 https://psl.noaa.gov/data/timeseries/month/for info
"""


def test_returns_aligned_year_month_value_arrays():
    years, months, sst = parse_nino34(SAMPLE)
    assert years.shape == months.shape == sst.shape
    assert len(years) > 0


def test_drops_sentinel_values():
    years, months, sst = parse_nino34(SAMPLE)
    assert not np.any(sst < -90)
    assert 1948 not in years          # the all-sentinel year vanishes
    assert (years == 1952).sum() == 10  # two trailing sentinels dropped


def test_reads_known_values():
    years, months, sst = parse_nino34(SAMPLE)
    jan50 = (years == 1950) & (months == 1)
    dec50 = (years == 1950) & (months == 12)
    assert sst[jan50][0] == pytest.approx(24.55)
    assert sst[dec50][0] == pytest.approx(25.29)


def test_months_run_one_to_twelve():
    years, months, sst = parse_nino34(SAMPLE)
    assert months.min() == 1 and months.max() == 12


def test_is_chronologically_ordered():
    years, months, sst = parse_nino34(SAMPLE)
    order = years * 12 + months
    assert np.all(np.diff(order) > 0)


def test_ignores_header_and_trailer_lines():
    years, _, _ = parse_nino34(SAMPLE)
    assert years.min() >= 1950 and years.max() <= 1952


def test_values_are_absolute_temperatures_not_anomalies():
    # Sanity guard: Nino 3.4 SST sits near 25-28 C, not near zero.
    _, _, sst = parse_nino34(SAMPLE)
    assert sst.mean() > 20.0
