"""
Does windtools.stats.effective_n reproduce the inline helper used in the
published analyses? Extraction must not move any reported number.

The agent's implementation uses a classical autocovariance estimator
(single series mean, population variance). Every inline helper in the
analysis scripts used np.corrcoef on the two lagged sub-series. These are
different estimators, so the difference must be measured before any
script is refactored to import from windtools.stats.
"""
import numpy as np
import xarray as xr
from windtools.stats import effective_n, corr_neff, deseason
from windtools.sam import sam_index
from windtools.asl import asl_timeseries


def inline_neff(x, y):
    rx = np.corrcoef(x[:-1], x[1:])[0, 1]
    ry = np.corrcoef(y[:-1], y[1:])[0, 1]
    n = len(x)
    return min(max(3.0, n * (1 - rx * ry) / (1 + rx * ry)), float(n))


ds = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
la = "latitude" if "latitude" in ds.coords else "lat"
lo = "longitude" if "longitude" in ds.coords else "lon"
t = "valid_time" if "valid_time" in ds.coords else "time"
pvar = "msl" if "msl" in ds else list(ds.data_vars)[0]

slp = (ds[pvar] / 100.0).squeeze(drop=True)
slp = slp.assign_coords({lo: slp[lo] % 360.0}).sortby(lo).sortby(la)
slp = slp.transpose(t, la, lo)
lons, lats = slp[lo].values, slp[la].values
months = slp[t].dt.month.values

outside = (lons < 170.0) | (lons > 298.0)
sam = sam_index(
    slp.sel({la: -40.0}, method="nearest").values[:, outside].mean(axis=1),
    slp.sel({la: -65.0}, method="nearest").values[:, outside].mean(axis=1),
)
asl = asl_timeseries(lons, lats, slp.values)

x = np.asarray(deseason(sam, months))
y = np.asarray(deseason(asl["central_pressure"], months))

print(f"records: {len(x)} months\n")
print("Published result: sector-excluded SAM vs ASL absolute central")
print("pressure, r = -0.585, n_eff ~ 490\n")
print(f"  correlation                    {np.corrcoef(x, y)[0, 1]:+.6f}")
inl = inline_neff(x, y)
ext = effective_n(x, y)
print(f"  inline effective_n (corrcoef)  {inl:.3f}")
print(f"  windtools.stats.effective_n    {ext:.3f}")
print(f"  absolute difference            {abs(ext - inl):.3f}")

r, ne, p = corr_neff(x, y)
print(f"\n  corr_neff -> r {r:+.6f}, n_eff {ne:.3f}, p {p:.3e}")
print(f"\n  verdict: {'EQUIVALENT at this precision'
                     if abs(ext - inl) < 1.0
                     else 'DIFFERS - decide which estimator to keep'}")
