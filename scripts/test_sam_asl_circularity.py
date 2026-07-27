"""
Sensitivity test: how much of the SAM-ASL correlation is circular?

The ASL sector (170-298E) is 36% of the 65S latitude circle, so the ASL
contributes to the southern node of a zonal-mean SAM index. Recompute SAM
over the complementary longitudes only and compare.
"""
import numpy as np
from windtools.stats import deseason as _deseason, corr_neff
import xarray as xr
from scipy import stats
from windtools.sam import sam_index
from windtools.asl import asl_timeseries

ds = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
pvar = "msl" if "msl" in ds else list(ds.data_vars)[0]
latname = "latitude" if "latitude" in ds.coords else "lat"
lonname = "longitude" if "longitude" in ds.coords else "lon"
timename = "valid_time" if "valid_time" in ds.coords else "time"

slp = (ds[pvar] / 100.0).squeeze(drop=True)
slp = slp.assign_coords({lonname: slp[lonname] % 360.0}).sortby(lonname)
slp = slp.sortby(latname).transpose(timename, latname, lonname)
lons, lats = slp[lonname].values, slp[latname].values
months = slp[timename].dt.month.values

outside = (lons < 170.0) | (lons > 298.0)
print(f"complementary longitudes: {outside.sum()} of {len(lons)}")

def zmean(lat, mask=None):
    row = slp.sel({latname: lat}, method="nearest")
    return (row.values[:, mask] if mask is not None else row.values).mean(axis=1)

sam_full = sam_index(zmean(-40.0), zmean(-65.0))
sam_excl = sam_index(zmean(-40.0, outside), zmean(-65.0, outside))
asl = asl_timeseries(lons, lats, slp.values)

def deseason(x):
    return _deseason(x, months)

def corr(x, y):
    return corr_neff(deseason(x), deseason(y))

print(f"\nSAM(full) vs SAM(sector excluded): r = "
      f"{np.corrcoef(deseason(sam_full), deseason(sam_excl))[0,1]:+.3f}")

print("\n                          full-circle SAM      sector-excluded SAM")
for label, s in (("absolute central P", asl["central_pressure"]),
                 ("relative central P", asl["relative_central_pressure"])):
    r1, n1, p1 = corr(sam_full, s)
    r2, n2, p2 = corr(sam_excl, s)
    print(f"{label:22s}  r={r1:+.3f} p={p1:.1e}    r={r2:+.3f} p={p2:.1e}")

print("\n--- Seasonal stratification (sector-excluded SAM, absolute P) ---")
for name, mm in (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)),
                 ("JJA", (6, 7, 8)), ("SON", (9, 10, 11))):
    sel = np.isin(months, mm)
    r = np.corrcoef(sam_excl[sel] - sam_excl[sel].mean(),
                    asl["central_pressure"][sel] - asl["central_pressure"][sel].mean())[0, 1]
    print(f"  {name}: r = {r:+.3f}")
