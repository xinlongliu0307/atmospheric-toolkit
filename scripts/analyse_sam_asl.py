"""
Is the Amundsen Sea Low distinct from the Southern Annular Mode?

Method choices, stated explicitly:
  - SAM: Gong and Wang (1999) form, computed from ERA5 zonal-mean MSLP at
    40S and 65S rather than station data. Normalised over the full record.
  - Correlations use deseasonalised anomalies (calendar-month climatology
    removed) to isolate interannual variability from the seasonal cycle.
  - Significance uses an effective sample size adjusted for lag-1
    autocorrelation in both series.
  - Total fields throughout; no linear detrending (reported separately).
"""
import numpy as np
import xarray as xr
from scipy import stats
from windtools.sam import sam_index
from windtools.asl import asl_timeseries
from windtools.stats import deseason, corr_neff

ds = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
pvar = "msl" if "msl" in ds else list(ds.data_vars)[0]
latname = "latitude" if "latitude" in ds.coords else "lat"
lonname = "longitude" if "longitude" in ds.coords else "lon"
timename = "valid_time" if "valid_time" in ds.coords else "time"

slp = (ds[pvar] / 100.0).squeeze(drop=True)  # Pa -> hPa
slp = slp.assign_coords({lonname: slp[lonname] % 360.0}).sortby(lonname)
slp = slp.sortby(latname).transpose(timename, latname, lonname)

lons = slp[lonname].values
lats = slp[latname].values
months = slp[timename].dt.month.values
print(f"field {slp.shape}, lat {lats.min():.0f} to {lats.max():.0f}, "
      f"lon {lons.min():.0f} to {lons.max():.0f}")

# --- SAM ---
p40 = slp.sel({latname: -40.0}, method="nearest").mean(dim=lonname).values
p65 = slp.sel({latname: -65.0}, method="nearest").mean(dim=lonname).values
sam = sam_index(p40, p65)

# --- ASL ---
asl = asl_timeseries(lons, lats, slp.values)

print("\n--- ASL climatology (sanity check) ---")
print("month  central(hPa)  relative(hPa)   lon(E)   lat(S)")
for m in range(1, 13):
    s = months == m
    print(f"  {m:02d}   {asl['central_pressure'][s].mean():9.1f}"
          f"   {asl['relative_central_pressure'][s].mean():10.1f}"
          f"   {asl['lon'][s].mean():7.1f}  {-asl['lat'][s].mean():6.1f}")

print("\n--- SAM vs ASL (deseasonalised) ---")
for label, series in (("absolute central pressure", asl["central_pressure"]),
                      ("relative central pressure", asl["relative_central_pressure"]),
                      ("ASL longitude", asl["lon"])):
    r, neff, p = corr_neff(deseason(sam, months), deseason(series, months))
    print(f"{label:28s} r = {r:+.3f}   n_eff = {neff:6.1f}   p = {p:.2e}")

print("\nExpectation: strong negative r with absolute central pressure "
      "(positive SAM deepens the circumpolar trough), weak r with relative "
      "central pressure (the sector-mean background is removed by construction).")
