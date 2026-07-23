"""Validate zw3_index against ERA5: winter amplitude max, phase clustering."""
import numpy as np
import xarray as xr
from windtools.zw3 import zw3_index

G = 9.80665

ds = xr.open_dataset("data/era5_z500_monthly_40S70S.nc")
# Variable and coordinate names vary between CDS backends; resolve them.
zvar = "z" if "z" in ds else list(ds.data_vars)[0]
latname = "latitude" if "latitude" in ds.coords else "lat"
lonname = "longitude" if "longitude" in ds.coords else "lon"
timename = "valid_time" if "valid_time" in ds.coords else "time"

circle = ds[zvar].sel({latname: -49.0}) / G  # geopotential -> height (m)

# Roll longitudes into 0..360 ascending so the phase convention anchors at 0E.
lons = circle[lonname].values % 360.0
order = np.argsort(lons)
lons_sorted = lons[order]

months = circle[timename].dt.month.values
amp = np.empty(circle.sizes[timename])
phs = np.empty(circle.sizes[timename])
for t in range(circle.sizes[timename]):
    z = circle.isel({timename: t}).values[order]
    result = zw3_index(lons_sorted, z)
    amp[t], phs[t] = result["amplitude"], result["phase"]

# Verdict 1: seasonal cycle of amplitude, expecting an austral winter max.
clim = {m: amp[months == m].mean() for m in range(1, 13)}
print("Monthly mean ZW3 amplitude (m):")
for m in range(1, 13):
    print(f"  {m:02d}: {clim[m]:6.1f}")
jja = np.mean([clim[6], clim[7], clim[8]])
djf = np.mean([clim[12], clim[1], clim[2]])
print(f"\nJJA mean: {jja:.1f} m | DJF mean: {djf:.1f} m")
print("WINTER MAX:", "PASS" if jja > djf else "FAIL")

# Verdict 2: phase clustering (circular stats on the 120-degree cycle).
ang = np.deg2rad(phs * 3.0)  # map 120-degree cycle onto the full circle
R = np.abs(np.mean(np.exp(1j * ang)))
mean_phase = (np.rad2deg(np.angle(np.mean(np.exp(1j * ang)))) / 3.0) % 120.0
print(f"\nPhase clustering R = {R:.2f} (0 = uniform, 1 = fixed)")
print(f"Circular mean phase = {mean_phase:.1f} deg")
print("QUASI-STATIONARY:", "PASS" if R > 0.3 else "FAIL")
