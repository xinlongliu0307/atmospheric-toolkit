"""Validate zw3_index against ERA5: winter amplitude max, phase clustering."""
import numpy as np
import xarray as xr
from windtools.zw3 import zw3_index

G = 9.80665

ds = xr.open_dataset("data/era5_z500_monthly_40S70S.nc")
zvar = "z" if "z" in ds else list(ds.data_vars)[0]
latname = "latitude" if "latitude" in ds.coords else "lat"
lonname = "longitude" if "longitude" in ds.coords else "lon"
timename = "valid_time" if "valid_time" in ds.coords else "time"

# Select the 49S circle, convert geopotential to height, and drop singleton
# dimensions such as a single pressure_level.
circle = (ds[zvar].sel({latname: -49.0}, method="nearest") / G).squeeze(drop=True)

# Put longitudes on 0-360 ascending so the phase convention anchors at 0E.
circle = circle.assign_coords({lonname: circle[lonname] % 360.0}).sortby(lonname)

# Guard: drop a duplicated endpoint if both -180 and 180 were returned.
vals, first = np.unique(circle[lonname].values, return_index=True)
if len(vals) != circle.sizes[lonname]:
    circle = circle.isel({lonname: np.sort(first)})

circle = circle.transpose(timename, lonname)
lons = circle[lonname].values
print(f"field dims {circle.dims} shape {circle.shape}")
print(f"lon range {lons[0]:.1f} to {lons[-1]:.1f}, n = {len(lons)}")

# The FFT requires an evenly spaced grid closing the full circle.
step = np.diff(lons)
assert np.allclose(step, step[0]), "longitudes are not evenly spaced"
assert np.isclose(lons[0] + 360.0 - lons[-1], step[0]), "longitudes do not close the circle"

values = circle.values
months = circle[timename].dt.month.values
nt = values.shape[0]
amp = np.empty(nt)
phs = np.empty(nt)
for t in range(nt):
    result = zw3_index(lons, values[t])
    amp[t], phs[t] = result["amplitude"], result["phase"]

# Verdict 1: seasonal cycle of amplitude, expecting an austral winter maximum.
clim = {m: amp[months == m].mean() for m in range(1, 13)}
print("\nMonthly mean ZW3 amplitude (m):")
for m in range(1, 13):
    print(f"  {m:02d}: {clim[m]:6.1f}")
jja = np.mean([clim[6], clim[7], clim[8]])
djf = np.mean([clim[12], clim[1], clim[2]])
print(f"\nJJA mean: {jja:.1f} m | DJF mean: {djf:.1f} m")
print("WINTER MAX:", "PASS" if jja > djf else "FAIL")

# Verdict 2: phase clustering (circular statistics on the 120-degree cycle).
ang = np.deg2rad(phs * 3.0)
R = np.abs(np.mean(np.exp(1j * ang)))
mean_phase = (np.rad2deg(np.angle(np.mean(np.exp(1j * ang)))) / 3.0) % 120.0
print(f"\nPhase clustering R = {R:.2f} (0 = uniform, 1 = fixed)")
print(f"Circular mean phase = {mean_phase:.1f} deg")
print("QUASI-STATIONARY:", "PASS" if R > 0.3 else "FAIL")
