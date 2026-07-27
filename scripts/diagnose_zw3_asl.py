"""
Diagnostics for the ZW3-ASL relationship.

Resolves three interpretive ambiguities:
  1. Range restriction: correlation is attenuated by low predictor
     variance, regression slope is not. Compare both by season.
  2. Boundary contamination: a grid-point minimum can be falsely placed
     on the sector edge (a documented weakness of absolute-pressure
     ASL location). Count how often this happens.
  3. Regression dilution: noise in the phase estimate biases the slope
     toward zero, so report an amplitude-stratified slope as a check
     (phase is better determined when the wave is strong).
"""
import numpy as np
from windtools.stats import deseason as _deseason, circ_mean, phase_anomaly
import xarray as xr
from scipy import stats
from windtools.zw3 import zw3_index
from windtools.asl import asl_timeseries

G = 9.80665


def coords(ds):
    return ("latitude" if "latitude" in ds.coords else "lat",
            "longitude" if "longitude" in ds.coords else "lon",
            "valid_time" if "valid_time" in ds.coords else "time")


dz = xr.open_dataset("data/era5_z500_monthly_40S70S.nc")
zla, zlo, zt = coords(dz)
zvar = "z" if "z" in dz else list(dz.data_vars)[0]
circle = (dz[zvar].sel({zla: -49.0}, method="nearest") / G).squeeze(drop=True)
circle = circle.assign_coords({zlo: circle[zlo] % 360.0}).sortby(zlo)
vals, first = np.unique(circle[zlo].values, return_index=True)
if len(vals) != circle.sizes[zlo]:
    circle = circle.isel({zlo: np.sort(first)})
circle = circle.transpose(zt, zlo)
zl, zv = circle[zlo].values, circle.values
amp = np.array([zw3_index(zl, zv[t])["amplitude"] for t in range(zv.shape[0])])
phs = np.array([zw3_index(zl, zv[t])["phase"] for t in range(zv.shape[0])])

dm = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
mla, mlo, mt = coords(dm)
pvar = "msl" if "msl" in dm else list(dm.data_vars)[0]
slp = (dm[pvar] / 100.0).squeeze(drop=True)
slp = slp.assign_coords({mlo: slp[mlo] % 360.0}).sortby(mlo)
slp = slp.sortby(mla).transpose(mt, mla, mlo)
asl = asl_timeseries(slp[mlo].values, slp[mla].values, slp.values)
months = slp[mt].dt.month.values


def deseason(x):
    return _deseason(x, months)


def circ_mean_120(p):
    return circ_mean(p, 120.0)


def phase_anom(p):
    return phase_anomaly(p, months, 120.0)


pa, lona = phase_anom(phs), deseason(asl["lon"])

print("--- 1. Range restriction: correlation vs slope by season ---")
print("season   sd(phase anom)  sd(ASL lon anom)      r      slope +/- se")
for name, mm in (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)),
                 ("JJA", (6, 7, 8)), ("SON", (9, 10, 11))):
    s = np.isin(months, mm)
    sl, _, r, _, se = stats.linregress(pa[s], lona[s])
    print(f"  {name}      {pa[s].std():6.2f} deg      {lona[s].std():6.2f} deg   "
          f"{r:+.3f}   {sl:+.3f} +/- {se:.3f}")

print("\n--- 2. Boundary contamination ---")
lon = asl["lon"]
west = np.isclose(lon, lon.min())
east = np.isclose(lon, lon.max())
print(f"ASL longitude range: {lon.min():.0f} to {lon.max():.0f} E "
      f"(sector 170-298 E, midpoint 234 E)")
print(f"months on the western edge: {west.sum()} ({100*west.mean():.1f}%)")
print(f"months on the eastern edge: {east.sum()} ({100*east.mean():.1f}%)")
lat = asl["lat"]
print(f"months on a latitude edge:  "
      f"{(np.isclose(lat, -80.0) | np.isclose(lat, -60.0)).sum()}")
print(f"mean ASL longitude {lon.mean():.1f} E vs sector midpoint 234.0 E")

print("\n--- 3. Regression dilution: slope by ZW3 amplitude tercile ---")
q1, q2 = np.percentile(amp, [33.3, 66.7])
for name, s in (("weak  ", amp <= q1), ("medium", (amp > q1) & (amp <= q2)),
                ("strong", amp > q2)):
    sl, _, r, _, se = stats.linregress(pa[s], lona[s])
    print(f"  {name} wave (n={s.sum():3d}): slope = {sl:+.3f} +/- {se:.3f}, "
          f"r = {r:+.3f}")
