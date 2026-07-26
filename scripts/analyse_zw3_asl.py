"""
Does the zonal wave 3 phase set the longitude of the Amundsen Sea Low?

Method choices:
  - ZW3 amplitude and phase from ERA5 monthly Z500 at 49S (total field).
  - ASL from ERA5 monthly MSLP, sector 170-298E / 60-80S.
  - ZW3 phase anomalies use WRAPPED differences on the 120-degree cycle,
    relative to each calendar month's circular-mean phase.
  - ASL longitude anomalies are plain deseasonalised anomalies (the ASL
    stays well inside the sector, so no wrapping is required).
  - Hypothesis: the ASL sits in the wave-3 trough at phase + 180 deg, so
    the regression slope of ASL longitude on phase anomaly should be ~1.
"""
import numpy as np
import xarray as xr
from scipy import stats
from windtools.zw3 import zw3_index
from windtools.asl import asl_timeseries

G = 9.80665


def coords(ds):
    return ("latitude" if "latitude" in ds.coords else "lat",
            "longitude" if "longitude" in ds.coords else "lon",
            "valid_time" if "valid_time" in ds.coords else "time")


# --- ZW3 from Z500 ---
dz = xr.open_dataset("data/era5_z500_monthly_40S70S.nc")
zla, zlo, zt = coords(dz)
zvar = "z" if "z" in dz else list(dz.data_vars)[0]
circle = (dz[zvar].sel({zla: -49.0}, method="nearest") / G).squeeze(drop=True)
circle = circle.assign_coords({zlo: circle[zlo] % 360.0}).sortby(zlo)
vals, first = np.unique(circle[zlo].values, return_index=True)
if len(vals) != circle.sizes[zlo]:
    circle = circle.isel({zlo: np.sort(first)})
circle = circle.transpose(zt, zlo)
zlons = circle[zlo].values
zvals = circle.values
amp = np.empty(zvals.shape[0])
phs = np.empty(zvals.shape[0])
for t in range(zvals.shape[0]):
    r = zw3_index(zlons, zvals[t])
    amp[t], phs[t] = r["amplitude"], r["phase"]
ztime = circle[zt].values

# --- ASL from MSLP ---
dm = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
mla, mlo, mt = coords(dm)
pvar = "msl" if "msl" in dm else list(dm.data_vars)[0]
slp = (dm[pvar] / 100.0).squeeze(drop=True)
slp = slp.assign_coords({mlo: slp[mlo] % 360.0}).sortby(mlo)
slp = slp.sortby(mla).transpose(mt, mla, mlo)
asl = asl_timeseries(slp[mlo].values, slp[mla].values, slp.values)
mtime = slp[mt].values

# --- align ---
assert len(ztime) == len(mtime) and np.array_equal(ztime, mtime), \
    "time axes differ between the Z500 and MSLP files"
months = slp[mt].dt.month.values
print(f"aligned records: {len(months)} months")

# --- anomalies ---
def deseason(x):
    out = np.empty_like(x, dtype=float)
    for m in range(1, 13):
        s = months == m
        out[s] = x[s] - x[s].mean()
    return out


def circ_mean_120(p):
    a = np.deg2rad(p * 3.0)
    return (np.rad2deg(np.angle(np.mean(np.exp(1j * a)))) / 3.0) % 120.0


def phase_anom(p):
    out = np.empty_like(p, dtype=float)
    for m in range(1, 13):
        s = months == m
        ref = circ_mean_120(p[s])
        out[s] = (p[s] - ref + 60.0) % 120.0 - 60.0   # wrapped into (-60, 60]
    return out


def corr(x, y):
    r = np.corrcoef(x, y)[0, 1]
    rx = np.corrcoef(x[:-1], x[1:])[0, 1]
    ry = np.corrcoef(y[:-1], y[1:])[0, 1]
    n = len(x)
    neff = min(max(3.0, n * (1 - rx * ry) / (1 + rx * ry)), float(n))
    t = r * np.sqrt((neff - 2) / max(1e-12, 1 - r**2))
    return r, neff, 2 * stats.t.sf(abs(t), df=neff - 2)


pa = phase_anom(phs)
lona = deseason(asl["lon"])
ampa = deseason(amp)
rcpa = deseason(asl["relative_central_pressure"])

print("\n--- Climatological check ---")
print(f"ZW3 mean phase {circ_mean_120(phs):.1f} deg -> sector trough at "
      f"{(circ_mean_120(phs) + 180.0) % 360.0:.1f} E")
print(f"ASL mean longitude {asl['lon'].mean():.1f} E")

print("\n--- Does ZW3 phase set ASL longitude? ---")
r, n, p = corr(pa, lona)
slope, icpt, _, _, se = stats.linregress(pa, lona)
print(f"phase anomaly vs ASL longitude anomaly: r = {r:+.3f}  "
      f"n_eff = {n:.0f}  p = {p:.2e}")
print(f"regression slope = {slope:+.3f} +/- {se:.3f} (hypothesis: +1.0)")

print("\n--- Does ZW3 amplitude set ASL relative depth? ---")
r, n, p = corr(ampa, rcpa)
print(f"amplitude anomaly vs relative central pressure: r = {r:+.3f}  "
      f"n_eff = {n:.0f}  p = {p:.2e}")

print("\n--- Seasonal breakdown (phase vs longitude) ---")
for name, mm in (("DJF", (12, 1, 2)), ("MAM", (3, 4, 5)),
                 ("JJA", (6, 7, 8)), ("SON", (9, 10, 11))):
    s = np.isin(months, mm)
    print(f"  {name}: r = {np.corrcoef(pa[s], lona[s])[0,1]:+.3f}")
