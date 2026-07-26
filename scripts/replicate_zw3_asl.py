"""
Out-of-sample replication of the ZW3-ASL findings.

The key result (strong-wave DJF slope near unity) was found through
extensive subsetting of the full record. This splits 1979-2025 into two
independent halves and asks whether each finding replicates in both.
A result present in one half only is not a result.
"""
import numpy as np
import xarray as xr
from scipy import stats
from windtools.zw3 import zw3_index
from windtools.asl import asl_timeseries
from windtools.sam import sam_index

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
mlons, mlats = slp[mlo].values, slp[mla].values
asl = asl_timeseries(mlons, mlats, slp.values)
months = slp[mt].dt.month.values
years = slp[mt].dt.year.values

outside = (mlons < 170.0) | (mlons > 298.0)
p40 = slp.sel({mla: -40.0}, method="nearest").values[:, outside].mean(axis=1)
p65 = slp.sel({mla: -65.0}, method="nearest").values[:, outside].mean(axis=1)
sam = sam_index(p40, p65)

halves = {"1979-2001": years <= 2001, "2002-2025": years >= 2002}


def circ_mean_120(p):
    a = np.deg2rad(p * 3.0)
    return (np.rad2deg(np.angle(np.mean(np.exp(1j * a)))) / 3.0) % 120.0


def anoms(sub):
    mm = months[sub]
    pa = np.empty(sub.sum())
    lo = np.empty(sub.sum())
    cp = np.empty(sub.sum())
    sm = np.empty(sub.sum())
    for m in range(1, 13):
        s = mm == m
        pa[s] = (phs[sub][s] - circ_mean_120(phs[sub][s]) + 60.0) % 120.0 - 60.0
        lo[s] = asl["lon"][sub][s] - asl["lon"][sub][s].mean()
        cp[s] = asl["central_pressure"][sub][s] - asl["central_pressure"][sub][s].mean()
        sm[s] = sam[sub][s] - sam[sub][s].mean()
    return pa, lo, cp, sm, mm, amp[sub]


def report(label, x, y):
    if len(x) < 8:
        print(f"  {label:32s} n too small")
        return
    sl, _, r, p, se = stats.linregress(x, y)
    print(f"  {label:32s} n={len(x):3d}  slope={sl:+.3f}+/-{se:.3f}  "
          f"r={r:+.3f}  p={p:.1e}")


for name, sub in halves.items():
    pa, lo, cp, sm, mm, a = anoms(sub)
    print(f"\n=== {name} (n = {sub.sum()}) ===")
    report("all months, phase vs lon", pa, lo)
    for sname, mset in (("DJF", (12, 1, 2)), ("JJA", (6, 7, 8))):
        s = np.isin(mm, mset)
        report(f"{sname}, phase vs lon", pa[s], lo[s])
        t2 = np.percentile(a[s], 66.7)
        strong = s.copy()
        strong[s] = a[s] > t2
        report(f"{sname} strong wave, phase vs lon", pa[strong], lo[strong])
    report("SAM vs absolute central P", sm, cp)

print("\nReplication criterion: a finding is credible only if the two halves "
      "agree in sign and their slopes overlap within about two standard errors.")
