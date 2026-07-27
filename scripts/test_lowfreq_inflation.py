"""
Why did the pooled slope exceed both half-sample slopes?

Hypothesis: shared low-frequency covariance between ZW3 phase and ASL
longitude survives full-record calendar-month centring but is removed by
within-half centring. Test by re-estimating the slope after removing
progressively more low-frequency variance.
"""
import numpy as np
from windtools.stats import deseason, circ_mean, phase_anomaly
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
c = (dz[zvar].sel({zla: -49.0}, method="nearest") / G).squeeze(drop=True)
c = c.assign_coords({zlo: c[zlo] % 360.0}).sortby(zlo)
v, f = np.unique(c[zlo].values, return_index=True)
if len(v) != c.sizes[zlo]:
    c = c.isel({zlo: np.sort(f)})
c = c.transpose(zt, zlo)
zl, zv = c[zlo].values, c.values
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
years = slp[mt].dt.year.values


def circ_mean_120(p):
    return circ_mean(p, 120.0)


def pa_full(p):
    return phase_anomaly(p, months, 120.0)


def ds_full(x):
    return deseason(x, months)


def detrend(x):
    t = np.arange(len(x), dtype=float)
    return x - np.polyval(np.polyfit(t, x, 1), t)


def highpass(x, win=61):
    k = np.ones(win) / win
    return x - np.convolve(np.pad(x, win // 2, mode="edge"), k, mode="valid")[:len(x)]


pa, lo = pa_full(phs), ds_full(asl["lon"])

print("--- Do the two series share low-frequency variability? ---")
for name, x in (("ZW3 phase anomaly", pa), ("ASL longitude anomaly", lo)):
    t = np.arange(len(x), dtype=float)
    sl, _, r, p, se = stats.linregress(t, x)
    print(f"  {name:24s} trend = {sl*120:+.2f} deg/decade  p = {p:.3f}")
lf_pa, lf_lo = pa - highpass(pa), lo - highpass(lo)
print(f"  correlation of the two low-frequency components: "
      f"{np.corrcoef(lf_pa, lf_lo)[0,1]:+.3f}")

print("\n--- Slope under progressively stricter low-frequency removal ---")
treatments = (("raw calendar-month anomalies", pa, lo),
              ("linearly detrended", detrend(pa), detrend(lo)),
              ("high-pass (5-yr running mean removed)", highpass(pa), highpass(lo)))
for label, x, y in treatments:
    print(f"\n  {label}")
    for sname, mset in (("all months", None), ("DJF", (12, 1, 2)), ("JJA", (6, 7, 8))):
        s = np.ones(len(months), bool) if mset is None else np.isin(months, mset)
        sl, _, r, p, se = stats.linregress(x[s], y[s])
        line = f"    {sname:10s} n={s.sum():3d}  slope={sl:+.3f}+/-{se:.3f}  r={r:+.3f}"
        if mset is not None:
            t2 = np.percentile(amp[s], 66.7)
            st = s.copy()
            st[s] = amp[s] > t2
            s2, _, r2, _, se2 = stats.linregress(x[st], y[st])
            line += f"   | strong wave n={st.sum():2d} slope={s2:+.3f}+/-{se2:.3f}"
        print(line)
