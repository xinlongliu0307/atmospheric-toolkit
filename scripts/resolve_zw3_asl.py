"""
Resolve the tension between the seasonal and amplitude-tercile results,
and test whether ZW3-ASL coupling depends on their meridional separation.

  A. Amplitude terciles on DESEASONALISED amplitude, so terciles carry no
     seasonal information. Settles whether dilution or season dominates.
  B. Within-season terciles for DJF and JJA, the two extreme seasons.
  C. Robustness with sector-edge months excluded.
  D. ZW3 computed at several latitudes: does coupling peak where the wave
     is closest to the ASL centre, and does that explain the seasonality?
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
zfield = (dz[zvar] / G).squeeze(drop=True)
zfield = zfield.assign_coords({zlo: zfield[zlo] % 360.0}).sortby(zlo)
vals, first = np.unique(zfield[zlo].values, return_index=True)
if len(vals) != zfield.sizes[zlo]:
    zfield = zfield.isel({zlo: np.sort(first)})
zfield = zfield.transpose(zt, zla, zlo)
zlons = zfield[zlo].values


def zw3_at(lat):
    row = zfield.sel({zla: lat}, method="nearest").values
    a = np.empty(row.shape[0])
    p = np.empty(row.shape[0])
    for t in range(row.shape[0]):
        r = zw3_index(zlons, row[t])
        a[t], p[t] = r["amplitude"], r["phase"]
    return a, p


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


def phase_anom(p, mask=None):
    return phase_anomaly(p, months if mask is None else months[mask], 120.0)


amp49, phs49 = zw3_at(-49.0)
pa, lona = phase_anom(phs49), deseason(asl["lon"])
ampa = deseason(amp49)

print("--- A. Terciles on DESEASONALISED amplitude ---")
q1, q2 = np.percentile(ampa, [33.3, 66.7])
for name, s in (("weak  ", ampa <= q1), ("medium", (ampa > q1) & (ampa <= q2)),
                ("strong", ampa > q2)):
    sl, _, r, _, se = stats.linregress(pa[s], lona[s])
    frac_djf = np.isin(months[s], (12, 1, 2)).mean()
    frac_jja = np.isin(months[s], (6, 7, 8)).mean()
    print(f"  {name} (n={s.sum():3d}, DJF {100*frac_djf:.0f}%, JJA {100*frac_jja:.0f}%): "
          f"slope = {sl:+.3f} +/- {se:.3f}, r = {r:+.3f}")

print("\n--- B. Within-season terciles (raw amplitude) ---")
for sname, mm in (("DJF", (12, 1, 2)), ("JJA", (6, 7, 8))):
    sel = np.isin(months, mm)
    a_s, pa_s, lo_s = amp49[sel], pa[sel], lona[sel]
    t1, t2 = np.percentile(a_s, [33.3, 66.7])
    print(f"  {sname}:")
    for label, m2 in (("weak  ", a_s <= t1), ("strong", a_s > t2)):
        sl, _, r, _, se = stats.linregress(pa_s[m2], lo_s[m2])
        print(f"    {label} (n={m2.sum():3d}): slope = {sl:+.3f} +/- {se:.3f}, r = {r:+.3f}")

print("\n--- C. Robustness: sector-edge months excluded ---")
lon = asl["lon"]
keep = ~(np.isclose(lon, lon.min()) | np.isclose(lon, lon.max()))
sl, _, r, _, se = stats.linregress(pa[keep], lona[keep])
print(f"  n = {keep.sum()} of {len(keep)}: slope = {sl:+.3f} +/- {se:.3f}, r = {r:+.3f}")

print("\n--- D. Coupling vs ZW3 reference latitude ---")
print("  lat    all-season slope        DJF slope         JJA slope")
for lat in (-45.0, -49.0, -55.0, -60.0, -65.0, -70.0):
    a_l, p_l = zw3_at(lat)
    pal = phase_anom(p_l)
    out = [f"  {abs(lat):.0f}S "]
    for mm in (None, (12, 1, 2), (6, 7, 8)):
        s = np.ones(len(months), bool) if mm is None else np.isin(months, mm)
        sl, _, r, _, se = stats.linregress(pal[s], lona[s])
        out.append(f"  {sl:+.3f}+/-{se:.3f} (r{r:+.2f})")
    print("".join(out))
print("\n  ASL centre latitude: DJF %.1fS, JJA %.1fS"
      % (-asl["lat"][np.isin(months, (12, 1, 2))].mean(),
         -asl["lat"][np.isin(months, (6, 7, 8))].mean()))
