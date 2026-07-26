"""
Robustness of H2 to the sea ice dipole node estimator, and an empirical
power assessment for a possible ZW3-phase chain test.

Estimator A (as used in the pre-registered analysis): argmax over
candidate longitudes of the west-minus-east mean contrast.
Estimator B: sign-weighted centroid of the longitudinal SIC anomaly
profile, a smooth alternative with no discrete search.

The ZW3 direct effect reported here is a POWER DIAGNOSTIC, not a
hypothesis test. It was not pre-registered and must not be interpreted
as evidence either way.
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


def prep(ds, var, scale=1.0):
    la, lo, t = coords(ds)
    f = (ds[var] * scale).squeeze(drop=True)
    f = f.assign_coords({lo: f[lo] % 360.0}).sortby(lo)
    v, first = np.unique(f[lo].values, return_index=True)
    if len(v) != f.sizes[lo]:
        f = f.isel({lo: np.sort(first)})
    return f.sortby(la).transpose(t, la, lo), la, lo, t


dm = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
pvar = "msl" if "msl" in dm else list(dm.data_vars)[0]
slp, mla, mlo, mt = prep(dm, pvar, 1 / 100.0)
asl = asl_timeseries(slp[mlo].values, slp[mla].values, slp.values)
months = slp[mt].dt.month.values
years = slp[mt].dt.year.values

dz = xr.open_dataset("data/era5_z500_monthly_40S70S.nc")
zla, zlo, zt = coords(dz)
zvar = "z" if "z" in dz else list(dz.data_vars)[0]
c = (dz[zvar].sel({zla: -49.0}, method="nearest") / G).squeeze(drop=True)
c = c.assign_coords({zlo: c[zlo] % 360.0}).sortby(zlo)
v, f = np.unique(c[zlo].values, return_index=True)
if len(v) != c.sizes[zlo]:
    c = c.isel({zlo: np.sort(f)})
c = c.transpose(zt, zlo)
phs = np.array([zw3_index(c[zlo].values, c.values[t])["phase"]
                for t in range(c.shape[0])])

di = xr.open_dataset("data/era5_siconc_monthly_50S80S.nc")
ivar = "siconc" if "siconc" in di else ("ci" if "ci" in di else list(di.data_vars)[0])
sic, ila, ilo, it = prep(di, ivar)
ilats, ilons = sic[ila].values, sic[ilo].values
band = (ilats >= -75.0) & (ilats <= -60.0)
sic_lon = np.nanmean(sic.values[:, band, :], axis=1)
m = (ilons >= 160.0) & (ilons <= 290.0)
lons, sub = ilons[m], sic_lon[:, m]

anom = np.empty_like(sub)
for k in range(1, 13):
    s = months == k
    anom[s] = sub[s] - sub[s].mean(axis=0)


def node_argmax(a):
    cand = np.arange(10, len(lons) - 10)
    out = np.empty(a.shape[0])
    for t in range(a.shape[0]):
        contrast = [a[t, :j].mean() - a[t, j:].mean() for j in cand]
        out[t] = lons[cand[int(np.argmax(contrast))]]
    return out


def node_centroid(a):
    """Sign-weighted centroid: midpoint between positive and negative masses."""
    out = np.empty(a.shape[0])
    for t in range(a.shape[0]):
        w = a[t] - a[t].mean()
        pos, neg = np.clip(w, 0, None), np.clip(-w, 0, None)
        if pos.sum() <= 0 or neg.sum() <= 0:
            out[t] = lons.mean()
            continue
        out[t] = 0.5 * ((lons * pos).sum() / pos.sum() + (lons * neg).sum() / neg.sum())
    return out


def deseason(x, mask):
    mm, xs = months[mask], x[mask]
    out = np.empty(mask.sum())
    for k in range(1, 13):
        s = mm == k
        out[s] = xs[s] - xs[s].mean()
    return out


def circ_mean_120(p):
    a = np.deg2rad(p * 3.0)
    return (np.rad2deg(np.angle(np.mean(np.exp(1j * a)))) / 3.0) % 120.0


def phase_anom(mask):
    mm, ps = months[mask], phs[mask]
    out = np.empty(mask.sum())
    for k in range(1, 13):
        s = mm == k
        out[s] = (ps[s] - circ_mean_120(ps[s]) + 60.0) % 120.0 - 60.0
    return out


def neff(x, y):
    rx = np.corrcoef(x[:-1], x[1:])[0, 1]
    ry = np.corrcoef(y[:-1], y[1:])[0, 1]
    n = len(x)
    return min(max(3.0, n * (1 - rx * ry) / (1 + rx * ry)), float(n))


samples = {"FULL 1979-2025": np.ones(len(months), bool),
           "half 1979-2001": years <= 2001,
           "half 2002-2025": years >= 2002}
estimators = {"A argmax  ": node_argmax(anom), "B centroid": node_centroid(anom)}

print("--- H2 robustness: ASL longitude -> ice dipole node ---")
for ename, nodes in estimators.items():
    print(f"  {ename}  (sd of node anomaly = {deseason(nodes, samples['FULL 1979-2025']).std():.1f} deg)")
    for sname, mask in samples.items():
        x, y = deseason(asl["lon"], mask), deseason(nodes, mask)
        sl, _, r, _, se = stats.linregress(x, y)
        ne = neff(x, y)
        t = r * np.sqrt((ne - 2) / max(1e-12, 1 - r ** 2))
        p = 2 * stats.t.sf(abs(t), df=ne - 2)
        print(f"      {sname:16s} slope = {sl:+.3f} +/- {se:.3f}  "
              f"r = {r:+.3f}  p = {p:.2e}")

print("\n--- Power diagnostic for a ZW3 chain test (NOT a hypothesis test) ---")
full = samples["FULL 1979-2025"]
pa = phase_anom(full)
for ename, nodes in estimators.items():
    y = deseason(nodes, full)
    sl_link2, _, _, _, _ = stats.linregress(deseason(asl["lon"], full), y)
    sl_link1 = 0.30  # replicated ZW3 phase -> ASL longitude slope
    predicted = sl_link1 * sl_link2
    exp_r = predicted * pa.std() / y.std()
    ne = neff(pa, y)
    detectable_r = 1.96 / np.sqrt(ne)
    obs_sl, _, obs_r, _, obs_se = stats.linregress(pa, y)
    print(f"  {ename}: predicted chain slope = {predicted:.3f} "
          f"-> expected r = {exp_r:.3f}")
    print(f"              detectable |r| at n_eff {ne:.0f} = {detectable_r:.3f}"
          f"   -> {'ADEQUATE' if abs(exp_r) > detectable_r else 'UNDERPOWERED'}")
    print(f"              observed: slope {obs_sl:+.3f} +/- {obs_se:.3f}, r {obs_r:+.3f}")
