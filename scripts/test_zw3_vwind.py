"""
Pre-registered sensitivity test (PREREGISTRATION_zw3_vwind.md).
Order is deliberate: analytic check, then degeneracy diagnostic, then the
single primary test.
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


def circle_at(path, var, lat, scale=1.0):
    ds = xr.open_dataset(path)
    la, lo, t = coords(ds)
    v = (ds[var] * scale).squeeze(drop=True) if var in ds else None
    if v is None:
        var = list(ds.data_vars)[0]
        v = (ds[var] * scale).squeeze(drop=True)
    c = v.sel({la: lat}, method="nearest")
    c = c.assign_coords({lo: c[lo] % 360.0}).sortby(lo)
    u, first = np.unique(c[lo].values, return_index=True)
    if len(u) != c.sizes[lo]:
        c = c.isel({lo: np.sort(first)})
    c = c.transpose(t, lo)
    return c[lo].values, c.values, c[t]


zl, zv, ztime = circle_at("data/era5_z500_monthly_40S70S.nc", "z", -49.0, 1 / G)
vl, vv, vtime = circle_at("data/era5_v500_monthly_40S70S.nc", "v", -49.0)
assert np.array_equal(ztime.values, vtime.values), "time axes differ"

phz = np.array([zw3_index(zl, zv[t])["phase"] for t in range(zv.shape[0])])
phv = np.array([zw3_index(vl, vv[t])["phase"] for t in range(vv.shape[0])])

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
    a = np.deg2rad(p * 3.0)
    return (np.rad2deg(np.angle(np.mean(np.exp(1j * a)))) / 3.0) % 120.0


def phase_anom(p, mask):
    mm, ps = months[mask], p[mask]
    out = np.empty(mask.sum())
    for k in range(1, 13):
        s = mm == k
        out[s] = (ps[s] - circ_mean_120(ps[s]) + 60.0) % 120.0 - 60.0
    return out


def deseason(x, mask):
    mm, xs = months[mask], x[mask]
    out = np.empty(mask.sum())
    for k in range(1, 13):
        s = mm == k
        out[s] = xs[s] - xs[s].mean()
    return out


def neff(x, y):
    rx = np.corrcoef(x[:-1], x[1:])[0, 1]
    ry = np.corrcoef(y[:-1], y[1:])[0, 1]
    n = len(x)
    return min(max(3.0, n * (1 - rx * ry) / (1 + rx * ry)), float(n))


print("=== P1. Analytic check: quadrature between v and Z phases ===")
diff = (phv - phz) % 120.0
off = circ_mean_120(diff)
R = np.abs(np.mean(np.exp(1j * np.deg2rad(diff * 3.0))))
print(f"  circular-mean offset = {off:.1f} deg  (predicted +30; "
      f"30 or 90 = quadrature)")
print(f"  concentration R = {R:.3f}")
print(f"  verdict: {'QUADRATURE' if min(abs(off-30), abs(off-90)) < 15 else 'NOT QUADRATURE'}")

full = np.ones(len(months), bool)
paz, pav = phase_anom(phz, full), phase_anom(phv, full)
print("\n=== D1. Degeneracy diagnostic (reported before the test) ===")
rd = np.corrcoef(paz, pav)[0, 1]
print(f"  corr(phase anomaly Z, phase anomaly v) = {rd:+.3f}")
print(f"  status: {'NEAR-DEGENERATE, test uninformative' if abs(rd) > 0.95 else 'informative'}")

print("\n=== H1. Primary test: ZW3(v) phase -> ASL longitude ===")
samples = {"FULL 1979-2025": full, "half 1979-2001": years <= 2001,
           "half 2002-2025": years >= 2002}
for name, mask in samples.items():
    x, y = phase_anom(phv, mask), deseason(asl["lon"], mask)
    sl, _, r, _, se = stats.linregress(x, y)
    ne = neff(x, y)
    t = r * np.sqrt((ne - 2) / max(1e-12, 1 - r ** 2))
    p = 2 * stats.t.sf(abs(t), df=ne - 2)
    lo95, hi95 = sl - 1.96 * se, sl + 1.96 * se
    print(f"  {name:16s} slope = {sl:+.3f} +/- {se:.3f}  "
          f"95% [{lo95:+.3f}, {hi95:+.3f}]  r = {r:+.3f}  p = {p:.2e}")
print("  (Z-based reference: +0.328 full, +0.308 and +0.207 in halves)")
