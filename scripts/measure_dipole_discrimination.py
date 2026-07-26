"""
Descriptive only. Two measurements, no relationship tested.

1. RED-NOISE NULL. White noise is the wrong null for a wavenumber-1
   criterion on a spatially autocorrelated field. Generate AR(1)-in-
   longitude profiles matched to the observed decorrelation and measure
   how often they pass r_squared >= 0.5.

2. ANTISYMMETRY. A dipole is antisymmetric about its node; a monopole is
   symmetric about its centre. Measure the antisymmetry statistic for the
   valid observed months, for monopoles, and for red noise, and see
   whether it separates them.

   A = -sum_d y(n+d) y(n-d) / mean(sum_d y(n+d)^2 + y(n-d)^2)
   A = +1 perfectly antisymmetric, -1 perfectly symmetric.
"""
import numpy as np
import xarray as xr
from windtools.dipole import dipole_node


def coords(ds):
    return ("latitude" if "latitude" in ds.coords else "lat",
            "longitude" if "longitude" in ds.coords else "lon",
            "valid_time" if "valid_time" in ds.coords else "time")


di = xr.open_dataset("data/era5_siconc_monthly_50S80S.nc")
ila, ilo, it = coords(di)
ivar = "siconc" if "siconc" in di else ("ci" if "ci" in di else list(di.data_vars)[0])
sic = di[ivar].squeeze(drop=True)
sic = sic.assign_coords({ilo: sic[ilo] % 360.0}).sortby(ilo)
v, first = np.unique(sic[ilo].values, return_index=True)
if len(v) != sic.sizes[ilo]:
    sic = sic.isel({ilo: np.sort(first)})
sic = sic.sortby(ila).transpose(it, ila, ilo)
lats, lons_all = sic[ila].values, sic[ilo].values
band = (lats >= -75.0) & (lats <= -60.0)
win = (lons_all >= 160.0) & (lons_all <= 290.0)
lons = lons_all[win]
prof = np.nanmean(sic.values[:, band, :], axis=1)[:, win]
months = sic[it].dt.month.values

anom = np.empty_like(prof)
for m in range(1, 13):
    s = months == m
    anom[s] = prof[s] - prof[s].mean(axis=0)

N = len(lons)


def antisymmetry(y, node_lon):
    n = int(np.argmin(np.abs(lons - node_lon)))
    d = np.arange(1, min(n, N - 1 - n) + 1)
    if len(d) < 5:
        return np.nan
    w = y - y.mean()
    a, b = w[n + d], w[n - d]
    denom = 0.5 * np.sum(a ** 2 + b ** 2)
    return float(-np.sum(a * b) / denom) if denom > 0 else np.nan


# observed lag-1 autocorrelation in longitude
lag1 = np.array([np.corrcoef(r[:-1], r[1:])[0, 1]
                 for r in anom if np.std(r) > 0])
rho = float(np.median(lag1))
print(f"observed profiles: median lag-1 autocorrelation in longitude = {rho:.3f}")
print(f"implied decorrelation length ~ {-1.0 / np.log(max(rho, 1e-6)):.1f} degrees\n")


def ar1_profile(rng, rho, n):
    e = rng.normal(0, 1, n)
    y = np.empty(n)
    y[0] = e[0]
    for i in range(1, n):
        y[i] = rho * y[i - 1] + np.sqrt(1 - rho ** 2) * e[i]
    return y


print("=== 1. Red-noise null ===")
rng = np.random.default_rng(0)
for r_test in (0.0, 0.5, 0.8, 0.9, 0.95, rho):
    res = [dipole_node(lons, ar1_profile(rng, r_test, N)) for _ in range(1000)]
    frac = np.mean([x["valid"] for x in res])
    med = np.median([x["r_squared"] for x in res])
    tag = "  <- matched to data" if abs(r_test - rho) < 1e-9 else ""
    print(f"  rho = {r_test:.3f}: valid {100 * frac:5.1f}%, "
          f"median r2 {med:.3f}{tag}")

print("\n=== 2. Antisymmetry ===")
res = [dipole_node(lons, anom[t]) for t in range(anom.shape[0])]
valid = np.array([r["valid"] for r in res])
obs_a = np.array([antisymmetry(anom[t], res[t]["node"])
                  if res[t]["valid"] else np.nan for t in range(len(res))])
obs_a = obs_a[~np.isnan(obs_a)]
print(f"  observed valid months (n = {len(obs_a)}):")
print(f"    mean {obs_a.mean():+.3f}, median {np.median(obs_a):+.3f}, "
      f"sd {obs_a.std():.3f}")
print(f"    deciles {np.round(np.percentile(obs_a, np.arange(10, 100, 10)), 2)}")
print(f"    fraction with A > 0.5: {100 * (obs_a > 0.5).mean():.1f}%")

C = 0.5 * (lons[0] + lons[-1])
print("\n  reference: monopoles (should be strongly NEGATIVE)")
for w in (20.0, 30.0, 40.0):
    y = np.exp(-((lons - C) / w) ** 2)
    r = dipole_node(lons, y)
    print(f"    width {w:4.0f}: r2 {r['r_squared']:.3f}, "
          f"A = {antisymmetry(y, r['node']):+.3f}")

print("\n  reference: ideal and smoothed step dipoles (should be near +1)")
for nl in (200.0, 225.0, 250.0):
    y = np.tanh((nl - lons) / 8.0)
    r = dipole_node(lons, y)
    print(f"    node {nl:5.0f}E: r2 {r['r_squared']:.3f}, "
          f"A = {antisymmetry(y, r['node']):+.3f}")

print("\n  reference: matched red noise")
rng = np.random.default_rng(1)
ra = []
for _ in range(1000):
    y = ar1_profile(rng, rho, N)
    r = dipole_node(lons, y)
    if r["valid"]:
        a = antisymmetry(y, r["node"])
        if not np.isnan(a):
            ra.append(a)
ra = np.array(ra)
if len(ra):
    print(f"    n = {len(ra)}, mean {ra.mean():+.3f}, "
          f"median {np.median(ra):+.3f}, "
          f"fraction A > 0.5: {100 * (ra > 0.5).mean():.1f}%")
