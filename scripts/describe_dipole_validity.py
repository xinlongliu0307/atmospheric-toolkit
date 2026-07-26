"""
Descriptive only: how often does a locatable sea ice dipole exist?

NO relationship is tested here. The purpose is to establish, before any
hypothesis is written, whether the construct assumed by the withdrawn H2
is present in the data at all, and how often.

Also measures the two dipole_node properties currently documented from
estimate rather than measurement: monopole sensitivity and the noise
margin on the r_squared threshold.
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
profile = np.nanmean(sic.values[:, band, :], axis=1)[:, win]
months = sic[it].dt.month.values

anom = np.empty_like(profile)
for m in range(1, 13):
    s = months == m
    anom[s] = profile[s] - profile[s].mean(axis=0)

res = [dipole_node(lons, anom[t]) for t in range(anom.shape[0])]
valid = np.array([r["valid"] for r in res])
r2 = np.array([r["r_squared"] for r in res])
node = np.array([r["node"] for r in res])
contrast = np.array([r["contrast"] for r in res])

print(f"window {lons[0]:.0f}-{lons[-1]:.0f}E, band 60-75S, "
      f"{len(months)} months\n")
print("=== How often is a dipole present? ===")
print(f"  valid months: {valid.sum()} of {len(valid)} "
      f"({100 * valid.mean():.1f}%)")
print(f"  r_squared: median {np.median(r2):.3f}, "
      f"25th {np.percentile(r2, 25):.3f}, 75th {np.percentile(r2, 75):.3f}, "
      f"max {r2.max():.3f}")

print("\n=== Valid fraction by month ===")
for m in range(1, 13):
    s = months == m
    print(f"  {m:02d}: {100 * valid[s].mean():5.1f}%   "
          f"median r2 {np.median(r2[s]):.3f}")

print("\n=== Node distribution, valid months only ===")
if valid.sum() >= 10:
    nv = node[valid]
    print(f"  n = {len(nv)}, mean {nv.mean():.1f}E, sd {nv.std():.1f} deg")
    print(f"  deciles: {np.round(np.percentile(nv, np.arange(10, 100, 10)), 0)}")
    edge = ((nv <= lons[0] + 5) | (nv >= lons[-1] - 5)).mean()
    print(f"  within 5 deg of a window edge: {100 * edge:.1f}% "
          "(argmax estimator was 43.4%)")
    print(f"  mean contrast: {contrast[valid].mean():.4f} "
          "(SIC fraction, peak to peak)")
else:
    print("  too few valid months to characterise")

print("\n=== Deferred measurements: dipole_node properties ===")
Lm = lons
Cm = 0.5 * (Lm[0] + Lm[-1])
print("  monopole sensitivity (centred Gaussian):")
for w in (8.0, 15.0, 20.0, 30.0, 40.0, 60.0):
    r = dipole_node(Lm, np.exp(-((Lm - Cm) / w) ** 2))
    print(f"    width {w:4.0f}: r2 = {r['r_squared']:.3f}  valid = {r['valid']}")

rng = np.random.default_rng(0)
nr2 = np.array([dipole_node(Lm, rng.normal(0, 1, len(Lm)))["r_squared"]
                for _ in range(3000)])
print(f"  noise r2 over 3000 draws: mean {nr2.mean():.4f}, "
      f"95th {np.percentile(nr2, 95):.3f}, max {nr2.max():.3f}, "
      f"fraction >= 0.5: {(nr2 >= 0.5).mean():.4f}")
