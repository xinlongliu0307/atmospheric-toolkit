"""
Is the argmax node estimator degenerate?

The pre-registered H2 result depends on estimator A, which B does not
reproduce. Test whether A is a genuine position estimate or a boundary-
saturated near-binary indicator, and compare against a third, principled
estimator: the phase of a single sinusoid fitted to the sector anomaly
profile, which has no discrete search and no centroid shrinkage.
"""
import numpy as np
from windtools.stats import deseason as _deseason
import xarray as xr
from scipy import stats
from windtools.asl import asl_timeseries


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

cand = np.arange(10, len(lons) - 10)
node_a = np.empty(anom.shape[0])
for t in range(anom.shape[0]):
    node_a[t] = lons[cand[int(np.argmax(
        [anom[t, :j].mean() - anom[t, j:].mean() for j in cand]))]]

print("--- Is estimator A boundary-saturated? ---")
lo_edge, hi_edge = lons[cand[0]], lons[cand[-1]]
at_lo = np.isclose(node_a, lo_edge)
at_hi = np.isclose(node_a, hi_edge)
print(f"  candidate range {lo_edge:.0f}-{hi_edge:.0f}E (width {hi_edge-lo_edge:.0f} deg)")
print(f"  sd of raw node = {node_a.std():.1f} deg; uniform would be "
      f"{(hi_edge-lo_edge)/np.sqrt(12):.1f} deg")
print(f"  at LOWER boundary: {at_lo.sum():3d} months ({100*at_lo.mean():.1f}%)")
print(f"  at UPPER boundary: {at_hi.sum():3d} months ({100*at_hi.mean():.1f}%)")
print(f"  at either boundary: {100*(at_lo|at_hi).mean():.1f}%")
print("  deciles:", np.round(np.percentile(node_a, np.arange(10, 100, 10)), 0))

# Estimator C: phase of a single sinusoid fitted across the sector window.
span = lons[-1] - lons[0]
theta = 2 * np.pi * (lons - lons[0]) / span
cos_t, sin_t = np.cos(theta), np.sin(theta)
node_c = np.empty(anom.shape[0])
for t in range(anom.shape[0]):
    w = anom[t] - anom[t].mean()
    a1 = 2 * np.dot(w, cos_t) / len(lons)
    b1 = 2 * np.dot(w, sin_t) / len(lons)
    ph = np.arctan2(b1, a1) % (2 * np.pi)
    node_c[t] = lons[0] + span * ph / (2 * np.pi)


def node_b(a):
    out = np.empty(a.shape[0])
    for t in range(a.shape[0]):
        w = a[t] - a[t].mean()
        pos, neg = np.clip(w, 0, None), np.clip(-w, 0, None)
        out[t] = (lons.mean() if pos.sum() <= 0 or neg.sum() <= 0
                  else 0.5 * ((lons * pos).sum() / pos.sum()
                              + (lons * neg).sum() / neg.sum()))
    return out


nb = node_b(anom)

print("\n--- Do the three estimators agree with each other? ---")
for n1, n2, lab in ((node_a, nb, "A vs B"), (node_a, node_c, "A vs C"),
                    (nb, node_c, "B vs C")):
    print(f"  {lab}: r = {np.corrcoef(n1, n2)[0,1]:+.3f}")


def deseason(x, mask):
    return _deseason(x[mask], months[mask])


print("\n--- H2 under all three estimators ---")
samples = {"FULL": np.ones(len(months), bool),
           "1979-2001": years <= 2001, "2002-2025": years >= 2002}
for lab, nodes in (("A argmax  ", node_a), ("B centroid", nb),
                   ("C sinusoid", node_c)):
    line = f"  {lab} (sd {deseason(nodes, samples['FULL']).std():5.1f} deg): "
    for sname, mask in samples.items():
        sl, _, r, _, se = stats.linregress(deseason(asl["lon"], mask),
                                           deseason(nodes, mask))
        line += f" {sname} {sl:+.3f}+/-{se:.3f} |"
    print(line)

print("\n--- Is A behaving as a binary west/east flag? ---")
flag = (node_a > np.median(node_a)).astype(float)
r_pb = np.corrcoef(deseason(asl["lon"], samples["FULL"]),
                   deseason(flag, samples["FULL"]))[0, 1]
print(f"  corr(ASL lon anomaly, median-split flag of A) = {r_pb:+.3f}")
print("  (comparable to the full A correlation implies A carries little "
      "more than a west/east dichotomy)")
