"""
Pre-registered analysis: ASL, ZW3, and Antarctic sea ice.
Specification: PREREGISTRATION_seaice.md (committed 7a15867, before any
sea ice data was examined).

H1  deeper ASL relative central pressure -> Ross ice gain, Bellingshausen
    ice loss. Test: corr(RCP anomaly, Ross-minus-Bellingshausen dipole).
    Directional prediction: NEGATIVE (RCP more negative = deeper).
H2  eastward ASL longitude -> eastward shift of the dipole node.
    Test: regress node-longitude anomaly on ASL-longitude anomaly.
    Directional prediction: POSITIVE slope, comparable to 1.
H3  SAM's ice signature is circumpolar, not sector-confined.
    Test: corr(SAM, circumpolar-mean SIC anomaly). Prediction: significant.
    Supporting (descriptive, not a primary test): pattern correlation
    between the SAM and RCP longitude-regression profiles should be low.

Caveat on provenance: ERA5 sea ice cover is prescribed from external
observational products, so it is observationally grounded but not
independent of the reanalysis framework. NSIDC passive-microwave
concentration would be the stricter choice for publication.
"""
import numpy as np
from windtools.stats import deseason as _deseason, corr_neff
import xarray as xr
from scipy import stats
from windtools.zw3 import zw3_index
from windtools.asl import asl_timeseries
from windtools.sam import sam_index

G = 9.80665
ALPHA = 0.017  # Bonferroni, three primary tests


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


# ---------- predictors ----------
dm = xr.open_dataset("data/era5_mslp_monthly_40S80S.nc")
pvar = "msl" if "msl" in dm else list(dm.data_vars)[0]
slp, mla, mlo, mt = prep(dm, pvar, 1 / 100.0)
mlons, mlats = slp[mlo].values, slp[mla].values
asl = asl_timeseries(mlons, mlats, slp.values)
outside = (mlons < 170.0) | (mlons > 298.0)
sam = sam_index(slp.sel({mla: -40.0}, method="nearest").values[:, outside].mean(axis=1),
                slp.sel({mla: -65.0}, method="nearest").values[:, outside].mean(axis=1))
months = slp[mt].dt.month.values
years = slp[mt].dt.year.values

# ---------- sea ice ----------
di = xr.open_dataset("data/era5_siconc_monthly_50S80S.nc")
ivar = "siconc" if "siconc" in di else ("ci" if "ci" in di else list(di.data_vars)[0])
sic, ila, ilo, it = prep(di, ivar)
assert sic.sizes[it] == len(months), "time axes differ"
ilats, ilons = sic[ila].values, sic[ilo].values
band = (ilats >= -75.0) & (ilats <= -60.0)
sic_lon = np.nanmean(sic.values[:, band, :], axis=1)   # (time, lon)
print(f"sea ice variable '{ivar}', band 60-75S: {band.sum()} lats, "
      f"{len(ilons)} lons, {sic.sizes[it]} months")


def sector(lo0, lo1):
    m = (ilons >= lo0) & (ilons <= lo1)
    return np.nanmean(sic_lon[:, m], axis=1)


ross = sector(160.0, 230.0)
amundsen = sector(230.0, 260.0)
bell = sector(260.0, 290.0)
circumpolar = np.nanmean(sic_lon, axis=1)


def deseason(x, mask):
    return _deseason(x[mask], months[mask])


def node_longitude(mask):
    """Longitude maximising (mean SIC anomaly west) - (mean east), 160-290E."""
    m = (ilons >= 160.0) & (ilons <= 290.0)
    lons = ilons[m]
    mm = months[mask]
    sub = sic_lon[:, m][mask]
    anom = np.empty_like(sub)
    for k in range(1, 13):
        s = mm == k
        anom[s] = sub[s] - sub[s].mean(axis=0)
    nodes = np.empty(anom.shape[0])
    cand = np.arange(10, len(lons) - 10)
    for t in range(anom.shape[0]):
        contrast = [anom[t, :c].mean() - anom[t, c:].mean() for c in cand]
        nodes[t] = lons[cand[int(np.argmax(contrast))]]
    return nodes - nodes.mean()


def corr(x, y):
    return corr_neff(x, y)


samples = {"FULL 1979-2025": np.ones(len(months), bool),
           "half 1979-2001": years <= 2001,
           "half 2002-2025": years >= 2002}

print("\n" + "=" * 66)
print("H1  deeper RCP -> Ross gain / Bellingshausen loss   (predict r < 0)")
print("=" * 66)
for name, mask in samples.items():
    d = deseason(ross, mask) - deseason(bell, mask)
    r, n, p = corr(deseason(asl["relative_central_pressure"], mask), d)
    verdict = "SUPPORTED" if (r < 0 and p < ALPHA) else "not supported"
    print(f"  {name:16s} r = {r:+.3f}  n_eff = {n:5.0f}  p = {p:.2e}   {verdict}")

print("\n" + "=" * 66)
print("H2  eastward ASL -> eastward dipole node   (predict slope > 0, ~1)")
print("=" * 66)
for name, mask in samples.items():
    x = deseason(asl["lon"], mask)
    y = node_longitude(mask)
    sl, _, r, p, se = stats.linregress(x, y)
    _, n, pa = corr(x, y)
    verdict = "SUPPORTED" if (sl > 0 and pa < ALPHA) else "not supported"
    print(f"  {name:16s} slope = {sl:+.3f} +/- {se:.3f}  r = {r:+.3f}  "
          f"p = {pa:.2e}   {verdict}")

print("\n" + "=" * 66)
print("H3  SAM signature is circumpolar   (predict significant)")
print("=" * 66)
for name, mask in samples.items():
    r, n, p = corr(deseason(sam, mask), deseason(circumpolar, mask))
    verdict = "SUPPORTED" if p < ALPHA else "not supported"
    print(f"  {name:16s} r = {r:+.3f}  n_eff = {n:5.0f}  p = {p:.2e}   {verdict}")

# ---------- descriptive support for H3 ----------
full = samples["FULL 1979-2025"]
sam_a = deseason(sam, full)
rcp_a = deseason(asl["relative_central_pressure"], full)
prof_sam, prof_rcp = np.empty(len(ilons)), np.empty(len(ilons))
for j in range(len(ilons)):
    col = np.empty(len(months))
    for k in range(1, 13):
        s = months == k
        col[s] = sic_lon[s, j] - np.nanmean(sic_lon[s, j])
    col = np.nan_to_num(col)
    prof_sam[j] = np.corrcoef(sam_a, col)[0, 1] if np.std(col) > 0 else 0.0
    prof_rcp[j] = np.corrcoef(rcp_a, col)[0, 1] if np.std(col) > 0 else 0.0

print("\n--- Descriptive (not a primary test): longitude profiles ---")
print(f"  pattern correlation SAM vs RCP profiles: "
      f"{np.corrcoef(prof_sam, prof_rcp)[0,1]:+.3f}  (low = distinguishable)")
print(f"  SAM profile: mean {prof_sam.mean():+.3f}, sd {prof_sam.std():.3f}")
print(f"  RCP profile: mean {prof_rcp.mean():+.3f}, sd {prof_rcp.std():.3f}")
print("  (circumpolar signature = large |mean|, small sd; "
      "sector-confined = small |mean|, large sd)")
print("\n  correlation with SIC by 30-degree longitude bin:")
print("     bin        SAM     RCP")
for lo0 in range(150, 300, 30):
    m = (ilons >= lo0) & (ilons < lo0 + 30)
    if m.sum():
        print(f"   {lo0:3d}-{lo0+30:3d}E   {prof_sam[m].mean():+.3f}  {prof_rcp[m].mean():+.3f}")
