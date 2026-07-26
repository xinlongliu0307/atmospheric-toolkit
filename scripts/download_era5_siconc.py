"""Download ERA5 monthly-mean sea ice cover, 50-80S, 1 degree."""
import cdsapi

c = cdsapi.Client()
c.retrieve(
    "reanalysis-era5-single-levels-monthly-means",
    {
        "product_type": "monthly_averaged_reanalysis",
        "variable": "sea_ice_cover",
        "year": [str(y) for y in range(1979, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": "00:00",
        "grid": [1.0, 1.0],
        "area": [-50, -180, -80, 180],
        "format": "netcdf",
    },
    "data/era5_siconc_monthly_50S80S.nc",
)
print("Done: data/era5_siconc_monthly_50S80S.nc")
