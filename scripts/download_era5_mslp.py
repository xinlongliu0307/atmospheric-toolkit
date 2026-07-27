"""Download ERA5 monthly-mean sea level pressure, 40-80S, 1 degree."""
import cdsapi

c = cdsapi.Client()
c.retrieve(
    "reanalysis-era5-single-levels-monthly-means",
    {
        "product_type": "monthly_averaged_reanalysis",
        "variable": "mean_sea_level_pressure",
        "year": [str(y) for y in range(1979, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": "00:00",
        "grid": [1.0, 1.0],
        "area": [-40, -180, -80, 180],
        "format": "netcdf",
    },
    "data/era5_mslp_monthly_40S80S.nc",
)
print("Done: data/era5_mslp_monthly_40S80S.nc")
