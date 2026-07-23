"""Download ERA5 monthly-mean 500 hPa geopotential, 40-70S, 1 degree."""
import cdsapi

c = cdsapi.Client()
c.retrieve(
    "reanalysis-era5-pressure-levels-monthly-means",
    {
        "product_type": "monthly_averaged_reanalysis",
        "variable": "geopotential",
        "pressure_level": "500",
        "year": [str(y) for y in range(1979, 2026)],
        "month": [f"{m:02d}" for m in range(1, 13)],
        "time": "00:00",
        "grid": [1.0, 1.0],
        "area": [-40, -180, -70, 180],  # N, W, S, E
        "format": "netcdf",
    },
    "data/era5_z500_monthly_40S70S.nc",
)
print("Done: data/era5_z500_monthly_40S70S.nc")
