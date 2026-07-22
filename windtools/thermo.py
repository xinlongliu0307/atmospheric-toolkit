def theta(pressure_hpa, temperature_k):
    return temperature_k * (1000 / pressure_hpa) ** 0.2854
