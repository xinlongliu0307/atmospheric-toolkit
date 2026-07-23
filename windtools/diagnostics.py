from windtools.speed import wind_speed, wind_direction
from windtools.thermo import theta
from windtools.asl import find_asl
from windtools.zw3 import zw3_index

DIAGNOSTICS = {
    "wind_speed": wind_speed,
    "wind_direction": wind_direction,
    "theta": theta,
    "asl": find_asl,
    "zw3": zw3_index,
}
