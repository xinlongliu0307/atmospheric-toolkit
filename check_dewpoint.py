import sys
from windtools.humidity import dewpoint

# For T=25 C, RH=50 %, the Magnus dewpoint is approximately 13.87 C.
result = dewpoint(25.0, 50.0)
if abs(result - 13.87) > 0.1:
    print(f"FAIL: dewpoint(25, 50) = {result:.2f}, expected ~13.87")
    sys.exit(1)
print(f"PASS: dewpoint(25, 50) = {result:.2f}")
sys.exit(0)
