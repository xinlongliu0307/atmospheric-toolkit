import numpy as np


def parse_nino34(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse the NOAA PSL Nino 3.4 monthly text format.

    The first line gives the start and end year and is ignored, as are the
    trailing sentinel and description lines. A data row is a four-digit
    year followed by exactly twelve numeric values. Entries equal to the
    missing-value sentinel are dropped.

    Note that the values are absolute sea surface temperatures in degrees
    Celsius, not anomalies, so they carry a strong annual cycle.

    Parameters:
        text: full contents of the downloaded file.

    Returns:
        tuple of (years, months, sst), aligned elementwise and in
        chronological order, with years and months as int and sst as float.

    Raises:
        ValueError: if no data rows are found.
    """
    years, months, sst = [], [], []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) != 13:
            continue
        head = parts[0]
        if not (len(head) == 4 and head.isdigit()):
            continue
        try:
            values = [float(v) for v in parts[1:]]
        except ValueError:
            continue
        year = int(head)
        for month, value in enumerate(values, start=1):
            if value <= -99.0:
                continue
            years.append(year)
            months.append(month)
            sst.append(value)

    if not years:
        raise ValueError("no data rows found in the Nino 3.4 text")

    years = np.array(years, dtype=int)
    months = np.array(months, dtype=int)
    sst = np.array(sst, dtype=float)
    order = np.argsort(years * 12 + months, kind="stable")
    return years[order], months[order], sst[order]
