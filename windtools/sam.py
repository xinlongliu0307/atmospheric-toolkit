import numpy as np
def sam_index(p40: np.ndarray, p65: np.ndarray) -> np.ndarray:
    """
    Calculate the Southern Annular Mode (SAM) index based on Gong and Wang (1999).

    Parameters:
        p40 (np.ndarray): Zonal-mean sea-level pressure at 40S.
        p65 (np.ndarray): Zonal-mean sea-level pressure at 65S.

    Returns:
        np.ndarray: SAM index as a numpy array.

    References:
        Gong, D. Y., & Wang, S. L. (1999). What drives the Southern Oscillation?
        Journal of Climate, 12(7), 2065-2085.
    """
    if len(p40) != len(p65):
        raise ValueError("The input arrays must have the same length.")

    p40_norm = (p40 - np.mean(p40)) / np.std(p40)
    p65_norm = (p65 - np.mean(p65)) / np.std(p65)

    return p40_norm - p65_norm
