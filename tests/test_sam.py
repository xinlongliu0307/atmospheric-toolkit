import numpy as np
import pytest
from windtools.sam import sam_index


def test_sam_index_normalisation_properties():
    # Any series produces a SAM index with ~zero mean by construction,
    # because both inputs are normalised before differencing.
    rng = np.random.default_rng(42)
    p40 = 1015.0 + rng.normal(0, 3, 120)
    p65 = 990.0 + rng.normal(0, 5, 120)
    sam = sam_index(p40, p65)
    assert sam.shape == (120,)
    assert abs(float(np.mean(sam))) < 1e-10


def test_sam_index_known_hand_case():
    # Three-point series chosen so the normalised values are exact.
    # p40 = [1014, 1015, 1016]: anomalies [-1, 0, 1], population std = sqrt(2/3).
    # p65 = [992, 990, 988]:    anomalies [2, 0, -2], population std = sqrt(8/3).
    # Normalised p40 = [-1.2247, 0, 1.2247]; normalised p65 = [1.2247, 0, -1.2247].
    # SAM = normalised p40 - normalised p65 = [-2.4495, 0, 2.4495].
    p40 = np.array([1014.0, 1015.0, 1016.0])
    p65 = np.array([992.0, 990.0, 988.0])
    sam = sam_index(p40, p65)
    expected = np.array([-2.44949, 0.0, 2.44949])
    assert sam == pytest.approx(expected, abs=1e-4)


def test_sam_positive_phase_sign_convention():
    # Positive SAM = anomalously high pressure at 40S and low at 65S.
    # The final point of the hand case has p40 above its mean and p65 below:
    # the index must be positive there.
    p40 = np.array([1014.0, 1015.0, 1016.0])
    p65 = np.array([992.0, 990.0, 988.0])
    assert sam_index(p40, p65)[-1] > 0


def test_sam_index_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        sam_index(np.zeros(10), np.zeros(9))
