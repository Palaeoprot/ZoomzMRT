"""
Unit tests for Waters MRT resolution-aware peak centroiding.
"""

import numpy as np
import pytest

from zoomzmrt.mrt_reader import _high_res_centroid, _high_res_profile_bin


def test_high_res_centroid_synthetic_peaks():
    """Test resolution-aware centroiding on two baseline resolved peaks."""
    # Peak 1 at 1570.6774, Peak 2 at 1570.6930
    mz1 = np.linspace(1570.670, 1570.684, 15)
    int1 = np.exp(-0.5 * ((mz1 - 1570.6774) / 0.002) ** 2) * 10000.0

    mz2 = np.linspace(1570.686, 1570.700, 15)
    int2 = np.exp(-0.5 * ((mz2 - 1570.6930) / 0.002) ** 2) * 8000.0

    mzs = np.concatenate([mz1, mz2])
    ints = np.concatenate([int1, int2])

    c_mz, c_int = _high_res_centroid(mzs, ints, resolving_power=250000.0, sn_threshold=2.0)
    assert len(c_mz) == 2
    assert abs(c_mz[0] - 1570.6774) < 0.0005
    assert abs(c_mz[1] - 1570.6930) < 0.0005


def test_high_res_profile_binning():
    """Test uniform high-resolution grid binning."""
    mzs = np.array([1000.001, 1000.002, 1000.005])
    ints = np.array([100.0, 200.0, 300.0])

    bin_mzs, bin_ints = _high_res_profile_bin(mzs, ints, bin_size_da=0.002)
    assert len(bin_mzs) >= 2
    assert np.sum(bin_ints) == np.sum(ints)
