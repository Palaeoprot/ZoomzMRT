"""
Unit tests for MALDI-FTICR fine isotopic preservation and multi-scan co-addition.
"""

import numpy as np
import pytest

from zoomzmrt.fticr_reader import _fticr_centroid, _fticr_coadd


def test_fticr_fine_isotopic_preservation():
    """Verify FTICR centroiding preserves fine isotopic splits separated by < 10 mDa."""
    mz_arr = np.array([1200.500, 1200.501, 1200.502, 1200.508, 1200.509, 1200.510])
    int_arr = np.array([5000.0, 10000.0, 5000.0, 4000.0, 8000.0, 4000.0])

    c_mz, c_int = _fticr_centroid(mz_arr, int_arr, resolving_power=500000.0, sn_threshold=2.0)
    assert len(c_mz) == 2
    assert abs(c_mz[0] - 1200.501) < 0.0005
    assert abs(c_mz[1] - 1200.509) < 0.0005


def test_fticr_multi_scan_coadd():
    """Verify FTICR multi-scan co-addition preserves intensities and centroid positions."""
    scan1_mz = np.array([1500.100, 1500.200])
    scan1_int = np.array([1000.0, 2000.0])

    scan2_mz = np.array([1500.101, 1500.201])
    scan2_int = np.array([1000.0, 2000.0])

    all_mz = np.concatenate([scan1_mz, scan2_mz])
    all_int = np.concatenate([scan1_int, scan2_int])
    sort_idx = np.argsort(all_mz)

    # Aggregation = mean (scale_factor = 0.5)
    mean_mz, mean_int = _fticr_coadd(all_mz[sort_idx], all_int[sort_idx], cluster_gap_da=0.005, scale_factor=0.5)
    assert len(mean_mz) == 2
    assert abs(mean_mz[0] - 1500.1005) < 0.0001
    assert abs(mean_int[0] - 1000.0) < 0.1
    assert abs(mean_int[1] - 2000.0) < 0.1
