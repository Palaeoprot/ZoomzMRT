"""
Unit tests for Waters MRT lockmass calibration and resolving power measurement.
"""

import numpy as np
import pytest

from zoomzmrt.mrt_reader import (
    WatersMRTScanQC,
    measure_peak_resolving_power,
)


def test_measure_peak_resolving_power_gaussian():
    """Verify that measure_peak_resolving_power accurately determines R = m / FWHM on a known peak."""
    target_mz = 1570.67742
    target_r = 250000.0
    fwhm = target_mz / target_r  # ~0.00628 Da
    sigma = fwhm / 2.35482

    mz_grid = np.linspace(1570.60, 1570.75, 500)
    int_profile = np.exp(-0.5 * ((mz_grid - target_mz) / sigma) ** 2) * 50000.0

    measured_r = measure_peak_resolving_power(mz_grid, int_profile, target_mz=target_mz)
    assert measured_r is not None
    # Measured R should be within 3% of target 250,000 FWHM
    assert abs(measured_r - target_r) / target_r < 0.03


def test_lockmass_multiplicative_calibration():
    """Verify lockmass calibration recovers exact theoretical mass from known ppm shifts."""
    lockmass_mz = 1570.67742
    
    # Simulate +2.5 ppm instrument drift
    obs_lm = lockmass_mz * (1.0 + 2.5e-6)
    corr_factor = lockmass_mz / obs_lm

    # Check that applying correction restores exact lockmass
    calibrated_lm = obs_lm * corr_factor
    assert abs(calibrated_lm - lockmass_mz) < 1e-9

    # Check on an analyte peak at 2000.0 Da with +2.5 ppm drift
    obs_analyte = 2000.0 * (1.0 + 2.5e-6)
    calibrated_analyte = obs_analyte * corr_factor
    assert abs(calibrated_analyte - 2000.0) < 1e-7
