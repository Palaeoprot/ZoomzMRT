"""
Unit tests for trapezoidal peak area integration and high-resolution deamidation / PQI calculation.
"""

import numpy as np
import pytest

from zoomzmrt.deamidation import (
    compute_high_res_deamidation,
    integrate_peak,
    peak_half_width,
)


def test_integrate_peak_gaussian_area():
    """Verify that integrate_peak integrates true peak area rather than merely returning apex height."""
    mz_grid = np.linspace(1105.50, 1105.65, 151)

    # Peak 1: Narrow Gaussian at 1105.54 (sigma = 0.002 Da, Height = 10000)
    # Area ~ Height * sigma * sqrt(2*pi) = 10000 * 0.002 * 2.5066 = 50.13
    p1 = np.exp(-0.5 * ((mz_grid - 1105.54) / 0.002) ** 2) * 10000.0

    # Peak 2: Broad Gaussian at 1105.60 (sigma = 0.004 Da, Height = 5000)
    # Area ~ Height * sigma * sqrt(2*pi) = 5000 * 0.004 * 2.5066 = 50.13 (Equal area, half apex height)
    p2 = np.exp(-0.5 * ((mz_grid - 1105.60) / 0.004) ** 2) * 5000.0

    total_int = p1 + p2

    area1, cent1, apex1 = integrate_peak(mz_grid, total_int, center_mz=1105.54, half_width_da=0.01)
    area2, cent2, apex2 = integrate_peak(mz_grid, total_int, center_mz=1105.60, half_width_da=0.015)

    assert abs(cent1 - 1105.54) < 0.0005
    assert abs(cent2 - 1105.60) < 0.0005
    assert apex1 > 9900.0
    assert apex2 > 4900.0

    # Areas should be approximately equal despite 2:1 apex ratio
    area_ratio = area1 / area2
    assert abs(area_ratio - 1.0) < 0.08


def test_direct_pqi_estimation_synthetic():
    # Simulate P1105 undeamidated (1105.5807) and deamidated (1106.5647) and 13C (1106.5841)
    # Area ratio: Undeamidated = 10000, Deamidated = 5000 -> Deamidation = 33.3%, PQI = 66.7%
    mzs = [1105.5807, 1106.5647, 1106.5841]
    ints = [10000.0, 5000.0, 5174.0]  # Expected 13C ~ 51.7% of undeamidated (48 carbons * 1.078%)

    summary = compute_high_res_deamidation(mzs, ints, sample_id="test_sample", dataset_id="test_ds")
    assert summary.n_markers_detected >= 1

    p1105 = [r for r in summary.results if r.marker_name == "COL1A1_P1105"][0]
    assert p1105.pqi_fraction is not None
    assert abs(p1105.pqi_fraction - (10000.0 / 15000.0)) < 0.01
    assert abs(p1105.deamidation_fraction - (5000.0 / 15000.0)) < 0.01
    assert p1105.quality_flag == "OK"


@pytest.mark.parametrize("resolving_power, expected_status", [
    (100000.0, "resolved"),
    (250000.0, "resolved"),
    (500000.0, "resolved"),
])
def test_19mda_separation_benchmark_across_resolving_power(resolving_power, expected_status):
    """Benchmark recovering deamidation and 13C across resolving powers from 100k to 500k."""
    # Target marker: COL1A1_P1105 (1105.5807 Da)
    # Deamidated M0: 1106.564716 Da
    # Undeamidated 13C1: 1106.584055 Da (separated by 19.339 mDa)
    mz_und = 1105.5807
    mz_deam = 1106.564716
    mz_c13 = 1106.584055

    fwhm = 1105.58 / resolving_power
    sigma = fwhm / 2.35482

    mz_grid = np.linspace(1105.50, 1106.65, 3000)
    int_profile = np.zeros_like(mz_grid)

    # True areas: Und = 10000, Deam = 5000, C13 = 5200 (52% of Und)
    int_profile += (10000.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_und) / sigma) ** 2)
    int_profile += (5000.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_deam) / sigma) ** 2)
    int_profile += (5200.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_c13) / sigma) ** 2)

    summary = compute_high_res_deamidation(
        mz_arr=mz_grid,
        int_arr=int_profile,
        sample_id="benchmark_sample",
        dataset_id="benchmark_ds",
        resolving_power=resolving_power,
    )

    p1105 = [r for r in summary.results if r.marker_name == "COL1A1_P1105"][0]
    assert p1105.pqi_fraction is not None
    # Check recovered area ratio (true is 10000 / 15000 = 0.6667)
    assert abs(p1105.pqi_fraction - 0.6667) < 0.05
    assert p1105.resolution_status == expected_status


def test_c13_anomaly_flagging():
    """Verify that anomalous 13C abundance triggers C13_RATIO_ANOMALY flag."""
    # Simulate P1105 undeamidated (10000) and deamidated (5000) but severely distorted 13C (100 instead of ~5200)
    mzs = [1105.5807, 1106.5647, 1106.5841]
    ints = [10000.0, 5000.0, 100.0]

    summary = compute_high_res_deamidation(mzs, ints, sample_id="anomaly_sample", dataset_id="test_ds")
    p1105 = [r for r in summary.results if r.marker_name == "COL1A1_P1105"][0]
    assert p1105.quality_flag == "C13_RATIO_ANOMALY"


@pytest.mark.parametrize("resolving_power", [100000.0, 150000.0, 250000.0, 500000.0])
def test_high_mass_marker_2044_deamidation_no_crosstalk(resolving_power):
    """Rigorous regression test for high-mass marker (~2044 Da).

    Verifies that the 8 ppm search window does NOT overlap or cross-assign adjacent
    deamidated M0 (2045.0161) and natural 13C1 (2045.0355) peaks separated by 19.339 mDa.
    """
    from zoomzmrt.isotopes import DeamidationMarker, ElementalComposition

    marker_2044 = DeamidationMarker(
        name="COL1A1_P2044",
        gene="COL1A1",
        sequence="GVQGPPGPAGPRGEQGPPGPPGPQGPR",
        site="Gln",
        theoretical_mz=2044.0321,
    )

    mz_und = 2044.0321
    mz_deam = mz_und + 0.984016   # 2045.016116
    mz_c13 = mz_und + 1.003355    # 2045.035455
    separation = mz_c13 - mz_deam # 0.019339 Da

    fwhm = 2044.0321 / resolving_power
    sigma = fwhm / 2.35482

    mz_grid = np.linspace(2044.00, 2045.10, 5000)
    int_profile = np.zeros_like(mz_grid)

    # True areas: Und = 20000, Deam = 10000 (33.3% deamidation, PQI = 66.7%), 13C = 19000 (95% of Und)
    true_und_area = 20000.0
    true_deam_area = 10000.0
    true_c13_area = 19000.0

    int_profile += (true_und_area / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_und) / sigma) ** 2)
    int_profile += (true_deam_area / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_deam) / sigma) ** 2)
    int_profile += (true_c13_area / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((mz_grid - mz_c13) / sigma) ** 2)

    summary = compute_high_res_deamidation(
        mz_arr=mz_grid,
        int_arr=int_profile,
        sample_id="high_mass_test",
        dataset_id="benchmark_2044",
        markers=[marker_2044],
        resolving_power=resolving_power,
    )

    res = summary.results[0]
    assert res.obs_mz_deam is not None
    assert res.obs_mz_c13_undeam is not None
    
    # 1. Critical assertion: At R >= 150k (MRT and FTICR), Deam and 13C peak positions are distinct
    if resolving_power >= 150000.0:
        assert res.obs_mz_deam != res.obs_mz_c13_undeam
        assert abs(res.obs_mz_deam - mz_deam) < 0.003
        assert abs(res.obs_mz_c13_undeam - mz_c13) < 0.003
        assert res.resolution_status in ("resolved", "partially_resolved")
    else:
        # At R=100k, separation (19.3 mDa) < FWHM (20.4 mDa), confirming partially_resolved / unresolved status
        assert res.resolution_status in ("partially_resolved", "unresolved")

    # 2. Both areas must be positive and non-zero
    assert res.area_undeam > 0
    assert res.area_deam > 0
    assert res.area_c13_undeam > 0
    assert res.measurement_method == "trapezoidal_area"

    # 3. Quantitation accuracy: at R >= 250k (MRT/FTICR), recovered PQI is accurate within 5%
    if resolving_power >= 250000.0:
        true_pqi = true_und_area / (true_und_area + true_deam_area) # 0.6667
        assert abs(res.pqi_fraction - true_pqi) < 0.05


def test_sparse_centroid_apex_fallback():
    """Verify that sparse centroid peaks (where trapezoidal area is 0) trigger apex_fallback."""
    mzs = [1105.5807, 1106.5647]
    ints = [10000.0, 5000.0]

    summary = compute_high_res_deamidation(mzs, ints, sample_id="sparse_sample", dataset_id="test_ds")
    p1105 = summary.results[0]
    assert p1105.measurement_method == "apex_fallback"
    assert p1105.quality_flag == "OK"
    assert abs(p1105.pqi_fraction - 0.6667) < 0.01

