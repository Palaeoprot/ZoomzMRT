"""
High-Resolution Direct Deamidation (Q -> E / N -> D) & PQI Calculator for Waters MRT & FT-ICR.

On high-resolution mass spectrometers (R >= 200,000 FWHM, sub-ppm accuracy):
- The deamidated monoisotopic peak (M0_deam, +0.984016 Da) and the natural 13C1 isotopic peak
  (M1_undeam, +1.003355 Da) are separated by Delta_m = 0.019339 Da (19.339 mDa).
- Deamidation fractions and Parchment Glutamine Index (PQI) are calculated by true numerical
  trapezoidal peak-area integration (integrate_peak) with local baseline subtraction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Literal, Sequence

import numpy as np
import pandas as pd

from zoomzmrt.isotopes import (
    DEFAULT_COLLAGEN_MARKERS,
    DeamidationMarker,
    ElementalComposition,
    peptide_composition,
)

# Mass constants
MONO_DEAMIDATION_DELTA = 0.984016   # OH - NH2
C13_DELTA = 1.003355                # 13C - 12C
SEPARATION_DELTA = C13_DELTA - MONO_DEAMIDATION_DELTA  # 0.019339 Da


# Compatibility helper for numpy trapz vs trapezoid (NumPy 2.0+)
def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


@dataclass
class PeptideDeamidationResult:
    """Deamidation metrics and peak area integration for a single peptide marker."""
    marker_name: str
    target_mz_undeam: float
    target_mz_deam: float
    target_mz_c13: float
    obs_mz_undeam: float | None = None
    obs_mz_deam: float | None = None
    obs_mz_c13_undeam: float | None = None
    area_undeam: float = 0.0
    area_deam: float = 0.0
    area_c13_undeam: float = 0.0
    intensity_undeam: float = 0.0          # Peak apex height
    intensity_deam: float = 0.0            # Peak apex height
    intensity_c13_undeam: float = 0.0      # Peak apex height
    ppm_error_undeam: float | None = None
    ppm_error_deam: float | None = None
    deamidation_fraction: float | None = None   # % Deamidation = Area_deam / (Area_und + Area_deam)
    pqi_fraction: float | None = None           # % Undeamidated = Area_und / (Area_und + Area_deam)
    resolution_status: str = "not_assessed"     # 'resolved', 'partially_resolved', 'unresolved', 'not_detected'
    separation_over_fwhm: float | None = None
    c13_ratio_observed: float | None = None
    c13_ratio_expected: float | None = None
    quality_flag: str = "OK"


@dataclass
class SampleDeamidationSummary:
    """Aggregated deamidation and PQI estimates for a sample."""
    sample_id: str
    dataset_id: str
    n_markers_detected: int = 0
    pqi_mean: float | None = None
    pqi_median: float | None = None
    deamidation_mean: float | None = None
    deamidation_median: float | None = None
    measured_resolving_power: float | None = None
    results: list[PeptideDeamidationResult] = field(default_factory=list)


def peak_half_width(
    mz: float,
    resolving_power: float = 250000.0,
    width_factor: float = 1.5,
) -> float:
    """Derive peak integration half-width from resolving power (FWHM = mz / R).

    Args:
        mz: Target m/z value.
        resolving_power: Nominal or measured resolving power (FWHM).
        width_factor: Multiplier for FWHM to encompass full peak bounds (default: 1.5x).

    Returns:
        Half-width in Da.
    """
    if resolving_power <= 0:
        resolving_power = 250000.0
    fwhm = mz / resolving_power
    return float(0.5 * width_factor * fwhm)


def find_peak_in_window(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    target_mz: float,
    tolerance_da: float,
) -> tuple[float | None, float, float | None]:
    """Locate the apex of the strongest peak near target_mz within tolerance_da."""
    if len(mz_arr) == 0:
        return None, 0.0, None

    mask = (mz_arr >= target_mz - tolerance_da) & (mz_arr <= target_mz + tolerance_da)
    if not np.any(mask):
        return None, 0.0, None

    sub_mz = mz_arr[mask]
    sub_int = int_arr[mask]
    best_idx = int(np.argmax(sub_int))
    
    obs_m = float(sub_mz[best_idx])
    obs_i = float(sub_int[best_idx])
    ppm_err = float((obs_m - target_mz) / target_mz * 1e6)

    return obs_m, obs_i, ppm_err


def integrate_peak(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    center_mz: float,
    half_width_da: float,
    baseline_subtraction: bool = True,
) -> tuple[float, float | None, float]:
    """Compute background-corrected trapezoidal peak area and centroid m/z.

    Args:
        mz_arr: Sorted array of m/z values.
        int_arr: Intensity array matching mz_arr.
        center_mz: Center m/z around which to integrate.
        half_width_da: Half-width window (Da).
        baseline_subtraction: Whether to subtract local minimum baseline.

    Returns:
        tuple of (peak_area, centroid_mz, apex_intensity)
    """
    if len(mz_arr) == 0:
        return 0.0, None, 0.0

    mask = (mz_arr >= center_mz - half_width_da) & (mz_arr <= center_mz + half_width_da)
    if not np.any(mask):
        return 0.0, None, 0.0

    mz_win = mz_arr[mask]
    int_win = int_arr[mask]
    apex_i = float(np.max(int_win))

    if len(mz_win) < 2:
        return apex_i, float(mz_win[0]), apex_i

    if baseline_subtraction:
        baseline = float(np.min(int_win))
        corrected = np.clip(int_win - baseline, 0.0, None)
    else:
        corrected = int_win

    area = _trapz(corrected, mz_win)
    
    # Fallback to apex height if trapezoidal area evaluates to 0 (e.g. baseline equals apex or centroided sparse peak)
    if area <= 0.0:
        area = apex_i

    sum_corr = float(np.sum(corrected))
    if sum_corr > 0.0:
        centroid_mz = float(np.sum(mz_win * corrected) / sum_corr)
    else:
        centroid_mz = float(mz_win[np.argmax(int_win)])

    return float(area), centroid_mz, apex_i


def compute_high_res_deamidation(
    mz_arr: Sequence[float] | np.ndarray,
    int_arr: Sequence[float] | np.ndarray,
    sample_id: str = "sample",
    dataset_id: str = "dataset",
    markers: Sequence[DeamidationMarker] | None = None,
    resolving_power: float = 250000.0,
    tolerance_ppm: float = 8.0,
    min_intensity: float = 100.0,
) -> SampleDeamidationSummary:
    """Directly calculate high-resolution deamidation and PQI using trapezoidal peak areas.

    Args:
        mz_arr: Array of calibrated m/z values.
        int_arr: Array of peak intensities.
        sample_id: Sample identifier.
        dataset_id: Dataset identifier.
        markers: List/tuple of DeamidationMarker instances (defaults to DEFAULT_COLLAGEN_MARKERS).
        resolving_power: Nominal or measured resolving power of the acquisition (FWHM).
        tolerance_ppm: High-resolution search tolerance.
        min_intensity: Minimum apex intensity required for a peak to be evaluated.

    Returns:
        SampleDeamidationSummary containing per-peptide results and aggregated PQI metrics.
    """
    m_arr = np.array(mz_arr, dtype=np.float64)
    i_arr = np.array(int_arr, dtype=np.float64)

    marker_list = markers or DEFAULT_COLLAGEN_MARKERS
    summary = SampleDeamidationSummary(
        sample_id=sample_id,
        dataset_id=dataset_id,
        measured_resolving_power=resolving_power,
    )

    valid_pqis = []
    valid_deams = []

    for m in marker_list:
        target_und = m.theoretical_mz
        target_deam = m.target_deam_mz
        target_c13 = m.target_13c_mz
        expected_c13_ratio = m.expected_13c_ratio

        # Calculate resolution-aware search tolerance & integration half-width
        tol_da = max(target_und * tolerance_ppm / 1e6, peak_half_width(target_und, resolving_power, width_factor=1.0))
        half_width_da = peak_half_width(target_und, resolving_power, width_factor=1.5)

        # 1. Locate apexes
        obs_und_m, apex_und_i, err_und = find_peak_in_window(m_arr, i_arr, target_und, tol_da)
        obs_deam_m, apex_deam_i, err_deam = find_peak_in_window(m_arr, i_arr, target_deam, tol_da)
        obs_c13_m, apex_c13_i, err_c13 = find_peak_in_window(m_arr, i_arr, target_c13, tol_da)

        # 2. Integrate trapezoidal peak areas
        area_und, cent_und_m, _ = integrate_peak(
            m_arr, i_arr, obs_und_m if obs_und_m is not None else target_und, half_width_da
        )
        area_deam, cent_deam_m, _ = integrate_peak(
            m_arr, i_arr, obs_deam_m if obs_deam_m is not None else target_deam, half_width_da
        )
        area_c13, cent_c13_m, _ = integrate_peak(
            m_arr, i_arr, obs_c13_m if obs_c13_m is not None else target_c13, half_width_da
        )

        res = PeptideDeamidationResult(
            marker_name=m.name,
            target_mz_undeam=target_und,
            target_mz_deam=target_deam,
            target_mz_c13=target_c13,
            obs_mz_undeam=cent_und_m or obs_und_m,
            obs_mz_deam=cent_deam_m or obs_deam_m,
            obs_mz_c13_undeam=cent_c13_m or obs_c13_m,
            area_undeam=area_und,
            area_deam=area_deam,
            area_c13_undeam=area_c13,
            intensity_undeam=apex_und_i,
            intensity_deam=apex_deam_i,
            intensity_c13_undeam=apex_c13_i,
            ppm_error_undeam=err_und,
            ppm_error_deam=err_deam,
            c13_ratio_expected=expected_c13_ratio,
        )

        # 3. Assess physical baseline resolution between deamidated M0 and undeamidated 13C1
        fwhm = target_und / resolving_power
        sep_ratio = SEPARATION_DELTA / fwhm if fwhm > 0 else 0.0
        res.separation_over_fwhm = float(sep_ratio)

        if obs_deam_m is not None and obs_c13_m is not None:
            if sep_ratio >= 1.5:
                res.resolution_status = "resolved"
            elif sep_ratio >= 0.8:
                res.resolution_status = "partially_resolved"
            else:
                res.resolution_status = "unresolved"
        elif obs_und_m is not None or obs_deam_m is not None:
            res.resolution_status = "resolved" if sep_ratio >= 1.5 else "partially_resolved"
        else:
            res.resolution_status = "not_detected"

        # 4. Calculate deamidation fraction from integrated peak areas
        tot_area = area_und + area_deam
        tot_apex = apex_und_i + apex_deam_i

        if (tot_area > 0 or tot_apex > 0) and tot_apex >= min_intensity:
            if tot_area > 0:
                deam_frac = area_deam / tot_area
                pqi_frac = area_und / tot_area
            else:
                deam_frac = apex_deam_i / tot_apex
                pqi_frac = apex_und_i / tot_apex

            res.deamidation_fraction = float(deam_frac)
            res.pqi_fraction = float(pqi_frac)

            # 5. Isotopic consistency QC
            if area_und > 0 and area_c13 > 0:
                obs_c13_ratio = area_c13 / area_und
                res.c13_ratio_observed = float(obs_c13_ratio)
                if abs(obs_c13_ratio - expected_c13_ratio) / expected_c13_ratio > 0.40:
                    res.quality_flag = "C13_RATIO_ANOMALY"
            elif apex_und_i > min_intensity and apex_c13_i > 0:
                obs_c13_ratio = apex_c13_i / apex_und_i
                res.c13_ratio_observed = float(obs_c13_ratio)
                if abs(obs_c13_ratio - expected_c13_ratio) / expected_c13_ratio > 0.40:
                    res.quality_flag = "C13_RATIO_ANOMALY"

            valid_pqis.append(pqi_frac)
            valid_deams.append(deam_frac)
            summary.n_markers_detected += 1
        else:
            res.quality_flag = "BELOW_INTENSITY_THRESHOLD"

        summary.results.append(res)

    if valid_pqis:
        summary.pqi_mean = float(np.mean(valid_pqis))
        summary.pqi_median = float(np.median(valid_pqis))
        summary.deamidation_mean = float(np.mean(valid_deams))
        summary.deamidation_median = float(np.median(valid_deams))

    return summary


def deamidation_summary_to_dataframe(summaries: list[SampleDeamidationSummary]) -> pd.DataFrame:
    """Flatten a list of sample deamidation summaries into a tabular DataFrame."""
    rows = []
    for s in summaries:
        base = {
            "dataset_id": s.dataset_id,
            "sample_id": s.sample_id,
            "n_markers_detected": s.n_markers_detected,
            "pqi_mean": s.pqi_mean,
            "pqi_median": s.pqi_median,
            "deamidation_mean": s.deamidation_mean,
            "deamidation_median": s.deamidation_median,
            "measured_resolving_power": s.measured_resolving_power,
        }
        for r in s.results:
            m_prefix = r.marker_name
            base[f"{m_prefix}_pqi"] = r.pqi_fraction
            base[f"{m_prefix}_deam"] = r.deamidation_fraction
            base[f"{m_prefix}_area_und"] = r.area_undeam
            base[f"{m_prefix}_area_deam"] = r.area_deam
            base[f"{m_prefix}_i_und"] = r.intensity_undeam
            base[f"{m_prefix}_i_deam"] = r.intensity_deam
            base[f"{m_prefix}_res_status"] = r.resolution_status
            base[f"{m_prefix}_flag"] = r.quality_flag
        rows.append(base)
    return pd.DataFrame(rows)
