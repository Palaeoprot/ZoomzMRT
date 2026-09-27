"""
High-Resolution Direct Deamidation (Q -> E / N -> D) & PQI Calculator for Waters MRT & FT-ICR.

On high-resolution mass spectrometers (R >= 200,000 FWHM, sub-ppm accuracy):
- The deamidated monoisotopic peak (M0_deam, +0.984016 Da) and the natural 13C1 isotopic peak
  (M1_undeam, +1.003355 Da) are fully baseline resolved by Delta_m = 0.019339 Da (19.3 mDa).
- No mathematical matrix deconvolution or least-squares unmixing is needed.
- Deamidation fractions and Parchment Glutamine Index (PQI) are measured by direct peak integration.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np
import pandas as pd

# Mass constants
MONO_DEAMIDATION_DELTA = 0.984016   # OH - NH2
C13_DELTA = 1.003355                # 13C - 12C
SEPARATION_DELTA = C13_DELTA - MONO_DEAMIDATION_DELTA  # 0.019339 Da

# Standard Collagen Deamidation Markers (COL1A1 & COL1A2)
# Reference: van Doorn et al. 2012, Wilson et al. 2012, Welker et al. 2016, Nair & Bethencourt et al. 2023
DEFAULT_COLLAGEN_DEAMIDATION_MARKERS = [
    {
        "name": "COL1A1_P1105",
        "gene": "COL1A1",
        "seq": "GVQGPPGPAGPR",
        "site": "Q502",
        "hyp_count": 1,
        "n_carbons": 48,
        "mz_undeam": 1105.5807,
    },
    {
        "name": "COL1A1_P1180",
        "gene": "COL1A1",
        "seq": "GQAGVMGFPGPK",
        "site": "Q391",
        "hyp_count": 1,
        "n_carbons": 52,
        "mz_undeam": 1180.5908,
    },
    {
        "name": "COL1A1_P1427",
        "gene": "COL1A1",
        "seq": "GSEGPQGVRGEPGPAGPR",
        "site": "Q178",
        "hyp_count": 1,
        "n_carbons": 68,
        "mz_undeam": 1427.6972,
    },
    {
        "name": "COL1A1_P1580",
        "gene": "COL1A1",
        "seq": "GATGAPGIAGAPGFPGAR",
        "site": "Q_or_G",
        "hyp_count": 1,
        "n_carbons": 72,
        "mz_undeam": 1580.7932,
    },
    {
        "name": "COL1A2_P1706",
        "gene": "COL1A2",
        "seq": "GIPGEFGLPGPAGAR",
        "site": "E/Q_marker",
        "hyp_count": 1,
        "n_carbons": 74,
        "mz_undeam": 1706.8837,
    },
    {
        "name": "COL1A1_P2043",
        "gene": "COL1A1",
        "seq": "GAPGADGPAGAPGTPGPQGIAGQR",
        "site": "Q774",
        "hyp_count": 2,
        "n_carbons": 94,
        "mz_undeam": 2043.9806,
    },
]


@dataclass
class PeptideDeamidationResult:
    """Deamidation metrics for a single peptide marker."""
    marker_name: str
    target_mz_undeam: float
    target_mz_deam: float
    obs_mz_undeam: float | None = None
    obs_mz_deam: float | None = None
    obs_mz_c13_undeam: float | None = None
    intensity_undeam: float = 0.0
    intensity_deam: float = 0.0
    intensity_c13_undeam: float = 0.0
    ppm_error_undeam: float | None = None
    ppm_error_deam: float | None = None
    deamidation_fraction: float | None = None   # % Deamidation = I_deam / (I_und + I_deam)
    pqi_fraction: float | None = None           # % Undeamidated = I_und / (I_und + I_deam)
    is_baseline_resolved: bool = True
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
    results: list[PeptideDeamidationResult] = field(default_factory=list)


def extract_peak_in_window(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    target_mz: float,
    tolerance_ppm: float = 10.0,
    tolerance_da: float = 0.008,
) -> tuple[float | None, float, float | None]:
    """Find peak closest to target_mz within high-resolution tolerance."""
    if len(mz_arr) == 0:
        return None, 0.0, None

    tol = max(target_mz * tolerance_ppm / 1e6, tolerance_da)
    mask = (mz_arr >= target_mz - tol) & (mz_arr <= target_mz + tol)
    
    if not np.any(mask):
        return None, 0.0, None

    sub_mz = mz_arr[mask]
    sub_int = int_arr[mask]
    best_idx = np.argmax(sub_int)
    
    obs_m = float(sub_mz[best_idx])
    obs_i = float(sub_int[best_idx])
    ppm_err = float((obs_m - target_mz) / target_mz * 1e6)

    return obs_m, obs_i, ppm_err


def compute_high_res_deamidation(
    mz_arr: Sequence[float] | np.ndarray,
    int_arr: Sequence[float] | np.ndarray,
    sample_id: str = "sample",
    dataset_id: str = "dataset",
    markers: list[dict[str, Any]] | None = None,
    tolerance_ppm: float = 8.0,
    min_intensity: float = 100.0,
) -> SampleDeamidationSummary:
    """Directly calculate high-resolution deamidation and PQI from resolved peaks.

    Args:
        mz_arr: Array of m/z values (calibrated).
        int_arr: Array of peak intensities.
        sample_id: Sample identifier.
        dataset_id: Dataset identifier.
        markers: List of marker dicts (defaults to DEFAULT_COLLAGEN_DEAMIDATION_MARKERS).
        tolerance_ppm: High-resolution search tolerance (e.g. 5-8 ppm for Waters MRT / FT-ICR).
        min_intensity: Minimum intensity required for a peak to be considered valid.

    Returns:
        SampleDeamidationSummary with per-peptide and aggregated PQI metrics.
    """
    m_arr = np.array(mz_arr, dtype=np.float64)
    i_arr = np.array(int_arr, dtype=np.float64)

    marker_list = markers or DEFAULT_COLLAGEN_DEAMIDATION_MARKERS
    summary = SampleDeamidationSummary(sample_id=sample_id, dataset_id=dataset_id)

    valid_pqis = []
    valid_deams = []

    for m in marker_list:
        m_name = m["name"]
        target_und = float(m["mz_undeam"])
        target_deam = target_und + MONO_DEAMIDATION_DELTA
        target_c13 = target_und + C13_DELTA
        n_c = m.get("n_carbons", 50)
        expected_c13_ratio = n_c * 0.0108  # ~1.08% per carbon

        # Extract undeamidated peak (M0_und)
        obs_und_m, obs_und_i, err_und = extract_peak_in_window(
            m_arr, i_arr, target_und, tolerance_ppm=tolerance_ppm
        )

        # Extract deamidated peak (M0_deam, +0.9840 Da)
        obs_deam_m, obs_deam_i, err_deam = extract_peak_in_window(
            m_arr, i_arr, target_deam, tolerance_ppm=tolerance_ppm
        )

        # Extract undeamidated 13C1 isotope (M1_und, +1.0034 Da)
        obs_c13_m, obs_c13_i, err_c13 = extract_peak_in_window(
            m_arr, i_arr, target_c13, tolerance_ppm=tolerance_ppm
        )

        res = PeptideDeamidationResult(
            marker_name=m_name,
            target_mz_undeam=target_und,
            target_mz_deam=target_deam,
            obs_mz_undeam=obs_und_m,
            obs_mz_deam=obs_deam_m,
            obs_mz_c13_undeam=obs_c13_m,
            intensity_undeam=obs_und_i,
            intensity_deam=obs_deam_i,
            intensity_c13_undeam=obs_c13_i,
            ppm_error_undeam=err_und,
            ppm_error_deam=err_deam,
            c13_ratio_expected=expected_c13_ratio,
        )

        # Evaluate detection and baseline separation
        tot_signal = obs_und_i + obs_deam_i
        if tot_signal >= min_intensity and (obs_und_i > 0 or obs_deam_i > 0):
            deam_frac = obs_deam_i / tot_signal
            pqi_frac = obs_und_i / tot_signal

            res.deamidation_fraction = float(deam_frac)
            res.pqi_fraction = float(pqi_frac)

            # C13 isotopic QC
            if obs_und_i > min_intensity and obs_c13_i > 0:
                obs_c13_ratio = obs_c13_i / obs_und_i
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
        }
        for r in s.results:
            m_prefix = r.marker_name
            base[f"{m_prefix}_pqi"] = r.pqi_fraction
            base[f"{m_prefix}_deam"] = r.deamidation_fraction
            base[f"{m_prefix}_i_und"] = r.intensity_undeam
            base[f"{m_prefix}_i_deam"] = r.intensity_deam
            base[f"{m_prefix}_flag"] = r.quality_flag
        rows.append(base)
    return pd.DataFrame(rows)
