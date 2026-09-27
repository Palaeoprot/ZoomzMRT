"""
Waters SELECT SERIES MRT (Multi-Reflecting Time-of-Flight) Reader Module.

High-resolution processing for Waters MRT MALDI-MS1 data:
  1. Scan-level laser shot QC and intensity filtering.
  2. Per-scan lockmass tracking, multiplicative alignment, and robust error metrics (Median/MAD).
  3. Empirical resolving power estimation (R = m / FWHM) on lockmass reference peak.
  4. Resolution-aware centroiding and high-resolution co-addition.
"""

from __future__ import annotations

import io
import math
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Sequence

import numpy as np
import pandas as pd
from pyteomics import mzml
from scipy.signal import find_peaks

# Standard lock mass references
KNOWN_LOCK_MASSES = {
    "glufib": 1570.67742,          # [Glu1]-Fibrinopeptide B [M+H]+ (EGVNDNEEGFFSAR)
    "leuenk": 556.27710,           # [Leu]-Enkephalin [M+H]+ (YGGFL)
    "bradykinin": 1060.56920,      # Bradykinin 1-9 [M+H]+ (RPPGFSPFR)
    "acth_18_39": 2465.19890,      # ACTH 18-39 [M+H]+
}


@dataclass
class WatersMRTScanQC:
    """Comprehensive quality control and calibration metrics for a Waters MRT acquisition."""
    file_id: str
    total_scans: int = 0
    included_scans: int = 0
    lockmass_detected_scans: int = 0
    lockmass_failed_scans: int = 0
    lockmass_theoretical: float = 1570.67742
    lockmass_median_ppm_error: float | None = None
    lockmass_mad_ppm: float | None = None
    lockmass_max_abs_ppm: float | None = None
    mean_ppm_shift: float = 0.0
    std_ppm_shift: float = 0.0
    measured_resolving_power: float | None = None
    nominal_resolving_power: float = 250000.0
    base_peak_mz: float = 0.0
    base_peak_intensity: float = 0.0
    total_raw_points: int = 0
    coadded_points: int = 0
    centroid_peaks: int = 0


def measure_peak_resolving_power(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    target_mz: float,
    window_da: float = 0.25,
) -> float | None:
    """Measure empirical resolving power (R = m / FWHM) around target_mz.

    Args:
        mz_arr: Sorted array of m/z values.
        int_arr: Intensity values matching mz_arr.
        target_mz: Theoretical or target peak m/z.
        window_da: Search window width in Da.

    Returns:
        Resolving power R (FWHM) or None if peak is not well-formed.
    """
    if len(mz_arr) < 5:
        return None

    mask = (mz_arr >= target_mz - window_da) & (mz_arr <= target_mz + window_da)
    if not np.any(mask):
        return None

    sub_mz = mz_arr[mask]
    sub_int = int_arr[mask]

    apex_idx = int(np.argmax(sub_int))
    apex_mz = sub_mz[apex_idx]
    apex_val = sub_int[apex_idx]

    if apex_val <= 0:
        return None

    half_max = apex_val / 2.0

    # Locate left half-max crossing
    left_mz = None
    for i in range(apex_idx, -1, -1):
        if sub_int[i] <= half_max:
            # Linear interpolation
            if i < apex_idx and (sub_int[i + 1] - sub_int[i]) > 0:
                frac = (half_max - sub_int[i]) / (sub_int[i + 1] - sub_int[i])
                left_mz = sub_mz[i] + frac * (sub_mz[i + 1] - sub_mz[i])
            else:
                left_mz = sub_mz[i]
            break

    # Locate right half-max crossing
    right_mz = None
    for i in range(apex_idx, len(sub_int)):
        if sub_int[i] <= half_max:
            # Linear interpolation
            if i > apex_idx and (sub_int[i - 1] - sub_int[i]) > 0:
                frac = (half_max - sub_int[i]) / (sub_int[i - 1] - sub_int[i])
                right_mz = sub_mz[i] - frac * (sub_mz[i] - sub_mz[i - 1])
            else:
                right_mz = sub_mz[i]
            break

    if left_mz is not None and right_mz is not None:
        fwhm = right_mz - left_mz
        if fwhm > 0:
            return float(apex_mz / fwhm)

    return None


def parse_waters_mrt_mzml(
    file_path: Path | str,
    dataset_id: str,
    lockmass_mz: float = 1570.67742,
    lockmass_window_da: float = 0.3,
    base_peak_threshold: float = 10000.0,
    sub_scan_min_intensity: float = 100.0,
    mz_min: float = 100.0,
    mz_max: float = 4000.0,
    bin_size_da: float = 0.002,
    centroid: bool = True,
    sn_threshold: float = 3.0,
    source_type: str = "external",
    instrument_name: str = "Waters SELECT SERIES MRT",
    nominal_resolving_power: float = 250000.0,
) -> tuple[dict[str, Any], WatersMRTScanQC]:
    """Parse a Waters MRT mzML file, performing scan-level lockmass calibration and co-addition.

    Returns:
        tuple of (record_dict, qc_metrics)
    """
    path_obj = Path(file_path)
    fname = path_obj.name
    sample_id = re.sub(r"\.(mzml|raw|txt)$", "", fname, flags=re.IGNORECASE)

    qc = WatersMRTScanQC(
        file_id=fname,
        lockmass_theoretical=lockmass_mz,
        nominal_resolving_power=nominal_resolving_power,
    )

    scans_mz = []
    scans_intensity = []
    scan_ppm_shifts = []

    with mzml.read(str(path_obj), huge_tree=True) as reader:
        for scan in reader:
            if scan.get("ms level") != 1:
                continue
            qc.total_scans += 1
            
            mz_arr = scan.get("m/z array")
            int_arr = scan.get("intensity array")
            if mz_arr is None or len(mz_arr) == 0:
                continue
            
            mz_arr = np.array(mz_arr, dtype=np.float64)
            int_arr = np.array(int_arr, dtype=np.float64)
            
            max_i = np.max(int_arr)
            if max_i < base_peak_threshold:
                continue
            
            # Lockmass evaluation in this scan
            lm_mask = (mz_arr >= lockmass_mz - lockmass_window_da) & (mz_arr <= lockmass_mz + lockmass_window_da)
            if np.any(lm_mask):
                lm_sub_i = int_arr[lm_mask]
                lm_sub_mz = mz_arr[lm_mask]
                best_lm_idx = np.argmax(lm_sub_i)
                obs_lm = lm_sub_mz[best_lm_idx]
                ppm_err = (obs_lm - lockmass_mz) / lockmass_mz * 1e6
                scan_ppm_shifts.append(ppm_err)
                qc.lockmass_detected_scans += 1
                
                # Multiplicative lockmass correction
                corr_factor = lockmass_mz / obs_lm
                mz_calibrated = mz_arr * corr_factor
            else:
                qc.lockmass_failed_scans += 1
                mz_calibrated = mz_arr
            
            # Scan filtering
            filter_mask = (mz_calibrated >= mz_min) & (mz_calibrated <= mz_max) & (int_arr >= sub_scan_min_intensity)
            if np.any(filter_mask):
                scans_mz.append(mz_calibrated[filter_mask])
                scans_intensity.append(int_arr[filter_mask])
                qc.included_scans += 1

    if not scans_mz:
        return {
            "file_id": fname,
            "dataset_id": dataset_id,
            "source_type": source_type,
            "raw_path": str(file_path),
            "sample_id": sample_id,
            "scan_number": qc.total_scans,
            "rt": None,
            "instrument": instrument_name,
            "mz": [],
            "intensity": [],
            "n_peaks": 0,
            "is_centroided": centroid,
            "extraction_strategy": "mrt_summed_lockmass",
            "instrument_status": "from_metadata",
            "instrument_serial": None,
            "schema_version": "0.1.0",
        }, qc

    if scan_ppm_shifts:
        shifts_arr = np.array(scan_ppm_shifts)
        qc.mean_ppm_shift = float(np.mean(shifts_arr))
        qc.std_ppm_shift = float(np.std(shifts_arr))
        qc.lockmass_median_ppm_error = float(np.median(shifts_arr))
        qc.lockmass_mad_ppm = float(np.median(np.abs(shifts_arr - np.median(shifts_arr))))
        qc.lockmass_max_abs_ppm = float(np.max(np.abs(shifts_arr)))
    else:
        qc.mean_ppm_shift = 0.0
        qc.std_ppm_shift = 0.0

    all_mzs = np.concatenate(scans_mz)
    all_intens = np.concatenate(scans_intensity)
    qc.total_raw_points = len(all_mzs)

    sort_idx = np.argsort(all_mzs)
    sorted_mz = all_mzs[sort_idx]
    sorted_int = all_intens[sort_idx]

    # Measure actual resolving power on lockmass peak or base peak
    measured_r = measure_peak_resolving_power(sorted_mz, sorted_int, target_mz=lockmass_mz)
    if measured_r is not None:
        qc.measured_resolving_power = measured_r
    else:
        # Fallback to base peak
        b_idx = int(np.argmax(sorted_int))
        b_mz = sorted_mz[b_idx]
        qc.measured_resolving_power = measure_peak_resolving_power(sorted_mz, sorted_int, target_mz=b_mz)

    eff_res = qc.measured_resolving_power or nominal_resolving_power

    if centroid:
        final_mz, final_int = _high_res_centroid(
            sorted_mz,
            sorted_int,
            resolving_power=eff_res,
            sn_threshold=sn_threshold,
        )
    else:
        final_mz, final_int = _high_res_profile_bin(sorted_mz, sorted_int, bin_size_da=bin_size_da)

    qc.coadded_points = len(final_mz)
    qc.centroid_peaks = len(final_mz)

    if len(final_int) > 0:
        max_idx = np.argmax(final_int)
        qc.base_peak_mz = float(final_mz[max_idx])
        qc.base_peak_intensity = float(final_int[max_idx])

    record = {
        "file_id": fname,
        "dataset_id": dataset_id,
        "source_type": source_type,
        "raw_path": str(file_path),
        "sample_id": sample_id,
        "scan_number": qc.included_scans,
        "rt": None,
        "instrument": instrument_name,
        "mz": [float(m) for m in final_mz],
        "intensity": [float(i) for i in final_int],
        "n_peaks": len(final_mz),
        "is_centroided": centroid,
        "extraction_strategy": "mrt_summed_lockmass",
        "instrument_status": "from_metadata",
        "instrument_serial": None,
        "schema_version": "0.1.0",
    }

    return record, qc


def _high_res_centroid(
    sorted_mz: np.ndarray,
    sorted_int: np.ndarray,
    resolving_power: float = 250000.0,
    sn_threshold: float = 3.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Perform resolution-aware intensity-weighted peak centroiding.

    Extraction windows scale dynamically with m/z and resolving power (FWHM = mz / R).
    """
    if len(sorted_mz) == 0:
        return np.array([]), np.array([])

    non_zero = sorted_int[sorted_int > 0]
    if len(non_zero) > 50:
        noise_level = float(np.percentile(non_zero, 25))
    elif len(non_zero) > 0:
        noise_level = float(np.min(non_zero))
    else:
        noise_level = 1.0

    if noise_level <= 0:
        noise_level = 1.0

    peaks, _ = find_peaks(sorted_int, height=noise_level * sn_threshold, distance=2)

    if len(peaks) == 0:
        # Fallback cluster centroiding
        median_mz = float(np.median(sorted_mz))
        adaptive_gap = (median_mz / resolving_power) * 1.5
        diffs = np.diff(sorted_mz)
        split_indices = np.where(diffs > adaptive_gap)[0] + 1
        clusters_mz = np.split(sorted_mz, split_indices)
        clusters_int = np.split(sorted_int, split_indices)

        centroid_mzs = []
        centroid_ints = []
        for cmz, cint in zip(clusters_mz, clusters_int):
            if len(cmz) == 0:
                continue
            max_i = np.max(cint)
            if max_i / noise_level >= sn_threshold:
                sum_i = np.sum(cint)
                c_mz = np.sum(cmz * cint) / sum_i
                centroid_mzs.append(c_mz)
                centroid_ints.append(float(max_i))
        return np.array(centroid_mzs, dtype=np.float64), np.array(centroid_ints, dtype=np.float64)

    centroid_mzs = []
    centroid_ints = []

    for p in peaks:
        peak_mz = sorted_mz[p]
        local_fwhm = peak_mz / resolving_power
        half_win = local_fwhm * 1.0

        # Find window boundaries by m/z range rather than fixed index count
        left = p
        while left > 0 and sorted_mz[p] - sorted_mz[left - 1] <= half_win:
            left -= 1

        right = p + 1
        while right < len(sorted_mz) and sorted_mz[right] - sorted_mz[p] <= half_win:
            right += 1

        w_mz = sorted_mz[left:right]
        w_int = sorted_int[left:right]
        sum_i = float(np.sum(w_int))
        if sum_i > 0:
            c_mz = float(np.sum(w_mz * w_int) / sum_i)
            centroid_mzs.append(c_mz)
            centroid_ints.append(float(np.max(w_int)))

    return np.array(centroid_mzs, dtype=np.float64), np.array(centroid_ints, dtype=np.float64)


def _high_res_profile_bin(
    sorted_mz: np.ndarray,
    sorted_int: np.ndarray,
    bin_size_da: float = 0.002,
) -> tuple[np.ndarray, np.ndarray]:
    """Co-add profile points onto high-resolution uniform grid."""
    if len(sorted_mz) == 0:
        return np.array([]), np.array([])
    
    min_mz = sorted_mz[0]
    max_mz = sorted_mz[-1]
    n_bins = int(math.ceil((max_mz - min_mz) / bin_size_da)) + 1
    
    bin_indices = np.floor((sorted_mz - min_mz) / bin_size_da).astype(np.int64)
    binned_intensities = np.zeros(n_bins, dtype=np.float64)
    np.add.at(binned_intensities, bin_indices, sorted_int)
    
    bin_mzs = min_mz + np.arange(n_bins) * bin_size_da
    mask = binned_intensities > 0
    return bin_mzs[mask], binned_intensities[mask]
