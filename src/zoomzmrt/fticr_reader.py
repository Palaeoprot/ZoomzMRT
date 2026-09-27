"""
MALDI-FTICR (Fourier Transform Ion Cyclotron Resonance) Reader Module.

Ultra-high-resolution processing (500,000 - 1,000,000+ FWHM, sub-ppm accuracy) for MALDI-FTICR
data (Bruker solariX 7T/9.4T/12T/15T, Thermo FTMS):
  1. Preserves fine isotopic structure (e.g. 13C vs 15N vs 34S isobars, deamidation splits).
  2. Multi-scan co-addition with explicit aggregation semantics ('mean', 'sum', 'none').
  3. Measured resolving power estimation (R = m / FWHM) on dominant base peak.
  4. Resolving-power scaled peak centroiding without multiplet smearing.
"""

from __future__ import annotations

import io
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Literal, Sequence

import numpy as np
import pandas as pd
from pyteomics import mzml

from zoomzmrt.mrt_reader import measure_peak_resolving_power


@dataclass
class FTICRScanQC:
    """Quality control metrics for a MALDI-FTICR acquisition."""
    file_id: str
    total_scans: int = 0
    included_scans: int = 0
    n_peaks: int = 0
    base_peak_mz: float = 0.0
    base_peak_intensity: float = 0.0
    measured_resolving_power: float | None = None
    nominal_resolving_power: float = 500000.0
    magnetic_field_tesla: float | None = None
    n_mz_gaps_3_25_mda: int = 0
    instrument_name: str = "Bruker solariX FT-ICR"


def parse_fticr_mzml(
    file_path: Path | str,
    dataset_id: str,
    mz_min: float = 200.0,
    mz_max: float = 5000.0,
    centroid: bool = True,
    sn_threshold: float = 2.5,
    aggregation: Literal["mean", "sum", "first", "none"] = "mean",
    source_type: str = "external",
    instrument_name: str = "Bruker solariX FT-ICR",
    nominal_resolving_power: float = 500000.0,
    magnetic_field_tesla: float | None = None,
) -> tuple[dict[str, Any], FTICRScanQC]:
    """Parse a MALDI-FTICR mzML file, preserving fine isotopic structure and resolving power."""
    path_obj = Path(file_path)
    fname = path_obj.name
    sample_id = re.sub(r"\.(mzml|mzxml|txt|raw)$", "", fname, flags=re.IGNORECASE)

    qc = FTICRScanQC(
        file_id=fname,
        instrument_name=instrument_name,
        nominal_resolving_power=nominal_resolving_power,
        magnetic_field_tesla=magnetic_field_tesla,
    )

    all_mz_scans = []
    all_int_scans = []
    inst_header = None
    inst_status = "from_metadata"

    with mzml.read(str(path_obj), huge_tree=True) as reader:
        for scan in reader:
            if scan.get("ms level") != 1:
                continue
            qc.total_scans += 1
            
            if not inst_header:
                cfg = scan.get("instrumentConfigurationRef")
                if cfg:
                    inst_header = str(cfg)
                    inst_status = "from_header"

            mz_arr = scan.get("m/z array")
            int_arr = scan.get("intensity array")
            if mz_arr is None or len(mz_arr) == 0:
                continue

            mz_arr = np.array(mz_arr, dtype=np.float64)
            int_arr = np.array(int_arr, dtype=np.float64)

            mask = (mz_arr >= mz_min) & (mz_arr <= mz_max)
            if np.any(mask):
                all_mz_scans.append(mz_arr[mask])
                all_int_scans.append(int_arr[mask])
                qc.included_scans += 1

    if not all_mz_scans:
        return {
            "file_id": fname,
            "dataset_id": dataset_id,
            "source_type": source_type,
            "raw_path": str(file_path),
            "sample_id": sample_id,
            "scan_number": qc.total_scans,
            "rt": None,
            "instrument": inst_header or instrument_name,
            "mz": [],
            "intensity": [],
            "n_peaks": 0,
            "is_centroided": centroid,
            "extraction_strategy": "fticr_centroid" if centroid else "fticr_profile",
            "instrument_status": inst_status,
            "instrument_serial": None,
            "schema_version": "0.1.0",
        }, qc

    if len(all_mz_scans) == 1 or aggregation in ("first", "none"):
        raw_mz = all_mz_scans[0]
        raw_int = all_int_scans[0]
    else:
        # Multi-scan co-addition with aggregation mode ('mean' or 'sum')
        concat_mz = np.concatenate(all_mz_scans)
        concat_int = np.concatenate(all_int_scans)
        s_idx = np.argsort(concat_mz)
        sorted_mz = concat_mz[s_idx]
        sorted_int = concat_int[s_idx]

        if aggregation == "mean":
            # Rescale summed intensity by number of scans
            raw_mz, raw_int = _fticr_coadd(
                sorted_mz,
                sorted_int,
                resolving_power=nominal_resolving_power,
                scale_factor=1.0 / len(all_mz_scans),
            )
        else:  # 'sum'
            raw_mz, raw_int = _fticr_coadd(
                sorted_mz,
                sorted_int,
                resolving_power=nominal_resolving_power,
                scale_factor=1.0,
            )

    # Measure actual empirical resolving power on base peak
    if len(raw_int) > 0:
        b_idx = int(np.argmax(raw_int))
        b_mz = raw_mz[b_idx]
        qc.measured_resolving_power = measure_peak_resolving_power(raw_mz, raw_int, target_mz=b_mz)

    eff_res = qc.measured_resolving_power or nominal_resolving_power

    if centroid:
        final_mz, final_int = _fticr_centroid(raw_mz, raw_int, resolving_power=eff_res, sn_threshold=sn_threshold)
        extraction_strat = "fticr_centroid"
    else:
        final_mz, final_int = raw_mz, raw_int
        extraction_strat = "fticr_profile"

    qc.n_peaks = len(final_mz)
    if len(final_int) > 0:
        max_idx = int(np.argmax(final_int))
        qc.base_peak_mz = float(final_mz[max_idx])
        qc.base_peak_intensity = float(final_int[max_idx])

    if len(final_mz) > 1:
        diffs = np.diff(final_mz)
        qc.n_mz_gaps_3_25_mda = int(np.sum((diffs > 0.003) & (diffs < 0.025)))

    record = {
        "file_id": fname,
        "dataset_id": dataset_id,
        "source_type": source_type,
        "raw_path": str(file_path),
        "sample_id": sample_id,
        "scan_number": qc.included_scans,
        "rt": None,
        "instrument": inst_header or instrument_name,
        "mz": [float(m) for m in final_mz],
        "intensity": [float(i) for i in final_int],
        "n_peaks": len(final_mz),
        "is_centroided": centroid,
        "extraction_strategy": extraction_strat,
        "instrument_status": inst_status,
        "instrument_serial": None,
        "schema_version": "0.1.0",
    }

    return record, qc


def _fticr_coadd(
    sorted_mz: np.ndarray,
    sorted_int: np.ndarray,
    cluster_gap_da: float | None = None,
    resolving_power: float = 500000.0,
    scale_factor: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Co-add closely spaced multi-scan FT-ICR points preserving fine resolution."""
    if len(sorted_mz) == 0:
        return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

    if cluster_gap_da is None:
        median_mz = float(np.median(sorted_mz))
        cluster_gap_da = max(0.0002, (median_mz / resolving_power) * 0.4)

    diffs = np.diff(sorted_mz)
    split_indices = np.where(diffs > cluster_gap_da)[0] + 1
    clusters_mz = np.split(sorted_mz, split_indices)
    clusters_int = np.split(sorted_int, split_indices)

    coadd_mz = []
    coadd_int = []

    for cmz, cint in zip(clusters_mz, clusters_int):
        if len(cmz) == 0:
            continue
        sum_i = float(np.sum(cint))
        if sum_i > 0:
            c_mz = float(np.sum(cmz * cint) / sum_i)
            coadd_mz.append(c_mz)
            coadd_int.append(sum_i * scale_factor)

    return np.array(coadd_mz, dtype=np.float64), np.array(coadd_int, dtype=np.float64)


def _fticr_centroid(
    mz_arr: np.ndarray,
    int_arr: np.ndarray,
    resolving_power: float = 500000.0,
    sn_threshold: float = 2.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Centroid ultra-high resolution FT-ICR peaks using adaptive narrow windows."""
    if len(mz_arr) == 0:
        return np.array([]), np.array([])

    median_mz = float(np.median(mz_arr))
    gap_da = max(0.0005, (median_mz / resolving_power) * 1.2)

    diffs = np.diff(mz_arr)
    split_indices = np.where(diffs > gap_da)[0] + 1
    clusters_mz = np.split(mz_arr, split_indices)
    clusters_int = np.split(int_arr, split_indices)

    all_maxs = [np.max(ci) for ci in clusters_int if len(ci) > 0]
    if len(all_maxs) > 20:
        noise_level = float(np.percentile(all_maxs, 25))
    elif len(all_maxs) > 5:
        noise_level = float(np.min(all_maxs))
    else:
        noise_level = 1.0

    if noise_level <= 0:
        noise_level = 1.0

    cent_mz = []
    cent_int = []

    for cmz, cint in zip(clusters_mz, clusters_int):
        if len(cmz) == 0:
            continue
        max_i = float(np.max(cint))
        if max_i / noise_level < sn_threshold:
            continue
        
        sum_i = float(np.sum(cint))
        c_mz = float(np.sum(cmz * cint) / sum_i)
        cent_mz.append(c_mz)
        cent_int.append(max_i)

    return np.array(cent_mz, dtype=np.float64), np.array(cent_int, dtype=np.float64)
