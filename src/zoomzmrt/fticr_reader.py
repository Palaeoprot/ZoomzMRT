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
import sqlite3
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


def parse_bruker_fticr_d(
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
    """Parse a Bruker raw .d directory (solariX FT-ICR BAF/SQLite).

    Extracts high-resolution centroid peaks or continuous profile spectra directly
    from analysis.baf using pyBaf2Sql C-API wrappers, avoiding external conversion artifacts.
    """
    try:
        from pyBaf2Sql.baf import close_storage, open_storage, read_double
        from pyBaf2Sql.init_baf2sql import init_baf2sql_api
    except ImportError as exc:
        raise ImportError(
            "pyBaf2Sql is required to read native Bruker .d directories. "
            "Please install it using: pip install git+https://github.com/gtluu/pyBaf2Sql.git"
        ) from exc

    path_obj = Path(file_path)
    if path_obj.is_file() and path_obj.name.lower() in ("analysis.baf", "analysis.sqlite"):
        d_dir = path_obj.parent
    else:
        d_dir = path_obj

    fname = d_dir.name
    sample_id = d_dir.stem if d_dir.suffix.lower() == ".d" else d_dir.name

    qc = FTICRScanQC(
        file_id=fname,
        instrument_name=instrument_name,
        nominal_resolving_power=nominal_resolving_power,
        magnetic_field_tesla=magnetic_field_tesla,
    )

    dll = init_baf2sql_api()
    handle = open_storage(dll, str(d_dir), raw_calibration=False)
    if handle == 0:
        handle = open_storage(dll, str(d_dir), raw_calibration=True)
    if handle == 0:
        return {
            "file_id": fname,
            "dataset_id": dataset_id,
            "source_type": source_type,
            "raw_path": str(file_path),
            "sample_id": sample_id,
            "scan_number": 0,
            "rt": None,
            "instrument": instrument_name,
            "mz": [],
            "intensity": [],
            "n_peaks": 0,
            "is_centroided": centroid,
            "extraction_strategy": "fticr_bruker_baf_error",
            "instrument_status": "unreadable_baf_storage",
            "instrument_serial": None,
            "schema_version": "0.1.0",
        }, qc


    sqlite_path = d_dir / "analysis.sqlite"
    conn = None
    all_mz_scans = []
    all_int_scans = []
    prof_mz_sample = None
    prof_int_sample = None
    inst_header = None
    inst_status = "from_metadata"

    try:
        conn = sqlite3.connect(sqlite_path)
        cur = conn.cursor()

        # Check instrument metadata properties
        try:
            cur.execute("SELECT Key, Value FROM Properties;")
            props = dict(cur.fetchall())
            sw = props.get("AcquisitionSoftware")
            vendor = props.get("AcquisitionSoftwareVendor")
            if sw and vendor:
                inst_header = f"{vendor} {sw}"
                inst_status = "from_header"
        except sqlite3.OperationalError:
            pass

        cur.execute("SELECT Id, Rt, LineMzId, LineIntensityId, ProfileMzId, ProfileIntensityId FROM Spectra;")
        rows = cur.fetchall()
        qc.total_scans = len(rows)

        for row in rows:
            spec_id, rt, line_mz_id, line_int_id, prof_mz_id, prof_int_id = row

            if centroid:
                if line_mz_id != 0:
                    lmz = np.array(read_double(dll, handle, line_mz_id), dtype=np.float64)
                    lint = np.array(read_double(dll, handle, line_int_id), dtype=np.float64)
                    mask = (lmz >= mz_min) & (lmz <= mz_max)
                    if np.any(mask):
                        all_mz_scans.append(lmz[mask])
                        all_int_scans.append(lint[mask])
                        qc.included_scans += 1
                elif prof_mz_id != 0:
                    pmz = np.array(read_double(dll, handle, prof_mz_id), dtype=np.float64)
                    pint = np.array(read_double(dll, handle, prof_int_id), dtype=np.float64)
                    mask = (pmz >= mz_min) & (pmz <= mz_max)
                    if np.any(mask):
                        cmz, cint = _fticr_centroid(
                            pmz[mask], pint[mask], resolving_power=nominal_resolving_power, sn_threshold=sn_threshold
                        )
                        all_mz_scans.append(cmz)
                        all_int_scans.append(cint)
                        qc.included_scans += 1
            else:  # Profile requested
                if prof_mz_id != 0:
                    pmz = np.array(read_double(dll, handle, prof_mz_id), dtype=np.float64)
                    pint = np.array(read_double(dll, handle, prof_int_id), dtype=np.float64)
                    mask = (pmz >= mz_min) & (pmz <= mz_max)
                    if np.any(mask):
                        all_mz_scans.append(pmz[mask])
                        all_int_scans.append(pint[mask])
                        qc.included_scans += 1
                elif line_mz_id != 0:
                    lmz = np.array(read_double(dll, handle, line_mz_id), dtype=np.float64)
                    lint = np.array(read_double(dll, handle, line_int_id), dtype=np.float64)
                    mask = (lmz >= mz_min) & (lmz <= mz_max)
                    if np.any(mask):
                        all_mz_scans.append(lmz[mask])
                        all_int_scans.append(lint[mask])
                        qc.included_scans += 1

            # Retain first profile scan for empirical resolution calculation if present
            if prof_mz_sample is None and prof_mz_id != 0:
                prof_mz_sample = np.array(read_double(dll, handle, prof_mz_id), dtype=np.float64)
                prof_int_sample = np.array(read_double(dll, handle, prof_int_id), dtype=np.float64)

    finally:
        close_storage(dll, handle, conn)

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
            "extraction_strategy": "fticr_bruker_baf_centroid" if centroid else "fticr_bruker_baf_profile",
            "instrument_status": inst_status,
            "instrument_serial": None,
            "schema_version": "0.1.0",
        }, qc

    if len(all_mz_scans) == 1 or aggregation in ("first", "none"):
        raw_mz = all_mz_scans[0]
        raw_int = all_int_scans[0]
    else:
        concat_mz = np.concatenate(all_mz_scans)
        concat_int = np.concatenate(all_int_scans)
        s_idx = np.argsort(concat_mz)
        sorted_mz = concat_mz[s_idx]
        sorted_int = concat_int[s_idx]

        scale = (1.0 / len(all_mz_scans)) if aggregation == "mean" else 1.0
        raw_mz, raw_int = _fticr_coadd(
            sorted_mz,
            sorted_int,
            resolving_power=nominal_resolving_power,
            scale_factor=scale,
        )

    # Measure empirical resolving power
    if len(raw_int) > 0:
        b_idx = int(np.argmax(raw_int))
        b_mz = raw_mz[b_idx]
        if prof_mz_sample is not None and prof_int_sample is not None:
            qc.measured_resolving_power = measure_peak_resolving_power(
                prof_mz_sample, prof_int_sample, target_mz=b_mz, window_da=0.1
            )
        if qc.measured_resolving_power is None:
            qc.measured_resolving_power = measure_peak_resolving_power(raw_mz, raw_int, target_mz=b_mz)

    final_mz = raw_mz
    final_int = raw_int
    extraction_strat = "fticr_bruker_baf_centroid" if centroid else "fticr_bruker_baf_profile"

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
