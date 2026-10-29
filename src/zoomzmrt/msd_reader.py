"""
Reader module for mMass .msd (Mass Spectrometry Data) XML files.

Parses zlib-compressed, base64-encoded profile spectrum arrays and peaklists.
"""

from __future__ import annotations

import base64
import re
import xml.etree.ElementTree as ET
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class MSDSpectrumQC:
    """Quality control metrics for ingested .msd file."""
    file_id: str
    sample_title: str = ""
    profile_points: int = 0
    peaklist_peaks: int = 0
    min_mz: float = 0.0
    max_mz: float = 0.0
    base_peak_mz: float = 0.0
    base_peak_intensity: float = 0.0
    measured_resolving_power: float | None = None
    nominal_resolving_power: float = 250000.0


def parse_mmass_msd(
    file_path: Path | str,
    dataset_id: str,
    instrument_name: str = "Waters SELECT SERIES MRT (mMass MSD)",
    nominal_resolving_power: float = 250000.0,
    prefer_peaklist: bool = False,
    source_type: str = "external",
) -> tuple[dict[str, Any], MSDSpectrumQC]:
    """Parse an mMass .msd XML file, decompressing binary spectrum arrays and peaklists.

    Args:
        file_path: Path to .msd file.
        dataset_id: Dataset identifier.
        instrument_name: Instrument display label.
        nominal_resolving_power: Nominal resolving power (FWHM).
        prefer_peaklist: If True, uses discrete centroided peaklist rather than raw profile.
        source_type: Source category.

    Returns:
        tuple of (record_dict, qc_metrics)
    """
    path_obj = Path(file_path)
    fname = path_obj.name
    sample_id = re.sub(r"\.msd$", "", fname, flags=re.IGNORECASE)

    tree = ET.parse(path_obj)
    root = tree.getroot()

    title_elem = root.find("description/title")
    sample_title = title_elem.text.strip() if title_elem is not None and title_elem.text else sample_id

    qc = MSDSpectrumQC(
        file_id=fname,
        sample_title=sample_title,
        nominal_resolving_power=nominal_resolving_power,
    )

    profile_mz = np.array([], dtype=np.float64)
    profile_int = np.array([], dtype=np.float64)

    spec_elem = root.find("spectrum")
    if spec_elem is not None:
        mz_elem = spec_elem.find("mzArray")
        int_elem = spec_elem.find("intArray")
        if mz_elem is not None and int_elem is not None and mz_elem.text and int_elem.text:
            mz_bytes = zlib.decompress(base64.b64decode(mz_elem.text.strip()))
            int_bytes = zlib.decompress(base64.b64decode(int_elem.text.strip()))
            prec = spec_elem.get("precision", "32")
            dtype_code = "<f8" if prec == "64" else "<f4"
            profile_mz = np.frombuffer(mz_bytes, dtype=dtype_code).astype(np.float64)
            profile_int = np.frombuffer(int_bytes, dtype=dtype_code).astype(np.float64)
            qc.profile_points = len(profile_mz)

    # Parse peaklist if present
    pl_mzs = []
    pl_ints = []
    pl_elem = root.find("peaklist")
    if pl_elem is not None:
        for p in pl_elem.findall("peak"):
            m = p.get("mz")
            i = p.get("intensity")
            if m is not None and i is not None:
                pl_mzs.append(float(m))
                pl_ints.append(float(i))
        qc.peaklist_peaks = len(pl_mzs)

    # Determine whether to use profile array or peaklist
    use_peaks = prefer_peaklist or (len(profile_mz) == 0 and len(pl_mzs) > 0)

    if use_peaks and len(pl_mzs) > 0:
        final_m = np.array(pl_mzs, dtype=np.float64)
        final_i = np.array(pl_ints, dtype=np.float64)
        is_centroided = True
        strategy = "msd_peaklist_centroid"
    elif len(profile_mz) > 0:
        final_m = profile_mz
        final_i = profile_int
        is_centroided = False
        strategy = "msd_binary_profile"
    else:
        final_m = np.array([], dtype=np.float64)
        final_i = np.array([], dtype=np.float64)
        is_centroided = False
        strategy = "msd_empty"

    if len(final_m) > 0:
        sort_idx = np.argsort(final_m)
        final_m = final_m[sort_idx]
        final_i = final_i[sort_idx]

        qc.min_mz = float(final_m[0])
        qc.max_mz = float(final_m[-1])
        b_idx = int(np.argmax(final_i))
        qc.base_peak_mz = float(final_m[b_idx])
        qc.base_peak_intensity = float(final_i[b_idx])

    record = {
        "file_id": fname,
        "dataset_id": dataset_id,
        "source_type": source_type,
        "raw_path": str(file_path),
        "sample_id": sample_title or sample_id,
        "scan_number": 1,
        "rt": None,
        "instrument": instrument_name,
        "mz": [float(m) for m in final_m],
        "intensity": [float(i) for i in final_i],
        "n_peaks": len(final_m),
        "is_centroided": is_centroided,
        "extraction_strategy": strategy,
        "instrument_status": "msd_parsed",
        "instrument_serial": None,
        "schema_version": "0.1.0",
    }

    return record, qc
