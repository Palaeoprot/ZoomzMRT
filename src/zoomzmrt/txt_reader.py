"""
Reader module for pre-summed text/TSV spectra (e.g. from MZSummer or mMass export).

Supports 2-column tab-, whitespace-, or comma-separated tables (m/z, intensity).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class TextSpectrumQC:
    """Quality control metrics for ingested 2-column text spectrum."""
    file_id: str
    total_points: int = 0
    min_mz: float = 0.0
    max_mz: float = 0.0
    base_peak_mz: float = 0.0
    base_peak_intensity: float = 0.0
    is_precalibrated: bool = False
    measured_resolving_power: float | None = None
    nominal_resolving_power: float = 250000.0


def parse_summed_text_spectrum(
    file_path: Path | str,
    dataset_id: str,
    instrument_name: str = "Waters SELECT SERIES MRT (Summed)",
    nominal_resolving_power: float = 250000.0,
    mz_min: float = 100.0,
    mz_max: float = 4000.0,
    source_type: str = "external",
) -> tuple[dict[str, Any], TextSpectrumQC]:
    """Parse a 2-column text or TSV spectrum file into standardized ZoomzMRT record.

    Args:
        file_path: Path to .txt or .tsv spectrum file.
        dataset_id: Dataset identifier.
        instrument_name: Instrument display label.
        nominal_resolving_power: Nominal resolving power (FWHM).
        mz_min: Minimum m/z cutoff.
        mz_max: Maximum m/z cutoff.
        source_type: Source category.

    Returns:
        tuple of (record_dict, qc_metrics)
    """
    path_obj = Path(file_path)
    fname = path_obj.name
    sample_id = re.sub(r"\.(txt|tsv|csv|mzml)$", "", fname, flags=re.IGNORECASE)

    is_lmc = bool(re.search(r"(_LMC|_lmc|corrected)", fname, flags=re.IGNORECASE))

    mzs = []
    intensities = []

    with open(path_obj, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", ";", "//")):
                continue
            # Split on tab, comma, or whitespace
            parts = re.split(r"[\t, ]+", stripped)
            if len(parts) < 2:
                continue
            try:
                m_val = float(parts[0])
                i_val = float(parts[1])
            except ValueError:
                # Header row (e.g. 'm/z\tIntensity')
                continue

            if mz_min <= m_val <= mz_max:
                mzs.append(m_val)
                intensities.append(i_val)

    qc = TextSpectrumQC(
        file_id=fname,
        total_points=len(mzs),
        is_precalibrated=is_lmc,
        nominal_resolving_power=nominal_resolving_power,
    )

    if not mzs:
        record = {
            "file_id": fname,
            "dataset_id": dataset_id,
            "source_type": source_type,
            "raw_path": str(file_path),
            "sample_id": sample_id,
            "scan_number": 1,
            "rt": None,
            "instrument": instrument_name,
            "mz": [],
            "intensity": [],
            "n_peaks": 0,
            "is_centroided": False,
            "extraction_strategy": "text_summed_profile",
            "instrument_status": "pre_summed_external",
            "instrument_serial": None,
            "schema_version": "0.1.0",
        }
        return record, qc

    arr_m = np.array(mzs, dtype=np.float64)
    arr_i = np.array(intensities, dtype=np.float64)

    # Sort by m/z
    sort_idx = np.argsort(arr_m)
    sorted_m = arr_m[sort_idx]
    sorted_i = arr_i[sort_idx]

    qc.min_mz = float(sorted_m[0])
    qc.max_mz = float(sorted_m[-1])

    b_idx = int(np.argmax(sorted_i))
    qc.base_peak_mz = float(sorted_m[b_idx])
    qc.base_peak_intensity = float(sorted_i[b_idx])

    record = {
        "file_id": fname,
        "dataset_id": dataset_id,
        "source_type": source_type,
        "raw_path": str(file_path),
        "sample_id": sample_id,
        "scan_number": 1,
        "rt": None,
        "instrument": instrument_name,
        "mz": [float(m) for m in sorted_m],
        "intensity": [float(i) for i in sorted_i],
        "n_peaks": len(sorted_m),
        "is_centroided": False,
        "extraction_strategy": "text_summed_profile",
        "instrument_status": "pre_summed_external",
        "instrument_serial": None,
        "schema_version": "0.1.0",
    }

    return record, qc
