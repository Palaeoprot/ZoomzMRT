"""
ZoomzMRT: High-Resolution MS1 Ingestion and Deamidation Engine for ZooMS Palaeoproteomics.
"""

from __future__ import annotations

from zoomzmrt.cli import run_pipeline
from zoomzmrt.deamidation import (
    DEFAULT_COLLAGEN_MARKERS,
    PeptideDeamidationResult,
    SampleDeamidationSummary,
    compute_high_res_deamidation,
    deamidation_summary_to_dataframe,
    find_peak_in_window,
    integrate_peak,
    peak_half_width,
)
from zoomzmrt.fticr_reader import (
    FTICRScanQC,
    parse_fticr_mzml,
)
from zoomzmrt.isotopes import (
    DeamidationMarker,
    ElementalComposition,
    peptide_composition,
)
from zoomzmrt.mrt_reader import (
    KNOWN_LOCK_MASSES,
    WatersMRTScanQC,
    measure_peak_resolving_power,
    parse_waters_mrt_mzml,
)
from zoomzmrt.parquet_writer import (
    SCHEMA_VERSION,
    SPEC_VERSION,
    ZOOMS_SPECTRA,
    dataset_id_for,
    validate_spectrum_record,
    write_sidecar,
    write_zooms_dataset,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "run_pipeline",
    "parse_waters_mrt_mzml",
    "WatersMRTScanQC",
    "measure_peak_resolving_power",
    "KNOWN_LOCK_MASSES",
    "parse_fticr_mzml",
    "FTICRScanQC",
    "compute_high_res_deamidation",
    "integrate_peak",
    "find_peak_in_window",
    "peak_half_width",
    "deamidation_summary_to_dataframe",
    "PeptideDeamidationResult",
    "SampleDeamidationSummary",
    "peptide_composition",
    "ElementalComposition",
    "DeamidationMarker",
    "DEFAULT_COLLAGEN_MARKERS",
    "write_zooms_dataset",
    "write_sidecar",
    "validate_spectrum_record",
    "dataset_id_for",
    "ZOOMS_SPECTRA",
    "SCHEMA_VERSION",
    "SPEC_VERSION",
]
