"""
ZoomzMRT: High-Resolution MS1 Ingestion Engine for Waters SELECT SERIES MRT and MALDI-FTICR.
"""

from zoomzmrt.deamidation import (
    DEFAULT_COLLAGEN_DEAMIDATION_MARKERS,
    PeptideDeamidationResult,
    SampleDeamidationSummary,
    compute_high_res_deamidation,
    deamidation_summary_to_dataframe,
)
from zoomzmrt.fticr_reader import FTICRScanQC, parse_fticr_mzml
from zoomzmrt.mrt_reader import (
    KNOWN_LOCK_MASSES,
    WatersMRTScanQC,
    parse_waters_mrt_mzml,
)
from zoomzmrt.parquet_writer import write_sidecar, write_zooms_dataset

__version__ = "0.1.0"

__all__ = [
    "parse_waters_mrt_mzml",
    "WatersMRTScanQC",
    "KNOWN_LOCK_MASSES",
    "parse_fticr_mzml",
    "FTICRScanQC",
    "compute_high_res_deamidation",
    "deamidation_summary_to_dataframe",
    "DEFAULT_COLLAGEN_DEAMIDATION_MARKERS",
    "SampleDeamidationSummary",
    "PeptideDeamidationResult",
    "write_zooms_dataset",
    "write_sidecar",
    "__version__",
]
