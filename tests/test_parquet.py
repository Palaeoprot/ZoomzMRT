"""
Unit tests for Parquet writing, schema validation, and sidecar metadata generation.
"""

import json
from pathlib import Path
import pyarrow.parquet as pq
import pytest

from zoomzmrt.parquet_writer import (
    SCHEMA_VERSION,
    SPEC_VERSION,
    ZOOMS_SPECTRA,
    validate_spectrum_record,
    write_sidecar,
    write_zooms_dataset,
)


def test_validate_spectrum_record_valid():
    rec = {
        "file_id": "test_1.mzML",
        "dataset_id": "Test_DS",
        "source_type": "external",
        "raw_path": "/path/to/test_1.mzML",
        "sample_id": "test_1",
        "scan_number": 1,
        "rt": None,
        "instrument": "Waters SELECT SERIES MRT",
        "mz": [1000.0, 1500.0],
        "intensity": [5000.0, 10000.0],
        "n_peaks": 2,
        "is_centroided": True,
        "extraction_strategy": "mrt_summed_lockmass",
        "instrument_status": "from_metadata",
        "instrument_serial": None,
        "schema_version": SCHEMA_VERSION,
    }
    # Should not raise
    validate_spectrum_record(rec)


def test_validate_spectrum_record_missing_field():
    rec = {
        "file_id": "test_1.mzML",
        "dataset_id": "Test_DS",
        "mz": [1000.0],
        "intensity": [5000.0],
        "n_peaks": 1,
    }
    with pytest.raises(ValueError, match="Record missing required schema fields"):
        validate_spectrum_record(rec)


def test_validate_spectrum_record_length_mismatch():
    rec = {
        "file_id": "test_1.mzML",
        "dataset_id": "Test_DS",
        "source_type": "external",
        "raw_path": "/path/to/test_1.mzML",
        "sample_id": "test_1",
        "scan_number": 1,
        "rt": None,
        "instrument": "Waters SELECT SERIES MRT",
        "mz": [1000.0, 1500.0],
        "intensity": [5000.0],  # Length 1 vs n_peaks 2
        "n_peaks": 2,
        "is_centroided": True,
        "extraction_strategy": "mrt_summed_lockmass",
        "instrument_status": "from_metadata",
        "instrument_serial": None,
        "schema_version": SCHEMA_VERSION,
    }
    with pytest.raises(ValueError, match="does not match intensity length"):
        validate_spectrum_record(rec)


def test_write_zooms_dataset_and_sidecar(tmp_path: Path):
    rec = {
        "file_id": "sample_a.mzML",
        "dataset_id": "Test_Dataset_2026",
        "source_type": "external",
        "raw_path": str(tmp_path / "sample_a.mzML"),
        "sample_id": "sample_a",
        "scan_number": 1,
        "rt": None,
        "instrument": "Waters SELECT SERIES MRT",
        "mz": [1105.5807, 1180.5908],
        "intensity": [15000.0, 8000.0],
        "n_peaks": 2,
        "is_centroided": True,
        "extraction_strategy": "mrt_summed_lockmass",
        "instrument_status": "from_metadata",
        "instrument_serial": None,
        "schema_version": SCHEMA_VERSION,
    }

    parquet_path = tmp_path / "zooms_ms1_maldi" / "dataset_id=Test_Dataset_2026" / "spectra.parquet"
    sidecar_dir = tmp_path / "experiments_metadata"

    total = write_zooms_dataset([rec], parquet_path, dataset_id="Test_Dataset_2026")
    assert total == 1
    assert parquet_path.exists()

    # Verify Parquet table
    table = pq.read_table(parquet_path)
    assert len(table) == 1
    assert table.num_columns == 16
    assert table.schema.metadata[b"zoomzpeak.spec_version"].decode("utf-8") == SPEC_VERSION

    # Generate sidecar
    sidecar_path = write_sidecar(
        dataset_id="Test_Dataset_2026",
        out_dir=sidecar_dir,
        parquet_path=parquet_path,
        metadata_overrides={"instrument": ["Waters SELECT SERIES MRT"]},
    )
    assert sidecar_path.exists()

    with open(sidecar_path, "r", encoding="utf-8") as f:
        meta = json.load(f)

    assert meta["dataset_id"] == "Test_Dataset_2026"
    assert meta["n_files"] == 1
    assert meta["software"]["name"] == "zoomzmrt"
    assert meta["software"]["version"] == "0.1.0"
