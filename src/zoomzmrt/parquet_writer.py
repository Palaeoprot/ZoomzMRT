"""
Parquet and Sidecar Metadata Writer for High-Resolution ZooMS.

Conforms to the 16-column ZOOMS_SPECTRA schema with ZSTD compression and PSI-MS metadata.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

SCHEMA_VERSION = "0.1.0"
SPEC_VERSION = "mzPeakMS-ZooMS/0.1-draft"

ZOOMS_SPECTRA = pa.schema(
    [
        ("file_id", pa.string()),
        ("dataset_id", pa.string()),
        ("source_type", pa.string()),
        ("raw_path", pa.string()),
        ("sample_id", pa.string()),
        ("scan_number", pa.int64()),
        ("rt", pa.float64()),
        ("instrument", pa.string()),
        ("mz", pa.list_(pa.float64())),
        ("intensity", pa.list_(pa.float64())),
        ("n_peaks", pa.int64()),
        ("is_centroided", pa.bool_()),
        ("extraction_strategy", pa.string()),
        ("instrument_status", pa.string()),
        ("instrument_serial", pa.string()),
        ("schema_version", pa.string()),
    ]
)


def dataset_id_for(path_or_name: Path | str) -> str:
    """Derive a canonical dataset_id from a folder name or string."""
    name = path_or_name.name if isinstance(path_or_name, Path) else path_or_name
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def write_zooms_dataset(
    rows: Iterable[dict[str, Any]],
    out_path: Path,
    dataset_id: str,
    extra_metadata: dict[str, str] | None = None,
    batch_size: int = 50,
) -> int:
    """Write an iterable of ZooMS spectrum records to Parquet with full ZoomzPeak metadata."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    schema = ZOOMS_SPECTRA

    meta = {
        "zoomzpeak.spec_version": SPEC_VERSION,
        "zoomzpeak.table_kind": "zooms_ms1_maldi_spectra",
        "zoomzpeak.builder": "zoomzmrt.parquet_writer",
        "zoomzpeak.cv_reference": "PSI-MS",
        "zoomzpeak.schema_version": SCHEMA_VERSION,
        "zoomzpeak.dataset_id": dataset_id,
        "zoomzpeak.created_utc": datetime.now(timezone.utc).isoformat(),
    }
    if extra_metadata:
        meta.update(extra_metadata)

    schema_with_meta = schema.with_metadata(meta)

    writer = None
    batch = []
    total = 0

    try:
        for row in rows:
            if row.get("dataset_id") != dataset_id:
                row["dataset_id"] = dataset_id
            if "schema_version" not in row or not row["schema_version"]:
                row["schema_version"] = SCHEMA_VERSION

            batch.append(row)
            if len(batch) >= batch_size:
                table = pa.Table.from_pylist(batch, schema=schema).replace_schema_metadata(meta)
                if writer is None:
                    writer = pq.ParquetWriter(
                        out_path,
                        schema_with_meta,
                        compression="zstd",
                        version="2.6",
                    )
                writer.write_table(table)
                total += len(batch)
                batch = []

        if batch:
            table = pa.Table.from_pylist(batch, schema=schema).replace_schema_metadata(meta)
            if writer is None:
                writer = pq.ParquetWriter(
                    out_path,
                    schema_with_meta,
                    compression="zstd",
                    version="2.6",
                )
            writer.write_table(table)
            total += len(batch)
    finally:
        if writer is not None:
            writer.close()

    return total


def write_sidecar(
    dataset_id: str,
    out_dir: Path,
    parquet_path: Path,
    metadata_overrides: dict[str, Any] | None = None,
) -> Path:
    """Generate or update the experiments_metadata/<dataset_id>.json sidecar."""
    out_dir.mkdir(parents=True, exist_ok=True)
    sidecar_path = out_dir / f"{dataset_id}.json"

    table = pq.read_table(parquet_path)
    df = table.select(
        ["file_id", "sample_id", "instrument", "source_type", "extraction_strategy"]
    ).to_pandas()

    instruments = sorted(df["instrument"].dropna().unique().tolist())
    source_types = sorted(df["source_type"].dropna().unique().tolist())
    strategies = sorted(df["extraction_strategy"].dropna().unique().tolist())

    sidecar_data = {
        "dataset_id": dataset_id,
        "n_files": int(df["file_id"].nunique()),
        "n_spectra_rows": int(len(df)),
        "n_unique_samples": int(df["sample_id"].nunique()),
        "instrument": instruments,
        "source_type": source_types,
        "extraction_strategy": strategies,
        "_metadata_source": {
            "n_files": "derived_from_parquet",
            "n_spectra_rows": "derived_from_parquet",
            "instrument": "from_metadata" if instruments else None,
        },
        "citation": None,
        "species": None,
        "sample_provenance": None,
        "_enrichment_note": (
            "citation/species/sample_provenance intentionally left blank -- these need "
            "the actual publication looked up per dataset_id, not guessed from the folder "
            "name. Fill in from the source paper before treating this sidecar as complete."
        ),
    }

    if metadata_overrides:
        sidecar_data.update(metadata_overrides)

    sidecar_path.write_text(json.dumps(sidecar_data, indent=2), encoding="utf-8")
    return sidecar_path
