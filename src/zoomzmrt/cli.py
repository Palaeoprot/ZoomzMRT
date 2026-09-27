"""
Command-line interface for ZoomzMRT.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from zoomzmrt.deamidation import compute_high_res_deamidation, deamidation_summary_to_dataframe
from zoomzmrt.fticr_reader import parse_fticr_mzml
from zoomzmrt.mrt_reader import KNOWN_LOCK_MASSES, parse_waters_mrt_mzml
from zoomzmrt.parquet_writer import (
    dataset_id_for,
    write_sidecar,
    write_zooms_dataset,
)

DEFAULT_OUTPUT_ROOT = Path.cwd() / "output"


def detect_instrument_type(files: Sequence[Path]) -> str:
    """Auto-detect instrument analyzer from file stems."""
    names_str = " ".join(f.name.lower() for f in files)
    if any(k in names_str for k in ["mrt", "select_series", "glufib", "efib"]):
        return "mrt"
    if any(k in names_str for k in ["fticr", "solarix", "solari_x", "7t", "9.4t", "12t", "15t", "ftms"]):
        return "fticr"
    raise ValueError(
        "Could not automatically detect instrument type from file names. "
        "Please specify --instrument-type explicitly (e.g. --instrument-type mrt or --instrument-type fticr)."
    )


def run_pipeline(
    input_path: Path | str,
    dataset_id: str | None = None,
    output_root: Path | str = DEFAULT_OUTPUT_ROOT,
    instrument_type: str = "auto",
    instrument_name: str | None = None,
    lockmass_name: str = "glufib",
    lockmass_mz: float | None = None,
    lockmass_window_da: float = 0.3,
    centroid: bool = True,
    sn_threshold: float = 3.0,
    base_peak_threshold: float = 10000.0,
    aggregation: str = "mean",
    source_type: str = "external",
    generate_qc_report: bool = True,
) -> dict[str, Any]:
    """Execute complete high-resolution ingestion pipeline."""
    in_path = Path(input_path)
    out_root = Path(output_root)

    if not in_path.exists():
        raise FileNotFoundError(f"Input path does not exist: {in_path}")

    ds_id = dataset_id or dataset_id_for(in_path)

    if in_path.is_file():
        files = [in_path]
    else:
        found = set()
        for ext in ("*.mzML", "*.mzml", "*.mzXML", "*.mzxml"):
            for p in in_path.rglob(ext):
                found.add(p.resolve())
        files = sorted(list(found))

    if not files:
        raise ValueError(f"No .mzML or .mzXML spectral files found in {in_path}")

    print(f"\n=======================================================")
    print(f" ZoomzMRT High-Resolution MS1 Ingestion Engine")
    print(f"=======================================================")
    print(f" Dataset ID:       {ds_id}")
    print(f" Input path:       {in_path}")
    print(f" Files discovered: {len(files)}")
    print(f" Target root:      {out_root}")

    if instrument_type == "auto":
        resolved_type = detect_instrument_type(files)
        print(f" Auto-detected type: {resolved_type.upper()}")
    else:
        resolved_type = instrument_type.lower()
        print(f" Specified type:     {resolved_type.upper()}")

    if lockmass_mz is None:
        lm_mz = KNOWN_LOCK_MASSES.get(lockmass_name.lower(), 1570.67742)
    else:
        lm_mz = float(lockmass_mz)

    if resolved_type == "mrt":
        inst_display = instrument_name or "Waters SELECT SERIES MRT"
        nominal_res = 250000.0
        analyzer_type = "MRT"
    elif resolved_type == "fticr":
        inst_display = instrument_name or "Bruker solariX FT-ICR"
        nominal_res = 500000.0
        analyzer_type = "FTICR"
    else:
        inst_display = instrument_name or "MALDI-TOF"
        nominal_res = 15000.0
        analyzer_type = "TOF"

    print(f" Instrument:       {inst_display} (Analyzer: {analyzer_type})")
    print(f" Nominal Res:      {nominal_res:,.0f} FWHM")
    if resolved_type == "mrt":
        print(f" Lockmass:         {lm_mz:.5f} Da ({lockmass_name})")

    records = []
    qc_data = []
    measured_res_list = []

    t0 = time.time()
    for idx, f in enumerate(files, 1):
        print(f" [{idx}/{len(files)}] Processing: {f.name}...", end="", flush=True)
        if resolved_type == "mrt":
            rec, qc = parse_waters_mrt_mzml(
                file_path=f,
                dataset_id=ds_id,
                lockmass_mz=lm_mz,
                lockmass_window_da=lockmass_window_da,
                base_peak_threshold=base_peak_threshold,
                centroid=centroid,
                sn_threshold=sn_threshold,
                source_type=source_type,
                instrument_name=inst_display,
                nominal_resolving_power=nominal_res,
            )
            records.append(rec)
            qc_data.append(asdict(qc))
            if qc.measured_resolving_power:
                measured_res_list.append(qc.measured_resolving_power)
            print(f" Done ({rec['n_peaks']} peaks, shift: {qc.mean_ppm_shift:+.2f} ppm)")
        elif resolved_type == "fticr":
            rec, qc = parse_fticr_mzml(
                file_path=f,
                dataset_id=ds_id,
                centroid=centroid,
                sn_threshold=sn_threshold,
                aggregation=aggregation,  # type: ignore
                source_type=source_type,
                instrument_name=inst_display,
                nominal_resolving_power=nominal_res,
            )
            records.append(rec)
            qc_data.append(asdict(qc))
            if qc.measured_resolving_power:
                measured_res_list.append(qc.measured_resolving_power)
            print(f" Done ({rec['n_peaks']} peaks, fine clusters: {qc.n_mz_gaps_3_25_mda})")

    elapsed = time.time() - t0
    print(f"\nCompleted spectral processing in {elapsed:.2f}s.")

    parquet_dir = out_root / "zooms_ms1_maldi" / f"dataset_id={ds_id}"
    parquet_path = parquet_dir / "spectra.parquet"
    metadata_dir = out_root / "experiments_metadata"
    sidecar_path = metadata_dir / f"{ds_id}.json"

    dataset_measured_res = float(np.median(measured_res_list)) if measured_res_list else None

    extra_meta = {
        "zoomzpeak.instrument_analyzer_type": analyzer_type,
        "zoomzpeak.nominal_resolving_power": str(nominal_res),
        "zoomzpeak.measured_resolving_power": str(dataset_measured_res) if dataset_measured_res else "unmeasured",
        "zoomzpeak.lockmass_standard": lockmass_name if resolved_type == "mrt" else "none",
        "zoomzpeak.lockmass_mz": str(lm_mz) if resolved_type == "mrt" else "none",
    }

    n_rows = write_zooms_dataset(
        rows=records,
        out_path=parquet_path,
        dataset_id=ds_id,
        extra_metadata=extra_meta,
    )

    # Compute high-resolution PQI and deamidation metrics with trapezoidal integration
    deam_summaries = []
    for rec, qc_dict in zip(records, qc_data):
        rec_res = qc_dict.get("measured_resolving_power") or nominal_res
        d_summary = compute_high_res_deamidation(
            mz_arr=rec["mz"],
            int_arr=rec["intensity"],
            sample_id=rec["sample_id"],
            dataset_id=ds_id,
            resolving_power=rec_res,
        )
        deam_summaries.append(d_summary)

    pqi_df = deamidation_summary_to_dataframe(deam_summaries)
    pqi_csv_path = parquet_dir / "pqi_report.csv"
    pqi_df.to_csv(pqi_csv_path, index=False)
    print(f" PQI Report saved: {pqi_csv_path}")

    # Compute high-res sidecar overrides with full calibration and resolution provenance
    median_shifts = [q.get("lockmass_median_ppm_error") for q in qc_data if q.get("lockmass_median_ppm_error") is not None]
    overall_median_shift = float(np.median(median_shifts)) if median_shifts else None
    overall_mad_shift = float(np.median([q.get("lockmass_mad_ppm") for q in qc_data if q.get("lockmass_mad_ppm") is not None])) if median_shifts else None

    valid_sample_pqis = [s.pqi_median for s in deam_summaries if s.pqi_median is not None]
    dataset_pqi_median = float(np.median(valid_sample_pqis)) if valid_sample_pqis else None

    sidecar_overrides = {
        "instrument": [inst_display],
        "mass_analyzer_type": analyzer_type,
        "resolution": {
            "nominal_resolving_power": nominal_res,
            "measured_resolving_power": dataset_measured_res,
            "resolution_source": "empirical_peak_fwhm" if dataset_measured_res else "nominal_specification",
        },
        "processing_configuration": {
            "instrument_type": resolved_type,
            "centroid": centroid,
            "sn_threshold": sn_threshold,
            "integration_method": "trapezoidal_area_baseline_subtracted",
            "lockmass": {
                "name": lockmass_name if resolved_type == "mrt" else None,
                "target_mz": lm_mz if resolved_type == "mrt" else None,
                "median_ppm_error": overall_median_shift,
                "mad_ppm_error": overall_mad_shift,
            } if resolved_type == "mrt" else None,
        },
        "pqi_glutamine_preservation": {
            "dataset_median_pqi": dataset_pqi_median,
            "samples_assessed": len(valid_sample_pqis),
            "quantitation_mode": "resolved_peak_area_ratio",
        },
        "ingestion_engine": "zoomzmrt",
        "ingestion_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }

    written_sidecar = write_sidecar(
        dataset_id=ds_id,
        out_dir=metadata_dir,
        parquet_path=parquet_path,
        metadata_overrides=sidecar_overrides,
    )

    if generate_qc_report:
        qc_df = pd.DataFrame(qc_data)
        qc_csv_path = parquet_dir / "qc_report.csv"
        qc_df.to_csv(qc_csv_path, index=False)
        print(f" QC Report saved:  {qc_csv_path}")

    print(f"\n=======================================================")
    print(f" Ingestion Summary for {ds_id}")
    print(f"=======================================================")
    print(f" Rows written:     {n_rows}")
    print(f" Parquet dataset:  {parquet_path}")
    print(f" Sidecar JSON:     {written_sidecar}")
    print(f"=======================================================\n")

    return {
        "dataset_id": ds_id,
        "n_rows": n_rows,
        "parquet_path": str(parquet_path),
        "sidecar_path": str(written_sidecar),
    }


def main():
    parser = argparse.ArgumentParser(description="ZoomzMRT: High-Resolution MS1 Ingestion Engine")
    parser.add_argument("input_path", type=str, help="Input directory containing mzML/mzXML files")
    parser.add_argument("--dataset-id", "-d", type=str, default=None, help="Canonical dataset_id")
    parser.add_argument("--output-root", "-o", type=str, default=str(DEFAULT_OUTPUT_ROOT), help="Output directory root")
    parser.add_argument("--instrument-type", "-t", choices=["mrt", "fticr", "tof", "auto"], default="auto", help="Instrument analyzer type")
    parser.add_argument("--instrument-name", type=str, default=None, help="Custom instrument display name")
    parser.add_argument("--lockmass-name", type=str, default="glufib", help="Standard lockmass name (glufib, leuenk, bradykinin)")
    parser.add_argument("--lockmass-mz", type=float, default=None, help="Custom theoretical lockmass m/z")
    parser.add_argument("--sn-threshold", type=float, default=3.0, help="Signal-to-noise peak threshold")
    parser.add_argument("--aggregation", choices=["mean", "sum", "first", "none"], default="mean", help="Multi-scan aggregation mode for FTICR (mean, sum, first)")
    parser.add_argument("--profile", action="store_true", help="Store continuous profile rather than centroided peaks")

    args = parser.parse_args()

    run_pipeline(
        input_path=args.input_path,
        dataset_id=args.dataset_id,
        output_root=args.output_root,
        instrument_type=args.instrument_type,
        instrument_name=args.instrument_name,
        lockmass_name=args.lockmass_name,
        lockmass_mz=args.lockmass_mz,
        centroid=not args.profile,
        sn_threshold=args.sn_threshold,
        aggregation=args.aggregation,
    )


if __name__ == "__main__":
    main()
