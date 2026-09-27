# ZooMS Parquet Schema Specification
**Date & Time:** 2026-09-27 10:05:00 (+02:00)

## 1. 16-Column ZooMS Spectrum Schema (`mzPeakMS-ZooMS/0.1-draft`)

All high-resolution and MALDI-TOF MS1 data in `parquet_master/ZooMS_parquet/zooms_ms1_maldi/` is stored in Apache Parquet format using ZSTD compression (v2.6) with the following 16 standard columns:

| Column Name | Apache Arrow Data Type | Nullable | Description |
| :--- | :--- | :--- | :--- |
| `file_id` | `pa.string()` | No | Original source spectral file name (e.g. `110924_h12_glufib.mzML`). |
| `dataset_id` | `pa.string()` | No | Master dataset identifier partition key (e.g. `Mitchell_2026_MRT_Collagen`). |
| `source_type` | `pa.string()` | No | Data repository provenance (`zenodo`, `pride`, `external`, `ourdata`). |
| `raw_path` | `pa.string()` | Yes | Full file system path to source file during ingestion. |
| `sample_id` | `pa.string()` | No | Unique biological/archaeological specimen identifier. |
| `scan_number` | `pa.int64()` | No | Total or included MS1 scan count. |
| `rt` | `pa.float64()` | Yes | Retention time in seconds (null for direct infusion / MALDI spot). |
| `instrument` | `pa.string()` | Yes | Instrument model display name (e.g. `Waters SELECT SERIES MRT`). |
| `mz` | `pa.list_(pa.float64())` | No | Sorted vector of calibrated m/z values. |
| `intensity` | `pa.list_(pa.float64())` | No | Vector of peak intensities matching the `mz` vector. |
| `n_peaks` | `pa.int64()` | No | Peak count (must strictly match `len(mz)` and `len(intensity)`). |
| `is_centroided` | `pa.bool_()` | No | `True` for peak-centroided data; `False` for continuous profile data. |
| `extraction_strategy` | `pa.string()` | No | Extraction method (`mrt_summed_lockmass`, `fticr_centroid`, `maldi_direct`). |
| `instrument_status` | `pa.string()` | Yes | Provenance of instrument metadata (`from_header`, `from_metadata`, `user_supplied`). |
| `instrument_serial` | `pa.string()` | Yes | Mass spectrometer hardware serial number if present. |
| `schema_version` | `pa.string()` | No | Schema release identifier (`0.1.0`). |

---

## 2. Sidecar Metadata JSON Schema

Accompanying each Parquet partition is a JSON sidecar in `experiments_metadata/<dataset_id>.json`:

```json
{
  "dataset_id": "Mitchell_2026_MRT_Collagen",
  "n_files": 12,
  "n_spectra_rows": 12,
  "n_unique_samples": 12,
  "instrument": ["Waters SELECT SERIES MRT"],
  "mass_analyzer_type": "MRT",
  "resolution": {
    "nominal_resolving_power": 250000.0,
    "measured_resolving_power": 274812.5,
    "resolution_source": "empirical_peak_fwhm"
  },
  "processing_configuration": {
    "instrument_type": "mrt",
    "centroid": true,
    "sn_threshold": 3.0,
    "integration_method": "trapezoidal_area_baseline_subtracted",
    "lockmass": {
      "name": "glufib",
      "target_mz": 1570.67742,
      "median_ppm_error": 0.42,
      "mad_ppm_error": 0.18
    }
  },
  "pqi_glutamine_preservation": {
    "dataset_median_pqi": 0.784,
    "samples_assessed": 12,
    "quantitation_mode": "resolved_peak_area_ratio"
  },
  "software": {
    "name": "zoomzmrt",
    "version": "0.1.0",
    "spec_version": "mzPeakMS-ZooMS/0.1-draft"
  },
  "citation": {
    "doi": "10.xxxx/xxxxx",
    "title": "...",
    "authors": ["..."],
    "year": 2026
  }
}
```
