# ZoomzMRT

**Date & Time:** 2026-09-27 09:25:00 (+02:00)

**High-Resolution MS1 Ingestion Engine for Waters SELECT SERIES MRT and MALDI-FTICR Data in ZooMS Palaeoproteomics.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Format: Parquet](https://img.shields.io/badge/Store-ZooMS__Parquet-success)](https://github.com/Palaeoprot/ZoomzPeak)
[![Standard: PSI-MS](https://img.shields.io/badge/Standard-PSI--MS%20mzPeak-orange)](https://www.psidev.info/)

---

## 1. Overview

`ZoomzMRT` is a specialized ingestion and processing pipeline designed to handle ultra-high-resolution mass spectrometry deposits for ZooMS (Zooarchaeology by Mass Spectrometry) and MALDI peptide mass fingerprinting.

Unlike conventional linear and reflectron MALDI-TOF instruments (1,000–15,000 FWHM), high-resolution instruments generate data with orders-of-magnitude higher resolving power and part-per-billion (ppb) mass accuracy:

1. **Waters SELECT SERIES MRT (Multi-Reflecting Time-of-Flight):**
   - Employs gridless electrostatic ion mirrors extending ion flight path to 47–50 meters.
   - 200,000 to 300,000+ FWHM resolution (in Resolution Enhancement Mode / REM).
   - Stable sub-ppm and ppb mass accuracy (< 200–500 ppb).
   - Acquired across multiple sub-scans / laser shot rasters per well.
2. **MALDI-FTICR (Fourier Transform Ion Cyclotron Resonance):**
   - Ultra-high resolving power (500,000 to 1,000,000+ FWHM at m/z 1000).
   - Frequency-domain FID Fourier transform acquisitions.
   - Full baseline resolution of fine isotopic multiplets (¹³C vs ¹⁵N vs ³⁴S isobars, deamidation +0.9840 Da vs ¹³C₁ +1.0034 Da).

---

## 2. Core Processing Pipeline

```
                               ┌─────────────────────────────────────────────────┐
                               │           Raw / mzML High-Res Deposit           │
                               └────────────────────────┬────────────────────────┘
                                                        │
                         ┌──────────────────────────────┴──────────────────────────────┐
                         ▼                                                             ▼
        ┌───────────────────────────────────┐                       ┌───────────────────────────────────┐
        │     Waters SELECT SERIES MRT      │                       │            MALDI-FTICR            │
        │   Multi-Reflecting TOF Scans      │                       │     Cyclotron FID Transients      │
        └─────────────────┬─────────────────┘                       └─────────────────┬─────────────────┘
                          │                                                           │
                          ▼                                                           ▼
        ┌───────────────────────────────────┐                       ┌───────────────────────────────────┐
        │ • Laser shot signal QC filtering  │                       │ • Fine multiplet preservation     │
        │ • Per-scan lockmass tracking      │                       │ • Narrow adaptive 1/m clustering  │
        │   ([Glu1]-Fib 1570.67742 Da)      │                       │ • Magnet field metadata recording │
        │ • Sub-100 ppb scan alignment      │                       └─────────────────┬─────────────────┘
        │ • High-res continuous co-addition │                                         │
        │ • Resolving-power centroiding     │                                         │
        └─────────────────┬─────────────────┘                                         │
                          │                                                           │
                          └─────────────────────────────┬─────────────────────────────┘
                                                        │
                                                        ▼
                                    ┌───────────────────────────────────────┐
                                    │    Standardized Parquet & Sidecar     │
                                    │    (mzPeakMS-ZooMS/0.1-draft)         │
                                    └───────────────────────────────────────┘
```

---

## 3. Installation

```bash
git clone https://github.com/Palaeoprot/ZoomzMRT.git
cd ZoomzMRT
pip install -e .
```

Dependencies: `numpy`, `scipy`, `pandas`, `pyarrow`, `pyteomics`, `pymzml`.

---

## 4. Usage

### Command Line Interface (CLI)

```bash
# Ingest Waters MRT mzML files with lock mass correction (Glu-Fib)
zoomzmrt "path/to/waters_mrt_dir" \
    --dataset-id Mitchell_2026_Waters_MRT \
    --instrument-type mrt \
    --lockmass-name glufib

# Ingest MALDI-FTICR files preserving fine isotopic structure
zoomzmrt "path/to/fticr_dir" \
    --dataset-id Study_2026_SolariX_FTICR \
    --instrument-type fticr \
    --instrument-name "Bruker solariX 7T FT-ICR"
```

### Python API

```python
from zoomzmrt import parse_waters_mrt_mzml, write_zooms_dataset, write_sidecar
from pathlib import Path

# Parse and calibrate single Waters MRT mzML file
record, qc = parse_waters_mrt_mzml(
    file_path="110924_h12_glufib.mzML",
    dataset_id="Mitchell_2026_Waters_MRT",
    lockmass_mz=1570.67742,
    centroid=True,
)

print(f"Co-added peaks: {qc.centroid_peaks}, Mean ppm shift: {qc.mean_ppm_shift:+.2f} ppm")

# Write directly to master Parquet store
parquet_path = Path("parquet_master/ZooMS_parquet/zooms_ms1_maldi/dataset_id=Mitchell_2026_Waters_MRT/spectra.parquet")
write_zooms_dataset([record], parquet_path, dataset_id="Mitchell_2026_Waters_MRT")
```

---

## 5. Metadata Schema & Conformance

`ZoomzMRT` outputs 100% byte-compatible Parquet tables following the 16-column `ZOOMS_SPECTRA` schema:

| Column | Type | Description |
| :--- | :--- | :--- |
| `file_id` | `string` | Source filename |
| `dataset_id` | `string` | Partition dataset identifier |
| `source_type` | `string` | `external` or `internal` |
| `raw_path` | `string` | Origin file system path |
| `sample_id` | `string` | Normalized sample name |
| `scan_number` | `int64` | Sub-scans co-added / acquisition index |
| `rt` | `float64` | Retained for format compatibility (null for MALDI) |
| `instrument` | `string` | E.g. `Waters SELECT SERIES MRT`, `Bruker solariX FT-ICR` |
| `mz` | `list<float64>` | Calibrated m/z values |
| `intensity` | `list<float64>` | Peak / profile intensities |
| `n_peaks` | `int64` | Number of data points in array |
| `is_centroided` | `bool` | `True` for centroided, `False` for profile |
| `extraction_strategy` | `string` | `mrt_summed_lockmass`, `fticr_centroid`, `fticr_profile` |
| `instrument_status` | `string` | `from_metadata` or `from_header` |
| `instrument_serial` | `string` | Serial number if reported |
| `schema_version` | `string` | `0.1.0` |

---

## 6. License

GNU General Public License v3.0 or later ([LICENSE](LICENSE)).
Part of the [PAASTA](https://paasta-community.github.io/) and [Palaeoprot](https://github.com/Palaeoprot) open-science community.
