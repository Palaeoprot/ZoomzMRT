# ZoomzMRT
**Date & Time:** 2026-09-27 10:10:00 (+02:00)

[![Tests & Quality Gates](https://github.com/Palaeoprot/ZoomzMRT/actions/workflows/tests.yml/badge.svg)](https://github.com/Palaeoprot/ZoomzMRT/actions/workflows/tests.yml)

High-resolution MS1 ingestion and deamidation analysis for ZooMS palaeoproteomics.

ZoomzMRT converts mzML/mzXML mass spectrometry data from high-resolution instruments into a standardized ZooMS Parquet representation, with instrument-specific lockmass calibration, resolution-aware peak area extraction, quality-control metrics, and direct resolved-peak deamidation measurements.

> **Status:** Beta / active development  
> **Python:** ≥ 3.10  
> **License:** GPL-3.0-or-later  
> **Repository:** [https://github.com/Palaeoprot/ZoomzMRT](https://github.com/Palaeoprot/ZoomzMRT)

---

## What ZoomzMRT Does

ZoomzMRT provides specialized processing pipelines for:

- **Waters SELECT SERIES MRT (Multi-Reflecting TOF)**
  - *Nominal capability:* 200,000–300,000+ FWHM.
  - MS1 mzML ingestion with scan-level laser shot filtering.
  - Per-scan lockmass tracking ([Glu1]-Fib, Leu-Enk, etc.) and multiplicative alignment.
  - Empirical resolving power estimation (*R* = *m* / FWHM) on reference peaks.
  - High-resolution co-addition and resolution-scaled peak centroiding.
- **MALDI-FTICR (Bruker solariX, Thermo FTMS)**
  - *Nominal capability:* 500,000–1,000,000+ FWHM.
  - Explicit multi-scan aggregation (`mean`, `sum`, or `none`).
  - Fine isotopic multiplet preservation (e.g. ¹³C vs ¹⁵N vs ³⁴S).
  - Accurate instrument and magnetic field metadata extraction.
- **Collagen Deamidation & PQI Analysis**
  - Direct background-subtracted trapezoidal peak area integration (`integrate_peak`).
  - At a nominal resolving power of 250,000 FWHM, the 19.339 mDa separation between deamidated *M*₀ (+0.9840 Da) and natural ¹³C₁ *M*₁ (+1.0034 Da) corresponds to ~3.2–4.4 FWHM across the ZooMS *m*/*z* range (1100–1600 Da); actual separation is evaluated per spectrum.
  - Exact elemental stoichiometry isotope QC (*C*, *H*, *N*, *O*, *S*).
  - Parchment Glutamine Index (PQI) and deamidation fraction calculations across diagnostic COL1A1/COL1A2 markers.
- **Standardized ZooMS Parquet Output**
  - Writes the current 16-column `mzPeakMS-ZooMS/0.1-draft` spectrum schema with ZSTD compression.
  - Generates comprehensive JSON metadata sidecars in `experiments_metadata/` logging both nominal and measured resolving power.
  - Exports sample-level `pqi_report.csv` and acquisition-level `qc_report.csv`.

---

## Installation

Clone the repository and install in editable mode:

```bash
git clone https://github.com/Palaeoprot/ZoomzMRT.git
cd ZoomzMRT
pip install -e .
```

For development and running tests:

```bash
pip install -e ".[dev]"
```

---

## Quick Start

### 1. Waters SELECT SERIES MRT
```bash
zoomzmrt path/to/mrt_mzml_dir \
    --dataset-id Mitchell_2026_MRT_Collagen \
    --instrument-type mrt \
    --lockmass-name glufib \
    --output-root ./output
```

### 2. MALDI-FTICR
```bash
zoomzmrt path/to/fticr_mzml_dir \
    --dataset-id Raymond_2024_FTICR_Bone \
    --instrument-type fticr \
    --instrument-name "Bruker solariX FT-ICR" \
    --aggregation mean \
    --output-root ./output
```

### 3. Continuous Profile Output
```bash
zoomzmrt path/to/mzml_dir \
    --dataset-id Profile_Acquisition \
    --instrument-type mrt \
    --profile \
    --output-root ./output
```

---

## Output Architecture

A processing run produces:

```
output/
├── zooms_ms1_maldi/
│   └── dataset_id=My_Dataset/
│       ├── spectra.parquet       # Standardized 16-column ZooMS Parquet
│       ├── pqi_report.csv        # Peptide-level & sample-level PQI metrics
│       └── qc_report.csv         # Lockmass shifts, scan counts, resolving power
│
└── experiments_metadata/
    └── My_Dataset.json           # JSON sidecar with full processing configuration
```

### `spectra.parquet`
Stores the current 16-column ZooMS spectrum schema (`file_id`, `dataset_id`, `sample_id`, `mz`, `intensity`, `n_peaks`, `is_centroided`, `extraction_strategy`, etc.) with ZSTD compression.

### `pqi_report.csv`
Contains sample-by-sample and peptide-by-peptide resolved deamidation ratios, PQI values, resolution status (`resolved`, `partially_resolved`, `unresolved`), and quality flags.

### `qc_report.csv`
Logs scan-level metrics including lockmass detection rate, median ppm error, Median Absolute Deviation (MAD), base peak stability, and empirical resolving power.

---

## Python API

```python
from pathlib import Path
from zoomzmrt import (
    parse_waters_mrt_mzml,
    parse_fticr_mzml,
    compute_high_res_deamidation,
    write_zooms_dataset,
    write_sidecar,
)

# 1. Parse and calibrate Waters MRT acquisition
record, qc = parse_waters_mrt_mzml(
    file_path=Path("sample_data") / "110924_h12_glufib.mzML",
    dataset_id="Example_MRT",
    lockmass_mz=1570.67742,
)

print(f"Measured Resolving Power: {qc.measured_resolving_power:,.0f} FWHM")
print(f"Lockmass Median Shift:    {qc.lockmass_median_ppm_error:+.2f} ppm (MAD: {qc.lockmass_mad_ppm:.2f} ppm)")

# 2. Compute resolved-peak deamidation & PQI using trapezoidal peak areas
summary = compute_high_res_deamidation(
    mz_arr=record["mz"],
    int_arr=record["intensity"],
    sample_id=record["sample_id"],
    dataset_id=record["dataset_id"],
    resolving_power=qc.measured_resolving_power or 250000.0,
)

print(f"Sample Median PQI: {summary.pqi_median:.3f}")
```

---

## Testing

Run the test suite with:

```bash
pytest -v
```

The test suite covers:
- True trapezoidal peak-area integration vs apex heights.
- Physical 19.339 mDa deamidation/¹³C separation benchmarks across resolving powers from 100k to 500k.
- Lockmass multiplicative calibration and empirical FWHM resolving power measurement.
- Resolution-aware peak centroiding.
- FT-ICR fine multiplet preservation and multi-scan co-addition.
- Peptide elemental stoichiometry and theoretical isotope distributions.
- 16-column Parquet schema validation and sidecar generation.

---

## Scientific Method & Documentation

- [Methodology & Instrument Physics](docs/methodology.md)
- [Deamidation & PQI Model](docs/deamidation.md)
- [ZooMS Parquet Schema](docs/parquet-schema.md)
- [Validation Benchmarks](docs/validation.md)

---

## Limitations

- High-resolution peak integration depends on accurate lockmass calibration and spectral acquisition quality.
- Resolving power varies across acquisitions and is measured empirically from data rather than assumed from instrument specifications.
- Parquet tables conform to the current draft `mzPeakMS-ZooMS/0.1-draft` specification.

---

## References

1. **Nair et al. (2023)** — *Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data.* Peer Community Journal, 3: e20.
2. **van Doorn et al. (2012)** — *Site-specific glutamine deamidation in ancient bone collagen.* Journal of Archaeological Science, 39(7): 2154–2163.
3. **Wilson et al. (2012)** — *Quantifying collagen degradation and glutamine deamidation in historical and archaeological parchment.* Analytical Chemistry, 84(21): 9051–9058.
4. **Welker et al. (2016)** — *Variations in glutamine deamidation for a prehistoric bone assemblage.* Journal of Proteomics, 152: 122–132.
5. **Yang et al. (2026)** — *MALDI Deamidation Score (MDS): A fast and flexible method for assessing deamidation in ZooMS data.* Journal of Proteomics, 314: 105406.

---

## License

ZoomzMRT is open source under the **GNU General Public License v3.0 or later** (GPL-3.0-or-later).
