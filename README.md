# ZoomzMRT

**Date & Time:** 2026-09-27 09:45:00 (+02:00)

**High-Resolution MS1 Ingestion, Calibration, and Direct Deamidation (PQI) Engine for Waters SELECT SERIES MRT and MALDI-FTICR Mass Spectrometry in Palaeoproteomics.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Format: Parquet](https://img.shields.io/badge/Store-ZooMS__Parquet-success)](https://github.com/Palaeoprot/ZoomzPeak)
[![Standard: PSI-MS](https://img.shields.io/badge/Standard-PSI--MS%20mzPeak-orange)](https://www.psidev.info/)
[![Tests: Passing](https://img.shields.io/badge/Tests-4%20Passed-brightgreen.svg)]()

---

## 1. Introduction & Scientific Context

In palaeoproteomics and biocodicology, **ZooMS** (Zooarchaeology by Mass Spectrometry) has traditionally relied on low- to medium-resolution linear or reflectron MALDI-TOF mass spectrometers (e.g. Bruker autoflex, ultraflex, ABI 4800). These instruments provide resolving powers of *R* ≈ 1,000–15,000 FWHM and mass accuracies of 10–50 ppm, requiring quadratic polynomial recalibration and mathematical isotopic deconvolution to estimate post-translational modifications.

**`ZoomzMRT`** is the reference ingestion and quantitative analysis engine built specifically for **ultra-high-resolution MS1 mass spectrometry**, bridging modern multi-reflecting and Fourier-transform mass spectrometers with the standardized [ZooMS Parquet master repository](https://github.com/Palaeoprot/ZoomzPeak) (`mzPeakMS-ZooMS/0.1-draft`).

```
                                  ┌──────────────────────────────────────────────────────────┐
                                  │               High-Resolution Raw Deposits               │
                                  └────────────────────────────┬─────────────────────────────┘
                                                               │
                              ┌────────────────────────────────┴────────────────────────────────┐
                              ▼                                                                 ▼
          ┌───────────────────────────────────────┐                         ┌───────────────────────────────────────┐
          │       Waters SELECT SERIES MRT        │                         │              MALDI-FTICR              │
          │  Multi-Reflecting TOF (47-50m path)   │                         │    Cyclotron Frequency / Transients   │
          │  200,000 - 300,000+ FWHM (REM mode)   │                         │       500,000 - 1,000,000+ FWHM       │
          │      Part-Per-Billion Accuracy        │                         │       Fine Isotopic Multiplets        │
          └───────────────────┬───────────────────┘                         └───────────────────┬───────────────────┘
                              │                                                                 │
                              ▼                                                                 ▼
          ┌───────────────────────────────────────┐                         ┌───────────────────────────────────────┐
          │ • Sub-scan laser shot QC filtering    │                         │ • Fine multiplet preservation         │
          │ • Per-scan lockmass tracking (GluFib) │                         │ • Adaptive 1/m resolving power window │
          │ • Sub-100 ppb scan alignment          │                         │ • Narrow local-maxima extraction      │
          │ • High-res continuous co-addition     │                         │ • Magnet field metadata recording     │
          │ • Resolving-power aware centroiding   │                         └───────────────────┬───────────────────┘
          └───────────────────┬───────────────────┘                                             │
                              │                                                                 │
                              └────────────────────────────────┬────────────────────────────────┘
                                                               │
                                                               ▼
                                          ┌─────────────────────────────────────────┐
                                          │     Direct Baseline Deamidation (PQI)   │
                                          │   • Baseline resolved M0(Deam) vs M1    │
                                          │   • Direct peak area ratio integration  │
                                          │   • Isotopic 13C consistency QC         │
                                          └────────────────────┬────────────────────┘
                                                               │
                              ┌────────────────────────────────┴────────────────────────────────┐
                              ▼                                                                 ▼
          ┌───────────────────────────────────────┐                         ┌───────────────────────────────────────┐
          │            spectra.parquet            │                         │             Sidecar JSON              │
          │  16-column ZOOMS_SPECTRA schema       │                         │  experiments_metadata/<dataset_id>   │
          │  ZSTD compression, PSI-MS CV metadata │                         │  mass_analyzer_type: "MRT" / "FTICR"  │
          │  extraction_strategy:                 │                         │  nominal_resolving_power_fwhm: 250k+  │
          │  "mrt_summed_lockmass" /              │                         │  pqi_glutamine_preservation metrics   │
          │  "fticr_centroid" / "fticr_profile"   │                         │  lockmass_calibration statistics      │
          └───────────────────────────────────────┘                         └───────────────────────────────────────┘
```

---

## 2. Core Instrument Technologies & Physics

### 2.1. Waters SELECT SERIES MRT (Multi-Reflecting Time-of-Flight)
The Waters SELECT SERIES MRT uses dual opposing gridless electrostatic ion mirrors with periodic central electrostatic lenses:
- **Extended Flight Path:** Reflects ion packets back and forth along a zig-zag trajectory, extending the effective flight path up to **47–50 meters** within a benchtop instrument footprint.
- **Ultra-High Resolving Power:** Standard operations deliver **200,000 FWHM**. In **Resolution Enhancement Mode (REM)**, ions undergo a second pass through the mirror array, exceeding **300,000 FWHM**.
- **Part-Per-Billion Accuracy:** Achieves < 200–500 ppb mass accuracy across the entire mass range independent of scan speed.
- **Sub-Scan Multi-Shot Rasters:** Acquired across tens to hundreds of individual laser shots. `ZoomzMRT` inspects each sub-scan, tracks the lock-mass peak ([Glu1]-Fibrinopeptide B, monoisotopic [M+H]+ = 1570.67742 Da), aligns all scans into exact sub-100 ppb space, and co-adds intensity arrays without float splitting or peak smearing.

### 2.2. MALDI-FTICR (Fourier Transform Ion Cyclotron Resonance)
Bruker solariX (7T, 9.4T, 12T, 15T) and Thermo FTMS systems measure the cyclotron frequency of orbiting ions trapped in a high-field superconducting magnet:
- **Ultra-Ultra-High Resolving Power:** Resolving power *R* scales inversely with mass (*R* ∝ 1/*m*), regularly surpassing **500,000 to 1,000,000+ FWHM** in the peptide fingerprinting window (*m/z* 800–2,500).
- **Fine Isotopic Resolution:** Directly separates isobaric multiplets that differ by millidaltons (e.g. splitting ¹³C vs ¹⁵N vs ³⁴S).
- **Transient Preservation:** `ZoomzMRT` processes profile transients and narrow centroids without aggressive smoothing windows that would otherwise merge fine isotopic peaks.

---

## 3. High-Resolution Glutamine Deamidation & Direct PQI

### 3.1. The Low-Resolution Deamidation Bottleneck
Glutamine deamidation (*Q* → *E*) converts an amide group into a carboxylic acid, adding **+0.984016 Da** (mass of OH − NH₂).

In natural peptide isotopic distributions, the undeamidated peptide possesses natural isotopic peaks (*M*₁, *M*₂, ...), where *M*₁ contains one ¹³C atom (+1.003355 Da). The mass difference between the deamidated monoisotopic peak (*M*₀<sup>deam</sup>) and the natural undeamidated ¹³C₁ isotope (*M*₁<sup>undeam</sup>) is:

$$\Delta m = 1.003355 - 0.984016 = 0.019339\text{ Da (19.3 mDa, }\approx 12.3\text{ ppm at }m/z\text{ 1570)}$$

In low-resolution MALDI-TOF (*R* ≈ 1,000–15,000 FWHM, peak width FWHM ≈ 0.2–0.5 Da), these two peaks **completely merge** into an unresolved envelope. Classic low-resolution algorithms (van Doorn et al. 2012, Wilson et al. 2012, Nair/Bethencourt et al. 2023 PQI, Yang/Brown et al. 2026 MDS) must perform complex mathematical deconvolution (weighted least-squares unmixing of theoretical isotopic envelopes). This mathematical unmixing is vulnerable to baseline noise, chemical noise, and overlapping adducts.

### 3.2. Direct Baseline Measurement on Waters MRT & FT-ICR
On the Waters SELECT SERIES MRT (*R* ≥ 250,000–300,000 FWHM):
- The peak width at *m/z* 1500 is:
  $$\text{FWHM} = \frac{1500}{250,000} = 0.0060\text{ Da}$$
- The 0.01934 Da mass split represents **> 3.2 × FWHM** (more than 3 full peak widths apart).
- **The deamidated peak and the undeamidated natural ¹³C₁ isotope are fully baseline-resolved discrete peaks.**

```
Low-Resolution MALDI-TOF (R ≈ 5,000 FWHM):          Waters MRT (R ≥ 250,000 FWHM):
        Unresolved Composite Envelope                         Baseline Resolved Peaks
                ┌─────────┐                                      ┌───┐       ┌───┐
                │ Composite│                                     │   │       │   │
                │ M1(13C) +│                                     │M0 │       │M1 │
                │ M0(Deam) │                                     │Deam       │13C│
     ┌───┐      │ Overlap  │                               ┌───┐ │   │       │   │
     │   │      │          │                               │   │ │   │       │   │
     │M0 │      │          │                               │M0 │ │   │       │   │
     │Und│      │          │                               │Und│ │   │       │   │
  ───┴───┴──────┴──────────┴────                        ───┴───┴─┴───┴───────┴───┴───
     m/z        m/z + 1.0 Da                               m/z   +0.9840 Da  +1.0034 Da
                                                                 └── Δm = 19.3 mDa ──┘
```

`ZoomzMRT` measures deamidation directly by integrating the discrete peak areas:

$$\text{PQI (Glutamine Preservation Index)} = \frac{I(M_0^{\text{undeam}})}{I(M_0^{\text{undeam}}) + I(M_0^{\text{deam}})}$$

$$\% \text{Deamidation} = \frac{I(M_0^{\text{deam}})}{I(M_0^{\text{undeam}}) + I(M_0^{\text{deam}})} \times 100\%$$

### 3.3. Standard Diagnostic Collagen Markers
`ZoomzMRT` evaluates a curated panel of diagnostic Type I collagen (COL1A1 and COL1A2) markers:

| Marker Name | Gene | Target Sequence | Deamidation Site | Monoisotopic *m/z* (Undeamidated) | Monoisotopic *m/z* (Deamidated) | First Isotope *m/z* (¹³C₁ Undeam) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`COL1A1_P1105`** | COL1A1 | `GVQGPPGPAGPR` | Q502 | 1105.5807 | 1106.5647 | 1106.5841 |
| **`COL1A1_P1180`** | COL1A1 | `GQAGVMGFPGPK` | Q391 | 1180.5908 | 1181.5748 | 1181.5942 |
| **`COL1A1_P1427`** | COL1A1 | `GSEGPQGVRGEPGPAGPR` | Q178 | 1427.6972 | 1428.6812 | 1428.7006 |
| **`COL1A1_P1580`** | COL1A1 | `GATGAPGIAGAPGFPGAR` | Q / Gln-variant | 1580.7932 | 1581.7772 | 1581.7966 |
| **`COL1A2_P1706`** | COL1A2 | `GIPGEFGLPGPAGAR` | Diagnostic Marker | 1706.8837 | 1707.8677 | 1707.8871 |
| **`COL1A1_P2043`** | COL1A1 | `GAPGADGPAGAPGTPGPQGIAGQR` | Q774 | 2043.9806 | 2044.9646 | 2044.9840 |

---

## 4. Master Parquet Table Schema (`ZOOMS_SPECTRA`)

All datasets written by `ZoomzMRT` conform 100% to the 16-column specification (`mzPeakMS-ZooMS/0.1-draft`) with ZSTD compression:

| Field | Type | Description |
| :--- | :--- | :--- |
| `file_id` | `string` | Source filename (e.g. `110924_h12_glufib.mzML`) |
| `dataset_id` | `string` | Canonical dataset identifier (partition key) |
| `source_type` | `string` | `external` or `internal` |
| `raw_path` | `string` | Origin file system path |
| `sample_id` | `string` | Sanitized biological/artefact sample identifier |
| `scan_number` | `int64` | Sub-scans co-added or scan index |
| `rt` | `float64` | Retained for schema compatibility (null for MALDI) |
| `instrument` | `string` | Declared model (e.g. `Waters SELECT SERIES MRT`, `Bruker solariX FT-ICR`) |
| `mz` | `list<float64>` | Calibrated *m/z* array |
| `intensity` | `list<float64>` | Peak / profile intensity array |
| `n_peaks` | `int64` | Number of points in *m/z* array |
| `is_centroided` | `bool` | `True` for centroided peak lists, `False` for profile |
| `extraction_strategy` | `string` | `mrt_summed_lockmass`, `fticr_centroid`, `fticr_profile` |
| `instrument_status` | `string` | `from_metadata` or `from_header` |
| `instrument_serial` | `string` | Instrument serial number if present |
| `schema_version` | `string` | `0.1.0` |

---

## 5. Installation

```bash
# Clone the public repository
git clone https://github.com/Palaeoprot/ZoomzMRT.git
cd ZoomzMRT

# Install dependencies and package in editable mode
pip install -e .
```

### Dependencies
- Python >= 3.10
- `numpy >= 1.24`
- `scipy >= 1.10`
- `pandas >= 2.0`
- `pyarrow >= 14.0`
- `pyteomics >= 4.6`
- `pymzml >= 2.5`
- `pytest >= 7.4` (for test suite)

---

## 6. Command Line Interface (CLI)

`ZoomzMRT` installs the `zoomzmrt` executable entrypoint:

```bash
# Ingest Waters SELECT SERIES MRT files with [Glu1]-Fib lockmass correction
zoomzmrt "C:\data\waters_mrt_run" \
    --dataset-id Mitchell_2026_Waters_MRT \
    --instrument-type mrt \
    --lockmass-name glufib

# Ingest MALDI-FTICR files with custom instrument name
zoomzmrt "C:\data\fticr_deposit" \
    --dataset-id Buckley_2026_SolariX_FTICR \
    --instrument-type fticr \
    --instrument-name "Bruker solariX 7T FT-ICR" \
    --sn-threshold 3.0

# Store continuous profile arrays rather than centroided peaks
zoomzmrt "C:\data\highres_profiles" \
    --dataset-id Study_2026_Profile \
    --instrument-type mrt \
    --profile
```

### CLI Output Artifacts
Each run automatically generates:
1. `spectra.parquet`: The compressed, conformant ZooMS Parquet table.
2. `pqi_report.csv`: Peptide-by-peptide and sample-by-sample table of raw intensities, PQI scores, and isotopic consistency flags.
3. `qc_report.csv`: Sub-scan laser shot statistics, lockmass detection rate, and ppm shift distribution.
4. `experiments_metadata/<dataset_id>.json`: Truthful sidecar record capturing analyzer type, nominal resolution, lock mass statistics, and dataset median PQI.

---

## 7. Python API Reference

```python
from pathlib import Path
from zoomzmrt import (
    parse_waters_mrt_mzml,
    parse_fticr_mzml,
    compute_high_res_deamidation,
    deamidation_summary_to_dataframe,
    write_zooms_dataset,
    write_sidecar,
)

# 1. Parse and calibrate Waters MRT mzML with scan-level lockmass tracking
record, qc = parse_waters_mrt_mzml(
    file_path="110924_h12_glufib.mzML",
    dataset_id="Mitchell_2026_Waters_MRT",
    lockmass_mz=1570.67742,
    centroid=True,
)

print(f"File: {qc.file_id}")
print(f"Total sub-scans: {qc.total_scans}, Included: {qc.included_scans}")
print(f"Lockmass ppm shift: {qc.mean_ppm_shift:+.2f} +/- {qc.std_ppm_shift:.2f} ppm")
print(f"Centroid peaks: {qc.centroid_peaks}")

# 2. Compute direct high-resolution deamidation and PQI
pqi_summary = compute_high_res_deamidation(
    mz_arr=record["mz"],
    int_arr=record["intensity"],
    sample_id=record["sample_id"],
    dataset_id=record["dataset_id"],
)

print(f"Sample PQI Median: {pqi_summary.pqi_median:.3f}")
for r in pqi_summary.results:
    if r.pqi_fraction is not None:
        print(f"  {r.marker_name}: PQI = {r.pqi_fraction:.3f}, % Deam = {r.deamidation_fraction*100:.1f}%, Status = {r.quality_flag}")

# 3. Write directly to ZooMS Parquet master store
parquet_path = Path("C:/Users/matth/Documents/parquet_master/ZooMS_parquet/zooms_ms1_maldi/dataset_id=Mitchell_2026_Waters_MRT/spectra.parquet")
write_zooms_dataset([record], parquet_path, dataset_id="Mitchell_2026_Waters_MRT")
```

---

## 8. Verification & Test Suite

The package includes an automated test suite verifying synthetic multiplet resolution, FT-ICR isotopic preservation, PQI calculations, and full pipeline execution:

```bash
pytest -v
```

```
============================= test session starts =============================
platform win32 -- Python 3.14.1, pytest-9.1.1, pluggy-1.6.0
collected 4 items

tests/test_zoomzmrt.py::test_high_res_centroid_synthetic PASSED          [ 25%]
tests/test_zoomzmrt.py::test_fticr_fine_isotopic_preservation PASSED     [ 50%]
tests/test_zoomzmrt.py::test_direct_pqi_estimation PASSED                [ 75%]
tests/test_zoomzmrt.py::test_zoomzmrt_pipeline PASSED                    [100%]

======================== 4 passed, 1 warning in 1.46s =========================
```

---

## 9. Literature & Scientific References

- **Waters SELECT SERIES MRT:**
  - Waters Corporation. *SELECT SERIES MRT Mass Spectrometer Specification Sheet* (720007621EN).
  - Giles et al. (2020) *High-resolution multi-reflecting time-of-flight mass spectrometry*.
- **Parchment Glutamine Index (PQI):**
  - Nair, B., Rodríguez Palomo, I., Markussen, B., Wiuf, C., Fiddyment, S., Collins, M. J. (2023). *Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data*. **Peer Community Journal**, 3, e10. [doi:10.24072/pcjournal.230](https://doi.org/10.24072/pcjournal.230).
- **Glutamine Deamidation in Bone Collagen & ZooMS:**
  - van Doorn, N. L., Hollund, H., Collins, M. J. (2012). *Site-specific deamidation of glutamine: a new marker of bone collagen deterioration*. **Rapid Communications in Mass Spectrometry**, 26(19), 2319–2327. [doi:10.1002/rcm.6351](https://doi.org/10.1002/rcm.6351).
  - Wilson, J., van Doorn, N. L., Collins, M. J. (2012). *Assessing the Extent of Bone Degradation Using Glutamine Deamidation in Collagen*. **Analytical Chemistry**, 84(21), 9041–9048. [doi:10.1021/ac301333t](https://doi.org/10.1021/ac301333t).
  - Welker, F., et al. (2016). *Variations in glutamine deamidation for a Châtelperronian bone assemblage*. **PNAS**, 113(40), 11162–11167. [doi:10.1073/pnas.1605834113](https://doi.org/10.1073/pnas.1605834113).
  - Yang, F., Rodríguez Palomo, I., Nair, B. A. B., Brown, S. (2026). *MALDI Deamidation Score (MDS): A fast and flexible method for assessing deamidation in ZooMS data and its application to the Denisova Cave bone assemblage*. **Journal of Proteomics**, 324, 105577. [doi:10.1016/j.jprot.2025.105577](https://doi.org/10.1016/j.jprot.2025.105577).

---

## 10. License

GNU General Public License v3.0 or later ([LICENSE](LICENSE)).
Developed as part of the [PAASTA](https://paasta-community.github.io/) and [Palaeoprot](https://github.com/Palaeoprot) open-science community.
