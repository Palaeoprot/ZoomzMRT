# ZoomzMRT Validation & Benchmarks
**Date & Time:** 2026-09-27 10:05:00 (+02:00)

## 1. 19.339 mDa Separation Benchmark

To establish quantitative fidelity, synthetic spectra simulating the 19.339 mDa separation between monoisotopic deamidated *M*₀ (1106.564716 Da) and intact ¹³C₁ *M*₁ (1106.584055 Da) were evaluated across varying resolving powers:

| Resolving Power (*R* = *m* / FWHM) | Peak FWHM (at m/z 1105.58) | Separation / FWHM Ratio | Resolution Status | True PQI | Recovered PQI | Recovery Error |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **100,000 FWHM** | 11.05 mDa | 1.75 | `resolved` | 0.6667 | 0.6642 | -0.37% |
| **150,000 FWHM** | 7.37 mDa | 2.62 | `resolved` | 0.6667 | 0.6661 | -0.09% |
| **250,000 FWHM (MRT)** | 4.42 mDa | 4.37 | `resolved` | 0.6667 | 0.6667 | < 0.01% |
| **500,000 FWHM (FTICR)** | 2.21 mDa | 8.75 | `resolved` | 0.6667 | 0.6667 | < 0.01% |

At resolving power *R* ≥ 100,000 FWHM, peak-area integration recovers the true deamidation and PQI fraction within < 0.5% absolute error.

---

## 2. Lockmass Calibration Robustness

Lockmass multiplicative alignment evaluated on synthetic and experimental datasets:
- **Zero-drift baseline:** Calibrated m/z equals raw m/z with 0.00 ppm error.
- **Instrument drift recovery:** Known ±2.5 ppm linear mass drifts are recovered to < 10⁻⁷ relative error.
- **Robust outlier metrics:** Quality reports log median error and Median Absolute Deviation (MAD) to prevent single-scan outliers from skewing summary statistics.
