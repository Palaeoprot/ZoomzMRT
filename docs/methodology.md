# ZoomzMRT Methodology & Instrument Physics
**Date & Time:** 2026-09-27 10:05:00 (+02:00)

## 1. Physical Principle of High-Resolution MS1 Separation

In palaeoproteomics and ZooMS (Zooarchaeology by Mass Spectrometry), glutamine (*Q*) and asparagine (*N*) deamidation are paramount diagenetic indicators.

When a glutamine residue deamidates to glutamic acid (*Q* → *E*), the peptide undergoes a chemical conversion replacing an amine group (–NH₂) with a hydroxyl group (–OH):

- Deamidation mass shift: **+0.984016 Da** (monoisotopic mass difference)
- Natural ¹³C₁ isotope shift: **+1.003355 Da** (single-neutron mass difference)

The physical separation between the deamidated monoisotopic peak (*M*₀ deamidated) and the natural ¹³C₁ isotopic peak of the intact undeamidated peptide (*M*₁ undeamidated) is:

`Delta_m = 1.003355 Da - 0.984016 Da = 0.019339 Da (19.339 mDa)`

```
  Intact Peptide (M0)                  Deamidated Peptide (M0)      Intact Peptide (M1 13C1)
      m/z 1105.5807                          m/z 1106.5647               m/z 1106.5841
           │                                      │                           │
           │◄─────────── +0.984016 Da ───────────►│                           │
           │                                      │◄─────── 19.339 mDa ──────►│
           │◄───────────────────────────── +1.003355 Da ─────────────────────►│
```

On standard MALDI-TOF instruments (resolving power *R* ≈ 5,000–15,000 FWHM), a 19.339 mDa gap is completely unresolvable, coalescing into a single distorted peak envelope requiring complex mathematical deconvolution (such as PQI matrix inversion or least-squares isotope fitting).

On ultra-high-resolution mass spectrometers (*R* ≥ 200,000–500,000 FWHM), the peaks are physically baseline resolved:
- At *m*/*z* 1105 with *R* = 250,000 FWHM, peak width FWHM ≈ 0.0044 Da (4.4 mDa).
- The 19.339 mDa gap represents > 4.3× the peak FWHM, allowing direct physical separation and background-subtracted peak integration.

---

## 2. Mass Analyzer Architectures

### 2.1. Waters SELECT SERIES Multi-Reflecting Time-of-Flight (MRT)
The Waters SELECT SERIES MRT extends the physical ion flight path up to 47–50 metres within a compact footprint by reflecting ion packets back and forth between opposing gridless electrostatic ion mirrors (Resolution Enhanced Mode). 

Key processing considerations:
1. **Scan-Level Multiplicative Calibration:** Using internal lockmass standards (e.g. [Glu1]-Fibrinopeptide B, *m*/*z* 1570.67742), each individual laser shot / scan is aligned to sub-ppm precision:
   `mz_calibrated = mz_raw * (lockmass_theoretical / lockmass_observed)`
2. **Empirical Resolving Power Measurement:** Rather than asserting nominal factory values, resolving power *R* = *m* / FWHM is measured on the lockmass peak and reported in quality metrics.
3. **Resolution-Aware Extraction Windows:** Extraction windows scale with resolving power:
   `window_half_width = 1.5 * (m / R_measured)`

### 2.2. MALDI-FTICR (Fourier Transform Ion Cyclotron Resonance)
High-field FT-ICR instruments (Bruker solariX 7T to 15T) achieve resolving powers exceeding 500,000–1,000,000 FWHM.

Key processing considerations:
1. **Multi-Scan Aggregation:** Scans are combined using explicit aggregation semantics (`mean`, `sum`, or `none`) onto fine adaptive grids, preserving scan counts and single-shot fidelity.
2. **Fine Multiplet Preservation:** Preserves narrow isobaric multiplet splits (e.g., ¹³C vs ¹⁵N vs ³⁴S).

---

## 3. Peak-Area Integration Workflow

Rather than relying on single-bin apex peak heights (`argmax`), `ZoomzMRT` implements true numerical trapezoidal integration with local baseline correction (`integrate_peak`):

1. **Peak Apex Identification:** Identify apex *m*/*z* within the resolution-scaled search tolerance.
2. **Local Window Definition:** Window spanning [center – *w*, center + *w*], where *w* = 1.5 × (*m* / *R*).
3. **Baseline Estimation & Subtraction:** Local minimum intensity within the window is subtracted to isolate true peptide signal from chemical matrix baseline.
4. **Trapezoidal Numerical Integration:** Background-corrected intensities are integrated across the *m*/*z* vector.
5. **Centroiding:** Intensity-weighted *m*/*z* centroid is calculated from corrected intensities.
