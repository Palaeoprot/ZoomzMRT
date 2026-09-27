# High-Resolution Deamidation & Parchment Glutamine Index (PQI)
**Date & Time:** 2026-09-27 10:05:00 (+02:00)

## 1. Diagnostic Collagen Markers

In bone, antler, and parchment palaeoproteomics, type I collagen (COL1A1 and COL1A2) contains several well-characterized glutamine (*Q*) sites whose deamidation rates reflect thermal age, post-mortem taphonomy, and preservation state.

| Marker Name | Protein Gene | Sequence | Target Site | Hyp Count | Monoisotopic m/z (M0) | Deamidated m/z (+0.9840 Da) | Natural 13C1 m/z (+1.0034 Da) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **COL1A1_P1105** | COL1A1 | `GVQGPPGPAGPR` | Q502 | 1 | 1105.5807 | 1106.5647 | 1106.5841 |
| **COL1A1_P1180** | COL1A1 | `GQAGVMGFPGPK` | Q391 | 1 | 1180.5908 | 1181.5748 | 1181.5942 |
| **COL1A1_P1427** | COL1A1 | `GSEGPQGVRGEPGPAGPR` | Q178 | 1 | 1427.6972 | 1428.6812 | 1428.7006 |
| **COL1A1_P1580** | COL1A1 | `GATGAPGIAGAPGFPGAR` | Q_or_G | 1 | 1580.7932 | 1581.7772 | 1581.7966 |
| **COL1A2_P1706** | COL1A2 | `GIPGEFGLPGPAGAR` | E/Q_marker | 1 | 1706.8837 | 1707.8677 | 1707.8871 |
| **COL1A1_P2043** | COL1A1 | `GAPGADGPAGAPGTPGPQGIAGQR` | Q774 | 2 | 2043.9806 | 2044.9646 | 2044.9840 |

---

## 2. Direct Area Ratio Calculations

With resolved high-resolution peaks, the deamidation fraction and PQI are computed directly from background-subtracted integrated peak areas:

### 2.1. Deamidation Fraction
`Deamidation Fraction = Area_deamidated / (Area_undeamidated + Area_deamidated)`

### 2.2. Parchment Glutamine Index (PQI)
`PQI Fraction = Area_undeamidated / (Area_undeamidated + Area_deamidated)`

*(Note: `PQI Fraction = 1.0 - Deamidation Fraction`)*

Sample-level summaries compute the median and mean PQI across all detected diagnostic markers.

---

## 3. Elemental Isotope Quality Control

To prevent chemical noise, isobaric adducts, or baseline anomalies from distorting deamidation measurements, `ZoomzMRT` computes the theoretical ¹³C₁ / *M*₀ ratio from the exact elemental stoichiometry (*C*, *H*, *N*, *O*, *S*) of the peptide:

`Expected 13C1 Ratio = Number of Carbons * 0.01078`

The observed ¹³C₁ ratio is calculated:
`Observed 13C1 Ratio = Area_13C1 / Area_undeamidated`

If the relative discrepancy exceeds 40%:
`|Observed - Expected| / Expected > 0.40`
the result is flagged with `C13_RATIO_ANOMALY`.

---

## 4. Resolution Status Assessment

For each marker, the ratio of observed separation to measured peak FWHM is determined:
`Separation / FWHM = 0.019339 Da / (m / R_measured)`

- **`resolved`**: Separation / FWHM ≥ 1.5 (peaks are fully separated).
- **`partially_resolved`**: 0.8 ≤ Separation / FWHM < 1.5 (peaks overlap slightly at base).
- **`unresolved`**: Separation / FWHM < 0.8 (peaks coalesce).
- **`not_detected`**: Peak signal is below the minimum intensity threshold.
