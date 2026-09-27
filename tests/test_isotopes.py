"""
Unit tests for peptide elemental composition and isotope calculations.
"""

import pytest
from zoomzmrt.isotopes import (
    DEFAULT_COLLAGEN_MARKERS,
    ElementalComposition,
    peptide_composition,
)


def test_peptide_composition_simple():
    # Gly-Ala: Gly (C2 H3 N1 O1) + Ala (C3 H5 N1 O1) + H2O = C5 H10 N2 O3
    comp = peptide_composition("GA")
    assert comp.C == 5
    assert comp.H == 10
    assert comp.N == 2
    assert comp.O == 3
    assert comp.S == 0


def test_hydroxyproline_composition():
    # Gly-Pro with 0 Hyp vs 1 Hyp
    comp_p = peptide_composition("GP", hyp_count=0)
    comp_hyp = peptide_composition("GP", hyp_count=1)
    
    assert comp_hyp.C == comp_p.C
    assert comp_hyp.H == comp_p.H
    assert comp_hyp.N == comp_p.N
    assert comp_hyp.O == comp_p.O + 1  # +1 oxygen for hydroxylation


def test_expected_13c1_ratio():
    comp = ElementalComposition(C=50, H=80, N=14, O=15, S=0)
    expected_c13 = comp.expected_13c1_ratio()
    # 50 * 0.01078 = 0.539
    assert abs(expected_c13 - 0.539) < 0.001


def test_collagen_markers_invariants():
    assert len(DEFAULT_COLLAGEN_MARKERS) >= 6
    for marker in DEFAULT_COLLAGEN_MARKERS:
        assert marker.theoretical_mz > 500.0
        assert marker.composition.C > 20
        assert marker.expected_13c_ratio > 0.2
        assert marker.target_deam_mz > marker.theoretical_mz
        assert marker.target_13c_mz > marker.target_deam_mz
        # Separation delta between deamidated M0 and undeamidated 13C1 must be ~19.339 mDa
        sep = marker.target_13c_mz - marker.target_deam_mz
        assert abs(sep - 0.019339) < 0.0001
