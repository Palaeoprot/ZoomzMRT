"""
Peptide Elemental Composition & Theoretical Isotope Calculator for High-Resolution MS.

Computes exact elemental stoichiometry (C, H, N, O, S) and theoretical M+1 / 13C1
isotopic abundances for collagen peptides and diagnostic ZooMS markers.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Mapping

# Standard monoisotopic elemental masses (Da)
ELEMENT_MASSES: Mapping[str, float] = {
    "C": 12.000000,
    "H": 1.007825032,
    "N": 14.003074004,
    "O": 15.994914620,
    "S": 31.972071174,
}

# Standard amino acid elemental compositions (residue, minus H2O)
AMINO_ACID_COMPOSITION: Mapping[str, Dict[str, int]] = {
    "A": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 0},   # Ala
    "R": {"C": 6, "H": 12, "N": 4, "O": 1, "S": 0},  # Arg
    "N": {"C": 4, "H": 6, "N": 2, "O": 2, "S": 0},   # Asn
    "D": {"C": 4, "H": 5, "N": 1, "O": 3, "S": 0},   # Asp
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},   # Cys
    "E": {"C": 5, "H": 7, "N": 1, "O": 3, "S": 0},   # Glu
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2, "S": 0},   # Gln
    "G": {"C": 2, "H": 3, "N": 1, "O": 1, "S": 0},   # Gly
    "H": {"C": 6, "H": 7, "N": 3, "O": 1, "S": 0},   # His
    "I": {"C": 6, "H": 11, "N": 1, "O": 1, "S": 0},  # Ile
    "L": {"C": 6, "H": 11, "N": 1, "O": 1, "S": 0},  # Leu
    "K": {"C": 6, "H": 12, "N": 2, "O": 1, "S": 0},  # Lys
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},   # Met
    "F": {"C": 9, "H": 9, "N": 1, "O": 1, "S": 0},   # Phe
    "P": {"C": 5, "H": 7, "N": 1, "O": 1, "S": 0},   # Pro
    "S": {"C": 3, "H": 5, "N": 1, "O": 2, "S": 0},   # Ser
    "T": {"C": 4, "H": 7, "N": 1, "O": 2, "S": 0},   # Thr
    "W": {"C": 11, "H": 10, "N": 2, "O": 1, "S": 0}, # Trp
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2, "S": 0},   # Tyr
    "V": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 0},   # Val
    "O": {"C": 5, "H": 7, "N": 1, "O": 2, "S": 0},   # Hyp (Hydroxyproline, P + O)
}

# Terrestrial natural isotopic relative abundances (ratio to primary isotope)
# 13C/12C ~ 0.01078, 15N/14N ~ 0.00368, 18O/16O ~ 0.00205, 2H/1H ~ 0.000115, 33S/32S ~ 0.0075
ABUNDANCE_13C = 0.01078
ABUNDANCE_15N = 0.00368
ABUNDANCE_18O = 0.00205
ABUNDANCE_17O = 0.00038
ABUNDANCE_2H = 0.000115
ABUNDANCE_33S = 0.00750


@dataclass(frozen=True)
class ElementalComposition:
    """Elemental formula counts for a molecule."""
    C: int
    H: int
    N: int
    O: int
    S: int = 0

    @property
    def monoisotopic_mass(self) -> float:
        """Calculate exact monoisotopic neutral mass."""
        return (
            self.C * ELEMENT_MASSES["C"]
            + self.H * ELEMENT_MASSES["H"]
            + self.N * ELEMENT_MASSES["N"]
            + self.O * ELEMENT_MASSES["O"]
            + self.S * ELEMENT_MASSES["S"]
        )

    def expected_13c1_ratio(self) -> float:
        """Expected relative intensity of the single 13C1 isotope peak to M0."""
        return self.C * ABUNDANCE_13C

    def expected_total_m1_ratio(self) -> float:
        """Expected relative intensity of the total M+1 isobaric envelope to M0."""
        return (
            self.C * ABUNDANCE_13C
            + self.N * ABUNDANCE_15N
            + self.H * ABUNDANCE_2H
            + self.O * ABUNDANCE_17O
            + self.S * ABUNDANCE_33S
        )


def peptide_composition(sequence: str, hyp_count: int = 0) -> ElementalComposition:
    """Calculate elemental composition of a peptide sequence (with terminal H and OH).

    Args:
        sequence: Standard 1-letter amino acid sequence.
        hyp_count: Number of Prolines that are hydroxylated to Hydroxyproline (+O per Hyp).

    Returns:
        ElementalComposition instance.
    """
    clean_seq = re.sub(r"[^A-Za-z]", "", sequence).upper()
    
    # Terminal H (N-terminus) + OH (C-terminus) = H2O
    c_total = 0
    h_total = 2
    n_total = 0
    o_total = 1
    s_total = 0

    for aa in clean_seq:
        if aa in AMINO_ACID_COMPOSITION:
            comp = AMINO_ACID_COMPOSITION[aa]
            c_total += comp["C"]
            h_total += comp["H"]
            n_total += comp["N"]
            o_total += comp["O"]
            s_total += comp["S"]
        else:
            raise ValueError(f"Unrecognized amino acid code in sequence: {aa}")

    # Add oxygen atoms for hydroxyprolines if sequence used 'P' rather than 'O'
    if hyp_count > 0:
        o_total += hyp_count

    return ElementalComposition(
        C=c_total,
        H=h_total,
        N=n_total,
        O=o_total,
        S=s_total,
    )


@dataclass(frozen=True)
class DeamidationMarker:
    """Immutable definition of a diagnostic collagen deamidation marker peptide."""
    name: str
    gene: str
    sequence: str
    site: str
    theoretical_mz: float
    hyp_count: int = 1
    charge: int = 1
    composition: ElementalComposition = None  # type: ignore

    def __post_init__(self):
        if self.composition is None:
            comp = peptide_composition(self.sequence, hyp_count=self.hyp_count)
            object.__setattr__(self, "composition", comp)

    @property
    def expected_13c_ratio(self) -> float:
        """Expected 13C1 / M0 ratio."""
        return self.composition.expected_13c1_ratio()

    @property
    def target_deam_mz(self) -> float:
        """Theoretical m/z for single deamidation (+0.984016 Da / z)."""
        return self.theoretical_mz + (0.984016 / self.charge)

    @property
    def target_13c_mz(self) -> float:
        """Theoretical m/z for single 13C isotope (+1.003355 Da / z)."""
        return self.theoretical_mz + (1.003355 / self.charge)


# Standard Collagen Deamidation Markers (COL1A1 & COL1A2)
# Verified references: van Doorn et al. 2012, Wilson et al. 2012, Welker et al. 2016, Nair & Bethencourt et al. 2023
DEFAULT_COLLAGEN_MARKERS: tuple[DeamidationMarker, ...] = (
    DeamidationMarker(
        name="COL1A1_P1105",
        gene="COL1A1",
        sequence="GVQGPPGPAGPR",
        site="Q502",
        theoretical_mz=1105.5807,
        hyp_count=1,
    ),
    DeamidationMarker(
        name="COL1A1_P1180",
        gene="COL1A1",
        sequence="GQAGVMGFPGPK",
        site="Q391",
        theoretical_mz=1180.5908,
        hyp_count=1,
    ),
    DeamidationMarker(
        name="COL1A1_P1427",
        gene="COL1A1",
        sequence="GSEGPQGVRGEPGPAGPR",
        site="Q178",
        theoretical_mz=1427.6972,
        hyp_count=1,
    ),
    DeamidationMarker(
        name="COL1A1_P1580",
        gene="COL1A1",
        sequence="GATGAPGIAGAPGFPGAR",
        site="Q_or_G",
        theoretical_mz=1580.7932,
        hyp_count=1,
    ),
    DeamidationMarker(
        name="COL1A2_P1706",
        gene="COL1A2",
        sequence="GIPGEFGLPGPAGAR",
        site="E/Q_marker",
        theoretical_mz=1706.8837,
        hyp_count=1,
    ),
    DeamidationMarker(
        name="COL1A1_P2043",
        gene="COL1A1",
        sequence="GAPGADGPAGAPGTPGPQGIAGQR",
        site="Q774",
        theoretical_mz=2043.9806,
        hyp_count=2,
    ),
)
