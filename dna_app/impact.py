from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .dna import find_orfs, longest_orf


Impact = Literal["silent", "missense", "nonsense", "frameshift", "unknown"]


@dataclass(frozen=True)
class ProteinImpact:
    impact: Impact
    ref_protein: str
    alt_protein: str
    ref_orf_len_nt: int | None
    alt_orf_len_nt: int | None
    notes: str


def analyze_protein_impact(ref_dna: str, alt_dna: str) -> ProteinImpact:
    """Best-effort impact analysis using the longest ORF on the + strand.

    This is an educational heuristic:
    - if indel length difference not multiple of 3 => frameshift (likely)
    - else compare translated proteins of longest ORF before/after
    """
    ref_orfs = find_orfs(ref_dna, include_reverse_complement=False)
    alt_orfs = find_orfs(alt_dna, include_reverse_complement=False)

    ref_best = longest_orf(ref_orfs)
    alt_best = longest_orf(alt_orfs)

    if ref_best is None or alt_best is None:
        return ProteinImpact(
            impact="unknown",
            ref_protein=ref_best.protein if ref_best else "",
            alt_protein=alt_best.protein if alt_best else "",
            ref_orf_len_nt=(len(ref_best.dna) if ref_best else None),
            alt_orf_len_nt=(len(alt_best.dna) if alt_best else None),
            notes="ORF non détectée dans au moins une séquence (analyse impact limitée).",
        )

    ref_len = len(ref_best.dna)
    alt_len = len(alt_best.dna)

    if (alt_len - ref_len) % 3 != 0:
        return ProteinImpact(
            impact="frameshift",
            ref_protein=ref_best.protein,
            alt_protein=alt_best.protein,
            ref_orf_len_nt=ref_len,
            alt_orf_len_nt=alt_len,
            notes="Différence de longueur ORF non multiple de 3 (décalage de cadre probable).",
        )

    if ref_best.protein == alt_best.protein:
        return ProteinImpact(
            impact="silent",
            ref_protein=ref_best.protein,
            alt_protein=alt_best.protein,
            ref_orf_len_nt=ref_len,
            alt_orf_len_nt=alt_len,
            notes="Protéines identiques (aucun changement au niveau AA).",
        )

    # Nonsense: alt protein ends earlier (often due to premature stop)
    if len(alt_best.protein) < len(ref_best.protein):
        return ProteinImpact(
            impact="nonsense",
            ref_protein=ref_best.protein,
            alt_protein=alt_best.protein,
            ref_orf_len_nt=ref_len,
            alt_orf_len_nt=alt_len,
            notes="Protéine mutée plus courte (stop prématuré possible).",
        )

    return ProteinImpact(
        impact="missense",
        ref_protein=ref_best.protein,
        alt_protein=alt_best.protein,
        ref_orf_len_nt=ref_len,
        alt_orf_len_nt=alt_len,
        notes="Protéines différentes (changement AA).",
    )
