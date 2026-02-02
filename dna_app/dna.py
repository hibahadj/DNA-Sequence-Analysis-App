from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal


VALID_BASES = {"A", "T", "C", "G"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


DNA_CODON_TABLE: dict[str, str] = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",
    # Leucine
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Methionine
    "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "TAT": "Y", "TAC": "Y",
    # Histidine
    "CAT": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAT": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartic acid
    "GAT": "D", "GAC": "D",
    # Glutamic acid
    "GAA": "E", "GAG": "E",
    # Cysteine
    "TGT": "C", "TGC": "C",
    # Tryptophan
    "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stop
    "TAA": "*", "TAG": "*", "TGA": "*",
}


@dataclass(frozen=True)
class SequenceStats:
    length: int
    count_a: int
    count_t: int
    count_c: int
    count_g: int
    gc_percent: float
    at_percent: float


@dataclass(frozen=True)
class ORF:
    start: int  # 0-based inclusive
    end: int  # 0-based exclusive
    frame: int  # 0,1,2 relative to sequence start
    strand: Literal["+", "-"]
    dna: str
    protein: str


def normalize_sequence(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if not ch.isspace())


def validate_dna(seq: str) -> tuple[bool, str]:
    if not seq:
        return False, "Séquence vide."
    bad = sorted({ch for ch in seq if ch not in VALID_BASES})
    if bad:
        return False, f"Caractères invalides: {', '.join(bad)} (A, T, C, G uniquement)."
    return True, "OK"


def stats(seq: str) -> SequenceStats:
    seq = normalize_sequence(seq)
    length = len(seq)
    a = seq.count("A")
    t = seq.count("T")
    c = seq.count("C")
    g = seq.count("G")
    gc = ((g + c) / length * 100.0) if length else 0.0
    at = ((a + t) / length * 100.0) if length else 0.0
    return SequenceStats(length=length, count_a=a, count_t=t, count_c=c, count_g=g, gc_percent=gc, at_percent=at)


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    seq = normalize_sequence(seq)
    return "".join(comp[b] for b in reversed(seq))


def split_codons(seq: str, frame: int = 0) -> list[str]:
    seq = normalize_sequence(seq)
    frame = frame % 3
    codons: list[str] = []
    for i in range(frame, len(seq) - 2, 3):
        codons.append(seq[i : i + 3])
    return codons


def translate_dna(seq: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    seq = normalize_sequence(seq)
    aa: list[str] = []
    for codon in split_codons(seq, frame=frame):
        amino = DNA_CODON_TABLE.get(codon, "X")
        if amino == "*" and stop_at_stop:
            break
        aa.append(amino)
    return "".join(aa)


def find_orfs(seq: str, include_reverse_complement: bool = False) -> list[ORF]:
    """Detect ORFs in 3 frames (optionally 6 frames).

    ORF definition used here: starts at ATG and ends at the first in-frame stop codon.
    """
    seq = normalize_sequence(seq)
    results: list[ORF] = []

    def scan(s: str, strand: Literal["+", "-"]):
        for frame in range(3):
            i = frame
            while i <= len(s) - 3:
                codon = s[i : i + 3]
                if codon == START_CODON:
                    start = i
                    j = i + 3
                    while j <= len(s) - 3:
                        codon2 = s[j : j + 3]
                        if codon2 in STOP_CODONS:
                            end = j + 3
                            dna_orf = s[start:end]
                            prot = translate_dna(dna_orf, frame=0, stop_at_stop=True)
                            results.append(
                                ORF(
                                    start=start,
                                    end=end,
                                    frame=frame,
                                    strand=strand,
                                    dna=dna_orf,
                                    protein=prot,
                                )
                            )
                            i = end
                            break
                        j += 3
                    else:
                        i += 3
                else:
                    i += 3

    scan(seq, "+")
    if include_reverse_complement:
        scan(reverse_complement(seq), "-")

    results.sort(key=lambda o: (o.strand, o.start, o.end))
    return results


def longest_orf(orfs: Iterable[ORF]) -> ORF | None:
    best: ORF | None = None
    for o in orfs:
        if best is None or (o.end - o.start) > (best.end - best.start):
            best = o
    return best
