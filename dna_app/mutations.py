from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .alignment import AlignmentResult


MutationType = Literal["substitution", "insertion", "deletion"]
SubstitutionClass = Literal["transition", "transversion", "other"]


@dataclass(frozen=True)
class Mutation:
    kind: MutationType
    # 1-based positions in the original (ungapped) sequences
    ref_pos: int | None
    alt_pos: int | None
    ref_base: str
    alt_base: str
    substitution_class: SubstitutionClass | None = None


_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def classify_substitution(ref_base: str, alt_base: str) -> SubstitutionClass:
    rb = ref_base.upper()
    ab = alt_base.upper()
    if (rb, ab) in _TRANSITIONS:
        return "transition"
    if rb in {"A", "G", "C", "T"} and ab in {"A", "G", "C", "T"}:
        return "transversion"
    return "other"


def detect_mutations_from_alignment(aln: AlignmentResult) -> list[Mutation]:
    muts: list[Mutation] = []

    ref_i = 0
    alt_i = 0

    for r, a in zip(aln.aligned_ref, aln.aligned_alt, strict=True):
        if r != "-":
            ref_i += 1
        if a != "-":
            alt_i += 1

        if r == a:
            continue

        if r == "-":
            muts.append(
                Mutation(
                    kind="insertion",
                    ref_pos=ref_i,  # insertion after current ref index
                    alt_pos=alt_i,
                    ref_base="-",
                    alt_base=a,
                    substitution_class=None,
                )
            )
        elif a == "-":
            muts.append(
                Mutation(
                    kind="deletion",
                    ref_pos=ref_i,
                    alt_pos=alt_i,  # deletion corresponds to current alt index
                    ref_base=r,
                    alt_base="-",
                    substitution_class=None,
                )
            )
        else:
            muts.append(
                Mutation(
                    kind="substitution",
                    ref_pos=ref_i,
                    alt_pos=alt_i,
                    ref_base=r,
                    alt_base=a,
                    substitution_class=classify_substitution(r, a),
                )
            )

    return muts


def mutation_rate(mutations: list[Mutation], ref_length: int) -> float:
    if ref_length <= 0:
        return 0.0
    # Rate per reference base
    return len(mutations) / ref_length


def summarize_mutations(mutations: list[Mutation]) -> dict[str, int]:
    summary = {
        "substitution": 0,
        "transition": 0,
        "transversion": 0,
        "insertion": 0,
        "deletion": 0,
    }

    for m in mutations:
        summary[m.kind] += 1
        if m.kind == "substitution" and m.substitution_class in {"transition", "transversion"}:
            summary[m.substitution_class] += 1

    return summary
