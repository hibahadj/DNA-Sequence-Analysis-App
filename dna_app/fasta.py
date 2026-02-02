from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str


def _clean_sequence_lines(lines: Iterable[str]) -> str:
    seq_parts: list[str] = []
    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            continue
        seq_parts.append(line)
    return "".join(seq_parts).upper()


def parse_fasta_text(text: str) -> FastaRecord:
    """Parse FASTA-like text.

    Accepts:
    - True FASTA (with a header line starting with '>')
    - Plain sequence pasted directly

    Returns a single record (the first header is used if present).
    """
    lines = text.splitlines()
    header = ""
    for raw in lines:
        line = raw.strip()
        if line.startswith(">"):
            header = line[1:].strip()
            break

    sequence = _clean_sequence_lines(lines)
    return FastaRecord(header=header, sequence=sequence)


def read_fasta_file(path: str | Path) -> FastaRecord:
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="ignore")
    return parse_fasta_text(text)
