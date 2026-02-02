from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AlignmentResult:
    aligned_ref: str
    aligned_alt: str
    midline: str  # '|' match, '*' mismatch, ' ' gap
    score: int


def needleman_wunsch_global(
    ref: str,
    alt: str,
    *,
    match: int = 1,
    mismatch: int = -1,
    gap: int = -1,
) -> AlignmentResult:
    """Simple global alignment (Needleman–Wunsch).

    Returns aligned strings (with '-') and a midline visualization.
    """
    ref = ref.upper()
    alt = alt.upper()

    n = len(ref)
    m = len(alt)

    # DP score matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    # traceback: 0 diag, 1 up (gap in alt), 2 left (gap in ref)
    tb = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * gap
        tb[i][0] = 1
    for j in range(1, m + 1):
        dp[0][j] = j * gap
        tb[0][j] = 2

    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            aj = alt[j - 1]
            diag = dp[i - 1][j - 1] + (match if ri == aj else mismatch)
            up = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap

            best = diag
            move = 0
            if up > best:
                best = up
                move = 1
            if left > best:
                best = left
                move = 2

            dp[i][j] = best
            tb[i][j] = move

    # traceback
    i, j = n, m
    a_ref: list[str] = []
    a_alt: list[str] = []

    while i > 0 or j > 0:
        move = tb[i][j] if (i > 0 and j > 0) else (1 if i > 0 else 2)

        if move == 0:
            a_ref.append(ref[i - 1])
            a_alt.append(alt[j - 1])
            i -= 1
            j -= 1
        elif move == 1:
            a_ref.append(ref[i - 1])
            a_alt.append("-")
            i -= 1
        else:
            a_ref.append("-")
            a_alt.append(alt[j - 1])
            j -= 1

    a_ref.reverse()
    a_alt.reverse()

    aligned_ref = "".join(a_ref)
    aligned_alt = "".join(a_alt)

    mid: list[str] = []
    for r, a in zip(aligned_ref, aligned_alt, strict=True):
        if r == "-" or a == "-":
            mid.append(" ")
        elif r == a:
            mid.append("|")
        else:
            mid.append("*")

    return AlignmentResult(aligned_ref=aligned_ref, aligned_alt=aligned_alt, midline="".join(mid), score=dp[n][m])
