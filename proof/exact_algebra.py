"""Exact integer-polynomial checks for the ordered big-block proof track.

Polynomials in t are lists of integer coefficients in ascending order.
No floating-point arithmetic is used.
"""
from __future__ import annotations

from math import comb, factorial
from pathlib import Path
import json
import sys

PYTHON = Path(sys.executable)


def trim(a: list[int]) -> list[int]:
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a or [0]


def add(a: list[int], b: list[int], scale: int = 1) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += scale * x
    return trim(out)


def scale(a: list[int], c: int) -> list[int]:
    return trim([c * x for x in a])


def mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for k, y in enumerate(b):
                out[i + k] += x * y
    return trim(out)


def deriv(a: list[int]) -> list[int]:
    return trim([i * a[i] for i in range(1, len(a))])


def t_deriv(a: list[int]) -> list[int]:
    return mul([0, 1], deriv(a))


def t_pow(a: list[int], p: int = 1) -> list[int]:
    return [0] * p + list(a)


def falling(x: int, q: int) -> int:
    out = 1
    for s in range(q):
        out *= x - s
    return out


def ordered_rows(j: int, nmax: int) -> list[list[int]]:
    """Coefficients of G=1/(1-B_hat), where B_hat=B_literal-t."""
    if j < 2:
        raise ValueError("j must be at least 2")
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [0] * (n // j + 1)
        for r in range(1, min(n, j - 1) + 1):
            c = comb(n, r)
            for m, a in enumerate(rows[n - r]):
                row[m] += c * a
        for r in range(j, n + 1):
            c = comb(n, r)
            for m, a in enumerate(rows[n - r]):
                row[m + 1] += c * a
        rows.append(trim(row))
    return rows


def apply_P(j: int, a: list[int]) -> list[int]:
    """The recovered operator P_j(D)=1+(t-1)D, independent of j."""
    return add(a, mul([-1, 1], deriv(a)))


def recurrence_residual(j: int, rows: list[list[int]], n: int) -> list[int]:
    """Residual of the guessed short recurrence, RHS-LHS."""
    if n < 1:
        return [0]
    rhs = add(scale(rows[n - 1], n), scale(t_deriv(rows[n - 1]), 2))
    for r in range(2, j):
        if r <= n:
            rhs = add(rhs, scale(rows[n - r], comb(n, r)))
    if n >= j:
        rhs = add(
            rhs,
            scale(
                t_pow(apply_P(j, rows[n - j])),
                falling(n - 1, j - 1) // factorial(j - 1),
            ),
        )
    return add(rhs, rows[n], -1)


def main() -> None:
    result = {"python": str(PYTHON), "candidate": {}, "examples": {}}
    for j in range(2, 9):
        nmax = 30
        rows = ordered_rows(j, nmax)
        bad = []
        for n in range(1, nmax + 1):
            res = recurrence_residual(j, rows, n)
            if res != [0]:
                bad.append({"n": n, "residual": res})
                break
        result["candidate"][str(j)] = {"checked_through": nmax, "first_bad": bad}
    for j in (2, 3, 4, 5):
        result["examples"][str(j)] = {str(n): rows for n, rows in enumerate(ordered_rows(j, 8))}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
