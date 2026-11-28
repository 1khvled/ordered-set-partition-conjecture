"""Cramer/convolution candidate M_n=A_n^{-1}B_n from first-block recurrence.

No root finding and no external dependencies.  All arithmetic is exact.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb
from explore import ordered_rows, det


def cramer_matrices(j: int, n: int) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    """Return A,B with O_{n,j}(t)=det(A+tB).

    Rows are r=1,...,n and columns m=0,...,n-1.  The entry is binom(r,r-m)
    for block size r-m<j, and t*binom(r,r-m) for block size >=j.
    """
    a = [[Fraction(0)] * n for _ in range(n)]
    b = [[Fraction(0)] * n for _ in range(n)]
    for i, r in enumerate(range(1, n + 1)):
        for k, m in enumerate(range(n)):
            d = r - m
            if d == 0:
                a[i][k] = Fraction(-1)
            elif 0 < d < j:
                a[i][k] = Fraction(comb(r, d))
            elif d >= j:
                b[i][k] = Fraction(comb(r, d))
    return a, b


def inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    n = len(matrix)
    aug = [row[:] + [Fraction(i == k) for k in range(n)] for i, row in enumerate(matrix)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if aug[i][k]), None)
        if pivot is None:
            raise ZeroDivisionError("singular matrix")
        aug[k], aug[pivot] = aug[pivot], aug[k]
        p = aug[k][k]
        aug[k] = [x / p for x in aug[k]]
        for i in range(n):
            if i != k and aug[i][k]:
                f = aug[i][k]
                aug[i] = [x - f * y for x, y in zip(aug[i], aug[k])]
    return [row[n:] for row in aug]


def matmul(a: list[list[Fraction]], b: list[list[Fraction]]) -> list[list[Fraction]]:
    n, p, q = len(a), len(b), len(b[0])
    return [[sum((a[i][k] * b[k][j] for k in range(p)), Fraction(0)) for j in range(q)] for i in range(n)]


def add(a: list[list[Fraction]], b: list[list[Fraction]]) -> list[list[Fraction]]:
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def det_poly_cramer(a: list[list[Fraction]], b: list[list[Fraction]], degree: int) -> list[Fraction]:
    """det(A+tB) by interpolation at 0,...,degree (exact)."""
    vals = []
    for t in range(degree + 1):
        vals.append(det(add(a, [[t * x for x in row] for row in b])))
    # Finite differences give coefficients in the binomial basis.
    diffs = vals[:]
    binom_coeff = []
    for k in range(degree + 1):
        binom_coeff.append(diffs[0])
        diffs = [diffs[i + 1] - diffs[i] for i in range(len(diffs) - 1)]
    out = [Fraction(0)] * (degree + 1)
    for k, c in enumerate(binom_coeff):
        # binom(t,k)=t(t-1)...(t-k+1)/k!, converted to monomials.
        falling = [Fraction(1)]
        for i in range(k):
            nxt = [Fraction(0)] * (len(falling) + 1)
            for r, x in enumerate(falling):
                nxt[r] -= i * x
                nxt[r + 1] += x
            falling = nxt
        for i, x in enumerate(falling):
            out[i] += c * x / _factorial(k)
    return out


def _factorial(k: int) -> int:
    out = 1
    for i in range(2, k + 1):
        out *= i
    return out


def matrix_string(m: list[list[Fraction]], max_n: int = 6) -> str:
    return ";".join(",".join(str(x) for x in row) for row in m[:max_n])


def main() -> None:
    for j in range(2, 8):
        rows = ordered_rows(j, 18)
        print(f"j={j}")
        symmetry_fail = None
        for n, p in enumerate(rows):
            if n == 0:
                continue
            a, b = cramer_matrices(j, n)
            da = det(a)
            assert da == p[0], (j, n, da, p[0])
            q = len(p) - 1
            assert det_poly_cramer(a, b, q) == p, (j, n)
            m = matmul(inverse(a), b)
            asym = [(i, k, m[i][k], m[k][i]) for i in range(n) for k in range(n) if m[i][k] != m[k][i]]
            if asym and symmetry_fail is None:
                symmetry_fail = (n, asym[0])
            if n <= max(4, 2 * j):
                print(f"  n={n:2d} detA={da} M={matrix_string(m)} asym={asym[:1]}")
        print("  first symmetry failure:", symmetry_fail)


if __name__ == "__main__":
    main()
