"""Exact, dependency-free probes for coefficient-derived PSD realizations.

Polynomials use ascending coefficients.  All arithmetic is Fraction/integer.
This is exploratory code; NOTES.md states what is and is not proved.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb
def ordered_rows(j: int, nmax: int) -> list[list[int]]:
    if j < 2:
        raise ValueError("j >= 2 required")
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [0] * (n // j + 1)
        for r in range(1, min(n, j - 1) + 1):
            c = comb(n, r)
            for k, a in enumerate(rows[n - r]):
                row[k] += c * a
        for r in range(j, n + 1):
            c = comb(n, r)
            for k, a in enumerate(rows[n - r]):
                row[k + 1] += c * a
        rows.append(row)
    return rows


def power_sums(coeffs: list[int], count: int) -> list[Fraction]:
    """Newton power sums s_k=sum(lambda_i^k) for roots lambda_i of p.

    p(t)=prod(1+lambda_i t); s_0=degree.  Fraction arithmetic handles the
    normalized coefficients directly.
    """
    d = len(coeffs) - 1
    c = [Fraction(a, coeffs[0]) for a in coeffs]
    s = [Fraction(d)]
    for k in range(1, count + 1):
        # Newton: s_k + c_1 s_{k-1}+...+c_d s_{k-d} + k c_k=0,
        # where c_k=0 for k>d and the final term occurs only for k<=d.
        rhs = sum((c[k - i] * s[i] for i in range(1, k) if k - i <= d), Fraction(0))
        if k <= d:
            s.append(-(k * c[k] + rhs))
        else:
            s.append(-rhs)
    return s


def poly_shift(p: list[Fraction], a: Fraction) -> list[Fraction]:
    """(x-a)p(x), ascending coefficients."""
    out = [Fraction(0)] * (len(p) + 1)
    for i, c in enumerate(p):
        out[i + 1] += c
        out[i] -= a * c
    return out


def poly_sub(a: list[Fraction], b: list[Fraction], scale: Fraction = Fraction(1)) -> list[Fraction]:
    out = [Fraction(0)] * max(len(a), len(b))
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] -= scale * c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_inner(a: list[Fraction], b: list[Fraction], moments: list[Fraction]) -> Fraction:
    return sum(x * y * moments[i + k] for i, x in enumerate(a) for k, y in enumerate(b))


def jacobi_parameters(coeffs: list[int]) -> tuple[list[Fraction], list[Fraction]]:
    """Canonical symmetric Jacobi J with det(I+tJ)=p(t), if moments are PD.

    Diagonal alpha[0..q-1]; off-diagonal squares beta[1..q-1].  Entries are
    obtained only from coefficients through exact moment Gram--Schmidt.
    """
    q = len(coeffs) - 1
    if q == 0:
        return [], []
    s = power_sums(coeffs, 2 * q)
    moments = [(-1) ** k * s[k] for k in range(2 * q + 1)]
    v_prev: list[Fraction] = []
    v = [Fraction(1)]
    norm_prev = None
    norm = poly_inner(v, v, moments)
    alpha: list[Fraction] = []
    beta: list[Fraction] = [Fraction(0)] * q
    for n in range(q):
        if norm == 0:
            raise ZeroDivisionError(f"singular moment functional at Lanczos step {n}")
        xv = [Fraction(0)] + v
        a = poly_inner(xv, v, moments) / norm
        alpha.append(a)
        w = poly_sub(xv, v, a)
        if n:
            w = poly_sub(w, v_prev, beta[n])
        if n < q - 1:
            v_next = w
            norm_next = poly_inner(v_next, v_next, moments)
            beta[n + 1] = norm_next / norm
            v_prev, v = v, v_next
            norm_prev, norm = norm, norm_next
    return alpha, beta


def det_tridiagonal(alpha: list[Fraction], beta: list[Fraction]) -> list[int]:
    """Coefficients of det(I+tJ), ascending, for a Jacobi matrix."""
    if not alpha:
        return [1]
    prev2 = [1]
    prev1 = [1, alpha[0]]
    for i in range(1, len(alpha)):
        # D_i=(1+alpha_i t)D_{i-1}-beta_i t^2 D_{i-2}.
        out = [Fraction(0)] * (i + 2)
        for k, c in enumerate(prev1):
            out[k] += c
            out[k + 1] += alpha[i] * c
        for k, c in enumerate(prev2):
            out[k + 2] -= beta[i] * c
        prev2, prev1 = prev1, out
    return [Fraction(c) for c in prev1]


def det(matrix: list[list[Fraction]]) -> Fraction:
    """Exact Gaussian elimination determinant."""
    a = [row[:] for row in matrix]
    n = len(a)
    out = Fraction(1)
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            out = -out
        p = a[k][k]
        out *= p
        for i in range(k + 1, n):
            f = a[i][k] / p
            for j in range(k + 1, n):
                a[i][j] -= f * a[k][j]
    return out


def hankel_leading_minors(coeffs: list[int]) -> list[Fraction]:
    q = len(coeffs) - 1
    s = power_sums(coeffs, 2 * q)
    h = [[(-1) ** (i + k) * s[i + k] for k in range(q + 1)] for i in range(q + 1)]
    return [det([row[:r] for row in h[:r]]) for r in range(1, q + 2)]


def matrix_rank_difference(a: list[list[Fraction]], b: list[list[Fraction]]) -> int:
    if len(a) != len(b) or any(len(x) != len(y) for x, y in zip(a, b)):
        raise ValueError("same shapes required")
    c = [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]
    m = [row[:] for row in c]
    rows, cols = len(m), len(m[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if m[i][col]), None)
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        p = m[rank][col]
        for i in range(rank + 1, rows):
            f = m[i][col] / p
            for j in range(col, cols):
                m[i][j] -= f * m[rank][j]
        rank += 1
        if rank == rows:
            break
    return rank


def jacobi_matrix(alpha: list[Fraction], beta: list[Fraction]) -> list[list[Fraction]]:
    q = len(alpha)
    out = [[Fraction(0)] * q for _ in range(q)]
    for i, a in enumerate(alpha):
        out[i][i] = a
    for i in range(1, q):
        # beta is a square; sqrt need not be rational.  Keep squared entries
        # in reports rather than constructing an irrational matrix.
        out[i][i - 1] = out[i - 1][i] = beta[i]
    return out


def main() -> None:
    verified = 0
    for j in range(2, 9):
        rows = ordered_rows(j, 18)
        print(f"j={j}")
        for n, p in enumerate(rows):
            if len(p) <= 1:
                continue
            alpha, beta = jacobi_parameters(p)
            normalized = [Fraction(a, p[0]) for a in p]
            assert det_tridiagonal(alpha, beta) == normalized
            minors = hankel_leading_minors(p)
            assert minors[-1] == 0
            assert all(x >= 0 for x in minors[:-1])
            assert all(x > 0 for x in alpha) and all(x > 0 for x in beta[1:])
            verified += 1
            if n <= 4 * j:
                print(
                    f"  n={n:2d} p={p} alpha={[str(x) for x in alpha]} "
                    f"beta2={[str(x) for x in beta[1:]]} H={[str(x) for x in minors[:-1]]}"
                )
        same_degree_ranks = []
        for n in range(1, len(rows)):
            if len(rows[n]) == len(rows[n - 1]) and len(rows[n]) > 1:
                a1, b1 = jacobi_parameters(rows[n - 1])
                a2, b2 = jacobi_parameters(rows[n])
                same_degree_ranks.append((n, matrix_rank_difference(jacobi_matrix(a1, b1), jacobi_matrix(a2, b2))))
        print("  same-degree squared-J difference ranks:", same_degree_ranks)
    print("verified Jacobi characteristic identities:", verified)


if __name__ == "__main__":
    main()
