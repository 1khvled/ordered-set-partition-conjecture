"""Consolidated exact verification for proof-psd-realization.

Run from the repository root:
    py -3.14 proof-psd-realization/verify.py
Only the Python standard library is used.  No root finding is performed.
"""
from __future__ import annotations

from fractions import Fraction
import json

# The imports below are kept explicit to make failures easy to locate.
from explore import (
    det_tridiagonal,
    hankel_leading_minors,
    jacobi_parameters,
    matrix_rank_difference,
    ordered_rows,
)
from cramer_probe import (
    det,
    det_poly_cramer,
    inverse,
    matmul,
    cramer_matrices,
)
from kill_test import (
    lanczos_from_moments,
    moments_of_ratio,
    primitive,
    ratio_polynomials,
)
from even_border_probe import even_border_data
from jacobi_chain_probe import fit_step


def zero(n: int) -> list[list[Fraction]]:
    return [[Fraction(0)] * n for _ in range(n)]


def verify_cramer(jmax: int = 8, nmax: int = 20) -> dict:
    pairs = updates = 0
    for j in range(2, jmax + 1):
        rows = ordered_rows(j, nmax)
        matrices = {}
        for n in range(1, nmax + 1):
            a, b = cramer_matrices(j, n)
            assert det(a) == rows[n][0]
            assert det_poly_cramer(a, b, len(rows[n]) - 1) == rows[n]
            matrices[n] = (a, b, matmul(inverse(a), b))
            pairs += 1
        for n in range(1, nmax):
            a, b, m = matrices[n]
            an, bn, mn = matrices[n + 1]
            ext = [row[:] for row in mn]
            for i in range(n):
                for k in range(n):
                    ext[i][k] -= m[i][k]
            assert matrix_rank_difference(ext, zero(n + 1)) <= 1
            # D=A_{n+1}^{-1} e_last (e^T-b^T M_n,0).
            f = [Fraction(0)] * n + [Fraction(1)]
            w = [bn[n][k] - sum(an[n][i] * m[i][k] for i in range(n)) for k in range(n)] + [Fraction(0)]
            inv_an = inverse(an)
            u = [inv_an[i][n] for i in range(n + 1)]
            assert matmul(an, ext) == [[f[i] * w[k] for k in range(n + 1)] for i in range(n + 1)]
            # Explicit positive vector u_r=O_r(0)/O_{n+1}(0).
            explicit_u = [Fraction(rows[r][0], rows[n + 1][0]) for r in range(n + 1)]
            assert u == explicit_u
            updates += 1
    return {"cramer_identities": pairs, "rank_one_nonsymmetric_updates": updates}


def verify_hankel_jacobi(jmax: int = 8, nmax: int = 18) -> dict:
    count = 0
    for j in range(2, jmax + 1):
        for n, p in enumerate(ordered_rows(j, nmax)):
            if len(p) <= 1:
                continue
            alpha, beta = jacobi_parameters(p)
            assert det_tridiagonal(alpha, beta) == [Fraction(x, p[0]) for x in p]
            minors = hankel_leading_minors(p)
            assert minors[-1] == 0 and all(x >= 0 for x in minors[:-1])
            assert all(x > 0 for x in alpha) and all(x > 0 for x in beta[1:])
            count += 1
    return {"hankel_jacobi_rows": count}


def verify_kill_tests() -> dict:
    rows = ordered_rows(2, 11)
    odd_digits = []
    even_digits = []
    for m in range(1, 6):
        even, odd = rows[2 * m], rows[2 * m + 1]
        num, den = ratio_polynomials(even, odd)
        num, den = primitive(num, den)
        size = len(den) - 1
        moments = moments_of_ratio(num, den, 2 * size + 2)
        alpha, beta, last_norm = lanczos_from_moments(moments, size)
        assert last_norm == 0
        assert det_tridiagonal(alpha, beta) == [Fraction(x, even[0]) for x in even]
        assert all(x > 0 for x in alpha) and all(x > 0 for x in beta[1:size])
        odd_digits.append(max(len(str(x.denominator)) for x in alpha))

        border_alpha, _, _, a, b, weight = even_border_data(rows[2 * m - 1], even)
        assert border_alpha > 0 and weight >= 0
        if a:
            assert all(x > 0 for x in a) and all(x > 0 for x in b[1:len(a)])
            even_digits.append(max(len(str(x.denominator)) for x in a))
        else:
            even_digits.append(0)
    assert odd_digits == [1, 3, 11, 35, 73]
    assert even_digits == [0, 1, 7, 24, 56]
    return {"odd_alpha_denominator_digits": odd_digits, "even_alpha_denominator_digits": even_digits}


def verify_obstructions() -> dict:
    rows = ordered_rows(2, 8)
    a, b = cramer_matrices(2, 4)
    m = matmul(inverse(a), b)
    # det(I+tM)=1+15/8 t+1/4 t^2, so zero has algebraic multiplicity 2.
    assert det_poly_cramer(a, b, 2) == rows[4]
    assert matrix_rank_difference(m, zero(4)) == 3  # geometric multiplicity 1
    assert all(m[i][3] == 0 for i in range(4))
    # Rational symmetric 2x2 would make 161/64 a sum of two rational squares.
    assert Fraction(15, 8) ** 2 - 4 * Fraction(1, 4) == Fraction(161, 64)
    assert 161 == 7 * 23

    # No three-term Favard chain even for the first degree-grouped j=2 rows.
    chain = [rows[0], rows[2], rows[4], rows[6]]
    _, _, residual = fit_step(chain[1], chain[2], chain[3])
    assert residual == [Fraction(0), Fraction(0), Fraction(0), Fraction(-2323, 11520)]
    return {"jordan_obstruction": "j=2,n=4", "rational_obstruction": "161/64", "favard_residual": str(residual[-1])}


def main() -> None:
    report = {
        "scope": "proof-psd-realization only",
        "arithmetic": "exact Fraction/integer; no roots",
        **verify_cramer(),
        **verify_hankel_jacobi(),
        **verify_kill_tests(),
        **verify_obstructions(),
        "status": "all exact checks passed; no all-index proof",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
