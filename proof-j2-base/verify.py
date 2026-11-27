"""Exact finite checks for the j=2 base-case note.

These checks are regression certificates and evidence only.  They are not an
all-index proof.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, factorial

import sympy as sp


x, z, t, u = sp.symbols("x z t u")


def rows(N: int, with_x: bool = False) -> list[list[int]]:
    """Ascending coefficients of O_{n,2}(t), or O_n(x,t) if with_x."""
    if not with_x:
        out = [[1]]
        for n in range(1, N + 1):
            row = [n * c for c in out[n - 1]] + [0]
            for r in range(2, n + 1):
                for k, c in enumerate(out[n - r]):
                    row[k + 1] += comb(n, r) * c
            out.append(row)
        return out

    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(comb(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def poly(a: list[int], var=t) -> sp.Expr:
    return sum(c * var**k for k, c in enumerate(a))


def conv(a: list[int], b: list[int]) -> list[int]:
    c = [0] * (len(a) + len(b) - 1)
    for i, p in enumerate(a):
        for j, q in enumerate(b):
            c[i + j] += p * q
    return c


def deriv(a: list[int]) -> list[int]:
    return [k * a[k] for k in range(1, len(a))]


def check_residue_identities() -> None:
    R = rows(9)
    U: list[list[int]] = []
    for n in range(10):
        v = [0] * (n + 1)
        for k in range(n + 1):
            for j, c in enumerate(conv(R[k], R[n - k])):
                v[j] += comb(n, k) * c
        U.append(v)

    for n in range(1, 9):
        p, dp = poly(R[n]), sp.diff(poly(R[n]), t)
        un, unm1, on = poly(U[n]), poly(U[n - 1]), poly(R[n])
        assert sp.expand(t * dp - (un - n * unm1 - on)) == 0
        assert sp.expand(
            poly(R[n + 1])
            - (2 * un + (t - 1) * n * unm1 - on)
        ) == 0
        assert sp.expand(
            poly(R[n + 1])
            - (on + 2 * t * dp + n * (1 + t) * unm1)
        ) == 0


def check_deformation_and_limit() -> None:
    O = rows(7, with_x=True)
    for n in range(1, 7):
        assert sp.expand(sp.diff(O[n], x) - n * sum(
            comb(n - 1, k) * O[k] * O[n - 1 - k]
            for k in range(n)
        )) == 0

    phi = [sp.Integer(1), sp.Integer(1)]
    for n in range(2, 21):
        phi.append(sp.expand(phi[n - 1] + u * phi[n - 2] / 2))
    for n in range(1, 20):
        formula = sum(
            factorial(n) * factorial(n - k)
            / (factorial(n - 2 * k) * 2**k * factorial(k)) * u**k
            for k in range(n // 2 + 1)
        )
        assert sp.expand(factorial(n) * phi[n] - formula) == 0

    v = sp.symbols("v")
    P = [sp.Integer(1), sp.Integer(1)]
    for n in range(2, 21):
        P.append(sp.expand(P[n - 1] - v * P[n - 2]))
    for n in range(2, 15):
        explicit = sum(
            (-1)**j * comb(n - j, j) * v**j
            for j in range(n // 2 + 1)
        )
        assert sp.expand(P[n] - explicit) == 0
        if n >= 4:
            roots = sp.nroots(P[n], n=40, maxsteps=300)
            assert all(abs(sp.im(r)) < sp.Float("1e-30") for r in roots)
            assert all(sp.re(r) > 0 for r in roots)


def associated_stirling(N: int, K: int) -> list[list[int]]:
    A = [[0] * (K + 2) for _ in range(N + 1)]
    A[0][0] = 1
    for m in range(1, N + 1):
        for k in range(1, min(K + 1, m // 2) + 1):
            if m >= 2:
                A[m][k] = k * A[m - 1][k] + (m - 1) * A[m - 2][k - 1]
    return A


def check_associated_lemmas() -> None:
    A = associated_stirling(300, 100)
    # Unweighted ratio monotonicity (the proved lemma).
    for m in range(2, 300):
        for k in range(1, 100):
            if A[m][k + 1] == 0 or A[m + 1][k + 1] == 0:
                continue
            assert A[m][k] * A[m + 1][k + 1] >= A[m][k + 1] * A[m + 1][k]
    # Stronger inequality (H), checked but not proved in the note.
    for m in range(2, 300):
        for k in range(1, 100):
            if A[m][k + 1] == 0 or A[m + 1][k] == 0:
                continue
            assert (
                k * A[m + 1][k + 1] * A[m][k]
                >= (k + 1) * A[m + 1][k] * A[m][k + 1]
            )


def check_resultants_and_discriminants() -> None:
    O = rows(6, with_x=True)
    for n in range(1, 6):
        res = sp.factor(sp.resultant(O[n], O[n + 1], t))
        assert res != 0
        # The observed residual cofactor has positive coefficients.
        q = sp.Poly(res / x**n, x)
        sg = sp.sign(q.LC())
        assert all(sg * c > 0 for c in q.all_coeffs())

    expected_disc4 = (2 * x + 1) * (6 * x + 1) * (60 * x**2 + 8 * x + 1)
    assert sp.factor(sp.discriminant(O[4], t) - expected_disc4) == 0
    for n in range(2, 7):
        disc = sp.Poly(sp.discriminant(O[n], t), x)
        assert all(c > 0 for c in disc.all_coeffs())


def check_obstructions() -> None:
    # Falling-basis ratio monotonicity is not sufficient.
    P = 10 + 11 * t + t**2
    Q = 10 + 22 * t + 100 * t**2
    assert sp.discriminant(Q, t) == -3516
    # P = 10 + 12(t)_1 + (t)_2, Q = 10 + 122(t)_1 + 100(t)_2.
    assert all(
        sp.Rational(122, 12) > 1 and sp.Rational(100, 1) > sp.Rational(122, 12)
        for _ in [0]
    )

    # Associated-Bell interlacing alone does not imply (H).
    P = (t + 1) * (t + 2) * (t + 3)
    Q = (t + sp.Rational(3, 2)) * (t + sp.Rational(5, 2)) * (t + 4)
    p = sp.Poly(P, t).all_coeffs()[::-1]
    q = sp.Poly(Q, t).all_coeffs()[::-1]
    assert q[2] * p[1] - 2 * q[1] * p[2] == -149

    # G_4 and O_4 do not interlace.
    assert sp.solve(24 + 45 * t + 6 * t**2, t) == [
        (-45 - sp.sqrt(1449)) / 12,
        (-45 + sp.sqrt(1449)) / 12,
    ]
    assert (-45 + sp.sqrt(1449)) / 12 < -sp.Rational(1, 2)

    # Fixed disjoint quadratic grouping for W_8(-s) fails.
    R = rows(8)
    W = [
        a - b
        for a, b in zip(
            conv(deriv(R[8]), R[6]),
            conv(R[8], deriv(R[6])),
        )
    ]
    while W and W[-1] == 0:
        W.pop()
    assert W[3 : 6] == [323112180, 106099140, 7509600]
    assert W[4] ** 2 > 4 * W[3] * W[5]


def main() -> None:
    check_residue_identities()
    print("residue identities: OK")
    check_deformation_and_limit()
    print("deformation and Fibonacci limit: OK")
    check_associated_lemmas()
    print("associated-Stirling finite checks: OK")
    check_resultants_and_discriminants()
    print("resultant/discriminant finite checks: OK")
    check_obstructions()
    print("exact obstructions: OK")


if __name__ == "__main__":
    main()
