#!/usr/bin/env python3
"""Exact checks accompanying NOTES.md.

These are finite certificates/evidence, not an all-index proof.  Root signs
are certified by rational interval arithmetic; coefficient checks use Python
integers only.
"""
from __future__ import annotations

from math import comb

from audit_corrected import (
    audit_range,
    derivative,
    exact_eval,
    ordered_rows,
    retraction_certificate,
)
from fractions import Fraction


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for k, y in enumerate(b):
            out[i + k] += x * y
    return out


def add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(n)]


def scale(a, c):
    return [c * x for x in a]


def times_t(a):
    return [0] + a


def wronskian(p, q):
    """p' q - p q' (newer/older orientation)."""
    return [x - y for x, y in zip(mul(derivative(p), q), mul(p, derivative(q)))]


def tQ(q):
    d = derivative(q)
    # tQ(q) = t q + t(t-1)q'.
    return add(times_t(q), times_t(add(scale(d, -1), times_t(d))))


def inner_row(rows, j, n):
    """H_n = t sum_{r=j}^n binom(n,r) O_{n-r}."""
    out = [0]
    for r in range(j, n + 1):
        out = add(out, scale(times_t(rows[n-r]), comb(n, r)))
    return out


def coefficient_checks(j, nmax):
    rows = ordered_rows(j, nmax)
    pair_checks = 0
    cross_checks = 0
    for newer in range(1, nmax + 1):
        for older in range(newer):
            w = wronskian(rows[newer], rows[older])
            assert all(c >= 0 for c in w), (j, newer, older, w)
            pair_checks += 1
    for n in range(2 * j, nmax + 1):
        h = inner_row(rows, j, n)
        q = rows[n-j]
        w = wronskian(h, q)
        assert all(c >= 0 for c in w), (j, n, w)
        cross_checks += 1
    return {"j": j, "nmax": nmax, "wronskian_pairs": pair_checks,
            "cross_H_O_pairs": cross_checks}


def algebra_checks():
    # The product term in W(tp,q) is essential.
    q = [24, 45, 6]  # O_{4,2}
    p = [1]          # O_{0,2}; tp=t
    lhs = wronskian(times_t(p), q)
    rhs = add(mul(p, q), times_t(wronskian(p, q)))
    assert lhs == rhs == [24, 0, -6]

    # Correct tQ coefficient arithmetic (the retracted script used x*q').
    q_linear = [6, 7]
    mu_linear = Fraction(-6, 7)
    assert tQ(q_linear) == [0, -1, 14]
    assert exact_eval(tQ(q_linear), mu_linear) == (
        mu_linear * (mu_linear - 1) * exact_eval(derivative(q_linear), mu_linear)
    )

    rows = ordered_rows(2, 6)
    h6 = inner_row(rows, 2, 6)
    a = [h6[k+1] if k+1 < len(h6) else 0 for k in range(len(rows[4]))]
    expected_a = [sum(comb(6, m) * (rows[m][k] if k < len(rows[m]) else 0)
                      for m in range(5)) for k in range(3)]
    assert a == expected_a == [517, 830, 90]
    assert wronskian(h6, rows[4]) == [12408, 39840, 40728, 8100, 540]

    # The natural all-lag strengthening G_{n,d} with weight C(n-1,d-1)
    # already fails at an actual consecutive pair.  This is NOT a WRS failure:
    # here d=1, while the target for j=2 has d=j=2.
    mu = Fraction(-6, 7)  # root of O_{3,2}=6+7t
    pval = exact_eval([24, 45, 6], mu)
    qp = exact_eval(derivative([6, 7]), mu)
    residue = -pval / qp
    weighted_rhs = mu * (mu - 1)
    assert residue == Fraction(498, 343)
    assert weighted_rhs == Fraction(78, 49)
    assert residue < weighted_rhs

    # Exact obstruction to the simple product-of-TN outer transform for j>=3.
    # K_{n,r}=[z^(n-r)]D_3; rows {1,2}, cols {0,1} give -1/2.
    assert Fraction(1 * 1 - 1 * Fraction(3, 2)) == Fraction(-1, 2)

    # Natural inner production matrix for j=3.
    g3, g4, g5 = Fraction(1, 6), Fraction(5, 24), Fraction(3, 10)
    assert g4 * g4 - g3 * g5 == Fraction(-19, 2880)


def main():
    algebra_checks()
    retraction_certificate()
    for j, nmax in [(2, 40), (3, 36), (4, 32), (5, 30), (6, 30)]:
        print("WRS range", audit_range(j, nmax))
    for j, nmax in [(2, 60), (3, 45), (4, 36), (5, 30), (6, 30)]:
        print("coefficient invariants", coefficient_checks(j, nmax))
    print("all exact checks passed")


if __name__ == "__main__":
    main()
