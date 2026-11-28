"""Exact audits for proof-falling-to-realroots/NOTES.md.

Standard library only.  All arithmetic is integral or Fraction arithmetic.
Finite checks audit symbolic identities and exact witnesses; they are not an
all-index enumeration proof.
"""

from fractions import Fraction
from itertools import combinations
from math import comb, factorial


# ---------------------------------------------------------------------------
# Polynomials in ascending coefficient order.


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def scale(a, c):
    return trim([c * x for x in a])


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def evaluate(a, x):
    out = 0
    for c in reversed(a):
        out = out * x + c
    return out


def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))])


def discriminant2(a):
    assert len(a) == 3 and a[2] != 0
    return a[1] * a[1] - 4 * a[0] * a[2]


def falling_basis(degree):
    out = [(1,)]
    for i in range(degree):
        out.append(mul(out[-1], (-i, 1)))
    return out


def stirling2_table(n):
    s = [[0] * (n + 1) for _ in range(n + 1)]
    s[0][0] = 1
    for k in range(1, n + 1):
        for i in range(1, k + 1):
            s[k][i] = s[k - 1][i - 1] + i * s[k - 1][i]
    return s


def to_falling(monomial):
    d = len(monomial) - 1
    s = stirling2_table(d)
    return tuple(sum(s[k][i] * monomial[k] for k in range(i, d + 1))
                 for i in range(d + 1))


def to_monomial(falling):
    basis = falling_basis(len(falling) - 1)
    out = (0,)
    for c, f in zip(falling, basis):
        out = add(out, scale(f, c))
    return out


# ---------------------------------------------------------------------------
# Bezout kernel (q(x)p(y)-q(y)p(x))/(x-y), by exact synthetic division.


def bezout_matrix(p, q):
    """Return the full coefficient matrix in 1,x,...,x^(dx-1)."""
    numerator = {}
    for i, qi in enumerate(q):
        for j, pj in enumerate(p):
            numerator[(i, j)] = numerator.get((i, j), 0) + qi * pj
            numerator[(j, i)] = numerator.get((j, i), 0) - qi * pj
    numerator = {ij: c for ij, c in numerator.items() if c}
    if not numerator:
        return tuple()

    dx = max(i for i, _ in numerator)
    if dx == 0:
        return tuple()

    by_x = [dict() for _ in range(dx + 1)]
    for (i, j), c in numerator.items():
        by_x[i][j] = by_x[i].get(j, 0) + c

    # If N(x,y)=(x-y)Q(x,y), then
    # Q_{dx-1}=N_dx and Q_{e-1}=N_e + y Q_e.
    quotient = [dict() for _ in range(dx)]
    quotient[dx - 1] = by_x[dx]
    for e in range(dx - 2, -1, -1):
        row = dict(by_x[e + 1])
        for ypow, c in quotient[e + 1].items():
            row[ypow + 1] = row.get(ypow + 1, 0) + c
        quotient[e] = {j: c for j, c in row.items() if c}

    remainder = dict(by_x[0])
    for ypow, c in quotient[0].items():
        remainder[ypow + 1] = remainder.get(ypow + 1, 0) + c
    remainder = {j: c for j, c in remainder.items() if c}
    assert not remainder, ("nonzero Bezout division remainder", remainder)

    return tuple(tuple(quotient[a].get(b, 0) for b in range(dx))
                 for a in range(dx))


def embed(matrix, size):
    out = [[0] * size for _ in range(size)]
    for i, row in enumerate(matrix):
        for j, c in enumerate(row):
            out[i][j] = c
    return out


def matrix_add_weighted(out, matrix, weight):
    for i, row in enumerate(matrix):
        for j, c in enumerate(row):
            out[i][j] += weight * c


def det2(m):
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def minors2(row_a, row_b):
    return tuple(row_a[i] * row_b[j] - row_a[j] * row_b[i]
                 for i, j in combinations(range(len(row_a)), 2))


def falling_bezout_decomposition(p, q):
    """Check M(p,q)=sum Delta_ij K_ij in the falling basis."""
    cp = to_falling(p)
    cq = to_falling(q)
    direct = bezout_matrix(p, q)
    size = len(direct)
    out = [[0] * size for _ in range(size)]
    degree = max(len(cp), len(cq)) - 1
    basis = falling_basis(degree)
    cpi = cp + (0,) * (degree + 1 - len(cp))
    cqi = cq + (0,) * (degree + 1 - len(cq))
    for i in range(degree):
        for j in range(i + 1, degree + 1):
            delta = cpi[i] * cqi[j] - cpi[j] * cqi[i]
            kij = bezout_matrix(basis[i], basis[j])
            matrix_add_weighted(out, embed(kij, size), delta)
    reconstructed = tuple(tuple(row) for row in out)
    assert reconstructed == direct
    return cp, cq, direct, reconstructed


# ---------------------------------------------------------------------------
# Exact O_{n,2} rows from the first-block recurrence.


def o_rows_j2(nmax):
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [0] * (n // 2 + 1)
        for r in range(1, n + 1):
            shift = 0 if r == 1 else 1
            for k, a in enumerate(rows[n - r]):
                if k + shift < len(row):
                    row[k + shift] += comb(n, r) * a
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Audits.


def main():
    print("proof-falling-to-realroots exact audit")

    # Basis conversion and Bezout division residuals.
    basis = falling_basis(4)
    for i, f in enumerate(basis):
        unit_falling = (0,) * i + (1,)
        assert to_monomial(unit_falling) == f
        assert to_falling(f) == unit_falling
    print("basis conversion residuals: 0")

    # Strict obstruction from NOTES section 5.
    p = (6, 5, 1)
    q = (4, 9, 2)
    cp, cq, mp, reconstructed = falling_bezout_decomposition(p, q)
    assert cp == (6, 6, 1)
    assert cq == (4, 11, 2)
    assert mp == ((34, 8), (8, 1))
    assert reconstructed == mp
    assert minors2(p, q) == (34, 8, 1)
    assert minors2(cp, cq) == (42, 8, 1)
    assert discriminant2(p) == 1
    assert discriminant2(q) == 49
    assert discriminant2(cp) == 12
    assert discriminant2(cq) == 89
    assert evaluate(p, -3) == 0
    assert Fraction(evaluate(q, -3), evaluate(derivative(p), -3)) == 5
    assert det2(mp) == -30

    delta01, delta02, delta12 = 42, 8, 1
    residual = (delta01 * delta12
                - delta02 * (delta12 + delta02))
    assert residual == det2(mp) == -30
    print("strict witness monomial minors:", minors2(p, q))
    print("strict witness falling minors:", minors2(cp, cq))
    print("strict witness row discriminants:",
          discriminant2(p), discriminant2(q),
          discriminant2(cp), discriminant2(cq))
    print("strict witness Bezout:", mp, "det residual:", residual)
    print("strict witness failed residue Q(-3)/P'(-3):",
          Fraction(evaluate(q, -3), evaluate(derivative(p), -3)))

    # The indefinite lag-two basis kernel.
    k02 = bezout_matrix(basis[0], basis[2])
    assert k02 == ((-1, 1), (1, 0))
    assert det2(k02) == -1
    print("K_02:", k02, "determinant:", det2(k02))

    # General decomposition residuals on several exact pairs.
    test_pairs = [
        ((2, 1), (1, 2, 1)),
        ((6, 5, 1), (4, 9, 2)),
        ((1, 3, 3, 1), (2, 7, 6, 1)),
        ((24, 45, 6), (120, 311, 110)),
    ]
    for pp, qq in test_pairs:
        _, _, direct, recon = falling_bezout_decomposition(pp, qq)
        assert direct == recon
    print("full falling-basis Bezout decomposition residuals: 0")

    # Actual O_4,O_5 base certificate.
    rows = o_rows_j2(5)
    assert rows[4] == [24, 45, 6]
    assert rows[5] == [120, 311, 110]
    c4 = to_falling(tuple(rows[4]))
    c5 = to_falling(tuple(rows[5]))
    assert c4 == (24, 51, 6)
    assert c5 == (120, 421, 110)
    _, _, m45, _ = falling_bezout_decomposition(tuple(rows[4]), tuple(rows[5]))
    assert m45 == ((2064, 1920), (1920, 3084))
    assert det2(m45) == 2678976
    d01 = c4[0] * c5[1] - c4[1] * c5[0]
    d02 = c4[0] * c5[2] - c4[2] * c5[0]
    d12 = c4[1] * c5[2] - c4[2] * c5[1]
    assert (d01, d02, d12) == (3984, 1920, 3084)
    assert d01 * d12 - d02 * (d12 + d02) == 2678976
    print("O_4 falling:", c4, "O_5 falling:", c5)
    print("O_4/O_5 Delta:", (d01, d02, d12))
    print("O_4/O_5 Bezout:", m45, "determinant:", det2(m45))

    # Finite global TP2 embedding with the floor(n/2) degree pattern.
    embedding = [
        (1, 0, 0, 0),
        (1, 0, 0, 0),
        (4, 1, 0, 0),
        (4, 1, 0, 0),
        (6, 6, 1, 0),
        (4, 11, 2, 0),
        (1, 4, 4, 4),
        (1, 4, 4, 4),
    ]
    for n, row in enumerate(embedding):
        expected_degree = n // 2
        actual_degree = max(i for i, c in enumerate(row) if c)
        assert actual_degree == expected_degree
    for i, j in combinations(range(len(embedding)), 2):
        assert all(x >= 0 for x in minors2(embedding[i], embedding[j]))
    assert embedding[4][:3] == cp
    assert embedding[5][:3] == cq
    print("global TP2 embedding: all 2x2 minors nonnegative")

    # Fisk/Brenti direction counterexample.
    transformed = to_monomial((6, 5, 1))
    assert transformed == (6, 4, 1)
    assert discriminant2(transformed) == -8
    print("falling-transform converse witness:", transformed,
          "discriminant:", discriminant2(transformed))

    # Polya-Schur test for a generic strictly increasing ratio sequence.
    lam = (Fraction(1), Fraction(2), Fraction(3))
    egf = tuple(lam[i] / lam[0] / factorial(i) for i in range(3))
    assert egf == (Fraction(1), Fraction(2), Fraction(3, 2))
    assert egf[1] * egf[1] - 4 * egf[0] * egf[2] == -2
    print("increasing-ratio Polya-Schur test discriminant:",
          egf[1] * egf[1] - 4 * egf[0] * egf[2])

    # Basis mismatch: monomial multiplication by the falling ratios is not Q.
    falling_ratios = tuple(Fraction(cq[i], cp[i]) for i in range(3))
    monomial_multiplier = tuple(falling_ratios[i] * p[i] for i in range(3))
    assert falling_ratios == (Fraction(2, 3), Fraction(11, 6), Fraction(2))
    assert monomial_multiplier == (Fraction(4), Fraction(55, 6), Fraction(2))
    assert monomial_multiplier != q
    print("falling ratios applied in monomial basis:", monomial_multiplier,
          "!= Q:", q)

    print("all exact audits passed")


if __name__ == "__main__":
    main()
