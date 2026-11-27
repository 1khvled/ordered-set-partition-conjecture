"""Exact checks for proof-weighted/NOTES.md.

All polynomial arithmetic is over integers or Fractions. No floating-point
root decision is used.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a or [0]


def add(a, b, scale=1):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += scale * x
    return trim(out)


def scale(a, c):
    return trim([c * x for x in a])


def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for k, y in enumerate(b):
            out[i + k] += x * y
    return trim(out)


def deriv(a):
    return trim([i * a[i] for i in range(1, len(a))])


def tmul(a):
    return mul([0, 1], a)


def qop(a):
    return add(a, mul([-1, 1], deriv(a)))


def eval_at(a, x):
    out = 0
    for c in reversed(a):
        out = out * x + c
    return out


def original_rows(j, nmax):
    rows = [[1]]
    for n in range(1, nmax + 1):
        row = [0] * (n // j + 1)
        for r in range(1, min(n, j - 1) + 1):
            for m, a in enumerate(rows[n - r]):
                row[m] += comb(n, r) * a
        for r in range(j, n + 1):
            for m, a in enumerate(rows[n - r]):
                row[m + 1] += comb(n, r) * a
        rows.append(trim(row))
    return rows


def short_recurrence_rhs(j, rows, n):
    rhs = add(scale(rows[n - 1], n), scale(tmul(deriv(rows[n - 1])), 2))
    for r in range(2, j):
        if r <= n:
            rhs = add(rhs, scale(rows[n - r], comb(n, r)))
    if n >= j:
        rhs = add(
            rhs,
            scale(tmul(qop(rows[n - j])), comb(n - 1, j - 1)),
        )
    return rhs


def check_original_recurrence():
    for j in range(2, 9):
        rows = original_rows(j, 35)
        for n in range(1, 36):
            assert short_recurrence_rhs(j, rows, n) == rows[n], (j, n)
    return "original recurrence checked for j=2..8, n=1..35"


def check_counterexample():
    # p=(t+2)(t+4), q=t+3, j=2, n=5, c=4.
    p = [8, 6, 1]
    q = [3, 1]
    mu = Fraction(-3)
    a5 = add(scale(p, 5), scale(tmul(deriv(p)), 2))
    assert a5 == [40, 42, 9]
    # A_5 has exact roots -10/3 and -4/3, so q is compatible with A_5.
    assert Fraction(-10, 3) < mu < Fraction(-4, 3)
    lhs_wrs = -eval_at(a5, mu) / eval_at(deriv(q), mu)
    rhs_wrs = 4 * mu * (mu - 1)
    assert lhs_wrs == 5
    assert rhs_wrs == 48
    assert lhs_wrs < rhs_wrs

    test_o5 = add(a5, scale(tmul(qop(q)), 4))
    assert test_o5 == [40, 50, 17]
    disc = test_o5[1] ** 2 - 4 * test_o5[2] * test_o5[0]
    assert disc == -220
    return {
        "A5": a5,
        "q_root": mu,
        "A5_compatible_with_q": True,
        "WRS_lhs": lhs_wrs,
        "WRS_rhs": rhs_wrs,
        "test_O5": test_o5,
        "discriminant": disc,
    }


def check_actual_margin():
    rows = original_rows(2, 5)
    p, q = rows[4], rows[3]
    mu = Fraction(-6, 7)
    a5 = add(scale(p, 5), scale(tmul(deriv(p)), 2))
    qp = eval_at(deriv(q), mu)
    last = 4 * mu * (mu - 1) * qp
    margin = -(eval_at(a5, mu) + last) / qp
    assert margin == Fraction(3222, 343)
    assert margin > 0
    assert eval_at(rows[5], mu) == -qp * margin
    return {"mu": mu, "WRS_margin": margin}


# Marker polynomials use keys (t_degree, u_degree, y1_degree, ..., yc_degree).
def madd(a, b, c=1):
    out = dict(a)
    for key, value in b.items():
        out[key] = out.get(key, 0) + c * value
        if out[key] == 0:
            del out[key]
    return out


def mscale(a, c):
    return {key: c * value for key, value in a.items() if c * value}


def mshift(a, index):
    out = {}
    for key, value in a.items():
        shifted = list(key)
        shifted[index] += 1
        out[tuple(shifted)] = value
    return out


def mK(a):
    return {key: key[1] * value for key, value in a.items() if key[1] * value}


def mdy(a, s):
    out = {}
    for key, value in a.items():
        exponent = key[2 + s]
        if exponent:
            shifted = list(key)
            shifted[2 + s] -= 1
            out[tuple(shifted)] = exponent * value
    return out


def check_marker_recurrence():
    for j in range(2, 7):
        c = j - 1
        nmax = 8
        rows = [{(0, 0) + (0,) * c: 1}]
        for n in range(1, nmax + 1):
            old = rows[n - 1]
            kp1 = madd(mK(old), old)
            # K*P + u*y1*(K+1)*P.
            new = madd(mK(old), mshift(mshift(kp1, 1), 2))
            for s in range(c):
                d = mdy(old, s)
                # -y_s*d + y_{s+1}*d, with y_c replaced by t.
                term = mscale(mshift(d, 2 + s), -1)
                term = madd(term, mshift(d, 0 if s + 1 == c else 3 + s))
                new = madd(new, term)
            rows.append(new)

            # Independent first-block recurrence, marking every block by u.
            direct = {}
            for size in range(1, n + 1):
                term = mscale(rows[n - size], comb(n, size))
                term = mshift(term, 1)
                term = mshift(term, 0 if size >= j else 2 + size - 1)
                direct = madd(direct, term)
            assert new == direct, (j, n)

        if j == 2:
            # R_{2,1}=y1*d/dy1 P_2 at y1=1 is 4u^2.
            r2 = {}
            for (mt, mu, y1), value in rows[2].items():
                r2[mu] = r2.get(mu, 0) + y1 * value
            assert {k: v for k, v in r2.items() if v} == {2: 4}
    return "marker recurrence checked for j=2..6, n=1..8; R_{2,1}=4u^2"


def check_refinement_obstructions():
    # P_4(x1,x2,x3,x4), specialized x2=x3=x4=1:
    # x^4+18x^2+8x+7. For x=-y it is y^4+18y^2-8y+7.
    quadratic_disc = (-8) ** 2 - 4 * 18 * 7
    assert quadratic_disc == -440
    # P_3(t,u)=tu+6tu^2+6u^3 at t=1/3 has quadratic factor 18u^2+6u+1.
    assert 6**2 - 4 * 18 * 1 == -36
    return {
        "size_weight_quartic": [7, 8, 18, 0, 1],
        "size_weight_quadratic_discriminant": quadratic_disc,
        "total_block_quadratic_discriminant": -36,
    }


def main():
    print(check_original_recurrence())
    print(check_counterexample())
    print(check_actual_margin())
    print(check_marker_recurrence())
    print(check_refinement_obstructions())


if __name__ == "__main__":
    main()
