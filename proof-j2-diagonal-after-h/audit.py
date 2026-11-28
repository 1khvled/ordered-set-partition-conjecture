"""Exact, standard-library-only audit for proof-j2-diagonal-after-h.

No numerical root finding is used.  Finite actual-family checks are explicitly
finite evidence; the rational countermodels are logical falsifications of the
stated H + factorial-row + component-interlacing implications.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import comb, factorial

ZERO = F(0)
ONE = F(1)


def trim(p):
    p = [F(x) for x in p]
    while p and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    n = max(len(p), len(q))
    return trim([
        (p[i] if i < len(p) else ZERO) + (q[i] if i < len(q) else ZERO)
        for i in range(n)
    ])


def sub(p, q):
    return add(p, [-x for x in q])


def scale(c, p):
    return trim([F(c) * x for x in p])


def mul(p, q):
    if not p or not q:
        return []
    r = [ZERO] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            r[i + j] += a * b
    return trim(r)


def der(p):
    return trim([i * p[i] for i in range(1, len(p))])


def U(p):
    """U p = 2 p + t p'."""
    tp = [i * p[i] for i in range(len(p))]
    return add(scale(2, p), tp)


def W(p, q):
    """W(p,q) = p'q - p q'."""
    return sub(mul(der(p), q), mul(p, der(q)))


def sturm_sequence(p):
    p = trim(p)
    if len(p) <= 1:
        return [p]
    seq = [p, der(p)]
    while seq[-1]:
        a, b = seq[-2], seq[-1]
        r = a[:]
        while len(r) >= len(b):
            c = r[-1] / b[-1]
            shift = len(r) - len(b)
            for i in range(len(b)):
                r[i + shift] -= c * b[i]
            r = trim(r)
        if not r:
            break
        seq.append(scale(-1, r))
    return seq


def variations(signs):
    signs = [s for s in signs if s]
    return sum(signs[i] != signs[i - 1] for i in range(1, len(signs)))


def real_root_count(p):
    """Number of distinct real roots, by exact Sturm counts at +/- infinity."""
    p = trim(p)
    if len(p) <= 1:
        return 0
    seq = sturm_sequence(p)
    at_minus = []
    at_plus = []
    for f in seq:
        if not f:
            at_minus.append(0)
            at_plus.append(0)
            continue
        lead_sign = 1 if f[-1] > 0 else -1
        at_plus.append(lead_sign)
        at_minus.append(lead_sign * (-1 if (len(f) - 1) % 2 else 1))
    return variations(at_minus) - variations(at_plus)


def positive_on_R(p):
    p = trim(p)
    return bool(p) and p[-1] > 0 and real_root_count(p) == 0


def rising(a, r):
    v = 1
    for i in range(r):
        v *= a + i
    return v


def associated_rows(N):
    A = [[0] * (N // 2 + 2) for _ in range(N + 1)]
    A[0][0] = 1
    for m in range(1, N + 1):
        for r in range(1, m // 2 + 1):
            A[m][r] = r * A[m - 1][r] + (m - 1) * A[m - 2][r - 1]
    return A


def G(m, A):
    """E_m/t = sum_r r! A_{m,r} t^(r-1)."""
    return trim([factorial(r) * A[m][r] for r in range(1, m // 2 + 1)])


def Q(m, a, A):
    """B_{m,a}/t = sum_r (a)_r A_{m,r} t^(r-1)."""
    return trim([rising(a, r) * A[m][r] for r in range(1, m // 2 + 1)])


def R(m, A):
    """Unweighted D_m/t = sum_r A_{m,r} t^(r-1)."""
    return trim([A[m][r] for r in range(1, m // 2 + 1)])


def quad_disc(p):
    p = trim(p)
    assert len(p) == 3
    return p[1] ** 2 - 4 * p[0] * p[2]


def cubic_disc(p):
    """Discriminant of a*t^3+b*t^2+c*t+d."""
    p = trim(p)
    assert len(p) == 4
    d, c, b, a = p
    return b * b * c * c - 4 * a * c ** 3 - 4 * b ** 3 * d - 27 * a * a * d * d + 18 * a * b * c * d


# ---- Bivariate exact check of the corrected component decomposition. ----

def badd(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, ZERO) + v
        if r[k] == 0:
            del r[k]
    return r


def bmul(p, q):
    r = {}
    for (i, j), a in p.items():
        for (k, l), b in q.items():
            key = (i + k, j + l)
            r[key] = r.get(key, ZERO) + a * b
    return {k: v for k, v in r.items() if v}


def bscale(c, p):
    return {k: F(c) * v for k, v in p.items() if F(c) * v}


def bmon(ix, it, c=1):
    return {(ix, it): F(c)}


def O_rows(N):
    rows = [bmon(0, 0)]
    x = bmon(1, 0)
    t = bmon(0, 1)
    for n in range(1, N + 1):
        p = bscale(n, bmul(x, rows[n - 1]))
        for r in range(2, n + 1):
            p = badd(p, bscale(comb(n, r), bmul(t, rows[n - r])))
        rows.append(p)
    return rows


def decomposition_audit(N=12):
    A = associated_rows(N)
    O = O_rows(N)
    bad = []
    for n in range(2, N + 1):
        rhs = bmon(n, 0, factorial(n))  # the m=0, r=0 term
        for k in range(0, n - 1):
            m, a = n - k, k + 1
            B = {}
            for r in range(1, m // 2 + 1):
                c = F(factorial(n), factorial(m)) * rising(a, r) * A[m][r]
                B = badd(B, bmon(k, r, c))
            rhs = badd(rhs, B)
        residual = badd(O[n], bscale(-1, rhs))
        if residual:
            bad.append((n, residual))
    print("corrected decomposition residual (n=2..12):", "ZERO" if not bad else bad)
    return not bad


# ---- Actual finite evidence (not a proof). ----

def actual_finite_audit(N=18):
    A = associated_rows(N)
    h_bad = []
    for m in range(2, N):
        for k in range(1, (m + 1) // 2 + 1):
            residual = k * A[m + 1][k + 1] * A[m][k] - (k + 1) * A[m + 1][k] * A[m][k + 1]
            if residual < 0:
                h_bad.append((m, k, residual))
    i1_bad = []
    i2_bad = []
    for M in range(5, N + 1):
        D = [None] + [Q(M - a, a, A) for a in range(1, M - 1)]
        for a in range(1, M - 3):  # m=M-a >=4; constant-constant pairs omitted
            if not positive_on_R(W(D[a], D[a + 1])):
                i1_bad.append((M, a))
        for a in range(1, M - 4):  # distance two; m>=5
            if not positive_on_R(W(D[a], D[a + 2])):
                i2_bad.append((M, a))
    print(f"actual H failures (m<= {N-1}):", len(h_bad))
    print(f"actual adjacent I1 failures (M<= {N}):", len(i1_bad), i1_bad[:5])
    print(f"actual distance-two I2 failures (M<= {N}):", len(i2_bad), i2_bad[:5])


# ---- Rational feasible countermodels. ----

def model_rows(b):
    # Only rows 2..6 are needed.  These are hypotheses for a logical
    # counterexample, not a claim that the associated-Stirling recurrence holds.
    A = [[0] * 5 for _ in range(7)]
    A[2][1] = 1
    A[3][1] = 1
    A[4][1], A[4][2] = 1, b
    A[5][1], A[5][2] = 1, F(3, 10)
    A[6][1], A[6][2], A[6][3] = 1, F(3, 4), F(1, 20)
    return A


def model_audit(b):
    A = model_rows(b)
    print("\ncountermodel A_{4,2} =", b)
    print("rows A4,A5,A6:", A[4][1:3], A[5][1:3], A[6][1:4])

    h_res = []
    for m in (4, 5):
        for k in (1, 2):
            h_res.append((
                m, k,
                k * A[m + 1][k + 1] * A[m][k]
                - (k + 1) * A[m + 1][k] * A[m][k + 1],
            ))
    print("H residuals (m,k,value):", h_res)

    g4, g5, g6 = G(4, A), G(5, A), G(6, A)
    factorial_pairs = [
        ("G5,G4", W(g5, g4)),
        ("G5,UG4", W(g5, U(g4))),
        ("G6,G5", W(g6, g5)),
        ("G6,UG5", W(g6, U(g5))),
    ]
    print("factorial/coupled W >0 on R:", {
        name: positive_on_R(w) for name, w in factorial_pairs
    })

    r4, r5, r6 = R(4, A), R(5, A), R(6, A)
    print("unweighted Bona-Mezo horizontal W >0:", {
        "R5,R4": positive_on_R(W(r5, r4)),
        "R6,R5": positive_on_R(W(r6, r5)),
    })

    horizontal = {}
    for a in (1, 2, 3):
        horizontal[a] = positive_on_R(W(Q(6, a, A), Q(5, a, A)))
    print("Q horizontal W >0 (a=1,2,3):", horizontal)

    vertical_bad = []
    for m in (4, 5, 6):
        for a in range(1, m):
            residual = sub(
                Q(m, a + 1, A),
                scale(F(1, a), add(U(Q(m, a, A)), scale(a - 1, Q(m, a, A)))),
            )
            if residual:
                vertical_bad.append((m, a, residual))
    print("vertical identity residuals:", "ZERO" if not vertical_bad else vertical_bad)

    d1, d2, d3 = Q(6, 1, A), Q(5, 2, A), Q(4, 3, A)
    w12, w23, w13 = W(d1, d2), W(d2, d3), W(d1, d3)
    print("D1,D2,D3:", d1, d2, d3)
    print("I1 W12,W23 >0 on R:", positive_on_R(w12), positive_on_R(w23))
    print("I2 W13:", w13, "real roots =", real_root_count(w13), ">0 on R =", positive_on_R(w13))
    print("Disc(D1+D3) =", quad_disc(add(d1, d3)))

    # M=7 (n=6), x weights after dividing O_6 by 6!: D1+6xD2+30x^2D3
    # +120x^3D4+360x^4D5, with D4=4 and D5=5; top term x^6/t.
    def weighted(x):
        s = add(d1, scale(6 * x, d2))
        s = add(s, scale(30 * x * x, d3))
        s = add(s, [480 * x ** 3])
        s = add(s, [1800 * x ** 4])
        o = add([x ** 6], [0] + s)  # x^6 + t*S
        return s, o

    s, o = weighted(ONE)
    print("x = 1, S =", s, "Disc_t S =", quad_disc(s))
    print("x = 1, O =", o, "Disc_t O =", cubic_disc(o))


def generic_boundary_audit():
    # S itself is real-rooted, but adding the independent top pole still fails.
    s = [F(2), F(3), F(1)]  # (t+1)(t+2)
    o = add([ONE], [ZERO] + s)  # 1 + t*S
    print("\ngeneric top-boundary witness:")
    print("S =", s, "Disc_t S =", quad_disc(s))
    print("1+tS =", o, "Disc_t(1+tS) =", cubic_disc(o))


def main():
    print("=== exact decomposition ===")
    decomposition_audit()
    print("\n=== actual finite evidence (NOT a proof) ===")
    actual_finite_audit()
    print("\n=== feasible countermodels ===")
    model_audit(F(1, 100))  # satisfies I1, violates I2
    model_audit(F(1, 10))   # satisfies I1 and I2, violates tail/top conditions
    generic_boundary_audit()


if __name__ == "__main__":
    main()
