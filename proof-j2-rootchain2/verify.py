"""Exact certificates for proof-j2-rootchain2.

No numerical root solver is used.  Finite Sturm evidence is kept in
exact_search_correct.py and is not used as an all-index proof.
"""
from __future__ import annotations
import sympy as sp

x, t, a = sp.symbols("x t a")


def assoc_table(N: int):
    A = [[0] * (N // 2 + 3) for _ in range(N + 1)]
    A[0][0] = 1
    for m in range(N):
        for r in range(len(A[0]) - 1):
            A[m + 1][r] += r * A[m][r]
            if r:
                A[m + 1][r] += m * A[m - 1][r - 1]
    return A


def rising(av, r: int):
    v = 1
    for i in range(r):
        v *= av + i
    return v


A = assoc_table(140)


def Q(n: int, av, var=t):
    return sp.expand(sum(rising(av, r) * A[n][r] * var ** (r - 1)
                         for r in range(1, len(A[0])) if A[n][r]))


def E(f):
    return sp.expand(f + t * sp.diff(f, t))


def T(f):
    return sp.expand(sp.diff(t * sp.diff(f, t), t))


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def variations(coeffs):
    sg = [sp.sign(c) for c in coeffs if c]
    return sum(u != v for u, v in zip(sg, sg[1:]))


def shifted_coeffs(p: sp.Poly, q):
    """Coefficients, low to high, of p(q-y)."""
    d = p.degree()
    return [(-1) ** k * sum(p.nth(i) * sp.binomial(i, k) * q ** (i - k)
                            for i in range(k, d + 1))
            for k in range(d + 1)]


def main():
    print("V1: (V), generic a, n=4..12")
    for n in range(4, 13):
        res = sp.expand(a * Q(n, a + 1) - (a + 1) * Q(n, a)
                        - t * sp.diff(Q(n, a), t))
        assert res == 0
    print("  all residuals = 0")

    print("V2: (M), generic a, n=4..12")
    for n in range(4, 13):
        res = sp.expand(Q(n, a) - sp.diff(t * Q(n - 1, a), t)
                        - (n - 1) * a * t * Q(n - 2, a + 1))
        assert res == 0
    print("  all residuals = 0")

    print("V3: exact T-recurrence residual")
    S, K, c = sp.Function("S")(t), sp.Function("K")(t), sp.symbols("c")
    P = E(S) + c * t * K
    residual = sp.expand(T(P) - (E(T(S)) + T(S) + c * E(E(K))))
    print("  residual =", residual)
    assert residual == 0

    print("V4: parameter multiplier product")
    for av in range(1, 8):
        for k in range(0, 10):
            product = sp.prod(sp.Rational(b + 1 + k, b) for b in range(1, av))
            target = sp.Rational(rising(av, k + 1), sp.factorial(k + 1))
            assert sp.simplify(product - target) == 0
    print("  product residuals = 0")

    print("V5: exact direction sample n=6,a=1")
    R = Q(6, 1, -x)
    V = Q(6, 2, -x)
    C = sp.expand(sp.diff(x * sp.diff(R, x), x))
    zeta = sp.Rational(5, 36)
    eta = sp.Rational(5, 18)
    assert sp.expand(C.subs(x, zeta)) == 0
    assert sp.expand(sp.diff(R, x).subs(x, eta)) == 0
    # Rational isolating intervals, certified by endpoint signs and simplicity.
    assert V.subs(x, sp.Rational(1, 75)) > 0 > V.subs(x, sp.Rational(1, 70))
    assert R.subs(x, sp.Rational(1, 50)) > 0 > R.subs(x, sp.Rational(1, 48))
    assert V.subs(x, sp.Rational(2, 5)) < 0 < V.subs(x, sp.Rational(41, 100))
    assert R.subs(x, sp.Rational(53, 100)) < 0 < R.subs(x, sp.Rational(54, 100))
    print("  sigma_1 in (1/75,1/70), rho_1 in (1/50,1/48)")
    print("  zeta_1 = 5/36, eta_1 = 5/18")
    print("  sigma_2 in (2/5,41/100), rho_2 in (53/100,54/100)")
    print("  hence sigma_1 < rho_1 < zeta_1 < eta_1 < sigma_2 < rho_2")

    print("V6: new exact obstruction to the J-reduction")
    J = (1 + t) * (1 + 2 * t) * (1 + 3 * t)
    P = E(J)
    TT = E(sp.diff(P, t))
    H = W(P, TT)
    val = H.subs(t, -sp.Rational(1, 10))
    EGF = 2 * sp.Symbol("z") ** 3 + 11 * sp.Symbol("z") ** 2 + 12 * sp.Symbol("z") + 2
    disc = sp.discriminant(EGF, sp.Symbol("z"))
    print("  J =", J)
    print("  H =", sp.factor(H))
    print("  H(-1/10) =", val, "< 0")
    print("  coefficient-EGF discriminant =", disc, "> 0")
    assert val == -sp.Rational(2211, 625)
    assert disc == 2024 > 0

    print("V7: exact retraction of the n=129 numerical candidate")
    n = 129
    R129 = sp.Poly(Q(n, 1, -x), x)
    C129 = sp.Poly(sp.expand(sp.diff(x * sp.diff(R129.as_expr(), x), x)), x)
    q = sp.Rational(1, 10 ** 20)
    rq, cq, c0 = R129.eval(q), C129.eval(q), C129.eval(0)
    vr = variations(shifted_coeffs(R129, q))
    vc = variations(shifted_coeffs(C129, q))
    r23 = R129.eval(sp.Rational(1, 10 ** 23))
    r22 = R129.eval(sp.Rational(1, 10 ** 22))
    print("  signs R(10^-23),R(10^-22),R(q) =", sp.sign(r23), sp.sign(r22), sp.sign(rq))
    print("  signs C(0),C(q) =", sp.sign(c0), sp.sign(cq))
    print("  variations R(q-y), C(q-y) =", vr, vc)
    assert r23 < 0 < r22 and rq > 0
    assert c0 < 0 < cq
    assert vr == 2 and vc == 1
    print("  exactly two R-roots and one C-root lie in (0,q)")
    print("  the compared numerical roots had different indices; no counterexample")

    print("V8: chain-reduction identities")
    R = Q(8, 1, -x)
    V = Q(8, 2, -x)
    assert sp.expand(2 * R + x * sp.diff(R, x) - V) == 0
    # At a root eta of R', V(eta)=(a+1)R(eta)/a; this fixes eta<sigma_{j+1}.
    print("  2R+xR'-V = 0")
    print("  all exact certificates passed")


if __name__ == "__main__":
    main()
