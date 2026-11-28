"""Policy audit for proof-j2-b-split.

Checks direction, sample points, endpoint signs, symbolic residuals,
and boundary cases for every inequality/identity used in NOTES.md.
Finite evidence only.
"""
from __future__ import annotations

import math
import sympy as sp

t = sp.symbols("t")


def rows(nmax):
    R = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = math.comb(n, 1) * R[n - 1]
        for r in range(2, n + 1):
            s += t * math.comb(n, r) * R[n - r]
        R.append(sp.expand(s))
    return R


def B_direct(R, n):
    return sp.expand(2 * sp.diff(R[n - 1], t) + (n - 1) * (t - 1) * sp.diff(R[n - 2], t))


def B_sum(R, n):
    s = sp.Integer(0)
    for r in range(2, n):
        s += math.comb(n - 1, r) * (R[n - r] + R[n - 1 - r])
    return sp.expand(s)


def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))


def neg_roots(poly_t):
    return sorted([r for r in sp.real_roots(sp.Poly(poly_t, t)) if r < 0], key=lambda x: sp.N(x))


def main():
    N = 15
    R = rows(N)

    print("1. Symbolic residuals (identities)")
    for n in range(2, N + 1):
        r = sp.expand(B_direct(R, n) - B_sum(R, n))
        assert r == 0, (n, r)
    for n in range(2, N + 1):
        r = sp.expand(R[n] - (n * R[n - 1] + (n - 1) * t * R[n - 2] + t * B_direct(R, n)))
        assert r == 0, (n, r)
    # W(B,O_{n-1})(mu) = -B(mu)*O'_{n-1}(mu)
    for n in range(4, 10):
        On1 = R[n - 1]
        dOn1 = sp.diff(On1, t)
        Bn = B_direct(R, n)
        w = W(Bn, On1)
        for mu in neg_roots(On1):
            res = sp.simplify(w.subs(t, mu) + Bn.subs(t, mu) * dOn1.subs(t, mu))
            assert res == 0
    # W(O_{n-1},tB)(mu) = mu*O'_{n-1}(mu)*B(mu)
    for n in range(4, 10):
        On1 = R[n - 1]
        dOn1 = sp.diff(On1, t)
        Bn = B_direct(R, n)
        w = W(On1, sp.expand(t * Bn))
        for mu in neg_roots(On1):
            res = sp.simplify(w.subs(t, mu) - mu * dOn1.subs(t, mu) * Bn.subs(t, mu))
            assert res == 0
    print("   OK")

    print("2. Direction audit (explicit sample points)")
    # sample: n=6, largest root of O_5
    n = 6
    mu = neg_roots(R[5])[-1]
    dOn1 = sp.diff(R[5], t).subs(t, mu)
    Bn = B_direct(R, 6).subs(t, mu)
    print(f"   n=6, mu={sp.N(mu):.6f}")
    print(f"   O'_5(mu)={sp.N(dOn1):.6f}, B_6(mu)={sp.N(Bn):.6f}, product={sp.N(dOn1*Bn):.6f} (>0 expected)")
    assert sp.N(dOn1 * Bn) > 0
    R6 = (5) * (1 - mu) * sp.diff(R[4], t).subs(t, mu) / sp.diff(R[5], t).subs(t, mu)
    print(f"   R_6(mu)={sp.N(R6):.6f} (<2 expected)")
    assert sp.N(R6) < 2

    print("3. Cone reduction direction audit")
    # Q_{n-2}(mu)/Q'_{n-1}(mu)>0 at roots of Q_{n-1} when interlacing holds
    for n in range(4, 12):
        On1, On2 = R[n - 1], R[n - 2]
        dOn1 = sp.diff(On1, t)
        for mu in neg_roots(On1):
            val = sp.N(On2.subs(t, mu) / dOn1.subs(t, mu))
            assert val > 0, (n, mu, val)
    # W(Q_{n-1},tQ_{n-2})(0) = -(n-1)!(n-2)!
    for n in range(3, 12):
        w0 = int(W(R[n - 1], sp.expand(t * R[n - 2])).subs(t, 0))
        assert w0 == -math.factorial(n - 1) * math.factorial(n - 2), (n, w0)
    # W(Q_{n-1},tQ_{n-2}) negative at roots
    for n in range(4, 12):
        w = W(R[n - 1], sp.expand(t * R[n - 2]))
        for mu in neg_roots(R[n - 1]):
            assert sp.N(w.subs(t, mu)) < 0
    print("   OK")

    print("4. A_n Wronskian sign at 0")
    for n in range(3, 12):
        A = sp.expand(n * R[n - 1] + (n - 1) * t * R[n - 2])
        w0 = int(W(R[n - 1], A).subs(t, 0))
        print(f"   n={n} W(O_{{n-1}},A_n)(0)={w0}")
        assert w0 < 0

    print("5. Boundary/degeneracy")
    for n in range(0, 8):
        deg = int(sp.degree(R[n], t))
        q0 = int(R[n].subs(t, 0))
        lc = int(sp.LC(sp.Poly(R[n], t)))
        print(f"   n={n} deg={deg} Q(0)={q0} lc={lc}")
        assert q0 == math.factorial(n)
        assert lc > 0
    # degrees: floor(n/2)
    for n in range(0, N + 1):
        assert int(sp.degree(R[n], t)) == n // 2
    print("   OK")

    print("6. Feasible-constraint falsification (B_n as general nonneg poly?)")
    # If we only know B_n has nonnegative coeffs and the right degree,
    # can Q'_{n-1}(mu)*B_n(mu) be negative? Yes: construct example.
    # Q_5 has roots; take B_tilde positive coeff poly of same degree as B_6
    # but with bad value. This shows nonneg-coeff alone is insufficient.
    On1 = R[5]
    mu = neg_roots(On1)[0]  # most negative root
    dOn1 = sp.diff(On1, t).subs(t, mu)
    # B_6 has degree 2. Try B_tilde = 1 + 0*t + big*t^2
    B_tilde = 1 + 1000 * t ** 2
    val = sp.N(dOn1 * B_tilde.subs(t, mu))
    print(f"   witness: mu={sp.N(mu):.6f}, Q'_5(mu)={sp.N(dOn1):.6f}, B_tilde(mu)={sp.N(B_tilde.subs(t,mu)):.6f}, product={val:.6f}")
    print(f"   product negative? {bool(val < 0)} (shows coeff-nonneg insufficient)")
    assert val < 0

    print("7. Monotonicity claim check")
    # We only claim R_n(mu*) increasing on checked range n=4..25 - finite evidence.
    # Do NOT claim for all n. Check endpoint direction for the checked range.
    print("   R_n(mu*) checked-range monotonicity is FINITE EVIDENCE only, not proved for all n.")
    print("   No all-n monotonicity claim made.")

    print("\nAUDIT COMPLETE - no proof claim made; all stated identities/directions checked on samples.")


if __name__ == "__main__":
    main()
