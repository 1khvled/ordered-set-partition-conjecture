"""Exploration for the FIRST D2C premise W(Q_{m,a}, Q_{m-2,a+1}) > 0.

NOT a proof. Exact symbolic/numeric exploration only.
"""
import sympy as sp

t, x = sp.symbols("t x", real=True)


def astirling_table(mmax):
    rmax = mmax // 2 + 2
    A = [[0] * rmax for _ in range(mmax + 1)]
    A[0][0] = 1
    for m in range(mmax):
        for r in range(rmax - 1):
            A[m + 1][r] += r * A[m][r]
            if r:
                A[m + 1][r] += m * A[m - 1][r - 1]
    return A


def rising(a, r):
    out = 1
    for k in range(r):
        out *= a + k
    return out


def Q(A, m, a, z):
    out = 0
    for r in range(1, len(A[0])):
        if A[m][r]:
            out += rising(a, r) * A[m][r] * z ** (r - 1)
    return sp.expand(out)


def W(f, g, z=t):
    return sp.expand(sp.diff(f, z) * g - f * sp.diff(g, z))


def B(f, g, z=x):
    return sp.expand(f * sp.diff(g, z) - sp.diff(f, z) * g)


A = astirling_table(24)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


sep("1. W(p,w) factorizations for small (m,a)")
for (m, a) in [(6, 1), (7, 1), (8, 1), (9, 1), (7, 2), (8, 2), (9, 2), (10, 3)]:
    p = Q(A, m, a, t)
    w = Q(A, m - 2, a + 1, t)
    Wpw = W(p, w)
    print(f"m={m} a={a}: W(p,w) = {sp.factor(Wpw)}")
    print(f"   deg p = {sp.Poly(p,t).degree()}, deg w = {sp.Poly(w,t).degree()}, "
          f"deg W = {sp.Poly(Wpw,t).degree()}, LC={sp.Poly(Wpw,t).LC()}")

sep("2. Compare W(p,w) with W(p,s), W(p,q), W(q,w), W(s,w)")
for (m, a) in [(7, 1), (8, 1), (8, 2), (9, 2)]:
    p = Q(A, m, a, t)
    s = Q(A, m - 1, a, t)
    q = Q(A, m - 1, a + 1, t)
    w = Q(A, m - 2, a + 1, t)
    print(f"m={m} a={a}")
    print("  W(p,w)=", sp.factor(W(p,w)))
    print("  W(p,s)=", sp.factor(W(p,s)))
    print("  W(p,q)=", sp.factor(W(p,q)))
    print("  W(q,w)=", sp.factor(W(q,w)))
    print("  W(s,w)=", sp.factor(W(s,w)))
