"""Check the two componentwise interlacing properties:
(A) B_{m,a} interlaces B_{m-1,a} (fixed a, m decreases)
(B) B_{m,a} interlaces B_{m,a+1} (fixed m, a increases)
Here B_{m,a}(t)=sum_r (a)_r {m brace r}_{>=2} t^r (no division by t).
"""
from __future__ import annotations
import sympy as sp

t = sp.symbols("t")


def assoc_stirling(N, K):
    A = [[0] * (K + 2) for _ in range(N + 1)]
    A[0][0] = 1
    for m in range(1, N + 1):
        for r in range(1, min(K + 1, m // 2) + 1):
            A[m][r] = r * A[m - 1][r] + (m - 1) * A[m - 2][r - 1]
    return A


def rising(a, r):
    v = sp.Integer(1)
    for i in range(r):
        v *= a + i
    return v


def B(m, a, A):
    return sp.expand(sum(rising(a, r) * A[m][r] * t**r for r in range(1, m // 2 + 1)))


def roots(p):
    rr = sp.nroots(p, n=30, maxsteps=200)
    if any(abs(sp.im(r)) > sp.Float("1e-16") for r in rr):
        return None
    rr = sorted(float(sp.re(r)) for r in rr)
    # remove the common zero root (r=0 coefficient is 0 for m>=1)
    return [r for r in rr if abs(r) > 1e-12]


def inter(a, b):
    a = sorted([v for v in a if abs(v) > 1e-12])
    b = sorted([v for v in b if abs(v) > 1e-12])
    m = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    for i in range(len(m) - 1):
        if abs(m[i][0] - m[i + 1][0]) < 1e-10 or m[i][1] == m[i + 1][1]:
            return False
    return True


def main():
    A = assoc_stirling(12, 6)
    print("=== (A) fixed a, m decreases ===")
    for a in [1, 2, 3]:
        for m in range(5, 12):
            p, q = B(m, a, A), B(m - 1, a, A)
            if p == 0 or q == 0:
                continue
            r1, r2 = roots(p), roots(q)
            print(f"a={a} m={m}->{m-1}: inter={inter(r1, r2)} degs={sp.Poly(p,t).degree()},{sp.Poly(q,t).degree()}")
    print("=== (B) fixed m, a increases ===")
    for m in range(5, 12):
        for a in [1, 2, 3]:
            p, q = B(m, a, A), B(m, a + 1, A)
            if p == 0 or q == 0:
                continue
            r1, r2 = roots(p), roots(q)
            print(f"m={m} a={a}->{a+1}: inter={inter(r1, r2)} degs={sp.Poly(p,t).degree()},{sp.Poly(q,t).degree()}")


if __name__ == "__main__":
    main()
