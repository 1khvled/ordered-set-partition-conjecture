"""Search for a low-degree recurrence linking B_{m,a} and B_{m-1,a+1}, where
B_{m,a}(t)=sum_r (a)_r {m brace r}_{>=2} t^r.
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


def main():
    A = assoc_stirling(14, 7)
    # fit B_{m-1,a+1} = (alpha0+alpha1 t) B_{m,a} + (beta0+beta1 t) B'_{m,a} + gamma t^2 B''_{m,a}
    for m in range(4, 11):
        for a in range(1, 5):
            P = B(m, a, A)
            Q = B(m - 1, a + 1, A)
            # try Q = c0*P + c1*t*P'  (constants c0,c1)
            c0, c1 = sp.symbols("c0 c1")
            expr = sp.expand(Q - (c0 * P + c1 * t * sp.diff(P, t)))
            eqs = [sp.Eq(sp.Poly(expr, t).coeff_monomial(t**k), 0) for k in range(sp.Poly(expr, t).degree() + 1)]
            sol = sp.solve(eqs, [c0, c1], dict=True)
            print(f"m={m} a={a}: Q=c0P+c1tP' -> {sol}")
            if not sol:
                c0, c1, c2 = sp.symbols("c0 c1 c2")
                expr = sp.expand(Q - (c0 * P + c1 * t * sp.diff(P, t) + c2 * t**2 * sp.diff(P, t, 2)))
                eqs = [sp.Eq(sp.Poly(expr, t).coeff_monomial(t**k), 0) for k in range(sp.Poly(expr, t).degree() + 1)]
                sol2 = sp.solve(eqs, [c0, c1, c2], dict=True)
                print(f"   with second derivative -> {sol2}")
            break


if __name__ == "__main__":
    main()
