"""Verify the operator identities for B_{m,a}(t)=sum_r (a)_r {m brace r}_{>=2} t^r:
(I)  B_{m,a+1} = B_{m,a} + (1/a) t B'_{m,a}    [predicted]
(II) fit B_{m,a} = c0 B_{m-1,a} + c1 t B'_{m-1,a} + c2 t^2 B''_{m-1,a}
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
    ok_I = True
    for m in range(2, 14):
        for a in range(1, 6):
            lhs = sp.expand(B(m, a + 1, A))
            rhs = sp.expand(B(m, a, A) + sp.Rational(1, a) * t * sp.diff(B(m, a, A), t))
            if sp.expand(lhs - rhs) != 0:
                ok_I = False
                print(f"(I) FAILS m={m} a={a}")
    print("(I) B_{m,a+1}=B_{m,a}+(1/a)t B'_{m,a}:", ok_I)

    c0, c1, c2 = sp.symbols("c0 c1 c2")
    for m in range(4, 11):
        for a in [1, 2, 3]:
            P = B(m, a, A)
            Q = B(m - 1, a, A)
            expr = sp.expand(Q - (c0 * P + c1 * t * sp.diff(P, t)))
            poly = sp.Poly(expr, t)
            eqs = [sp.Eq(poly.coeff_monomial(t**k), 0) for k in range(poly.degree() + 1)]
            sol = sp.solve(eqs, [c0, c1], dict=True)
            print(f"(II) m={m} a={a}: Q=c0P+c1tP' from P=B_{m,a},Q=B_{m-1,a}: {sol}")


if __name__ == "__main__":
    main()
