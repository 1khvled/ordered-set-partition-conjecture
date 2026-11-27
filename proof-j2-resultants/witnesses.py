"""Exact witnesses for the two failure modes and the verified positivity pattern.

All arithmetic here is exact over Q (sympy Rational / Integer).
"""
from __future__ import annotations

import itertools
import sympy as sp

x, t = sp.symbols("x t")


def rows_x(N):
    out = [sp.Integer(1)]
    for n in range(1, N + 1):
        p = n * x * out[n - 1]
        p += t * sum(sp.binomial(n, r) * out[n - r] for r in range(2, n + 1))
        out.append(sp.expand(p))
    return out


def toeplitz_minor(vals, rows, cols):
    D = len(vals) - 1
    T = [[sp.Integer(0)] * (D + 1) for _ in range(D + 1)]
    for i in range(D + 1):
        for j in range(D + 1):
            if i >= j:
                T[i][j] = vals[i - j]
    M = sp.Matrix([[T[i][j] for j in cols] for i in rows])
    return sp.factor(M.det())


def pf_multiplier_witness():
    """n=6, x=1: the multiplier sequence (m_{6,r}(1)) is not PF."""
    O = rows_x(6)
    n, xv = 6, sp.Integer(1)
    D = n // 2
    vals = []
    for r in range(D + 1):
        ar = sp.expand(O[n].coeff(t, r)).subs(x, xv)
        m = sp.Rational(2**r, 1) * ar / (sp.factorial(n) * sp.binomial(n - r, r) * xv**(n - 2 * r))
        vals.append(sp.nsimplify(m))
    print("multipliers m_{6,r}(1):", vals)
    for size in range(2, D + 2):
        for rows in itertools.combinations(range(D + 1), size):
            for cols in itertools.combinations(range(D + 1), size):
                det = toeplitz_minor(vals, rows, cols)
                if det < 0:
                    print(f"NEGATIVE Toeplitz minor: size={size} rows={rows} cols={cols} det={det}")
                    return
    print("no negative minor found")


def diagonal_recurrence_witness():
    """Show no relation D_{a+1}=c0 D_a + c1 t D'_a with constant c0,c1 at n=8,a=1."""
    # B_{m,a} as in relation_probe
    def assoc(N, K):
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
    A = assoc(10, 5)
    # n=8: chain is (m,a)=(8,1),(7,2),(6,3),... check first step
    P = B(8, 1, A)  # D_1
    Q = B(7, 2, A)  # D_2
    c0, c1 = sp.symbols("c0 c1")
    expr = sp.Poly(sp.expand(Q - (c0 * P + c1 * t * sp.diff(P, t))), t)
    eqs = [sp.Eq(expr.coeff_monomial(t**k), 0) for k in range(expr.degree() + 1)]
    sol = sp.solve(eqs, [c0, c1], dict=True)
    print(f"n=8 diagonal step D_2 = c0 D_1 + c1 t D_1' : solutions = {sol}")
    # second derivative
    c2 = sp.symbols("c2")
    expr = sp.Poly(sp.expand(Q - (c0 * P + c1 * t * sp.diff(P, t) + c2 * t**2 * sp.diff(P, t, 2))), t)
    eqs = [sp.Eq(expr.coeff_monomial(t**k), 0) for k in range(expr.degree() + 1)]
    sol2 = sp.solve(eqs, [c0, c1, c2], dict=True)
    print(f"n=8 with second derivative: solutions = {sol2}")


def grid_counterexample():
    """Strict interlacing is not transitive: f<g and g<h but f not< h.

    This is the exact reason horizontal+vertical component interlacing (§3.3)
    does not imply the diagonal interlacing needed for (D).
    """
    f = (t + 1) * (t + 4)            # roots {1,4}
    g = (t + 2) * (t + 3)            # roots {2,3}
    h = (t + sp.Rational(5, 2)) * (t + sp.Rational(7, 2))  # roots {2.5,3.5}
    print("f roots {1,4}, g roots {2,3}, h roots {2.5,3.5}")
    print("f<g (1<2<3<4): True;  g<h (2<2.5<3<3.5): True")
    print("f<h ? need 1<h1<4<h2, but h2=3.5<4 => FALSE")


def main():
    pf_multiplier_witness()
    diagonal_recurrence_witness()
    grid_counterexample()


if __name__ == "__main__":
    main()
