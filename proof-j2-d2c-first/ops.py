"""Test whether {Q_{m,a}}_{m} (fixed a) is an orthogonal-polynomial-type
sequence satisfying a 3-term recurrence, and explore Wronskian recurrences.
"""
import sympy as sp
from explore import A, Q, W, t

a, m = sp.symbols("a m")

sep = lambda s: (print("\n" + "=" * 60), print(s), print("=" * 60))

sep("1. Test 3-term recurrence in m for fixed a")
for a0 in [1, 2, 3]:
    print("a =", a0)
    for mm in range(4, 9):
        f0 = sp.Poly(Q(A, mm, a0, t), t)
        f1 = sp.Poly(Q(A, mm - 1, a0, t), t)
        f2 = sp.Poly(Q(A, mm - 2, a0, t), t)
        alpha, beta, gamma = sp.symbols("alpha beta gamma")
        expr = sp.Poly(sp.expand((alpha * t + beta) * f1.as_expr() + gamma * f2.as_expr() - f0.as_expr()), t)
        eqs = [c for c in expr.all_coeffs()]
        sol = sp.solve(eqs, [alpha, beta, gamma], dict=True)
        print(f"  m={mm}: sol={sol}, degs=({f0.degree()},{f1.degree()},{f2.degree()})")
