"""Search for a 3-term recurrence in r for F_r=T_r with polynomial coefficients."""
import sympy as sp
from explore import T, MAXR

xs = sp.symbols('x')

# Try  F_r = (a0 + a1 x + a2 x^2) F_{r-1} + (b0 + b1 x + b2 x^2) F_{r-2}
# solve for a_i, b_i as rational functions of r, per parity, using exact linear systems.
for parity in (0, 1):
    print('=== parity r mod 2 =', parity, '===')
    for r in [p for p in range(5, 12) if p % 2 == parity][:3]:
        n = T[r].degree()
        # number of unknowns: 3+3=6, equations n+1 (could be fewer unknowns)
        Arows, brhs = [], []
        for j in range(n + 1):
            row = []
            for i in range(3):
                row.append(sp.expand(xs ** i * T[r - 1].as_expr()).coeff(xs, j))
            for i in range(3):
                row.append(sp.expand(xs ** i * T[r - 2].as_expr()).coeff(xs, j))
            Arows.append(row)
            brhs.append(T[r].as_expr().coeff(xs, j))
        sol = sp.linsolve((sp.Matrix(Arows), sp.Matrix(brhs)))
        print('  r=%d n=%d solution:' % (r, n), sol)
