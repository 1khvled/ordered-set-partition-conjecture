"""Find empirical second-order ODE in u for psi_n(x,u)."""
import sympy as sp
from probe2 import psi

x, u = sp.symbols('x u')
for n in range(4, 9):
    f = psi(n).as_expr()
    # ansatz: (A2 u^2 + A1 u + A0) f'' + (B1 u + B0) f' + (C1 u + C0) f = 0
    A2, A1, A0, B1, B0, C1, C0 = sp.symbols('A2 A1 A0 B1 B0 C1 C0')
    expr = (A2*u**2 + A1*u + A0)*sp.diff(f,u,2) + (B1*u+B0)*sp.diff(f,u) + (C1*u+C0)*f
    eqs = sp.Poly(sp.expand(expr), u).all_coeffs()
    sol = sp.solve(eqs, [A2,A1,A0,B1,B0,C1,C0], dict=True)
    print('n=%d sol=%s' % (n, sol))
