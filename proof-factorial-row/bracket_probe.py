"""Check the two brackets in W_r expression for sign."""
import sympy as sp
from explore import T
xs = sp.symbols('x')

for r in range(5, 12):
    A = T[r-1].as_expr()
    C = T[r-2].as_expr()
    Ap = sp.diff(A, xs); App = sp.diff(Ap, xs)
    Cp = sp.diff(C, xs); Cpp = sp.diff(Cp, xs)
    B1 = xs**2*(A*App - Ap**2) - A**2
    B2 = xs**2*(3*A*Cp - 2*Ap*C) + xs**3*(A*Cpp - Cp*Ap)
    W = sp.expand((r-1)*B2 + B1)
    vals1 = [float(B1.subs(xs, v)) for v in [-100, -3, -1, -0.1, -0.001, 0, 0.1, 1]]
    vals2 = [float(B2.subs(xs, v)) for v in [-100, -3, -1, -0.1, -0.001, 0, 0.1, 1]]
    print('r=%d' % r)
    print('  B1 signs', [1 if v>0 else -1 if v<0 else 0 for v in vals1])
    print('  B2 signs', [1 if v>0 else -1 if v<0 else 0 for v in vals2], ['%.3g'%v for v in vals2])
