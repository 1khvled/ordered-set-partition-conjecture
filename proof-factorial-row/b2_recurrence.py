"""Symbolic derivation: B2(A,C) where A=(xC)'+(m)(x^2D)' ; express in terms of C,D."""
import sympy as sp
xs, m = sp.symbols('x m')
C = sp.Function('C')(xs)
D = sp.Function('D')(xs)
A = sp.diff(xs*C, xs) + m*sp.diff(xs**2*D, xs)


def B2(P, Q):
    Pp, Qp = sp.diff(P, xs), sp.diff(Q, xs)
    return sp.expand(xs**2*(3*P*Qp - 2*Pp*Q) + xs**3*(P*sp.diff(Qp, xs) - Qp*Pp))


expr = sp.expand(B2(A, C))
print('B2(A,C) expanded (show):')
print(sp.collect(expr, [sp.diff(C,xs), sp.diff(D,xs)]))
