"""Symbolic W(A,D) with A=(xC)'+m(x^2E)', D=(x^2C)'/x."""
import sympy as sp
xs, m = sp.symbols('x m')
C = sp.Function('C')(xs)
E = sp.Function('E')(xs)
A = sp.diff(xs*C, xs) + m*sp.diff(xs**2*E, xs)
D = 2*C + xs*sp.diff(C, xs)
W = sp.expand(sp.diff(A, xs)*D - A*sp.diff(D, xs))
print(sp.collect(W, [sp.diff(C,xs), sp.diff(E,xs)]))
