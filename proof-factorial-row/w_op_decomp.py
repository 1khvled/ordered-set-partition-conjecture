"""Decompose W(D_C, x D_E) using C=(xE)'+m(x^2 H)'."""
import sympy as sp
xs, m = sp.symbols('x m')
E = sp.Function('E')(xs)
H = sp.Function('H')(xs)
C = sp.diff(xs*E, xs) + m*sp.diff(xs**2*H, xs)
Dc = 2*C + xs*sp.diff(C, xs)
De = 2*E + xs*sp.diff(E, xs)
Q = xs*De
W = sp.expand(sp.diff(Dc,xs)*Q - Dc*sp.diff(Q,xs))
print('W(Dc, xDe) =')
print(sp.collect(W, [sp.diff(E,xs), sp.diff(H,xs)]))
