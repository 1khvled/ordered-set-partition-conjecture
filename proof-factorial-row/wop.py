"""Wronskian of P=(x^2C)'/x and Q=(x^2E)' vs W(C, xE)."""
import sympy as sp
xs = sp.symbols('x')
C = sp.Function('C')(xs)
E = sp.Function('E')(xs)
P = 2*C + xs*sp.diff(C, xs)
Q = 2*xs*E + xs**2*sp.diff(E, xs)
W = sp.expand(sp.diff(P,xs)*Q - P*sp.diff(Q,xs))
print('W(P,Q) =')
print(sp.collect(W, [sp.diff(C,xs), sp.diff(E,xs)]))
print()
W0 = sp.expand(sp.diff(C,xs)*(xs*E) - C*sp.diff(xs*E, xs))
print('W(C,xE) =', W0)
print()
# log-derivative bracket simplification
e, d, ep, dp = sp.symbols('e d ep dp')
c = e + d
cp = ep + dp
B = (3*c + xs*c**2 + xs*cp)*(2+xs*e) - (2+xs*c)*(2+4*xs*e+xs**2*e**2+xs*ep)
print('bracket simplified:', sp.factor(sp.expand(B)))
