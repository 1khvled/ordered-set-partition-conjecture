"""Probe: can d_x O_n be written as (a*t+b)*O_n' + c*O_n ?"""
import sympy as sp
from gen import On, x, t

a, b, c = sp.symbols('a b c')
for n in range(2, 11):
    P = sp.Poly(On(n), t)
    dP = sp.Poly(sp.diff(On(n), x), t)
    deg = P.degree()
    ans = sp.Poly(a * t * On(n) + 0, t)  # placeholder
    expr = sp.expand((a * t + b) * sp.diff(On(n), t) + c * On(n) - sp.diff(On(n), x))
    eqs = sp.Poly(expr, t).all_coeffs()
    sol = sp.solve(eqs, (a, b, c), dict=True)
    print('n=%d deg=%d  solution=%s' % (n, deg, sol))
