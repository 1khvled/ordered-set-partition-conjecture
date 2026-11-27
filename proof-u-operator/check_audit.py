"""Check the external audit's counterexample against the actual operator lemma."""
from __future__ import annotations

import sympy as sp

x = sp.symbols("x")


def W(f, g):
    return sp.expand(sp.diff(f, x) * g - f * sp.diff(g, x))


def U(f):
    return sp.expand(2 * f + x * sp.diff(f, x))


def roots(poly):
    return sorted(sp.Poly(poly, x).all_roots(), key=lambda r: float(sp.N(r)))


p = x + sp.Rational(9, 10)   # external counterexample p
q = x + 1                    # external counterexample q
Up, Uq = U(p), U(q)
Q = x * Uq

print("p root", roots(p), "q root", roots(q))
print("Up", Up, "root", roots(Up))
print("xUq roots", roots(Q))
print("p interlaces xUq?", roots(Q)[0] < roots(p)[0] < roots(Q)[1])
print("Up interlaces xUq?", roots(Q)[0] < roots(Up)[0] < roots(Q)[1])
w = sp.factor(W(Up, Q))
print("W(Up,xUq) =", w)
print("discriminant =", sp.discriminant(w, x))
print("W(0) =", w.subs(x, 0))

eps = sp.Rational(1, 100)
Ae = p * (x + 3 * eps) * (x + eps)
Be = q * (x + 4 * eps) * (x + 2 * eps)
print("Aeps roots", [sp.N(r) for r in roots(Ae)])
print("Beps roots", [sp.N(r) for r in roots(Be)])
ra, rb = roots(Ae), roots(Be)
print("interlaced B<A<B<A<B<A?", all(rb[i] < ra[i] < rb[i+1] for i in range(len(rb)-1)) and ra[-1] > rb[-1])
