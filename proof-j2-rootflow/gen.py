"""Exact generation of the j=2 deformation O_n(x,t) and root-flow probes.

This file is EXPLORATION only; no claim of proof is made here.
"""
import sympy as sp

x, t, z = sp.symbols('x t z')


def Cseries(n):
    return sum(z**r / sp.factorial(r) for r in range(2, n + 1))


def On(n):
    """O_n(x,t) = n! [z^n] 1/(1-xz-t C(z)), exact in Q[x,t]."""
    C = Cseries(n)
    total = 0
    for k in range(n + 1):
        term = sp.expand((x * z + t * C) ** k)
        poly = sp.Poly(term, z)
        if n <= poly.degree():
            total += poly.coeff_monomial(z**n)
    return sp.expand(sp.factorial(n) * total)


def poly_t(expr):
    return sp.Poly(sp.expand(expr), t)


if __name__ == '__main__':
    for n in range(0, 8):
        p = On(n)
        print(n, sp.factor(p))
