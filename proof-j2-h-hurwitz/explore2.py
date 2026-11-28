"""Explore interlacing landscape numerically (evidence only)."""
import numpy as np
from common import Family, t

F = Family(30)
E, T = Family.E, Family.T

def roots(poly):
    cs = [float(c) for c in poly.all_coeffs()]
    r = np.roots(cs)
    return sorted([x.real for x in r if abs(x.imag) < 1e-7 and x.real < 0])

for n in (8,10,12,14):
    print('==== n', n)
    for name, poly in [('P_{n-2}', F.Qpoly(n-2)), ('P_{n-1}', F.Qpoly(n-1)),
                       ('P_n', F.Qpoly(n)), ('T_n', T(F.Qpoly(n))),
                       ('U P_{n-1}', F.Qpoly(n-1)+E(F.Qpoly(n-1))),
                       ('T_{n-1}', T(F.Qpoly(n-1)))]:
        rr = roots(poly)
        print(f'  {name:11s}:', ['%.6g'%x for x in rr])
