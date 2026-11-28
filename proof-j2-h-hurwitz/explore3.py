"""Check candidate common interlacers numerically (evidence only)."""
import numpy as np
from common import Family
F = Family(30)
E, T = Family.E, Family.T

def roots(poly):
    cs=[float(c) for c in poly.all_coeffs()]
    return sorted([z.real for z in np.roots(cs) if abs(z.imag)<1e-7 and z.real<0])

for n in (10,12,14):
    print('==== n',n)
    for name,p in [('P_n',F.Qpoly(n)), ('T_n',T(F.Qpoly(n))),
                   ('U P_{n-2}', F.Qpoly(n-2)+E(F.Qpoly(n-2))),
                   ('E P_{n-1}', E(F.Qpoly(n-1))),
                   ('E^2 P_{n-2}', E(E(F.Qpoly(n-2))))]:
        print(' ',name, ['%.6g'%x for x in roots(p)])
