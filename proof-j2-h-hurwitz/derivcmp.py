"""Compare P_n' with P_{n-2}, Q_n, T_n."""
import numpy as np
from common import Family, t
F=Family(24); E=Family.E
def roots(p):
    cs=[float(c) for c in p.all_coeffs()]
    return sorted([z.real for z in np.roots(cs) if abs(z.imag)<1e-8 and z.real<0])
for n in (10,12,14):
    print('n',n)
    for name,p in [("P_n'",F.Qpoly(n).diff()),('P_{n-2}',F.Qpoly(n-2)),
                   ('Q_n=U P_{n-2}',F.Qpoly(n-2)+E(F.Qpoly(n-2))),
                   ('T_n',F.T(F.Qpoly(n)))]:
        print(' ',name, ['%.6g'%x for x in roots(p)])
