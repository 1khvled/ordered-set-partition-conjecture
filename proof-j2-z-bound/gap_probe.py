"""Finite evidence: smallest-root gap vs coefficient ratios."""
import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)

def coefs(m, av, neg=False):
    K = m//2
    c = np.array([rising(av, r)*A[m][r] for r in range(1, K+1)], dtype=float)
    if neg:
        c = c*np.array([(-1)**i for i in range(K)])
    return c

def smallest(m, av):
    c = coefs(m, av, True)
    rts = sorted(z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0)
    return rts[0] if rts else None

for n in [6,8,10,12,14,20,26]:
    for av in [2,3,5,10,50]:
        rho = smallest(n, av)
        sig = smallest(n, av+1)
        a0 = A[n][1]; a1 = rising(av+1,1)*A[n][2]
        b0 = A[n][1]; b1 = rising(av+2,1)*A[n][2]
        r0 = a0/a1; s0 = b0/b1
        print(f"n={n:2d} a={av:2d}: rho={rho:.6g} r0={r0:.6g} rho/r0={rho/r0:.6f} | sigma={sig:.6g} s0={s0:.6g} sig/s0={sig/s0:.6f} | sig/rho={sig/rho:.6f} (a+1)/(a+2)={(av+1)/(av+2):.6f}")
