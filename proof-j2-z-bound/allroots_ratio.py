"""Finite evidence: ratio sigma_j / rho_j for all roots and several n,a."""
import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)
def roots(m, av):
    c = np.array([((-1)**(r-1))*rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
    return np.sort([z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0])
for n in [6,8,10,14,20]:
    for av in [2,5,50]:
        rr = roots(n, av); ss = roots(n, av+1)
        print(f"n={n} a={av}: rho/sigma ratios = {[round(s/r,6) for r,s in zip(rr,ss)]}, a/(a+1)={av/(av+1):.6f}")
