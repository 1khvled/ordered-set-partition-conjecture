"""Compare roots zeta of (xR')' with rho_j, sigma_j (finite evidence)."""
import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)
def coefs(m, av):
    return np.array([rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
def roots_alt(m, av):
    c = coefs(m,av)*np.array([(-1)**(r-1) for r in range(1,m//2+1)])
    return np.sort([z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0])
def zeta(m,av):
    c = coefs(m,av)*np.array([(-1)**(r-1) for r in range(1,m//2+1)])
    d1=np.polyder(c[::-1])[::-1]; d2=np.polyder(d1[::-1])[::-1]
    g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1
    return np.sort([z.real for z in np.roots(g[::-1]) if abs(z.imag)<1e-9 and z.real>0])
for n in [6,8,10,14]:
    for av in [2,5]:
        rr=roots_alt(n,av); ss=roots_alt(n,av+1); zz=zeta(n,av)
        print(f"n={n} a={av}")
        print("  rho:",[round(r,6) for r in rr])
        print("  sig:",[round(s,6) for s in ss])
        print("  zet:",[round(z,6) for z in zz])
