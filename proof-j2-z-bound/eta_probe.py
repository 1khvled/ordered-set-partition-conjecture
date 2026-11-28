import numpy as np
from ops import astirling_table, rising
A=astirling_table(30)
def roots(m,av):
    c=np.array([(-1)**(r-1)*rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
    return np.sort([z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0])
def eta(m,av):
    c=np.array([(-1)**(r-1)*rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
    d=np.polyder(c[::-1])
    return np.sort([z.real for z in np.roots(d) if abs(z.imag)<1e-9 and z.real>0])
def zeta(m,av):
    c=np.array([(-1)**(r-1)*rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
    d1=np.polyder(c[::-1])[::-1]; d2=np.polyder(d1[::-1])[::-1]
    g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1
    return np.sort([z.real for z in np.roots(g[::-1]) if abs(z.imag)<1e-9 and z.real>0])
for n in [8,10,14]:
    for av in [2,5]:
        print(f"n={n} a={av}")
        print("  rho:",[round(r,5) for r in roots(n,av)])
        print("  sig:",[round(s,5) for s in roots(n,av+1)])
        print("  eta:",[round(e,5) for e in eta(n,av)])
        print("  zeta:",[round(z,5) for z in zeta(n,av)])
