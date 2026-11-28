"""Finite evidence: K=(xP0')'/P0' at roots of P_theta along multiplier path."""
import numpy as np
from ops import astirling_table, rising
A=astirling_table(30)
def alpha(n,av):
    v=1; al=[]
    for j in range(n//2):
        if j>0: v*=av+j
        al.append(float(v)*A[n][j+1])
    return np.array(al)
def K_at(n,av,xv):
    al=alpha(n,av)
    p=al*np.array([(-1)**j for j in range(len(al))])
    d1=np.polyder(p[::-1])[::-1]; d2=np.polyder(d1[::-1])[::-1]
    g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1
    return np.polyval(g[::-1],xv)/np.polyval(d1[::-1],xv)
def roots_ptheta(n,av,th):
    al=alpha(n,av)
    p=al*np.array([(1+th*j)*(-1)**j for j in range(len(al))])
    return np.sort([z.real for z in np.roots(p[::-1]) if abs(z.imag)<1e-8 and z.real>0])
for n in [8,12,20]:
    for av in [2,5]:
        worst=1e9
        for th in np.linspace(0,1/(av+1),51):
            for z in roots_ptheta(n,av,th):
                worst=min(worst,K_at(n,av,z))
        print(f"n={n} a={av}: min K over path = {worst:.5f}")
