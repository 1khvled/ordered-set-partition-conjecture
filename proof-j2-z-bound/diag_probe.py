"""Finite evidence: rho_j(Q_{n,a}) vs rho_1(Q_{n-2(j-1), a+j-1})."""
import numpy as np
from ops import astirling_table, rising
A=astirling_table(30)
def roots(m,av):
    if m<4: return []
    c=np.array([(-1)**(r-1)*rising(av,r)*A[m][r] for r in range(1,m//2+1)],float)
    return np.sort([z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0])
for n in [10,14,20]:
    for av in [2,5]:
        rr=roots(n,av); ss=roots(n,av+1)
        print(f"n={n} a={av}")
        for j,r in enumerate(rr):
            m=n-2*j; aa=av+j
            r1=roots(m,aa)
            r2=roots(m,aa+1)
            if len(r1) and len(r2):
                print(f"  j={j+1}: rho_j={r:.6g} rho1(Q_{m},{aa})={r1[0]:.6g} ratio={r/r1[0]:.4f} | sig_j={ss[j]:.6g} rho1(Q_{m},{aa+1})={r2[0]:.6g} ratio={ss[j]/r2[0]:.4f}")
