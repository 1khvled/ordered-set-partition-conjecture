"""Check decreasing-term condition for (xR')' alternating series."""
import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)
for n in [6,8,10,14,20,26]:
    for av in [1,2,5,50]:
        al = np.array([rising(av+1,j)*A[n][j+1] for j in range(n//2)], float)  # alpha_j, j=0..k-1
        k = len(al)
        c = al*np.array([(-1)**j for j in range(k)])
        roots = sorted(z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0)
        z = roots[0]
        ratios = al[:-1]/al[1:]
        factors = np.array([((j+1)/j)**2*z/ratios[j-1] for j in range(1,k)])
        print(f"n={n:2d} a={av:2d}: rho1={z:.4g} max_j ((j+1)/j)^2 rho1/r_j = {factors.max():.4f} (j={factors.argmax()+1})")
