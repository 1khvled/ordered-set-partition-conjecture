import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)
for n in [6,8,10,14,20,26]:
    for av in [1,2,5,50]:
        al = np.array([rising(av+1,j)*A[n][j+1] for j in range(n//2)], float)
        c = al*np.array([(-1)**j for j in range(n//2)])
        roots = sorted(z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0)
        z = roots[0]
        r = al[:-1]/al[1:]  # r[0]=alpha_0/alpha_1 etc.
        factors = [((j+1)/j)**2 * z / r[j] for j in range(1, len(r))]
        print(f"n={n:2d} a={av:2d}: max_j ((j+1)/j)^2 rho1/r_j = {max(factors):.4f} (j={np.argmax(factors)+1})")
