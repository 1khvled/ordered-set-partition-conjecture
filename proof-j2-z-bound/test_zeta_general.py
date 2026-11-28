"""Finite evidence: does rho_j < zeta_j hold for geometric-spacing real-rooted polys?"""
import numpy as np
rng = np.random.default_rng(3)
bad=0; tot=0
for trial in range(5000):
    k = rng.integers(2,7)
    base = 10**rng.uniform(-2,0)
    fac = rng.uniform(1.6, 8)
    rs = base*fac**np.arange(k)
    c = np.array([1.0])
    for r in rs: c = np.polymul(c, [-r,1.0])
    d1=np.polyder(c[::-1])[::-1]; d2=np.polyder(d1[::-1])[::-1]
    g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1
    rho = np.sort([z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0])
    zet = np.sort([z.real for z in np.roots(g[::-1]) if abs(z.imag)<1e-9 and z.real>0])
    for j in range(min(len(rho),len(zet))):
        tot+=1
        if zet[j] <= rho[j]: bad+=1
print("pairs checked:",tot,"violations rho_j<zeta_j:",bad)
