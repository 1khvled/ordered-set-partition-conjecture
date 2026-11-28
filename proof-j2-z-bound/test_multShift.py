"""Finite evidence: sigma_i >= e^{-theta} rho_i under multiplier 1+theta*j."""
import numpy as np
rng=np.random.default_rng(5)
def rts(co): return np.sort([z.real for z in np.roots(co[::-1]) if abs(z.imag)<1e-7 and z.real>0])
bad=0;tot=0;ex=None
for trial in range(5000):
    k=int(rng.integers(2,8))
    base=10**rng.uniform(-2,0); fac=rng.uniform(1.1,10)
    rs=base*fac**np.arange(k)
    p=np.array([c*(-1)**(k-i) for i,c in enumerate(np.poly(-rs))])[::-1]  # low-to-high alt
    th=float(rng.uniform(0.02,0.5))
    g=np.zeros(len(p))
    for j,c in enumerate(p): g[j]=c*(1+th*j)
    rho=rts(p); sig=rts(g)
    for i in range(min(len(rho),len(sig))):
        tot+=1
        if sig[i] < np.exp(-th)*rho[i]*(1-1e-9):
            bad+=1
            if ex is None: ex=(rs,th,i,sig[i]/rho[i],np.exp(-th))
print("checked",tot,"violations",bad)
if ex: print("example",ex)
