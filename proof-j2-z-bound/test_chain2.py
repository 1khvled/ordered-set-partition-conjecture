"""Finite evidence: chain holds for strongly separated roots?"""
import numpy as np
rng=np.random.default_rng(13)
def rts(co): return np.sort([z.real for z in np.roots(co[::-1]) if abs(z.imag)<1e-7 and z.real>0])
for faclo in [2,3,4,6,10]:
    bad=0;tot=0
    for trial in range(2000):
        k=int(rng.integers(3,8)); base=10**rng.uniform(-2,0)
        fac=rng.uniform(faclo,faclo*2)
        rs=base*fac**np.arange(k)
        p=np.array([c*(-1)**(k-i) for i,c in enumerate(np.poly(-rs))])[::-1]
        d1=np.polyder(p[::-1])[::-1]; d2=np.polyder(d1[::-1])[::-1]
        g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1
        av=2.
        v=np.zeros(len(p))
        for j,c in enumerate(p):
            v[j]=c*((av+1)+j)/av
        rho=rts(p);eta=rts(d1);zet=rts(g);sig=rts(v)
        m=min(len(rho),len(eta),len(zet),len(sig))
        for j in range(m-1):
            tot+=1
            if not (sig[j]<rho[j]<zet[j]<eta[j]<sig[j+1]): bad+=1
    print(f"fac>={faclo}: checked {tot}, violations {bad}")
