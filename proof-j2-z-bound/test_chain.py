"""Finite evidence: interlacing chain sigma<rho<zeta<eta<sigma_next for PF sequences."""
import numpy as np
rng=np.random.default_rng(11)
def rts(co):
    # co: low-to-high coefficients; positive roots
    return np.sort([z.real for z in np.roots(co[::-1]) if abs(z.imag)<1e-7 and z.real>0])
bad=0; tot=0; ex=None
for trial in range(3000):
    k=int(rng.integers(3,8))
    base=10**rng.uniform(-2,0); fac=rng.uniform(1.5,6)
    rs=base*fac**np.arange(k)
    # polynomial in t with roots -rs: high-to-low = np.poly(-rs)
    ch=np.poly(-rs)  # high-to-low
    # p(x) = poly(-x)*(-1)^deg? subst t=-x -> coefficients flip sign by (-1)^(deg-i)
    n=k
    p=np.array([ch[i]*(-1)**(n-i) for i in range(n+1)])  # high-to-low in x
    pl=p[::-1]  # low-to-high
    d1h=np.polyder(p); d1=d1h[::-1]
    d2h=np.polyder(d1h); d2=d2h[::-1]
    g=np.polymul(np.array([1.,0.]),d2[::-1])[::-1]+d1  # (x p')' low-to-high
    av=2.
    v=np.zeros(len(pl))
    for i,ci in enumerate(pl):
        v[i]+=(av+1)*ci
        if i>=1: v[i]+=i*ci
    v/=av
    rho=rts(pl); eta=rts(d1); zet=rts(g); sig=rts(v)
    m=min(len(rho),len(eta),len(zet),len(sig))
    for j in range(m-1):
        tot+=1
        ok = sig[j]<rho[j]<zet[j]<eta[j]<sig[j+1]
        if not ok:
            bad+=1
            if ex is None: ex=(rs,j,sig,rho,zet,eta)
print("checked",tot,"failures",bad)
if ex: print("example roots",ex[0],"j",ex[1]); print("sig",ex[2]); print("rho",ex[3]); print("zet",ex[4]); print("eta",ex[5])
