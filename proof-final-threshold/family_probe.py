import sympy as sp
from math import comb
t=sp.symbols('t')
def O(nmax,j):
    a=[sp.Integer(1)]
    for n in range(1,nmax+1):
        v=sum(comb(n,r)*a[n-r] for r in range(1,min(j-1,n)+1))
        v+=sum(comb(n,r)*t*a[n-r] for r in range(j,n+1))
        a.append(sp.expand(v))
    return a
def H(nmax,j):
    h=[sp.Integer(1)]
    for n in range(1,nmax+1):
        if n<j: h.append(sp.Integer(0)); continue
        v=sum(comb(n,k)*h[n-k] for k in range(1,min(j-1,n)+1))
        v+=sum(comb(n,q)*t*h[n-q] for q in range(j,n+1))
        h.append(sp.expand(v))
    return h
def roots(p):
    p=sp.Poly(sp.expand(p),t)
    if p.degree()==0: return []
    rs=sp.nroots(p,n=25,maxsteps=500)
    return sorted(float(sp.re(r)) for r in rs if abs(sp.im(r))<1e-14)
def interlace(p,q):
    if not p or not q: return True
    a,b=(p,q) if len(p)<=len(q) else (q,p)
    m=len(a)
    for i in range(m):
        lo=b[i] if i<len(b) else None
        hi=b[i+1] if i+1<len(b) else None
        if lo is not None and a[i] < lo-1e-9: return False
        if hi is not None and a[i] > hi+1e-9: return False
    return True
for j in [3,4,5]:
    N=16
    Oj=O(N,j); Hj=H(N,j)
    print('=== j',j)
    for n in range(j,N+1):
        bad=[]
        for r in range(1,j):
            if n-r<0: continue
            if not interlace(roots(Hj[n]),roots(Oj[n-r])):
                bad.append(r)
        if bad: print(' H_%d does NOT interlace O_{%d-r} for r='%(n,n),bad)
