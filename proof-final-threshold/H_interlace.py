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

def H(nmax,j,Oj):
    # H_n = sum_{m>=1} t^m m! S_{>=j}(n,m); via EGF H = 1/(1 - t C_j/(1-A_j))
    # use recurrence from notes: H_n = sum_{k=1}^{j-1} binom(n,k)H_{n-k} + t sum_{q=j}^n binom(n,q)H_{n-q}
    h=[sp.Integer(1)]
    for n in range(1,nmax+1):
        if n < j:
            h.append(sp.Integer(0)); continue
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
    # p,q lists sorted roots; check alternating weak
    if not p or not q: return True
    a,b=(p,q) if len(p)<=len(q) else (q,p)
    # smaller a inside larger b: b1<=a1<=b2<=...<=b_len_a<=a_last? handle
    m=len(a)
    # require a_i between b_i and b_{i+1} for i<m (with b_0=-inf, b_{m+1}=+inf convention)
    for i in range(m):
        lo = b[i] if i < len(b) else None
        hi = b[i+1] if i+1 < len(b) else None
        if lo is not None and a[i] < lo - 1e-9: return False
        if hi is not None and a[i] > hi + 1e-9: return False
    return True

for j in [3,4]:
    N=14
    Oj=O(N,j); Hj=H(N,j,Oj)
    print('=== j',j)
    for n in range(j, N+1):
        for s in range(0, 2*j):
            m=n-s
            if m<j: continue
            rh=roots(Hj[m])
            ro=roots(Oj[n])
            if not interlace(rh,ro):
                print(' FAIL H_%d vs O_%d'%(m,n), 'deg',len(rh),len(ro))
