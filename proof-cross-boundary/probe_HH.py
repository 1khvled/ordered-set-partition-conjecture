import sympy as sp
t=sp.symbols('t')
def rows(j,nmax):
    O=[sp.Integer(1)]
    for n in range(1,nmax+1):
        s=sp.Integer(0)
        for r in range(1,min(j-1,n)+1): s+=sp.binomial(n,r)*O[n-r]
        for r in range(j,n+1): s+=t*sp.binomial(n,r)*O[n-r]
        O.append(sp.expand(s))
    return O
def W(p,q): return sp.expand(sp.diff(p,t)*q-p*sp.diff(q,t))
def Hpoly(O,j,n): return sp.expand(t*sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))
for j in range(2,7):
    nmax=min(5*j,30)
    O=rows(j,nmax)
    bad=0; tested=0
    for n in range(2*j,nmax+1):
        w=W(Hpoly(O,j,n), Hpoly(O,j,n-j))
        if w==0: continue
        tested+=1
        rc=sp.count_roots(sp.Poly(w,t))
        if rc!=0: bad+=1
    print(f'j={j} tested={tested} bad={bad}')
