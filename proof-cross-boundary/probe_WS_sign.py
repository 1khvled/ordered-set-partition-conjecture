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
for j in range(2,5):
    nmax=min(4*j,20)
    O=rows(j,nmax)
    for n in range(j,nmax+1):
        q=O[n-j]
        S=sp.expand(sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))
        w=sp.expand(sp.diff(S,t)*q - S*sp.diff(q,t))
        if w==0: continue
        lc = sp.LC(sp.Poly(w,t))
        val0 = w.subs(t,0)
        print(f'j={j} n={n} lc={lc} w(0)={val0} deg={sp.degree(w,t)}')
