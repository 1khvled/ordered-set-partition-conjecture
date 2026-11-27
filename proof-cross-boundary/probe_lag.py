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

j=3; n=12
O=rows(j,n)
q=O[n-j] # O_9
for k in range(0, n-j):
    w = sp.expand(sp.diff(q,t)*O[k] - q*sp.diff(O[k],t)) # W(q,O_k)
    rc = sp.count_roots(sp.Poly(w,t)) if w!=0 else 0
    print(k, 'lag', n-j-k, 'degW', sp.degree(w,t), 'real_roots', rc)
