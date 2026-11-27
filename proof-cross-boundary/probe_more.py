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
def Hpoly(O,j,n):
    return sp.expand(t*sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))

# check W(S_n, O_{n-j}) all coefficients nonpositive
print('--- W(S_n,O_{n-j}) coeff signs ---')
for j in range(2,8):
    nmax=min(5*j,28)
    O=rows(j,nmax)
    bad=0; tested=0
    for n in range(j,nmax+1):
        q=O[n-j]
        S=sp.expand(sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))
        w=W(S,q)
        if w==0: continue
        tested+=1
        coeffs=sp.Poly(w,t).all_coeffs()
        if any(c>0 for c in coeffs):
            bad+=1
            print('POS',j,n)
    print(f'j={j} tested={tested} bad={bad}')

# check W(H_n,H_k) signs
print('--- W(H_n,H_k) real roots ---')
for j in range(2,6):
    nmax=min(4*j,22)
    O=rows(j,nmax)
    Hs=[sp.Integer(1)]+[Hpoly(O,j,n) for n in range(1,nmax+1)]
    bad=0; tested=0
    for n in range(j,nmax+1):
        for k in range(0,n):
            w=W(Hs[n],Hs[k])
            if w==0: continue
            tested+=1
            rc=sp.count_roots(sp.Poly(w,t))
            if rc!=0: bad+=1
    print(f'j={j} tested={tested} bad={bad}')
