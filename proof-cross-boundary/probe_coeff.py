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
for j in range(2,8):
    nmax=min(5*j,28)
    O=rows(j,nmax)
    bad=0; tested=0
    for n in range(j,nmax+1):
        S=sp.expand(sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))
        H=t*S
        for d in range(1,j+1):
            w=W(sp.expand(H),O[n-d])
            tested+=1
            coeffs=sp.Poly(w,t).all_coeffs()
            if any(c<0 for c in coeffs):
                bad+=1
                print('BAD',j,n,d)
    print(f'j={j} tested={tested} bad={bad}')
