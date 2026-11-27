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

for j in range(2,7):
    nmax=min(4*j,24)
    O=rows(j,nmax)
    bad=0; tested=0
    for n in range(j,nmax+1):
        S=sp.Integer(0)
        for r in range(j,n+1): S+=sp.binomial(n,r)*O[n-r]
        S=sp.expand(S)
        for d in range(1,j+1):
            q=O[n-d]
            if sp.degree(q,t)==0: continue
            # roots
            for mu in sp.real_roots(q):
                qp = sp.diff(q,t).subs(t,mu)
                val = S.subs(t,mu)
                tested+=1
                if val*qp <= 0:
                    bad+=1
                    print('BAD',j,n,d,mu)
    print(f'j={j} tested={tested} bad={bad}')
