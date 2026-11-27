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
    nmax=min(5*j,28)
    O=rows(j,nmax)
    bad=0; tested=0
    for n in range(2*j, nmax+1):
        q=O[n-j]
        qroots=sp.real_roots(q)
        for k in range(0,n-j):
            for mu in qroots:
                val = sp.expand(sp.diff(q,t)*O[k] - q*sp.diff(O[k],t)).subs(t,mu)
                # W(q,O_k)(mu) = q'(mu)*O_k(mu)
                tested+=1
                if val <= 0:
                    bad+=1
                    print('BAD',j,n,k, 'lag',n-j-k)
    print(f'j={j} tested={tested} bad={bad}')
