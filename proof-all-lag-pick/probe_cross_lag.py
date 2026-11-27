import sympy as sp
t=sp.Symbol('t')

def rows(j,nmax):
    O=[sp.Integer(1)]
    for n in range(1,nmax+1):
        s=sp.Integer(0)
        for r in range(1,min(j-1,n)+1):
            s+=sp.binomial(n,r)*O[n-r]
        for r in range(j,n+1):
            s+=t*sp.binomial(n,r)*O[n-r]
        O.append(sp.expand(s))
    H=[]
    for n in range(nmax+1):
        if n==0: H.append(sp.Integer(1))
        else:
            # H_n = t * sum_{r=j}^{n} binom(n,r) O_{n-r}
            s2=sp.Integer(0)
            for r in range(j,n+1):
                s2+=sp.binomial(n,r)*O[n-r]
            H.append(sp.expand(t*s2))
    return O,H

def W(p,q):
    return sp.expand(sp.diff(p,t)*q-p*sp.diff(q,t))

bad=[]
tested=0
for j in range(2,8):
    nmax=min(5*j,32)
    O,H=rows(j,nmax)
    for n in range(j,nmax+1):
        for d in range(1,j+1):
            w=W(H[n],O[n-d])
            if w==0: continue
            tested+=1
            cnt=sp.polys.polytools.count_roots(sp.Poly(w,t))
            if cnt!=0:
                bad.append((j,n,d,cnt))
                print('BAD',j,n,d,cnt,flush=True)
    print(f'j={j} done tested {tested}',flush=True)
print('TOTAL',tested,'bad',len(bad))
