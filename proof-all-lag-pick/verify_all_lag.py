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
    H=[sp.Integer(1)]
    for n in range(1,nmax+1):
        s2=sp.Integer(0)
        for r in range(j,n+1):
            s2+=sp.binomial(n,r)*O[n-r]
        H.append(sp.expand(t*s2))
    return O,H

def W(p,q):
    return sp.expand(sp.diff(p,t)*q-p*sp.diff(q,t))

def deg(p):
    p=sp.expand(p)
    if p==0: return -1
    return sp.Poly(p,t).degree()

# 1) All-lag Pick: W(O_n,O_{n-d}) has no real roots
bad=[]; tested=0
for j in range(2,10):
    nmax=min(5*j,34)
    O,_=rows(j,nmax)
    for n in range(j,nmax+1):
        for d in range(1,j+1):
            w=W(O[n],O[n-d])
            if w==0: continue
            tested+=1
            if sp.polys.polytools.count_roots(sp.Poly(w,t))!=0:
                bad.append(('pick',j,n,d))
print('all-lag pick tested',tested,'bad',len(bad))

# 2) Cross-family all-lag: W(H_n,O_{n-d}) has no real roots
bad2=[]; tested2=0
for j in range(2,9):
    nmax=min(5*j,30)
    O,H=rows(j,nmax)
    for n in range(j,nmax+1):
        for d in range(1,j+1):
            w=W(H[n],O[n-d])
            if w==0: continue
            tested2+=1
            if sp.polys.polytools.count_roots(sp.Poly(w,t))!=0:
                bad2.append(('cross',j,n,d))
print('cross all-lag tested',tested2,'bad',len(bad2))

# 3) Decomposition identity
# W(O_n,O_{n-d}) = W(H_n,O_{n-d}) + sum_{r<d} binom(n,r)*W(O_{n-r},O_{n-d}) - sum_{r>d} binom(n,r)*W(O_{n-d},O_{n-r})
id_bad=[]
for j in range(2,7):
    O,H=rows(j,4*j)
    for n in range(j,4*j+1):
        for d in range(1,j+1):
            lhs=W(O[n],O[n-d])
            rhs=W(H[n],O[n-d])
            for r in range(1,j):
                if r==d: continue
                if r<d:
                    rhs = sp.expand(rhs + sp.binomial(n,r)*W(O[n-r],O[n-d]))
                else:
                    rhs = sp.expand(rhs - sp.binomial(n,r)*W(O[n-d],O[n-r]))
            if sp.expand(lhs-rhs)!=0:
                id_bad.append((j,n,d))
print('decomposition identity bad',len(id_bad))

# 4) Degree-diff-0 cases exist
count0=0
for j in range(2,7):
    O,_=rows(j,4*j)
    for n in range(j,4*j+1):
        for d in range(1,j+1):
            if deg(O[n])-deg(O[n-d])==0:
                count0+=1
print('deg-diff-0 cases',count0)

# 5) Mixed-sign decomposition: show a negative term occurs for d<j
mixed=0
for j in range(3,7):
    O,H=rows(j,4*j)
    for n in range(j,4*j+1):
        for d in range(1,j): # d<j
            # find if any r>d term is present
            has_neg=False
            for r in range(d+1,j):
                w= W(O[n-d],O[n-r]) # should be Pick-positive by conjecture
                if w!=0:
                    has_neg=True
            if has_neg: mixed+=1
print('cases with negative decomposition terms (d<j)',mixed)

# 6) check Pick residues negative for a sample
res_bad=[]
for j in [2,3,4]:
    O,_=rows(j,4*j)
    for n in range(j,4*j+1):
        for d in range(1,j+1):
            p,q=O[n],O[n-d]
            # at roots of q, p(mu)/q'(mu) < 0
            qpoly=sp.Poly(q,t)
            if qpoly.degree()<=0: continue
            # use Sturm/count and sign via factorization is hard; check via W>0 at roots
            # W(mu)= -r_i * q'(mu)^2, so W(mu)>0 means r_i<0
            # already checked W has no real roots => W>0 everywhere (lc>0?)
            w=W(p,q)
            lc=sp.Poly(w,t).LC()
            if lc<=0:
                res_bad.append((j,n,d,'lc',lc))
print('W leading coeff <=0 cases',len(res_bad))
