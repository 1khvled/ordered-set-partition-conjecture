import numpy as np, sympy as sp
from math import comb
from scipy.optimize import linprog
t=sp.symbols('t')
def O(nmax,j):
    a=[sp.Integer(1)]
    for n in range(1,nmax+1):
        v=sum(comb(n,r)*a[n-r] for r in range(1,min(j-1,n)+1))
        v+=sum(comb(n,r)*t*a[n-r] for r in range(j,n+1))
        a.append(sp.expand(v))
    return a

def basis(P,k,maxr=2):
    out=[]
    for r in range(maxr+1):
        q=(t**r)*sp.diff(P[k],t,r)
        out.append(sp.Poly(sp.expand(q),t))
    return out

for j,n in [(4,8),(3,8),(4,9),(5,12)]:
    P=O(n,j-1); Q=O(n,j)
    targets=sp.Poly(sp.expand(Q[n]),t)
    Bs=[]; names=[]
    for k in range(n+1):
        for r,q in enumerate(basis(P,k)):
            Bs.append(q); names.append((k,r))
    deg=max(targets.degree(), *(q.degree() for q in Bs))
    A=np.array([[float(q.coeff_monomial(t**i)) for q in Bs] for i in range(deg+1)])
    b=np.array([float(targets.coeff_monomial(t**i)) for i in range(deg+1)])
    res=linprog(np.zeros(len(Bs)),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
    print('j',j,'n',n,'feasible',res.success)
    if res.success:
        print(' ',[(names[i],round(v,4)) for i,v in enumerate(res.x) if v>1e-8])
