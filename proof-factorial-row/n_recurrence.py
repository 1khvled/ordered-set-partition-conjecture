"""Try to express N_r (level-1 Wronskian) via previous N/M and squares, positive coeffs."""
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from explore import T
xs = sp.symbols('x')
Uop = lambda f: sp.expand(2*f + xs*sp.diff(f, xs))
def W(p,R): return sp.expand(sp.diff(p,xs)*R - p*sp.diff(R,xs))
def W0(r): return W(T[r].as_expr(), xs*T[r-1].as_expr())
def Nr(r):
    p,q = Uop(T[r-2].as_expr()), Uop(T[r-3].as_expr())
    return sp.Poly(sp.expand(-W(p,xs*q)/4), xs)

for r in range(8, 14):
    target = Nr(r)
    items=[]
    for j in range(0,4):
        items.append((xs**(2*j)*Nr(r-1).as_expr(), 'x%d Nprev'%(2*j)))
        items.append((xs**(2*j)*Nr(r-2).as_expr(), 'x%d Nprev2'%(2*j)))
        for t in (r-2,r-3):
            items.append((xs**(2*j)*(-W0(t)), 'x%d M%d'%(2*j,t)))
    for s in (r-1,r-2,r-3):
        for j in range(0,4):
            items.append((xs**(2*j)*T[s].as_expr()**2,'x%d F%s2'%(2*j,s)))
    N=target.degree()+1
    M=np.zeros((N,len(items))); rhs=np.zeros(N)
    for d in range(N):
        rhs[d]=int(target.coeff_monomial(xs**d))
        for i,(p,_) in enumerate(items):
            M[d,i]=int(sp.Poly(sp.expand(p),xs).coeff_monomial(xs**d))
    res=linprog(np.zeros(len(items)),A_eq=M,b_eq=rhs,bounds=(0,None),method='highs')
    print('r=%d feasible=%s'%(r,res.success),end=' ')
    if res.success:
        print([(nm,round(v,3)) for (_,nm),v in zip(items,res.x) if v>1e-8][:10])
    else: print(res.message)
