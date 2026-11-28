"""Test: if S surrounds T (negative simple roots), is W(E S, U T)>0 on R?"""
import sympy as sp, random
t=sp.symbols('t'); random.seed(7)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
bad=[]
for trial in range(120):
    m=random.randint(2,4)
    u=sorted(random.sample(range(1,40),m+1), reverse=True)  # S roots -u_i, u descending
    S=1
    for x in u: S*= (1+x*t)
    # T roots -v_j with u_j > v_j > u_{j+1}  (S surrounds T)
    v=[]
    for j in range(m):
        lo,hi=u[j+1],u[j]
        vj=random.randint(lo+1,hi-1) if hi-lo>1 else None
        if vj is None: break
        v.append(vj)
    if len(v)<m: continue
    T=1
    for x in v: T*=(1+x*t)
    assert sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo)==0
    H=sp.Poly(W(E(S),U(T)),t)
    if H.count_roots(-sp.oo,sp.oo): bad.append((u,v,sp.factor(H.as_expr())))
print('failures',len(bad))
for b in bad[:5]: print(b[0],b[1],b[2])
