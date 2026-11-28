"""Stress-test operator lemma: S surrounds T => W(ES, U T)>0 on R."""
import sympy as sp, random
t=sp.symbols('t'); random.seed(19)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
bad=[]
for trial in range(300):
    m=random.randint(1,5)
    # random decreasing positive u_i (S roots -u_i, u descending)
    u=sorted([random.uniform(0.2,50) for _ in range(m+1)], reverse=True)
    S=1
    for x in u: S*=(1+x*t)
    v=sorted([random.uniform(u[i+1],u[i]) for i in range(m)], reverse=True)
    T=1
    for x in v: T*=(1+x*t)
    assert sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo)==0
    H=sp.Poly(W(E(S),U(T)),t)
    if H.count_roots(-sp.oo,sp.oo) or sp.LC(H,t)<=0 or H.eval(0)<=0:
        bad.append((u,v,sp.factor(H.as_expr())))
print('stress failures',len(bad))
for b in bad[:3]: print(b)
