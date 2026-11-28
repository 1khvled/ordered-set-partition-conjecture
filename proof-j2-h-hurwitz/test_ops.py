"""Test operator-lemma candidates on random real-rooted K (evidence)."""
import sympy as sp, itertools
t=sp.symbols('t')
def E(f): return sp.expand(f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
tests=[]
for roots in itertools.product([1,2,3,5,8], repeat=3):
    K=1
    for r in roots: K*= (1+r*t)
    H=sp.Poly(W(t*K,E(E(K))),t)
    rr=H.count_roots(-sp.oo,sp.oo)
    if rr: tests.append((roots,rr))
print('W(tK,E^2K) real roots count failures', len(tests), tests[:10])
# also W(ES, U T) when S,T interlace with known orientation: take S=(1+t)(1+3t)(1+5t), T=(1+2t)(1+4t) (S surrounds T)
S=(1+t)*(1+3*t)*(1+5*t); T=(1+2*t)*(1+4*t)
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
print('S,T interlace:', sp.Poly(W(S,T),t).count_roots(-sp.oo,sp.oo))
print('W(ES,UT) roots', sp.Poly(W(E(S),U(T)),t).count_roots(-sp.oo,sp.oo), sp.factor(W(E(S),U(T))))
