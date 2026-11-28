"""Test whether 4h0h2-h1^2>0 holds generically under interlacing hypotheses."""
import sympy as sp, random
t=sp.symbols('t'); random.seed(3)
def E(f): return sp.expand(f+t*sp.diff(f,t))
def U(f): return sp.expand(2*f+t*sp.diff(f,t))
def W(f,g): return sp.expand(sp.diff(f,t)*g-f*sp.diff(g,t))
bad=[]; good=0
for trial in range(200):
    m=random.randint(1,3)
    us=sorted(random.sample(range(1,30),m+1),reverse=True)  # S roots
    S=1
    for x in us: S*=(1+x*t)
    # T strictly interlacing S with S surrounding: roots between us
    vs=[random.randint(us[i+1]+1,us[i]-1) for i in range(m) if us[i]-us[i+1]>1]
    if len(vs)<m: continue
    T=1
    for x in vs: T*=(1+x*t)
    # K negative-rooted, degree m (or m-1): choose roots
    ks=sorted(random.sample(range(1,30),m),reverse=True)
    K=1
    for x in ks: K*=(1+x*t)
    A=E(S); B=t*K
    alpha=U(T); beta=E(E(K))
    h0=W(A,alpha); h2=W(B,beta); h1=W(A,beta)+W(B,alpha)
    D=sp.Poly(sp.expand(4*h0*h2-h1**2),t)
    if D.count_roots(-sp.oo,sp.oo) or sp.LC(D,t)<0:
        bad.append((us,vs,ks,D.count_roots(-sp.oo,sp.oo),sp.factor(D.as_expr())))
    else: good+=1
print('good',good,'bad',len(bad))
for b in bad[:4]: print(b[:4],b[4])
