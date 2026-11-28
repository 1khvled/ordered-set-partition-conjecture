"""Test eta_i < beta_i for S surrounding T (sympy exact/numeric)."""
import sympy as sp, random
t=sp.symbols('t'); random.seed(23)
bad=[]
for trial in range(600):
    m=random.randint(1,4)
    u=sorted([random.uniform(0.2,40) for _ in range(m+1)],reverse=True)
    v=sorted([random.uniform(u[i+1],u[i]) for i in range(m)],reverse=True)
    S=1
    for x in u: S*= (1+sp.Float(x)*t)
    dS=sp.diff(S,t)
    eta=sorted([float(sp.re(r)) for r in sp.nroots(dS) if abs(sp.im(r))<1e-10])
    eta=[x for x in eta if x<0]
    beta=sorted([-x for x in v])
    if len(eta)==m and not all(eta[i]<beta[i] for i in range(m)):
        bad.append((u,v,eta,beta))
print('bad',len(bad))
for b in bad[:3]: print(b)
