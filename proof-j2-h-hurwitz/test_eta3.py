"""Correct test: eta_i (critical pts of S) vs beta_i (roots of T), sorted ascending."""
import sympy as sp, random
t=sp.symbols('t'); random.seed(31)
bad=[]
for trial in range(500):
    m=random.randint(1,5)
    u=sorted([random.uniform(0.2,50) for _ in range(m+1)],reverse=True)
    v=[random.uniform(u[i+1],u[i]) for i in range(m)]
    S=1
    for x in u: S*= (1+sp.Float(x)*t)
    T=1
    for x in v: T*= (1+sp.Float(x)*t)
    eta=sorted([float(sp.re(r)) for r in sp.nroots(sp.diff(S,t)) if abs(sp.im(r))<1e-10 and sp.re(r)<0])
    beta=sorted([float(sp.re(r)) for r in sp.nroots(T) if abs(sp.im(r))<1e-10 and sp.re(r)<0])
    rho=sorted([-1.0/x for x in u])
    # verify interlacing rho_i<beta_i<rho_{i+1}
    if not all(rho[i]<beta[i]<rho[i+1] for i in range(m)): 
        continue
    if not all(eta[i]<beta[i] for i in range(m)):
        bad.append((rho,beta,eta))
print('bad',len(bad))
for b in bad[:3]: print(b)
