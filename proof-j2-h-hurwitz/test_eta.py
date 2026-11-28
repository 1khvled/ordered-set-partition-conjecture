"""Test: S surrounds T => each critical point eta_i of S lies left of beta_i (T root)."""
import numpy as np, random
random.seed(23)
bad=0
for trial in range(3000):
    m=random.randint(1,5)
    u=sorted([random.uniform(0.1,50) for _ in range(m+1)],reverse=True)
    v=sorted([random.uniform(u[i+1],u[i]) for i in range(m)],reverse=True)
    # S roots -u_i, S' roots via np.roots of derivative of prod(1+u t)
    import numpy.polynomial.polynomial as P
    cp=np.array([1.0]);
    for x in u: cp=P.polymul(cp,np.array([1.0,x]))
    dc=P.polyder(cp)
    eta=np.sort(np.roots(list(reversed(dc.tolist()))).real)
    eta=eta[eta<0]
    if len(eta)!=m: continue
    beta=np.sort(-np.array(v))  # T roots -v_j ascending since v descending
    if not all(eta[i] < beta[i] for i in range(m)):
        bad+=1
        if bad<4: print('counterexample',u,v,eta,beta)
print('bad',bad)
