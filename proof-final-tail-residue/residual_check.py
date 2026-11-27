import numpy as np
from math import comb

def compute_O_coeffs(j, nmax):
    Os = [None]*(nmax+1)
    Os[0] = [1]
    for n in range(1, nmax+1):
        deg = n//j
        coeff = [0]*(deg+1)
        for r in range(1, min(j, n+1)):
            if n-r < 0: continue
            Onr = Os[n-r]
            c = comb(n,r)
            for k,val in enumerate(Onr):
                if k <= deg:
                    coeff[k] += c*val
        for r in range(j, n+1):
            if n-r < 0: continue
            Onr = Os[n-r]
            c = comb(n,r)
            for k,val in enumerate(Onr):
                if k+1 <= deg:
                    coeff[k+1] += c*val
        Os[n] = coeff
    return Os

def polyval(coeffs, x):
    res = 0.0+0.0j
    for k in range(len(coeffs)-1, -1, -1):
        res = res*x + coeffs[k]
    return res

j=2; n=63
Os = compute_O_coeffs(j,n)
q = Os[n-j]
scale = sum(abs(c) for c in q)
rts = np.roots(list(reversed(q)))
good=0; bad=0
for mu in rts:
    if abs(mu.imag)>1e-7: 
        continue
    mu=float(mu.real)
    qmu = polyval(q, mu)
    rel = abs(qmu)/scale
    if rel < 1e-8:
        good+=1
    else:
        bad+=1
        if bad<=5:
            print(f"spurious mu={mu:.6f} relres={rel:.3e}")
print(f"good real roots={good} spurious={bad} total degree={len(q)-1}")
