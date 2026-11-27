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

def polyderiv(coeffs):
    return [k*coeffs[k] for k in range(1,len(coeffs))]

for j in [2,3,4]:
    n = 30 if j==2 else 24
    Os = compute_O_coeffs(j,n)
    q = Os[n-j]
    rts = np.roots(list(reversed(q)))
    maxsumpos=0; worst=None
    for mu in rts:
        if abs(mu.imag)>1e-7: continue
        mu=float(mu.real)
        qp = polyval(polyderiv(q), mu).real
        recent=0.0; tail=0.0; sumpos=0.0
        for r in range(1,n+1):
            if n-r<0: continue
            val = polyval(Os[n-r], mu).real
            c = comb(n,r)
            if r<j:
                term = c*(val/qp); recent+=term
            elif r==j:
                term=0.0
            else:
                term = mu*c*(val/qp); tail+=term
                if term>0: sumpos+=term
        if sumpos>maxsumpos:
            maxsumpos=sumpos; worst=(mu,recent,tail,sumpos)
    print(f"j={j} n={n}: max positive tail sum={maxsumpos:.4e} worst={worst}")
