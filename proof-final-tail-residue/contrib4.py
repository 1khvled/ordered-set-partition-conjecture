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

worst_ratio=0; worst=None
for j in [2,3,4,5]:
    for n in [15,25,35]:
        Os = compute_O_coeffs(j,n)
        q = Os[n-j]
        rts = np.roots(list(reversed(q)))
        for mu in rts:
            if abs(mu.imag)>1e-7: continue
            mu=float(mu.real)
            qp = polyval(polyderiv(q), mu).real
            recent=0.0; sumpos=0.0
            for r in range(1,n+1):
                if n-r<0: continue
                val = polyval(Os[n-r], mu).real
                c = comb(n,r)
                if r<j:
                    recent += c*(val/qp)  # recent is negative sum; store signed
                elif r==j:
                    pass
                else:
                    term = mu*c*(val/qp)
                    if term>0: sumpos+=term
            # recent signed is negative; magnitude:
            rec_mag = -recent if recent<0 else recent
            ratio = sumpos/rec_mag if rec_mag>0 else float('inf')
            if ratio>worst_ratio:
                worst_ratio=ratio; worst=(j,n,mu,rec_mag,sumpos)
print("worst ratio sumpos/recent_mag:", worst_ratio)
print("worst:", worst)
