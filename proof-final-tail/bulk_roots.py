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

def roots_for(j,n):
    Os = compute_O_coeffs(j,n)
    On = Os[n]
    m = n//j; s=n%j
    cnm = On[m]
    coeffs_y = [0]*(m+1)
    for ell in range(m+1):
        coeffs_y[ell] = (On[m-ell]/cnm) * ((-1)**ell)
    roots = np.roots(coeffs_y[::-1])
    return roots,m,s

if __name__=="__main__":
    for j in [3]:
        for n in [30,60,120]:
            roots,m,s = roots_for(j,n)
            nonreal = sum(1 for r in roots if abs(r.imag)>1e-7)
            rr = sorted(roots, key=lambda r:r.real)
            print(f"j={j} n={n} m={m} s={s} nonreal={nonreal} max_imag={max(abs(r.imag) for r in roots):.2e}")
            print("  real min/max:", min(r.real for r in roots), max(r.real for r in roots))
            # histogram of real roots counts positive/negative
            pos = sum(1 for r in roots if r.real>0 and abs(r.imag)<1e-7)
            neg = sum(1 for r in roots if r.real<0 and abs(r.imag)<1e-7)
            print("  pos/neg real:", pos, neg)
