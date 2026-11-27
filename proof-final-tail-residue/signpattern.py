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

for j in [3]:
    for n in [20, 40]:
        Os = compute_O_coeffs(j, n)
        q = Os[n-j]
        rts = [float(r.real) for r in np.roots(list(reversed(q))) if abs(r.imag) < 1e-7]
        for mu in [rts[0], rts[len(rts)//2], rts[-1]]:
            qp = polyval([k*Os[n-j][k] for k in range(1, len(Os[n-j]))], mu).real
            neg_ms = [m for m in range(0, n-j) if (polyval(Os[m], mu).real / qp) < 0]
            print(f"j={j} n={n} mu={mu:.3f}: num older X<0 = {len(neg_ms)}; max such m = {max(neg_ms) if neg_ms else None}")
