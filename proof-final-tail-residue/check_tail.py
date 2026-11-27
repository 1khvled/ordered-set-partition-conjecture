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

def check(j, nmax):
    Os = compute_O_coeffs(j, nmax)
    rows = []
    for n in range(j+1, nmax+1):
        q = Os[n-j]
        mdeg = len(q)-1
        coeffs_high = list(reversed(q))
        rts = np.roots(coeffs_high)
        for mu in rts:
            if abs(mu.imag) > 1e-7:
                continue
            mu = float(mu.real)
            qp = polyval(polyderiv(q), mu)
            if abs(qp) < 1e-12:
                continue
            recent = 0.0
            for r in range(1, j):
                if n-r < 0: continue
                Onr = Os[n-r]
                val = polyval(Onr, mu)
                recent += comb(n,r)*(-val/qp)
            tail = 0.0
            for r in range(j+1, n+1):
                if n-r < 0: continue
                Onr = Os[n-r]
                val = polyval(Onr, mu)
                tail += comb(n,r)*(val/qp)
            tail = mu*tail
            margin = (recent - tail).real
            rows.append((j,n,mu,margin))
    bad = [x for x in rows if x[3] <= 0]
    print(f"j={j} n<= {nmax}: total real-root evals={len(rows)} bad={len(bad)}")
    if bad:
        print("  first bad:", bad[:5])
    return rows,bad

if __name__=="__main__":
    for j in [2,3,4,5]:
        for N in [20,40]:
            check(j,N)