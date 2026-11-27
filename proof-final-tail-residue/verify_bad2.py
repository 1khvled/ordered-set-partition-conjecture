import mpmath
mpmath.mp.dps = 200
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

def polyval_mpf(coeffs, x):
    res = mpmath.mpf(0)
    for k in range(len(coeffs)-1, -1, -1):
        res = res*x + mpmath.mpf(coeffs[k])
    return res

j=2; n=63
Os = compute_O_coeffs(j,n)
q = Os[n-j]
deriv_coeffs = [k*Os[n-j][k] for k in range(1,len(Os[n-j]))]

def f(x): return polyval_mpf(q, x)

mu0 = mpmath.mpf("-0.2979493543628123")
try:
    mu = mpmath.findroot(f, (mpmath.mpf("-0.30"), mpmath.mpf("-0.29")))
    print("mu =", mu)
    qp = polyval_mpf(deriv_coeffs, mu)
    print("qp =", qp)
    recent = mpmath.mpf(0)
    for r in range(1,j):
        if n-r<0: continue
        val = polyval_mpf(Os[n-r], mu)
        recent += comb(n,r)*(-val/qp)
    tail = mpmath.mpf(0)
    for r in range(j+1, n+1):
        if n-r<0: continue
        val = polyval_mpf(Os[n-r], mu)
        tail += mu*comb(n,r)*(val/qp)
    margin = recent - tail
    print("recent =", recent)
    print("tail =", tail)
    print("margin =", margin)
except Exception as e:
    print("ERR", e)
