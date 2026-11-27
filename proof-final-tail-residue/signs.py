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

j=3; n=20
Os = compute_O_coeffs(j,n)
q = Os[n-j]
rts = np.roots(list(reversed(q)))
mu = [float(r.real) for r in rts if abs(r.imag)<1e-7][3]  # pick middle-ish
qp = polyval([k*Os[n-j][k] for k in range(1,len(Os[n-j]))], mu).real
print(f"j={j} n={n} mu={mu:.4f} q'={qp:.4e}")
for r in range(0,n+1):
    if n-r<0: continue
    val = polyval(Os[n-r], mu).real
    X = val/qp
    print(f"  r={r:2d} n-r={n-r:2d} X={X: .4e} sign={'+' if X>0 else '-'}")
