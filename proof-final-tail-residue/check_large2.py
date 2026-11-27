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

j=2; N=80
Os = compute_O_coeffs(j,N)
bad_list=[]
for n in range(j+1, N+1):
    q = Os[n-j]
    rts = np.roots(list(reversed(q)))
    for mu in rts:
        if abs(mu.imag)>1e-7: continue
        mu=float(mu.real)
        qp = polyval(polyderiv(q), mu)
        if abs(qp)<1e-8: continue
        recent=0.0; tail=0.0
        for r in range(1,j):
            if n-r<0: continue
            val=polyval(Os[n-r], mu)
            recent += comb(n,r)*(-val/qp)
        for r in range(j+1, n+1):
            if n-r<0: continue
            val=polyval(Os[n-r], mu)
            tail += mu*comb(n,r)*(val/qp)
        margin=(recent-tail).real
        if margin<=0:
            bad_list.append((n,mu,margin,qp.real))
print("num bad:", len(bad_list))
for b in bad_list[:20]:
    print(b)
