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

j=3; n=20
Os = compute_O_coeffs(j,n)
q = Os[n-j]
rts = np.roots(list(reversed(q)))
# pick a root
for mu in rts:
    if abs(mu.imag)>1e-7: continue
    mu=float(mu.real)
    qp = polyval(polyderiv(q), mu)
    print(f"\nroot mu={mu:.6f} q'={qp.real:.6e}")
    # contributions by r (in O_n/q' decomposition)
    # term r: for r<j: C(n,r) O_{n-r}(mu)/q'
    # for r>=j+1: mu*C(n,r) O_{n-r}(mu)/q'  (r=j zero)
    contribs=[]
    for r in range(1,n+1):
        if n-r<0: continue
        val = polyval(Os[n-r], mu)
        c = comb(n,r)
        if r < j:
            term = c*(val/qp)
        elif r == j:
            term = 0.0+0.0j  # mu*C*q(mu)/qp but q(mu)=0
        else:
            term = mu*c*(val/qp)
        contribs.append((r, term.real))
    # print sorted by magnitude
    for r,t in contribs:
        print(f"  r={r:2d} contrib={t: .6e}")
    # cumulative
    cum=0.0
    print("  cumulative:")
    for r,t in contribs:
        cum+=t
        print(f"    after r={r:2d}: {cum:.6e}")
    break
