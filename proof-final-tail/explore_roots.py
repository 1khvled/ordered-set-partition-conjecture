import numpy as np
from fractions import Fraction

def binom(n,r):
    if r<0 or r>n: return 0
    # exact int
    from math import comb
    return comb(n,r)

def compute_O_coeffs(j, nmax):
    # O_n as list of ints (coeffs in t), length floor(n/j)+1
    # O_0 = [1]
    Os = [None]*(nmax+1)
    Os[0] = [1]
    for n in range(1, nmax+1):
        deg = n//j
        coeff = [0]*(deg+1)
        # small blocks r=1..j-1: weight 1
        for r in range(1, min(j, n+1)):
            if r >= j: break
            if n-r < 0: continue
            Onr = Os[n-r]
            c = binom(n,r)
            for k,val in enumerate(Onr):
                if k <= deg:
                    coeff[k] += c*val
        # large blocks r=j..n: weight t -> shift k+1
        for r in range(j, n+1):
            if n-r < 0: continue
            Onr = Os[n-r]
            c = binom(n,r)
            for k,val in enumerate(Onr):
                if k+1 <= deg:
                    coeff[k+1] += c*val
        Os[n] = coeff
    return Os

def reversed_Q_roots(j, n):
    m = n//j
    s = n % j
    Os = compute_O_coeffs(j, n)
    On = Os[n]
    deg = len(On)-1  # floor(n/j)=m
    # coefficients in t: On[k] for k=0..m
    # reversed: coefficient for ell (0..m): On[m-ell]
    # normalize by On[m]
    cnm = On[m]
    # polynomial in y: sum_{ell=0}^m (On[m-ell]/cnm) * (-y)^ell
    # = sum_{ell} a_ell y^ell with a_ell = (On[m-ell]/cnm)*(-1)^ell
    coeffs_y = [0]*(m+1)
    for ell in range(m+1):
        coeffs_y[ell] = (On[m-ell]/cnm) * ((-1)**ell)
    # numpy.roots expects highest degree first
    roots = np.roots(coeffs_y[::-1])
    return roots, m, s, On

if __name__ == "__main__":
    for j in [3,4]:
        for n in [12,24,48]:
            roots,m,s,On = reversed_Q_roots(j,n)
            nonreal = sum(1 for r in roots if abs(r.imag)>1e-9)
            print(f"j={j} n={n} m={m} s={s} deg={m} nonreal={nonreal} max_imag={max(abs(r.imag) for r in roots):.3e}")
            # print few smallest abs roots
            sr = sorted(roots, key=lambda r:(abs(r),))
            print("  smallest 5 abs roots:", [(round(r.real,4),round(r.imag,4)) for r in sr[:5]])
