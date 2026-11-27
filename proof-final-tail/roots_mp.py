import mpmath
mpmath.mp.dps = 80
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

def roots_mp(j,n):
    Os = compute_O_coeffs(j,n)
    On = Os[n]
    m = n//j
    cnm = On[m]
    # coeffs highest degree first: for ell=m..0: On[m-ell]/cnm * (-1)^ell
    coeffs = []
    for ell in range(m, -1, -1):
        val = mpmath.mpf(On[m-ell]) / mpmath.mpf(cnm) * ((-1)**ell)
        coeffs.append(val)
    roots = mpmath.polyroots(coeffs)
    nonreal = sum(1 for r in roots if abs(mpmath.im(r))>1e-30)
    maximag = max(abs(mpmath.im(r)) for r in roots)
    return roots, m, nonreal, maximag

if __name__=="__main__":
    for n in [30,60,120]:
        roots,m,nonreal,maximag = roots_mp(3,n)
        print(f"j=3 n={n} m={m} nonreal={nonreal} max_imag={maximag}")
