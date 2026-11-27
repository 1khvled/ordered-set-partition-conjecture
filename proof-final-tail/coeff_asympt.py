import mpmath
from math import comb, factorial
mpmath.mp.dps = 50

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

def D_j(j):
    from fractions import Fraction
    return (j) * (Fraction(j+2,j+1))**j  # j!?? wait D_j = j! * ((j+2)/(j+1))^j
# Actually D_j = factorial(j)*((j+2)/(j+1))**j
from math import factorial as fact
def Dval(j):
    return fact(j)*( (j+2)/(j+1) )**j

if __name__=="__main__":
    j=3
    for n in [30,60]:
        m=n//j; s=n%j
        Os=compute_O_coeffs(j,n)
        On=Os[n]
        cnm=On[m]
        Dm = Dval(j)*(m**j)
        print(f"j={j} n={n} m={m} s={s}")
        for ell in [0,1,2,m//4,m//2,3*m//4,m-1,m]:
            if ell>m: continue
            B = On[m-ell]/cnm
            asym = (Dm**ell) * (fact(s)/fact(s+j*ell))  # s!/(s+j ell)!
            rel = B/asym
            print(f"  ell={ell:3d} rel_err={float(rel):.6e}")
