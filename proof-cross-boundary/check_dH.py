import sympy as sp
t=sp.symbols('t')
def rows(j,nmax):
    O=[sp.Integer(1)]
    for n in range(1,nmax+1):
        s=sp.Integer(0)
        for r in range(1,min(j-1,n)+1): s+=sp.binomial(n,r)*O[n-r]
        for r in range(j,n+1): s+=t*sp.binomial(n,r)*O[n-r]
        O.append(sp.expand(s))
    return O
def W(p,q): return sp.expand(sp.diff(p,t)*q-p*sp.diff(q,t))
def Hpoly(O,j,n):
    return sp.expand(t*sum(sp.binomial(n,k)*O[k] for k in range(0,n-j+1)))

j=3; nmax=15
O=rows(j,nmax)
# precompute d_m = m![z^m]1/(1-A_j)
# A_j = z/1!+z^2/2! for j=3
from math import factorial
d = [0]*(nmax+1)
# 1/(1-A) = sum A^r. compute coefficients.
# A coeffs: a1=1, a2=1/2
Acoeff = [sp.Rational(0)]*(nmax+1)
Acoeff[1]=sp.Rational(1,1)
Acoeff[2]=sp.Rational(1,2)
# power series 1/(1-A)
invA = [sp.Rational(0)]*(nmax+1)
invA[0]=sp.Rational(1)
for m in range(1,nmax+1):
    # invA[m] = sum_{r=1}^{m} Acoeff[r]*invA[m-r]
    s=sp.Rational(0)
    for r in range(1,m+1):
        s += Acoeff[r]*invA[m-r]
    invA[m]=s
for m in range(nmax+1):
    d[m] = sp.factorial(m)*invA[m]

for n in range(j, nmax+1):
    q=O[n-j]
    Hn=Hpoly(O,j,n)
    # expand O_{n-j} in d*H basis
    rhs = sp.Integer(0)
    for k in range(0,n-j+1):
        Hk = sp.Integer(1) if k==0 else Hpoly(O,j,k)
        rhs += sp.binomial(n-j,k)*d[n-j-k]*Hk
    rhs=sp.expand(rhs)
    print(n, sp.expand(rhs-q)==0)
