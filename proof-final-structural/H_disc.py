import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction
import sympy as sp

def series_inv_one_minus_A(j,N):
    A=[Fraction(0)]*(N+1)
    for r in range(1,j):
        if r<=N: A[r]=Fraction(1, math.factorial(r))
    D=[Fraction(0)]*(N+1); D[0]=Fraction(1)
    for n in range(1,N+1):
        s=Fraction(0)
        for r in range(1,min(j,n+1)):
            if r<=j-1: s+=A[r]*D[n-r]
        D[n]=s
    return D
def series_C(j,N):
    exp=[Fraction(1,math.factorial(r)) for r in range(N+1)]
    C=[Fraction(0)]*(N+1)
    for r in range(N+1):
        if r>=j: C[r]=exp[r]
    return C
def mul_series(a,b,N):
    c=[Fraction(0)]*(N+1)
    for i in range(N+1):
        for k in range(N+1-i): c[i+k]+=a[i]*b[k]
    return c

t=sp.symbols('t')
for j in [2,3,4,5]:
    N=14
    C=series_C(j,N); D=series_inv_one_minus_A(j,N)
    g=mul_series(C,D,N)
    H=[None]*(N+1); H[0]=sp.Poly(1,t)
    for n in range(1,N+1):
        s=sp.Poly(0,t)
        for k in range(1,n+1):
            if k<len(g) and g[k]!=0:
                s = s + sp.Poly(t*g[k]*H[n-k].as_expr(), t)
        H[n]=sp.Poly(sp.simplify(s.as_expr()), t)
    # check discriminants for degree>=2
    bad=[]
    for r in range(N+1):
        p=H[r].as_expr()
        deg=sp.degree(p,t)
        if deg>=2:
            disc=sp.discriminant(sp.Poly(p,t), t)
            if disc<0:
                bad.append((r, disc))
    print('j',j,'H discriminant failures:', bad)
