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
for j in [2,3,4]:
    N=10
    C=series_C(j,N); D=series_inv_one_minus_A(j,N)
    g=mul_series(C,D,N)
    # H = 1/(1 - t*g): H_n = t * sum_{k=1}^n g_k H_{n-k}, H_0=1
    H=[None]*(N+1)
    H[0]=sp.Poly(1,t)
    for n in range(1,N+1):
        s=sp.Poly(0,t)
        for k in range(1,n+1):
            if k<len(g) and g[k]!=0:
                s = s + sp.Poly(t*g[k]*H[n-k].as_expr(), t)
        H[n]=sp.Poly(sp.simplify(s.as_expr()), t)
    nonreal_total=0
    for r in range(N+1):
        poly=H[r].as_expr()
        if H[r].degree()<=0:
            print('j',j,'r',r,'deg',H[r].degree()); continue
        roots=[complex(rt) for rt in sp.nroots(poly)]
        nr=[rt for rt in roots if abs(rt.imag)>1e-9]
        if nr:
            nonreal_total+=1
            print('j',j,'r',r,'NONREAL',len(nr), nr[:2])
    if nonreal_total==0:
        print('j',j,'H all real-rooted up to',N)
