# transform_lib.py
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
    C=[Fraction(0)]*(N+1)
    for r in range(j,N+1):
        C[r]=Fraction(1, math.factorial(r))
    return C

def mul_series(a,b,N):
    c=[Fraction(0)]*(N+1)
    for i in range(N+1):
        for k in range(N+1-i): c[i+k]+=a[i]*b[k]
    return c

def compute_H(j,N):
    C=series_C(j,N); D=series_inv_one_minus_A(j,N)
    g=mul_series(C,D,N)
    t=sp.symbols('t')
    H=[None]*(N+1)
    H[0]=sp.Poly(1,t)
    for n in range(1,N+1):
        s=sp.Poly(0,t)
        for k in range(1,n+1):
            if k<len(g) and g[k]!=0:
                s = s + sp.Poly(t*g[k]*H[n-k].as_expr(), t)
        H[n]=sp.Poly(sp.simplify(s.as_expr()), t)
    return H,D,g

def compute_O(j,N):
    H,D,g=compute_H(j,N)
    t=sp.symbols('t')
    O=[None]*(N+1)
    for n in range(N+1):
        s=sp.Poly(0,t)
        for r in range(n+1):
            d=D[n-r]
            s = s + sp.Poly(d*sp.factorial(n)*H[r].as_expr(), t)
        O[n]=sp.Poly(sp.expand(s.as_expr()), t)
    return O,H,D

def poly_roots(poly_expr):
    if poly_expr.is_constant(): return []
    try:
        rts=sp.nroots(poly_expr)
    except Exception:
        rts=sp.solve(sp.Poly(poly_expr, sp.symbols('t')), sp.symbols('t'))
    return [complex(x) for x in rts]

def is_real_rooted(poly_expr, tol=1e-9):
    rts=poly_roots(poly_expr)
    return all(abs(z.imag)<tol for z in rts)

def interlace_check(P,Q):
    # P degree m, Q degree m or m+1? Check standard: roots of P and Q interlace if sorted roots alternate.
    # We'll just check real-rootedness of both and sign of resultant? Skip for now.
    return None
