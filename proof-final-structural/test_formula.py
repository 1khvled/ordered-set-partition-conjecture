import sys, math
sys.path.insert(0,r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction
import sympy as sp
exec(open('H_check.py').read().split('t=sp.symbols')[0])
t=sp.symbols('t')
j=2; N=6
C=series_C(j,N); D=series_inv_one_minus_A(j,N)
g=mul_series(C,D,N)
H=[None]*(N+1)
H[0]=sp.Poly(1,t)
for n in range(1,N+1):
    s=sp.Poly(0,t)
    for k in range(1,n+1):
        if k<len(g) and g[k]!=0:
            s = s + sp.Poly(t*g[k]*H[n-k].as_expr(), t)
    H[n]=sp.Poly(sp.simplify(s.as_expr()), t)
# O_n = n! sum_r d_{n-r} H_r / r!  (correct)
for n in range(N+1):
    s=sp.Poly(0,t)
    for r in range(n+1):
        d=D[n-r]
        s = s + sp.Poly((d*sp.factorial(n)/sp.factorial(r))*H[r].as_expr(), t)
    print('n',n,'O_formula1=',sp.expand(s.as_expr()))
# structural formula: n! sum d H_r
for n in range(N+1):
    s=sp.Poly(0,t)
    for r in range(n+1):
        d=D[n-r]
        s = s + sp.Poly(d*sp.factorial(n)*H[r].as_expr(), t)
    print('n',n,'O_formula2=',sp.expand(s.as_expr()))
