import sys, math
sys.path.insert(0,r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from fractions import Fraction
import sympy as sp
t=sp.symbols('t')
def O_rec(j,N):
    O=[None]*(N+1); O[0]=sp.Poly(1,t)
    for n in range(1,N+1):
        s=sp.Poly(0,t)
        for r in range(1,j):
            if r<=n: s=s+sp.Poly(math.comb(n,r)*O[n-r].as_expr(),t)
        for r in range(j,n+1):
            s=s+sp.Poly(t*math.comb(n,r)*O[n-r].as_expr(),t)
        O[n]=sp.Poly(sp.expand(s.as_expr()),t)
    return O
O=O_rec(2,6)
for n in range(7): print('n',n,O[n].as_expr())
