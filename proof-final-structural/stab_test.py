import sys, math
sys.path.insert(0, r'C:\Users\Abdelli\Documents\ordered-big-blocks')
from sympy import Poly, Symbol, nroots, series, exp, factorial
z=Symbol('z'); t=Symbol('t')
# Test if B_j truncated has all real roots for t=1 (all blocks marked) -> should be e^z-1? no.
for j in [2,3]:
    for N in [8,12]:
        # B = A + t*C truncated to degree N
        B=0
        for r in range(1,j):
            B += z**r/math.factorial(r)
        for r in range(j,N+1):
            B += t*z**r/math.factorial(r)
        for tval in [0.5, 1.0, 2.0]:
            p = Poly(B.subs(t,tval), z)
            roots = nroots(p.as_expr())
            nonreal = [rt for rt in roots if abs(complex(rt).imag)>1e-9]
            print('j',j,'N',N,'t',tval,'deg',p.degree(),'nonreal',len(nonreal))
