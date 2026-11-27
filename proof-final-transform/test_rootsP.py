import sys, math
sys.path.insert(0,r'C:\Users\Abdelli\Documents\ordered-big-blocks')
import sympy as sp
z=sp.symbols('z')
for j in [2,3,4,5,6,7]:
    A=sum(z**r/math.factorial(r) for r in range(1,j))
    P=1-A
    rts=sp.nroots(sp.Poly(sp.simplify(P*z**0), z) if False else sp.Poly(sp.simplify(P), z))
    print('j',j,'degree',sp.degree(P,z))
    print('  roots:', [complex(x) for x in rts])
