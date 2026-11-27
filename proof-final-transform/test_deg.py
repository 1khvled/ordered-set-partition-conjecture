import transform_lib as tl
import sympy as sp
for j in [2,3,4,5]:
    N=15
    O,H,D=tl.compute_O(j,N)
    print('j',j)
    print('H degs:', [sp.degree(H[r].as_expr(),sp.symbols('t')) for r in range(N+1)])
    print('O degs:', [sp.degree(O[n].as_expr(),sp.symbols('t')) for n in range(N+1)])
