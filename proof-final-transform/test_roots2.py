import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
j=3; N=10
O,H,D=tl.compute_O(j,N)
def roots(expr):
    if sp.degree(expr,t)<=0: return []
    return sorted([complex(x) for x in sp.nroots(expr)], key=lambda z:z.real)
for n in [6,7,8]:
    print('n',n)
    print('  H', [round(z.real,4) for z in roots(H[n].as_expr())])
    for r in range(1,j):
        print('  O_{n-r}'.format(), n-r, [round(z.real,4) for z in roots(O[n-r].as_expr())])
