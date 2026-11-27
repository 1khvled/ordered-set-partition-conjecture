import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
j=3; N=10
O,H,D=tl.compute_O(j,N)
def roots(expr):
    if sp.degree(expr,t)<=0: return []
    return sorted([complex(x) for x in sp.nroots(expr)], key=lambda z:(z.real,z.imag))
for n in range(j,N+1):
    hr=roots(H[n].as_expr())
    for r in range(1,j):
        orr=roots(O[n-r].as_expr())
        print('j',j,'n',n,'r',r,'degH',len(hr),'degO',len(orr))
