import transform_lib as tl
import sympy as sp
t=sp.symbols('t')
j=2; N=4
O,H,D=tl.compute_O(j,N)
for n in range(N+1):
    print('n',n,'O',O[n].as_expr(),'H',H[n].as_expr())
print('check n=2:', sp.expand(O[2].as_expr() - (H[2].as_expr() + 2*O[1].as_expr())))
