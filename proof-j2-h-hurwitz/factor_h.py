"""Factor H_n and inspect structure."""
import sympy as sp
from common import Family
F=Family(16)
for n in range(6,13):
    H=F.H(n).as_expr()
    print('n=',n,'factor=',sp.factor(H))
