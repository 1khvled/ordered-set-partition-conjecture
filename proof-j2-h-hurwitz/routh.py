"""Factor Routh first-column entries of H_n."""
import sympy as sp
from common import Family, t
F=Family(14)
def routh(poly):
    c=list(reversed(poly.all_coeffs())); d=poly.degree()
    rows=[[c[j] for j in range(0,d+1,2)],[c[j] for j in range(1,d+1,2)]]
    while len(rows[-1])>1:
        p=rows[-2]; q=rows[-1]
        rows.append([sp.cancel(-(p[0]*q[i+1]-q[0]*p[i+1])/q[0]) for i in range(len(q)-1)])
    return rows
for n in (8,10,12):
    print('=== n',n)
    rows=routh(F.H(n))
    for i,r in enumerate(rows):
        print(' row',i,'first',sp.factor(r[0]))
