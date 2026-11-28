import sympy as sp
from ops import A, Q, x
for n in [4,5,6]:
    for av in [1,2]:
        q=Q(n,av,x); qm=Q(n-1,av,x); qv=Q(n-2,av+1,x)
        resid=sp.expand(q-(sp.diff(x*qm,x)+(n-1)*av*x*qv))
        print(n,av,'q=',q,'qm=',qm,'qv=',qv,'resid=',resid)
