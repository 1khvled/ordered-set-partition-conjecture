import sympy as sp
from math import comb

t = sp.symbols('t')

def O(nmax, j):
    arr = [sp.Integer(1)]
    for n in range(1, nmax+1):
        v = sum(comb(n,r)*arr[n-r] for r in range(1, min(j-1,n)+1))
        v += sum(comb(n,r)*t*arr[n-r] for r in range(j, n+1))
        arr.append(sp.expand(v))
    return arr

def realroots(p):
    rs = sp.nroots(p, n=20, maxsteps=200)
    return sorted([float(sp.re(r)) for r in rs if abs(sp.im(r))<1e-12])

for j in [2,3,4]:
    P = O(12, j-1)
    Q = O(12, j)
    print('=== j', j)
    for n in [3,4,5,6,7]:
        print(' n', n)
        print('  P', sp.factor(P[n]), realroots(P[n]))
        print('  Q', sp.factor(Q[n]), realroots(Q[n]))
        A = sum(comb(n,k)*P[k] for k in range(n+1))
        print('  A', sp.factor(A), realroots(A))
        # feedback term via Eq I
        Fb = (1-t)*comb(n,j-1)*sum(comb(n-j+1,a)*Q[a]*P[n-j+1-a] for a in range(n-j+2)) if n>=j-1 else 0
        print('  Fb', sp.factor(sp.expand(Fb)), realroots(sp.expand(Fb)) if Fb!=0 else [])
