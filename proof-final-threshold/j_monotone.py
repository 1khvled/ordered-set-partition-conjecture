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

def roots(p):
    p = sp.Poly(p, t)
    if p.degree() == 0: return []
    rs = sp.nroots(p, n=25, maxsteps=500)
    out = []
    for r in rs:
        if abs(sp.im(r)) > 1e-15:
            out.append('C')
        else:
            out.append(float(sp.re(r)))
    return out

N=14
for n in range(2, N+1):
    print('n=', n)
    prev = None
    for j in range(2, 8):
        p = O(n,j)[n]
        rs = roots(p)
        if any(isinstance(r,str) for r in rs):
            print(' j',j,'COMPLEX',rs); continue
        rs = sorted(rs)
        print(' j',j, ' '.join('%.4f'%r for r in rs))
        if prev is not None and rs and prev:
            # compare: all roots of j below all roots of j-1?
            print('   sep?', max(rs) < min(prev))
        prev = rs
