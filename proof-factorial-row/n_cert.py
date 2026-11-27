"""LP certificate for N = -W(U F_{r-2}, x U F_{r-3})/4 as positive combination."""
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from explore import T
xs = sp.symbols('x')
Uop = lambda f: sp.expand(2*f + xs*sp.diff(f, xs))


def W(p, R):
    return sp.expand(sp.diff(p, xs)*R - p*sp.diff(R, xs))


def W0(r):
    return W(T[r].as_expr(), xs*T[r-1].as_expr())


for r in range(6, 13):
    p, q = Uop(T[r-2].as_expr()), Uop(T[r-3].as_expr())
    Npoly = sp.Poly(sp.expand(-W(p, xs*q)/4), xs)
    items = []
    for j in range(0, 4):
        items.append((xs**(2*j)*(-W0(r-1)), 'x%d W0_%d' % (2*j, r-1)))
        items.append((xs**(2*j)*(-W0(r-2)), 'x%d W0_%d' % (2*j, r-2)))
    for s in (r-1, r-2, r-3):
        for j in range(0, 5):
            items.append((xs**(2*j)*T[s].as_expr()**2, 'x%d F%s2' % (2*j, s)))
    N = Npoly.degree()+1
    M = np.zeros((N, len(items))); rhs = np.zeros(N)
    for d in range(N):
        rhs[d] = int(Npoly.coeff_monomial(xs**d))
        for i,(p2,_) in enumerate(items):
            M[d,i] = int(sp.Poly(sp.expand(p2), xs).coeff_monomial(xs**d))
    res = linprog(np.zeros(len(items)), A_eq=M, b_eq=rhs, bounds=(0,None), method='highs')
    print('r=%d feasible=%s' % (r,res.success), end=' ')
    if res.success:
        print([(nm,round(v,4)) for (_,nm),v in zip(items,res.x) if v>1e-8][:8])
    else:
        print(res.message)
