"""LP search for a nonnegative-coefficient positivity certificate of -W_r."""
import numpy as np
import sympy as sp
from scipy.optimize import linprog
from explore import T
xs = sp.symbols('x')


def W(r):
    F = T[r].as_expr()
    G = xs * T[r - 1].as_expr()
    return sp.expand(sp.diff(F, xs) * G - F * sp.diff(G, xs))


for r in range(5, 13):
    target = sp.Poly(-W(r), xs)
    items = []
    for j in (0, 1, 2, 3):
        items.append((xs**(2*j)*(-W(r-1)), 'x%dWprev' % (2*j)))
        items.append((xs**(2*j)*(-W(r-2)), 'x%dWprev2' % (2*j)))
    for s in (r-1, r-2, r-3, r-4):
        if s < 2:
            continue
        for j in (0, 1, 2, 3, 4):
            items.append((xs**(2*j)*T[s].as_expr()**2, 'x%dF%s2' % (2*j, s)))
    N = target.degree() + 1
    M = np.zeros((N, len(items)))
    rhs = np.zeros(N)
    for d in range(N):
        rhs[d] = int(target.coeff_monomial(xs**d))
        for i, (p, _) in enumerate(items):
            M[d, i] = int(sp.Poly(sp.expand(p), xs).coeff_monomial(xs**d))
    res = linprog(np.zeros(len(items)), A_eq=M, b_eq=rhs, bounds=(0, None), method='highs')
    print('r=%d feasible=%s' % (r, res.success), end=' ')
    if res.success:
        print([(nm, round(v, 6)) for (p, nm), v in zip(items, res.x) if v > 1e-8])
    else:
        print(res.message)
