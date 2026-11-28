"""Exploration: structural identities and Hurwitz check (evidence, not proof)."""
import sympy as sp
from common import Family, t
E, T, W = Family.E, Family.T, Family.W

F = Family(40)
print('== recurrence identities ==')
for n in range(4, 16):
    assert F.Qpoly(n) == E(F.Qpoly(n-1)) + (n-1)*sp.Poly(t,t)*(F.Qpoly(n-2)+E(F.Qpoly(n-2)))
    for a in (1,2,3):
        P = F.Qpoly(n,a)
        lhs = sp.Poly(sp.expand(a*F.Q(n,a+1)),t)
        rhs = sp.Poly(sp.expand((a+1)*F.Q(n,a) + t*sp.diff(F.Q(n,a),t)),t)
        assert lhs == rhs
print('P recurrence and (V) hold exactly')

print('== H samples ==')
for n in (6,7,8,10,12):
    H = F.H(n)
    print(n, 'deg', H.degree(), 'coeffs', H.all_coeffs())

print('== exact Routh-Hurwitz test (small n) ==')
def routh_signs(poly):
    c = list(reversed(poly.all_coeffs()))  # c[0]=constant
    d = poly.degree()
    rows = []
    row0 = [c[j] for j in range(0, d+1, 2)]
    row1 = [c[j] for j in range(1, d+1, 2)]
    rows.append(row0); rows.append(row1)
    while len(row1) > 1:
        prev = rows[-2]; cur = rows[-1]
        nxt = [sp.simplify(-(prev[0]*cur[i+1]-cur[0]*prev[i+1])/cur[0]) for i in range(len(cur)-1)]
        if all(x == 0 for x in nxt):
            break
        rows.append(nxt)
        if len(nxt) <= 1: break
    return [r[0] for r in rows]

for n in range(6, 19):
    H = F.H(n)
    rs = routh_signs(H)
    ok = all(sp.simplify(x) > 0 for x in rs)
    print('n=',n, 'first-column signs ok', ok, 'len', len(rs))

print('== numeric roots of H (evidence) ==')
try:
    import numpy as np
    for n in (10, 14, 20, 30):
        H = F.H(n)
        cs = [float(c) for c in H.all_coeffs()]
        roots = np.roots(cs)
        print('n=',n,'max Re', max(r.real for r in roots), 'min |Im| among near-real', min(abs(r.imag) for r in roots))
except Exception as ex:
    print('numpy unavailable', ex)
