"""Probe xU_r roots vs G_r=xF_{r-1} roots; Wronskian signs."""
import numpy as np
import sympy as sp
from explore import T, U, MAXR

xs = sp.symbols('x')


def rr(p):
    if p.degree() == 0:
        return np.array([]), 0.0
    r = np.roots([float(v) for v in p.all_coeffs()])
    return np.sort(r.real), max(abs(r.imag))


print('pattern of roots of P_r:=xU_r (P) vs G_r:=xF_{r-1} (g)')
for r in range(4, 14):
    P = sp.Poly(sp.expand(xs * U[r].as_expr()), xs)
    G = sp.Poly(sp.expand(xs * T[r - 1].as_expr()), xs)
    rp, ep = rr(P)
    rg, eg = rr(G)
    assert ep < 1e-6 and eg < 1e-6, (r, ep, eg)
    allr = sorted([(v, 'P') for v in rp] + [(v, 'g') for v in rg])
    print(r, ''.join(s for _, s in allr), 'alt=', all(allr[i][1] != allr[i + 1][1] for i in range(len(allr) - 1)))

print()
print('Wronskian W(F_r,G_r)=F_r\' G_r - F_r G_r\' sign at sample points (expect constant):')
for r in range(4, 11):
    f = T[r].as_expr()
    g = xs * T[r - 1].as_expr()
    W = sp.Poly(sp.expand(sp.diff(f, xs) * g - f * sp.diff(g, xs)), xs)
    samples = [-100, -10, -1, -0.1, -0.001, 0, 1]
    vals = [float(W.as_expr().subs(xs, v)) for v in samples]
    print('r=%d degW=%d' % (r, W.degree()), ['%+.3g' % v for v in vals],
          'signs:', set(1 if v > 0 else (-1 if v < 0 else 0) for v in vals))

print()
print('Check: roots of F_r vs roots of P_r (should interlace by Rolle) and vs roots of G_r')
for r in range(4, 11):
    f = T[r]
    P = sp.Poly(sp.expand(xs * U[r].as_expr()), xs)
    G = sp.Poly(sp.expand(xs * T[r - 1].as_expr()), xs)
    rf, _ = rr(f)
    rp, _ = rr(P)
    rg, _ = rr(G)
    def pattern(a, b, sa, sb):
        allr = sorted([(v, sa) for v in a] + [(v, sb) for v in b])
        return ''.join(s for _, s in allr)
    print('r=%d' % r, 'F,P:', pattern(rf, rp, 'F', 'p'), ' F,G:', pattern(rf, rg, 'F', 'g'))
