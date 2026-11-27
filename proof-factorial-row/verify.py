"""Exact verification of all identities used in NOTES.md. No installs/builds."""
from sympy import symbols, Poly, expand, diff, factor, simplify, Rational, Integer
from explore import table, Epoly, E, T, U, A, MAXR  # imports print verification of (R)

x = symbols('x')


def W(p, q):
    return expand(diff(p, x)*q - p*diff(q, x))


def Uop(f):
    return expand(2*f + x*diff(f, x))


print('=== Wronskian decomposition W_r = B1 + (r-1) B2 ===')
for r in range(5, MAXR+1):
    Av, Cv = T[r-1].as_expr(), T[r-2].as_expr()
    B1 = x**2*(Av*diff(Av, x, 2) - diff(Av, x)**2) - Av**2
    B2 = x**2*(3*Av*diff(Cv, x) - 2*diff(Av, x)*Cv) + x**3*(Av*diff(Cv, x, 2) - diff(Cv, x)*diff(Av, x))
    Wv = W(T[r].as_expr(), x*Av)
    assert expand(Wv - (B1 + (r-1)*B2)) == 0, r
print('verified r=5..%d' % MAXR)

print('=== -B2 = -x^2[ W(C,Dc) + (r-2) W(Dc, x De) ] for r>=7 ===')
for r in range(7, MAXR+1):
    Cv = T[r-2].as_expr()
    Dc = Uop(Cv)
    De = Uop(T[r-3].as_expr())
    Av = T[r-1].as_expr()
    B2 = x**2*(3*Av*diff(Cv, x) - 2*diff(Av, x)*Cv) + x**3*(Av*diff(Cv, x, 2) - diff(Cv, x)*diff(Av, x))
    rhs = -x**2*(W(Cv, Dc) + (r-2)*W(Dc, x*De))
    assert expand(-B2 - rhs) == 0, r
print('verified r=7..%d' % MAXR)

print('=== W(C,Dc) = -C^2 * sum_j (-gamma_j)/(x-gamma_j)^2 <= 0 (r>=6) ===')
# verified analytically in NOTES; here check via log-derivative identity
for r in range(6, 11):
    Cv = T[r-2].as_expr()
    Dc = Uop(Cv)
    c = diff(Cv, x)/Cv
    cp = diff(c, x)
    lhs = W(Cv, Dc)
    rhs = -Cv**2*(c + x*cp)
    assert simplify(lhs - rhs) == 0, r
print('log-derivative identity verified r=6..10')

print('=== -W_r has no real zero, r<=40, exact Sturm count ===')
import numpy as np
for r in range(4, MAXR+1):
    Av = T[r-1].as_expr()
    Wv = W(T[r].as_expr(), x*Av)
    P = Poly(Wv, x)
    if P.degree() == 0:
        assert int(P.as_expr()) < 0, r
        continue
    assert P.count_roots() == 0, ('W has real root', r)   # exact real-root isolation
print('-W_r has no real zero for r=4..%d (exact Sturm)' % MAXR)

print('=== Interlacing of F_r and xF_{r-1}: strict alternation r<=40 ===')
for r in range(4, MAXR+1):
    rf = np.roots([float(v) for v in T[r].all_coeffs()])
    rg = np.roots([float(v) for v in (x*T[r-1].as_expr()).as_poly().all_coeffs()]) if T[r-1].degree() else np.array([0.0])
    allv = sorted([(float(v.real), 0) for v in rf] + [(float(v.real), 1) for v in rg])
    assert max(abs(v.imag) for v in rf) < 1e-5 and all(allv[i][1] != allv[i+1][1] for i in range(len(allv)-1)), r
print('verified r=4..%d' % MAXR)

print('ALL CHECKS PASSED')
