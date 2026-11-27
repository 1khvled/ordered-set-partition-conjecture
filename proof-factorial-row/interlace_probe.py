"""Numerically probe interlacing structure of T_r and x T_{r-1}."""
import numpy as np
from explore import T, U, MAXR


def roots(p):
    c = [int(v) for v in p.all_coeffs()]
    return np.sort(np.roots(c).real) if max(abs(np.roots(c).imag)) < 1e-6 else None


def realroots(p):
    c = [float(v) for v in p.all_coeffs()]
    if p.degree() == 0:
        return np.array([]), 0.0
    r = np.roots(c)
    return np.sort(r.real), max(abs(r.imag))


print('r | deg T_r | deg xT_{r-1} | interleaving pattern (T=uppercase, xTprev=lowercase)')
for r in range(4, 16):
    rt, et = realroots(T[r])
    rp, ep = realroots(T[r - 1])
    rp = np.sort(np.concatenate([rp, [0.0]]))
    assert et < 1e-6 and ep < 1e-6, (r, et, ep)
    allr = sorted([(a, 'T') for a in rt] + [(a, 't') for a in rp])
    # check strict alternation
    ok = all(allr[i][1] != allr[i + 1][1] for i in range(len(allr) - 1))
    print(r, len(rt), len(rp), ''.join(s for _, s in allr), 'alternating=', ok)

print()
print('Signs of T_{r-2} at roots of T_{r-1}, and of U_r at roots of T_{r-1}:')
for r in range(5, 13):
    rt, _ = realroots(T[r - 1])
    f2 = T[r - 2].as_expr()
    fU = U[r].as_expr()
    import sympy as sp
    xs = sp.symbols('x')
    s2 = [sp.N(f2.subs(xs, sp.Float(v, 30))) for v in rt]
    sU = [sp.N(fU.subs(xs, sp.Float(v, 30))) for v in rt]
    print('r=%d' % r, 'T_{r-2} signs:', [float(v) > 0 for v in s2],
          'U_r signs:', [float(v) > 0 for v in sU])
