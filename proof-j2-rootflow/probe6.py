"""Check monotonicity of f(u)=d_x psi_n/psi_n' and pole orientation."""
import sympy as sp
import mpmath as mp
from probe2 import psi

x, u = sp.symbols('x u')


def analyze(n, X):
    Ps = psi(n).as_expr().subs(x, sp.Rational(X))
    dPdu = sp.diff(Ps, u)
    dPdx = sp.diff(psi(n).as_expr(), x).subs(x, sp.Rational(X))
    d2Pdu2 = sp.diff(Ps, u, 2)
    d2Pdudx = sp.diff(dPdx, u)
    fprime_num = sp.expand(d2Pdudx*dPdu - dPdx*d2Pdu2)
    cps = sp.nroots(sp.Poly(dPdu, u), n=30, maxsteps=200)
    cps = sorted([float(sp.re(c)) for c in cps if abs(sp.im(c)) < 1e-20])
    print('n=%d x=%s critical pts=%s' % (n, X, [round(v,5) for v in cps]))
    for c in cps:
        if c >= 0: continue
        val = float(fprime_num.subs(u, sp.Rational(str(c))))
        jumpprod = float((dPdx*d2Pdu2).subs(u, sp.Rational(str(c))))
        print('   c=%.6f  N(c)=%.6g  dPdx*d2Pdu2=%.6g' % (c, val, jumpprod))
    # sample f' on a grid of negative u avoiding poles
    bad = []
    for k in range(1, 200):
        v = -0.01*k
        den = float(dPdu.subs(u, v))
        if den == 0 or abs(den)<1e-10: continue
        Nval = float(fprime_num.subs(u, v))
        if Nval < 0: bad.append((v, Nval))
    print('   N<0 samples:', len(bad), bad[:5])

for n in (4,6,8,10):
    analyze(n, 1)
