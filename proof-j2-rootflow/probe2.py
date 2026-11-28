"""Probe scaled root velocities and gap monotonicity (EXPLORATION only)."""
import sympy as sp
import mpmath as mp

x, u = sp.symbols('x u')


def psi(n, X=None):
    """psi_n(x,u) = x^{-n} O_n(x,x^2 u) as sympy Poly in u (rational x)."""
    vals = [sp.Integer(1)]  # psi_0
    for m in range(1, n + 1):
        expr = m * vals[m - 1]
        for r in range(2, m + 1):
            expr += u * sp.binomial(m, r) * x**(2 - r) * vals[m - r]
        vals.append(sp.expand(expr))
    if X is not None:
        vals = [sp.expand(v.subs(x, X)) for v in vals]
    return sp.Poly(vals[n], u)


def residues(n, X):
    # build numerator directly from symbolic psi then substitute
    Ps = psi(n)
    dPs = sp.Poly(sp.diff(Ps.as_expr(), x), u)
    Xu = sp.Rational(X) if isinstance(X, int) else X
    Psx = sp.Poly(Ps.as_expr().subs(x, Xu), u)
    dPsx = sp.Poly(dPs.as_expr().subs(x, Xu), u)
    rem = dPsx.rem(Psx)
    coeffs = [mp.mpf(str(c)) for c in Psx.all_coeffs()]
    roots = sorted([mp.re(r) for r in mp.polyroots(coeffs, maxsteps=200, error=False)])
    rho = []
    for rr in roots:
        coeffs_rem = [mp.mpf(str(c)) for c in rem.all_coeffs()]  # descending
        num = mp.polyval(coeffs_rem, rr)
        dPdu = sp.Poly(sp.diff(Psx.as_expr(), u), u)
        den = mp.polyval([mp.mpf(str(c)) for c in dPdu.all_coeffs()], rr)
        rho.append(float(num / den))
    return roots, rho


if __name__ == '__main__':
    mp.mp.dps = 40
    for n in (4, 6, 8, 10):
        for X in (1, 2, 10):
            roots, rho = residues(n, X)
            gaps = [roots[i+1] - roots[i] for i in range(len(roots)-1)]
            rv = [rho[i+1] - rho[i] for i in range(len(rho)-1)]
            print('n=%d x=%s roots=%s' % (n, X, [float(v) for v in roots]))
            print('   rho=%s' % [round(v, 6) for v in rho])
            print('   gaps=%s  dgap/dsign(rho_next-rho)=%s' % ([float(v) for v in gaps], [round(v, 6) for v in rv]))
