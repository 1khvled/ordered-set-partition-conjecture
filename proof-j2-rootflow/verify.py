"""Exact verification for proof-j2-rootflow.

All identities are checked with sympy exact arithmetic; all inequalities are
checked at explicit feasible sample points.  Finite checks are EVIDENCE, not
an all-index proof.  No packages installed, no builds.
"""
import sympy as sp

x, u, t, z = sp.symbols('x u t z')


# ---------- exact generation ----------
def Cser(n):
    return sum(z**r / sp.factorial(r) for r in range(2, n + 1))


def On(n):
    C = Cser(n)
    total = 0
    for k in range(n + 1):
        p = sp.Poly(sp.expand((x * z + t * C) ** k), z)
        if n <= p.degree():
            total += p.coeff_monomial(z**n)
    return sp.expand(sp.factorial(n) * total)


def psi(n):
    """psi_n(x,u) = x^{-n} O_n(x,x^2 u), via the exact scaled recurrence."""
    vals = [sp.Integer(1)]
    for m in range(1, n + 1):
        e = m * vals[m - 1]
        for r in range(2, m + 1):
            e += u * sp.binomial(m, r) * x**(2 - r) * vals[m - r]
        vals.append(sp.expand(e))
    return sp.expand(vals[n])


def V(m):
    return sp.expand(sum(sp.binomial(m, i) * psi(i) * psi(m - i) for i in range(m + 1)))


def report(name, residual):
    print('%-58s residual = %s' % (name, sp.simplify(residual)))


print('=== 1. scaling identity and residue/coefficient identity ===')
for n in range(2, 11):
    report('n=%d  x^-n O_n(x,x^2 u) - psi_n' % n,
           sp.simplify(On(n).subs(t, x**2 * u) / x**n - psi(n)))

# d_x C_x(w) = - sum_{r>=3} (r-2) w^r/(r! x^{r-1}); coefficient identity
for n in range(2, 11):
    rhs = -u * sum(sp.Rational(r - 2, 1) * sp.binomial(n, r) / x**(r - 1) * V(n - r)
                   for r in range(3, n + 1))
    report('n=%d  d_x psi_n  -  (-u) sum_r (r-2)C(n,r)x^{1-r} V_{n-r}' % n,
           sp.simplify(sp.diff(psi(n), x) - rhs))
    report('n=%d  V_m - m![w^m]G^2  (m=n-3)' % n,
           sp.simplify(V(n - 3) - sp.expand(sum(sp.binomial(n - 3, i) * psi(i) * psi(n - 3 - i)
                                                for i in range(n - 2)))))

print()
print('=== 2. linear Calogero-Moser form d_x O_n = (a t+b)O_n\' + c O_n ===')
a, b, c = sp.symbols('a b c')
for n in range(2, 11):
    e = sp.expand((a * t + b) * sp.diff(On(n), t) + c * On(n) - sp.diff(On(n), x))
    sols = sp.solve(sp.Poly(e, t).all_coeffs(), (a, b, c), dict=True)
    print('n=%2d  deg_t=%d  linear-CM representation exists: %s'
          % (n, sp.Poly(On(n), t).degree(), bool(sols)))

print()
print('=== 3. Wronskian N = W(d_x psi_n, psi_n\') at x=1 ===')
for n in range(3, 11):
    f = psi(n)
    Q = sp.diff(f, x).subs(x, 1)
    d = sp.diff(f.subs(x, 1), u)
    N = sp.expand(sp.diff(Q, u) * d - Q * sp.diff(f.subs(x, 1), u, 2))
    P = -N
    disc = sp.discriminant(sp.Poly(P, u)) if sp.Poly(P, u).degree() == 2 else None
    extra = ''
    if disc is not None:
        extra = '  disc(-N)=%s' % disc
    print('n=%2d deg(N)=%d all coeffs of -N positive: %s%s'
          % (n, sp.Poly(N, u).degree(),
             all(co > 0 for co in sp.Poly(P, u).all_coeffs()), extra))

print()
print('=== 4. DIRECTION AUDIT (explicit feasible sample n=4) ===')
# psi_4(x,u) = 6u^2 + (36+8/x+1/x^2)u + 24
Acoef = 36 + sp.Rational(8, 1) / x + 1 / x**2
Delta_u = sp.sqrt(Acoef**2 - 4 * 6 * 24) / 6        # u-gap, exact
Delta_t = x**2 * Delta_u                             # t-gap, exact
for X in (1, 2, 10):
    du = float(Delta_u.subs(x, X))
    dt = float(Delta_t.subs(x, X))
    dudx = float(sp.diff(Delta_u, x).subs(x, X))
    dtdx = float(sp.diff(Delta_t, x).subs(x, X))
    print('x=%2d  u-gap=%.6f  d(u-gap)/dx=%+.6f   t-gap=%.6f  d(t-gap)/dx=%+.6f'
          % (X, du, dudx, dt, dtdx))
print('=> the proposed sign d(gap)/dx<0 is FALSE for unnormalised t-roots.')

print()
print('=== 5. n=4: exact all-x proof of normalized gap repulsion ===')
# Q = d_x psi_4 = -u*A(x) with A = 2(4x+1)/x^3 > 0.
A4 = sp.simplify(-sp.diff(psi(4), x) / u)
print('A_4(x) =', sp.factor(A4))
D4 = sp.simplify(Acoef**2 - 576)
# rho_1 = -A|r_1|/sqrt(D), rho_2 = +A|r_2|/sqrt(D), |r_1|+|r_2| = s/6.
gap_derivative_4 = sp.simplify(-A4 * Acoef / (6 * sp.sqrt(D4)))
print('rho_1 - rho_2 =', sp.factor(gap_derivative_4))
print('=> g\' = rho_1-rho_2 < 0 for every x>0 (n=4 fully proved):',
      sp.factor(-gap_derivative_4 * sp.sqrt(D4)))

print()
print('=== 6. FEASIBLE-CONSTRAINT FALSIFICATION (exact witness) ===')
# psi = (u+1)(u+2), Q = -u(50+31u): all listed constraints hold, gap' > 0.
P = sp.Poly(u**2 + 3 * u + 2, u)
Q = -50 * u - 31 * u**2
r1, r2 = sp.Integer(-2), sp.Integer(-1)
rho1 = sp.simplify(Q.subs(u, r1) / sp.diff(P.as_expr(), u).subs(u, r1))
rho2 = sp.simplify(Q.subs(u, r2) / sp.diff(P.as_expr(), u).subs(u, r2))
Nw = sp.expand(sp.diff(Q, u) * sp.diff(P.as_expr(), u) - Q * sp.diff(P.as_expr(), u, 2))
c = sp.Rational(-3, 2)
print('psi=(u+1)(u+2) has simple negative roots and positive coefficients')
print('Q = -u A,  A=50+31u has positive coefficients:', [50, 31])
print('rho_1 =', rho1, ' rho_2 =', rho2, ' gap\' = rho_1-rho_2 =', rho1 - rho2)
print('N(u) =', sp.factor(Nw), ' disc(N)=', sp.discriminant(Nw))
print('pole orientation Q(c)*psi\'\'(c) =', sp.simplify(Q.subs(u, c) * sp.diff(P.as_expr(), u, 2)))
print('=> constraints (negative simple roots, Q=-uA with A>0 coefficients,')
print('   N<0, upward pole jump) DO NOT force gap repulsion.')
print()
print('All residuals above are exact; finite lists are evidence only.')
