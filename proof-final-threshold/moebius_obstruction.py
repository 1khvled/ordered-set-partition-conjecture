"""Rigorous algebraic check: no fixed Moebius map in t sends the pole of F_{j-1}
to the pole of F_j, for any j>=2.

Write x = e^z, P = E_{j-2}(z), Q = 2-P (old), p = z^{j-1}/(j-1)!.
Then tau_{j-1} = Q/(x-P), tau_j = (Q-p)/(x-P-p).
Assume tau_j = (a*tau_{j-1}+b)/(c*tau_{j-1}+d) with constants a,b,c,d.
After cross-multiplication the two sides are quadratic polynomials in x with
z-dependent coefficients; comparing the x^2 and x coefficients gives b=0 and
a*(Q-p)=d*Q, impossible for constants a,d since (Q-p)/Q is nonconstant.
"""
import sympy as sp

x, z, a, b, c, d = sp.symbols('x z a b c d')

for j in range(2, 9):
    P = sum(z**k / sp.factorial(k) for k in range(j - 1))
    Q = 2 - P
    p = z**(j - 1) / sp.factorial(j - 1)
    Pp, Qp = P + p, Q - p
    # tau_old = Q/(x-P), tau_new = Qp/(x-Pp)
    lhs = sp.expand(Qp * (c * Q + d * (x - P)))          # Qp*(c*tau_old+d)*(x-P)
    rhs = sp.expand((x - Pp) * (a * Q + b * (x - P)))    # (x-Pp)*(a*tau_old+b)*(x-P)
    poly = sp.Poly(sp.expand(lhs - rhs), x)
    x2 = sp.simplify(poly.coeff_monomial(x**2))   # = -b
    x1 = sp.simplify(poly.coeff_monomial(x))       # = d*Qp - (a*Q - b*(Pp+P))
    print('j=%d  x2 coeff = %s  (should be -b)' % (j, x2))
    assert sp.simplify(x2 + b) == 0
    # with b=0, compare x coefficient: must have d*Qp = a*Q identically in z
    x1_b0 = sp.simplify(x1.subs(b, 0))
    print('    x coeff with b=0: %s' % sp.factor(x1_b0))
    # Show d*Qp - a*Q cannot vanish identically unless a=d=0
    # solve by comparing two z values numerically -> a=d=0
    vals = [sp.Rational(1, 2), sp.Rational(3, 2)]
    sol = sp.solve([sp.Eq(x1_b0.subs(z, v).subs(a, 0), 0) for v in vals], [d], dict=True)
    # direct ratio Qp/Q is nonconstant:
    ratio = sp.simplify(Qp / Q)
    assert sp.diff(ratio, z) != 0, j
print('Moebius obstruction verified for j=2..8')
