"""Exact symbolic verification of the identities used in NOTES.md.

All residuals are simplified with SymPy over QQ(e^z) where possible; no
numeric approximation enters a certificate.
"""
import sympy as sp

z, w, t = sp.symbols('z w t')
z1, z2, t1, t2 = sp.symbols('z1 z2 t1 t2')


def E(k):
    return sum(z**r / sp.factorial(r) for r in range(k + 1))


def data(j):
    Ej = E(j - 1)
    A = Ej - 1
    C = sp.exp(z) - Ej
    B = 1 - A
    D = 1 - A - t * C
    return A, C, B, D


# 1. alpha'(z) numerator = e^z (p-2) + 2 E_{j-2}
for j in range(2, 7):
    A, C, B, D = data(j)
    alpha = B / C
    p = z**(j - 1) / sp.factorial(j - 1)
    Ep = sum(z**r / sp.factorial(r) for r in range(j - 1))
    claimed = sp.exp(z) * (p - 2) + 2 * Ep
    resid = sp.simplify(sp.diff(alpha, z) * C**2 - claimed)
    print("alpha' j=", j, "residual", resid)
    assert resid == 0

# 2. Phi(z,w) = C(z)-C(w)+C(w)A(z)-C(z)A(w) = C(z)B(w)-C(w)B(z)
#    and W = n!(n-j)![z^n w^{n-j}] Phi/(D(z)^2 D(w)^2)
for j in range(2, 6):
    A, C, B, D = data(j)
    Cw = C.subs(z, w)
    Aw = A.subs(z, w)
    Phi = C - Cw + Cw * A - C * Aw
    Phi2 = C * B.subs(z, w) - Cw * B
    print("Phi identity j=", j, sp.simplify(Phi - Phi2))
    assert sp.simplify(Phi - Phi2) == 0

# 3. j=2: Psi(z,z) = (2-z)e^z - 2 changes sign on z>0
j = 2
A, C, B, D = data(j)
Cw = C.subs(z, w)
Aw = A.subs(z, w)
Phi = sp.expand(C - Cw + Cw * A - C * Aw)
Psi_diag = sp.simplify(sp.diff(Phi, z).subs(w, z))
print("j=2 Psi(z,z) =", Psi_diag)
print("  at z=1:", sp.simplify(Psi_diag.subs(z, 1)), "=", sp.N(Psi_diag.subs(z, 1)))
print("  at z=3:", sp.simplify(Psi_diag.subs(z, 3)), "=", sp.N(Psi_diag.subs(z, 3)))

# 4. alpha' sign change at j=2: values at z=1 and z=2
alpha2 = B / C
d2 = sp.simplify(sp.diff(alpha2, z) * C**2)
print("j=2 alpha' numerator at 1:", sp.simplify(d2.subs(z, 1)), sp.N(d2.subs(z, 1)))
print("j=2 alpha' numerator at 2:", sp.simplify(d2.subs(z, 2)), sp.N(d2.subs(z, 2)))

# 5. Pick kernel check for the concrete ratio O_5/O_3 (j=2):
#    R(t) = (110t^2+311t+120)/(7t+6) = (110/7)t+1517/49 + (3222/343)/(mu-t), mu=-6/7
mu = sp.Rational(-6, 7)
p = 120 + 311 * t + 110 * t**2
q = 6 + 7 * t
R = sp.cancel(p / q)
check = sp.simplify(R - (sp.Rational(110, 7) * t + sp.Rational(1517, 49)
                         + sp.Rational(3222, 343) / (mu - t)))
print("Pick decomposition residual O_5/O_3:", check)
assert check == 0

# 6. two-variable kernel minor, general rank-one resolvent form
a1, a2, w1, w2 = sp.symbols('a1 a2 w1 w2')
F11 = w1 / (a1 - t1)
F12 = w1 / (a1 - t2)
F21 = w2 / (a2 - t1)
F22 = w2 / (a2 - t2)
minor = sp.expand(F11 * F22 - F12 * F21)
claimed = (w1 * w2 * (a1 - a2) * (t2 - t1)
           / ((a1 - t1) * (a1 - t2) * (a2 - t1) * (a2 - t2)))
print("kernel minor residual:", sp.simplify(minor - claimed))
assert sp.simplify(minor - claimed) == 0
