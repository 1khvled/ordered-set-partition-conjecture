"""Exact symbolic verification of the chapter-renewal Mobius-flow identities."""
import sympy as sp

z, t, t0 = sp.symbols('z t t0')

for j in [2, 3, 4, 5]:
    A = sum(z**r / sp.factorial(r) for r in range(1, j))
    C = sp.exp(z) - sum(z**r / sp.factorial(r) for r in range(0, j))
    B = 1 - A
    p = z**(j-1) / sp.factorial(j-1)
    g = C / B
    D = B - t*C
    H = B / D
    F = 1 / D

    phi = B*t0 / (C*t0 + sp.exp(z))

    # (i) H(z, phi_z(t0)) = 1 + exp(-z) C(z) t0
    r1 = sp.simplify(H.subs(t, phi) - (1 + sp.exp(-z)*C*t0))
    # (ii) F(z, phi_z(t0)) = (1+exp(-z) C t0)/B
    r2 = sp.simplify(F.subs(t, phi) - (1 + sp.exp(-z)*C*t0)/B)
    # (iii) g solves (1-A) g' = (2-p) g + p
    r3 = sp.simplify(B*sp.diff(g, z) - ((2-p)*g + p))
    # (iv) PDE N H_z = t[2+(t-1)p]H_t + t p H
    r4 = sp.simplify(B*sp.diff(H, z) - (t*(2+(t-1)*p)*sp.diff(H, t) + t*p*H))
    # (v) phi solves d/dz phi = -(2-p)/B phi - p/B phi^2
    r5 = sp.simplify(sp.diff(phi, z) - (-(2-p)/B*phi - p/B*phi**2))
    # (vi) Riccati for F: F' = (2+(t-1)p) F^2 - F
    r6 = sp.simplify(sp.diff(F, z) - ((2+(t-1)*p)*F**2 - F))
    print('j=%d residuals: %s' % (j, [r1, r2, r3, r4, r5, r6]))
