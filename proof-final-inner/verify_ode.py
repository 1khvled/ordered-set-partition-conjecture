"""Exploration/verification only: symbolic check of the g_j and H_j ODEs for j=2,3,4,5."""
import sympy as sp

z, t = sp.symbols('z t')

for j in [2, 3, 4, 5]:
    A = sum(z**r / sp.factorial(r) for r in range(1, j))
    E = sum(z**r / sp.factorial(r) for r in range(0, j))
    C = sp.exp(z) - E
    N = 1 - A
    p = z**(j - 1) / sp.factorial(j - 1)
    g = C / N
    H = 1 / (1 - t * g)

    # Lemma 2: N*g' = (2-p)*g + p
    e1 = sp.simplify(sp.diff(g, z) * N - ((2 - p) * g + p))
    # Lemma 3: N*H_z = t(2+(t-1)p) H_t + t p H
    e2 = sp.simplify(sp.diff(H, z) * N - (t * (2 + (t - 1) * p) * sp.diff(H, t) + t * p * H))
    # also: N*H_z = (2+(t-1)p)H^2 + (p-2)H
    e3 = sp.simplify(sp.diff(H, z) * N - ((2 + (t - 1) * p) * H**2 + (p - 2) * H))
    print('j=%d  g-ODE residual=%s  PDE residual=%s  Riccati residual=%s'
          % (j, sp.simplify(e1), sp.simplify(e2), sp.simplify(e3)))
