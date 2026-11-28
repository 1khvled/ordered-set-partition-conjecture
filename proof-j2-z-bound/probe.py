"""Finite numeric probe of quantities at roots of V (labelled finite evidence)."""
import numpy as np
from ops import astirling_table, rising

A = astirling_table(30)

def Qc(m, av):
    deg = m // 2 - 1
    c = np.zeros(deg + 1)
    for r in range(1, m // 2 + 1):
        c[r - 1] = rising(av, r) * A[m][r]
    return c  # low-to-high in t

def Qx(m, av):
    # Q(m,av,-x): coefficients in x
    c = Qc(m, av)
    return c * np.array([(-1.0)**i for i in range(len(c))])

def ev(c, z):
    return np.polyval(c[::-1], z)

dev = lambda c: np.polyder(c[::-1])[::-1]

def probe(n, av):
    R = Qx(n, av); V = Qx(n, av + 1); H = Qx(n - 1, av + 1)
    S = Qx(n + 1, av)
    Y = 2 * H + np.array([0.0] + list(dev(H)))
    roots = np.roots(V[::-1])
    out = []
    for s in roots:
        if s.imag ** 2 < 1e-12 and s.real > 0:
            s = s.real
            Rv = ev(R, s); Vp = ev(dev(V), s)
            Z = s*s*ev(dev(dev(R)), s)/Rv
            out.append(dict(s=s, Z=Z, w=(av+1)*(av+2)-Z, g=-s*Vp/Rv,
                            SR=ev(S, s)/Rv, HR=ev(H, s)/Rv, YR=ev(Y, s)/Rv))
    return sorted(out, key=lambda d: d['s'])

for (n, av) in [(6,2),(8,2),(10,2),(8,5),(12,2),(12,5),(14,3)]:
    print(f"\nn={n} a={av}: 3a+2={3*av+2}, (a+1)^2={(av+1)**2}")
    for d in probe(n, av):
        print("  s=%.5f Z=%.4f w=%.4f g=%.4f S/R=%.3f H/R=%.3f Y/R=%.3f"
              % (d['s'], d['Z'], d['w'], d['g'], d['SR'], d['HR'], d['YR']))
