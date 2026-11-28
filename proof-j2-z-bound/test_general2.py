"""Finite evidence: sign of Z at roots of v=((a+1)f+xf')/a for random real-rooted f."""
import numpy as np

def ev(c, z): return np.polyval(c[::-1], z)
dev = lambda c: np.polyder(c[::-1])[::-1]
def poly_from_roots(rs):
    c = np.array([1.0])
    for r in rs: c = np.polymul(c, [-r, 1.0])
    return c

rng = np.random.default_rng(1)
cnt = 0
bad = []
for trial in range(5000):
    k = rng.integers(2, 8)
    rs = np.sort(rng.uniform(0.01, 100, size=k))
    av = rng.uniform(1.1, 20)
    f = poly_from_roots(rs)
    n = len(f) - 1
    cf = np.zeros(n + 1)
    for i, ci in enumerate(f):
        cf[i] += (av + 1) * ci
        if i >= 1: cf[i] += i * ci
    cf /= av
    roots = sorted([z.real for z in np.roots(cf[::-1]) if abs(z.imag)<1e-7 and z.real>0])
    Zs = [s*s*ev(dev(dev(f)), s)/ev(f, s) for s in roots]
    for j, Z in enumerate(Zs):
        cnt += 1
        if j >= 1 and Z > 0:
            bad.append((av, rs, j, Z, roots))
print("total non-first roots:", cnt, "positive Z count:", len(bad))
if bad:
    print("example:", bad[0])
