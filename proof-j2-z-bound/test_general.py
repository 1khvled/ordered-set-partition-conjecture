"""Test: for arbitrary real-rooted f with positive coeffs, is -x v'/f > a at roots of v,
where a v = (a+1) f + x f' ?  (label finite evidence)"""
import numpy as np

def poly_from_roots(rs):
    c = np.array([1.0])
    for r in rs:
        c = np.polymul(c, [-r, 1.0])
    return c  # low-to-high

def evalp(c, x):
    return np.polyval(c[::-1], x)

def test(rs, av):
    c = poly_from_roots(rs)
    # v = ((a+1) f + x f')/a
    n = len(c) - 1
    cf = np.zeros(n + 1)
    for i, ci in enumerate(c):
        cf[i] += (av + 1) * ci
        if i >= 1:
            cf[i] += i * ci
    cf = cf / av
    roots = np.roots(cf[::-1])
    vals = []
    for z in roots:
        if abs(z.imag) < 1e-8 and z.real > 0:
            xv = evalp(np.polyder(cf[::-1])[::-1], z.real) * z.real
            fv = evalp(c, z.real)
            vals.append(-xv / fv)
    return vals

rng = np.random.default_rng(0)
worst = None
for trial in range(20000):
    k = rng.integers(2, 7)
    rs = np.sort(rng.uniform(0.01, 100, size=k))
    av = rng.uniform(1.2, 10)
    vals = test(rs, av)
    for v in vals:
        if worst is None or v < worst[0]:
            worst = (v, av, rs.copy())
print("min -xv'/f over random f:", worst[0], "a=", worst[1], "roots=", worst[2])
