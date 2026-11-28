"""Test lemma: (x f')' <= 0 on (0, first root) for real-rooted f with positive coeffs."""
import numpy as np
from ops import astirling_table, rising
A = astirling_table(30)
def coefs(m, av):
    return np.array([rising(av, r)*A[m][r] for r in range(1, m//2+1)], float)
def evp(c, xs): return np.polyval(c[::-1], xs)
def dev(c): return np.polyder(c[::-1])[::-1]

def worst_on_first(c):
    d1 = dev(c); d2 = dev(d1)
    g = np.polymul(np.array([1.,0.]), d2[::-1])[::-1] + d1  # x f'' + f'
    roots = sorted(z.real for z in np.roots(c[::-1]) if abs(z.imag)<1e-9 and z.real>0)
    z = roots[0]
    xs = np.linspace(1e-8*z, z*(1-1e-6), 400)
    return np.max(evp(g, xs)), z

print("Actual Q(n,a):")
for n in [6,8,10,14,20]:
    for av in [1,2,5,20]:
        m, z = worst_on_first(coefs(n, av)*np.array([(-1)**i for i in range(n//2)]))
        print(f"  n={n} a={av}: max(xR''+R') on (0,rho1) = {m:.6g} (rho1={z:.4g})")

print("Random real-rooted with positive coefficients:")
rng = np.random.default_rng(7)
cnt=0; bad=0; worst=(-1e18,None)
for trial in range(3000):
    k = rng.integers(2,8)
    rs = np.sort(rng.uniform(0.01,100,k))
    c = np.array([1.0])
    for r in rs: c = np.polymul(c, [-r,1.0])
    d1=dev(c); d2=dev(d1)
    g = np.polymul(np.array([1.,0.]), d2[::-1])[::-1] + d1
    m,z = worst_on_first(c)
    cnt+=1
    if m>1e-8: bad+=1
    if m>worst[0]: worst=(m,(rs,z))
print(f"  {cnt} trials, {bad} with max>0; worst={worst}")
