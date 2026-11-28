"""Empirically test: positive log-concave coefficients => Hurwitz stable?"""
import numpy as np, random
random.seed(11)
bad=[]
for trial in range(20000):
    d=random.randint(2,7)
    c=[1.0]
    # build log-concave sequence randomly: choose ratios decreasing
    ratios=sorted([random.uniform(0.05,20) for _ in range(d)], reverse=True)
    for r in ratios: c.append(c[-1]*r)
    roots=np.roots(list(reversed(c)))
    if max(z.real for z in roots) >= 0:
        bad.append((c,max(z.real for z in roots)))
print('counterexamples',len(bad))
for b in bad[:5]: print(b[0],b[1])
