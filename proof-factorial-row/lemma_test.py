"""Test candidate general lemma:
if f and x*g interlace (negative roots), then h=(x(f+lam*x*g))' interlaces x*f.
Uses descending-coefficient numpy.poly1d throughout.
"""
import numpy as np
import random

random.seed(0)


def coeffs_from_roots(roots):
    c = np.poly1d([1.0])
    for r in roots:
        c *= np.poly1d([1.0, -r])
    return c


def roots_desc(c):
    r = np.roots(c.c)
    return np.sort(r.real), max(abs(r.imag))


def alternates(a, b):
    allr = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    return all(allr[i][1] != allr[i + 1][1] for i in range(len(allr) - 1))


fails = 0
trials = 0
for n in range(1, 5):
    for ng in range(max(0, n - 1), n + 2):
        for _ in range(150):
            total = n + ng
            slots = sorted([-random.uniform(0.1, 50) for _ in range(total)])
            for first in ('f', 'g'):
                labels = [first if i % 2 == 0 else ('g' if first == 'f' else 'f') for i in range(total)]
                if labels.count('f') != n or labels.count('g') != ng:
                    continue
                fr = [s for s, l in zip(slots, labels) if l == 'f']
                gr = [s for s, l in zip(slots, labels) if l == 'g']
                f = coeffs_from_roots(fr)
                g = coeffs_from_roots(gr)
                lam = random.choice([0.1, 0.5, 1, 2, 7, 50, 500])
                # p = x( f + lam x g )
                inner = np.polyadd(f, lam * np.polymul(np.poly1d([1.0, 0.0]), g))
                p = np.polymul(inner, np.poly1d([1.0, 0.0]))
                h = np.polyder(p)
                xf = np.polymul(f, np.poly1d([1.0, 0.0]))
                rh, eh = roots_desc(h)
                rf, ef = roots_desc(xf)
                if eh > 1e-5 or ef > 1e-5:
                    continue
                trials += 1
                if not alternates(rh, rf):
                    fails += 1
                    if fails < 5:
                        print('FAIL n=%d ng=%d lam=%s' % (n, ng, lam))
                        print(' f roots', fr, ' xg roots', gr + [0.0])
                        print(' h roots', rh, ' xf roots', rf)
print('trials', trials, 'fails', fails)
