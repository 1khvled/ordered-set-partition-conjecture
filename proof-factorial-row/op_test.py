"""Test: if p and x*q interlace (negative roots), then (x^2 p)'/x and (x^2 q)' interlace."""
import numpy as np
import random
random.seed(1)


def mk(roots):
    c = np.poly1d([1.0])
    for r in roots:
        c *= np.poly1d([1.0, -r])
    return c


def roots(c):
    r = np.roots(c.c)
    return np.sort(r.real), max(abs(r.imag))


def alt(a, b):
    allr = sorted([(v,0) for v in a] + [(v,1) for v in b])
    return all(allr[i][1] != allr[i+1][1] for i in range(len(allr)-1))


fails = trials = 0
for np_ in range(1, 5):
    for nq in range(max(0,np_-1), np_+2):
        for _ in range(200):
            total = np_ + nq
            slots = sorted([-random.uniform(0.05, 50) for _ in range(total)])
            for first in ('p','q'):
                labels = [first if i%2==0 else ('q' if first=='p' else 'p') for i in range(total)]
                if labels.count('p')!=np_ or labels.count('q')!=nq: continue
                pr = [s for s,l in zip(slots,labels) if l=='p']
                qr = [s for s,l in zip(slots,labels) if l=='q']
                p, q = mk(pr), mk(qr)
                # P = 2p + x p'  ; Q = x(2q + x q')
                Pp = 2*p + np.polymul(np.poly1d([1.,0.]), np.polyder(p))
                Q = np.polymul(np.poly1d([1.,0.]), 2*q + np.polymul(np.poly1d([1.,0.]), np.polyder(q)))
                rp, ep = roots(Pp); rq, eq = roots(Q)
                if ep>1e-5 or eq>1e-5: continue
                trials += 1
                if not alt(rp, rq):
                    fails += 1
                    if fails < 4:
                        print('FAIL np=%d nq=%d first=%s' % (np_,nq,first))
                        print(' p roots',pr,' xq roots',qr+[0.0])
                        print(' P roots',rp,' Q roots',rq)
print('trials',trials,'fails',fails)
