"""Test operator lemma with the two correct orientations.

p ~ F_j, q ~ F_{j-1}; p and x q interlace.
Case A (j even): deg p = deg q + 1, pattern p, xq, p, ..., 0.
Case B (j odd):  deg p = deg q,     pattern xq, p, xq, ..., 0.
Claim: P=2p+xp' and Q=x(2q+xq') interlace in the corresponding pattern.
"""
import numpy as np
import random
random.seed(2)


def mk(roots):
    c = np.poly1d([1.0])
    for r in roots:
        c *= np.poly1d([1.0, -r])
    return c


def rts(c):
    r = np.roots(c.c) if c.order > 0 else np.array([])
    return np.sort(r.real), max(abs(r.imag)) if len(r) else 0.0


def alt(a, b):
    allr = sorted([(v,0) for v in a] + [(v,1) for v in b])
    return all(allr[i][1] != allr[i+1][1] for i in range(len(allr)-1))


fails = trials = 0
for nq in range(0, 5):
    for case in ('A','B'):
        np_ = nq+1 if case=='A' else nq
        for _ in range(300):
            total = np_ + nq
            slots = sorted([-random.uniform(0.05, 50) for _ in range(total)])
            if case == 'A':
                # p leftmost, then alternate; xq has 0 at right end -> pattern p,q,p,q,...,q,0
                labels = ['p' if i%2==0 else 'q' for i in range(total)]
            else:
                # xq leftmost: pattern q,p,q,p,...,p? xq has nq+1 roots incl 0
                labels = ['q' if i%2==0 else 'p' for i in range(total)]
            if labels.count('p')!=np_ or labels.count('q')!=nq: continue
            pr = [s for s,l in zip(slots,labels) if l=='p']
            qr = [s for s,l in zip(slots,labels) if l=='q']
            p, q = mk(pr), mk(qr)
            Pp = 2*p + np.polymul(np.poly1d([1.,0.]), np.polyder(p))
            Q = np.polymul(np.poly1d([1.,0.]), 2*q + np.polymul(np.poly1d([1.,0.]), np.polyder(q)))
            rp, ep = rts(Pp); rq, eq = rts(Q)
            if ep>1e-5 or eq>1e-5: continue
            trials += 1
            if not alt(rp, rq):
                fails += 1
                if fails < 5:
                    print('FAIL case',case,'np',np_,'nq',nq)
                    print(' p',pr,' xq',qr+[0])
                    print(' P',rp,' Q',rq)
print('trials',trials,'fails',fails)
