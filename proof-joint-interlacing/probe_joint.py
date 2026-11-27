"""Probe joint interlacing structure between H_{n,j} and O_{n-d,j}.

Checks:
- W(H_n, O_{n-d}) real-root-free (boundary Wronskian)
- root interlacing pattern of H_n vs O_{n-d}
- common interlacer candidates
- HKO/compatible: a*H_n + b*O_{n-d} real-rooted for sampled a,b
Evidence only.
"""
import sympy as sp
from itertools import combinations

t = sp.Symbol('t')

def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax+1):
        s = sp.Integer(0)
        for r in range(1, min(j-1, n)+1):
            s += sp.binomial(n, r) * O[n-r]
        for r in range(j, n+1):
            s += t * sp.binomial(n, r) * O[n-r]
        O.append(sp.expand(s))
    H = [sp.Integer(1)]  # H_0
    for n in range(1, nmax+1):
        # H_n = t * sum_{r=j}^n binom(n,r) O_{n-r}
        s = sp.Integer(0)
        for r in range(j, n+1):
            s += sp.binomial(n, r) * O[n-r]
        H.append(sp.expand(t*s))
    return O, H

def W(p, q):
    return sp.expand(sp.diff(p, t)*q - p*sp.diff(q, t))

def rr(poly):
    p = sp.Poly(poly, t)
    if p.degree() <= 0:
        return 0
    return sp.count_roots(p)

def real_roots(poly):
    p = sp.Poly(poly, t)
    if p.degree() <= 0:
        return []
    # isolate real roots
    intervals = sp.polys.polytools.real_roots(p)  # exact algebraic
    return sorted([sp.N(r, 30) for r in intervals])

def is_interlaced(roots_a, roots_b):
    # check if roots_a interlace roots_b (both sorted ascending, all real, simple)
    # returns 'a<b<a...' or 'b<a<b...' or None
    if len(roots_a) == 0 or len(roots_b) == 0:
        return 'trivial'
    # merge
    allr = [('a',r) for r in roots_a] + [('b',r) for r in roots_b]
    allr.sort(key=lambda x: x[1])
    labels = [x[0] for x in allr]
    # check alternating
    ok_ab = all(labels[i] != labels[i+1] for i in range(len(labels)-1))
    return ok_ab

def interlace_pattern(p, q, name_p='p', name_q='q'):
    rp = real_roots(p)
    rq = real_roots(q)
    dp = sp.Poly(p,t).degree() if p!=0 else -1
    dq = sp.Poly(q,t).degree() if q!=0 else -1
    if len(rp) != max(dp,0) or len(rq) != max(dq,0):
        return {'real_rooted': False, 'rp':rp,'rq':rq}
    # check interlacing
    # merge
    allr = [(name_p,r) for r in rp] + [(name_q,r) for r in rq]
    allr.sort(key=lambda x: float(x[1]))
    labels = [x[0] for x in allr]
    return {'real_rooted': True, 'labels': labels, 'rp':rp, 'rq':rq, 'dp':dp,'dq':dq}

def main():
    for j in range(2, 6):
        nmax = min(5*j, 24)
        O,H = rows(j, nmax)
        print(f"\n=== j={j} ===")
        for n in range(j, nmax+1):
            for d in range(1, j+1):
                q = O[n-d]
                p = H[n]
                if p == 0 or q == 0:
                    continue
                w = W(p,q)
                cnt = rr(w) if w != 0 else -1
                pat = interlace_pattern(p,q,'H','O')
                # also check W(H_n, H_{n-?}) or common interlacer with O_{n-d-1}
                # check HKO: sample linear combos
                hko_bad = False
                for a,b in [(1,1),(1,-1),(1,2),(2,-1),(1,3),(-1,3)]:
                    lc = sp.expand(a*p+b*q)
                    if lc == 0:
                        continue
                    if sp.Poly(lc,t).degree() > 0 and rr(lc) != sp.Poly(lc,t).degree():
                        hko_bad = True
                        break
                # check if H_n interlaces O_{n-d}
                inter = False
                if pat['real_rooted']:
                    labels = pat['labels']
                    inter = all(labels[i] != labels[i+1] for i in range(len(labels)-1))
                print(f" n={n} d={d} degH={sp.degree(p,t)} degO={sp.degree(q,t)} Wreal={cnt} interlace={inter} HKO_bad={hko_bad} labels={pat.get('labels')}")

if __name__ == '__main__':
    main()
