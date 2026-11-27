"""Finite evidence probe: is W(O_n,O_{n-d})(t) > 0 for ALL real t for the
actual ordered big-block family (i.e. is O_n/O_{n-d} a Pick/Nevanlinna
function)?  Exact integer polynomial arithmetic + exact Sturm root counts.

Evidence only -- no float is used as a certificate.
"""
import sympy as sp

t = sp.symbols('t')


def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sp.Integer(0)
        for r in range(1, min(j - 1, n) + 1):
            s += sp.binomial(n, r) * O[n - r]
        for r in range(j, n + 1):
            s += t * sp.binomial(n, r) * O[n - r]
        O.append(sp.expand(s))
    return O


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


def real_root_count_exact(poly):
    # Sturm sequence over QQ; count distinct real roots.
    return sp.polys.polytools.count_roots(sp.Poly(poly, t))


def main():
    bad = []
    tested = 0
    for j in range(2, 9):
        nmax = min(6 * j, 40)
        O = rows(j, nmax)
        for n in range(j, nmax + 1):
            for d in range(1, j + 1):
                w = W(O[n], O[n - d])
                if w == 0:
                    continue
                tested += 1
                cnt = real_root_count_exact(w)
                if cnt != 0:
                    bad.append((j, n, d, cnt, sp.factor(w)))
    print("tested", tested, "bad", len(bad))
    for b in bad[:10]:
        print("BAD", b[0], b[1], b[2], "real roots of W:", b[3])


if __name__ == '__main__':
    main()
