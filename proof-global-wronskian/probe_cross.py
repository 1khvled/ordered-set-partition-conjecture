"""Evidence probe: is W(H_{n,j}, O_{n-j,j}) sign-definite (no real zeros)?
H_{n,j} = O_n - sum_{r=1}^{j-1} binom(n,r) O_{n-r}.
Exact Sturm counts; evidence only.
"""
import sympy as sp

t = sp.symbols('t')


def rows(j, nmax):
    O = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        s = sum(sp.binomial(n, r) * O[n - r] for r in range(1, min(j - 1, n) + 1))
        s += sum(t * sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
        O.append(sp.expand(s))
    return O


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


def main():
    bad = []
    tested = 0
    for j in range(2, 9):
        nmax = min(6 * j, 40)
        O = rows(j, nmax)
        for n in range(j, nmax + 1):
            H = sp.expand(O[n] - sum(sp.binomial(n, r) * O[n - r]
                                     for r in range(1, min(j - 1, n) + 1)))
            w = W(H, O[n - j])
            if w == 0:
                continue
            tested += 1
            cnt = sp.polys.polytools.count_roots(sp.Poly(w, t))
            if cnt:
                bad.append((j, n, cnt))
    print('tested', tested, 'bad', len(bad), bad[:20])


if __name__ == '__main__':
    main()
