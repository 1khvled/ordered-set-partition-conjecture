"""Exact exploratory checks for the ordered-big-block Bezoutian proof track.

Run with a Python interpreter having SymPy installed.  All printed matrix
entries, factorizations, characteristic polynomials, and signs are exact.
"""
from __future__ import annotations

import sympy as sp


x, y, t = sp.symbols("x y t")


def bezout_kernel(q, p, a=x, b=y):
    """Root-positive Bezoutian (q(b)p(a)-q(a)p(b))/(a-b).

    Thus B(q,p)(mu,mu)=-q'(mu)p(mu) when q(mu)=0, which is the
    positive quantity in WRS.  This is the negative of the convention
    (q(a)p(b)-q(b)p(a))/(a-b).
    """
    return sp.cancel((q.subs(t, b) * p.subs(t, a) - q.subs(t, a) * p.subs(t, b)) / (a - b))


def bezout_matrix(q, p):
    k = sp.Poly(bezout_kernel(q, p), x, y)
    dx, dy = k.degree(x), k.degree(y)
    assert dx == dy, (p, q, k)
    return sp.Matrix(dx + 1, dy + 1, lambda i, j: k.coeff_monomial(x**i * y**j))


def qop(p):
    return sp.expand(p + (t - 1) * sp.diff(p, t))


def tqop(p):
    return sp.expand(t * qop(p))


def rows(j, nmax):
    out = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        r = sp.Integer(0)
        for s in range(1, min(n, j - 1) + 1):
            r += sp.binomial(n, s) * out[n - s]
        for s in range(j, n + 1):
            r += t * sp.binomial(n, s) * out[n - s]
        out.append(sp.expand(r))
    return out


def show_intrinsic():
    print("Intrinsic kernels B(q,tQ(q)) for arbitrary real-rooted q")
    examples = [
        t + 3,
        (t + 1) * (t + 4),
        (t + 1) * (t + 2) * (t + 5),
        (t + 1) * (t + 2) * (t + 4) * (t + 7),
    ]
    for q in examples:
        b = bezout_matrix(q, tqop(q))
        print("q =", sp.factor(q))
        print("B =", b.tolist())
        print("charpoly =", sp.factor(b.charpoly().as_expr()))
        print("principal minors =", [sp.factor(b[:k, :k].det()) for k in range(1, b.rows + 1)])


def check_decomposition_and_obstruction():
    q = t + 3
    p = (t + 2) * (t + 4)
    a5 = sp.expand(5 * p + 2 * t * sp.diff(p, t))
    h = sp.expand(t * (t - 1) * sp.diff(q, t))
    tq = sp.expand(t * q)
    tqop_q = sp.expand(t * qop(q))

    assert sp.expand(tqop_q - tq - h) == 0
    assert bezout_matrix(q, tq) == sp.Matrix([[9, 3], [3, 1]])
    assert bezout_matrix(q, tqop_q) == bezout_matrix(q, tq) + bezout_matrix(q, h)
    assert bezout_matrix(q, a5) == sp.Matrix([[86, 27], [27, 9]])
    assert bezout_matrix(q, a5).det() == 45
    assert bezout_matrix(q, tqop_q) == sp.Matrix([[6, 6], [6, 2]])
    assert bezout_matrix(q, tqop_q).det() == -24

    test = sp.expand(a5 + 4 * tqop_q)
    assert test == 17 * t**2 + 50 * t + 40
    assert sp.discriminant(test, t) == -220
    assert bezout_matrix(q, test) == sp.Matrix([[110, 51], [51, 17]])
    assert bezout_matrix(q, test).det() == -731
    return "all-index decomposition specialized exactly; compatible-window obstruction verified"


def show_original():
    print("\nOriginal-family Bezoutians")
    for j, nmax in [(2, 8), (3, 9), (4, 10)]:
        rr = rows(j, nmax)
        print("j =", j)
        for n in range(j, nmax + 1):
            q, p = rr[n - j], rr[n]
            b = bezout_matrix(q, p)
            minors = [sp.factor(b[:k, :k].det()) for k in range(1, b.rows + 1)]
            print(f"n={n}, deg q={sp.degree(q,t)}, B={b.tolist()}, minors={minors}")


if __name__ == "__main__":
    print(check_decomposition_and_obstruction())
    show_intrinsic()
    show_original()
