"""Exact polynomial operations for the Z-bound study. No installs; sympy only."""
import sympy as sp

a, x, t = sp.symbols("a x t", positive=True)


def astirling_table(mmax):
    rmax = mmax // 2 + 2
    A = [[0] * rmax for _ in range(mmax + 1)]
    A[0][0] = 1
    for m in range(mmax):
        for r in range(rmax - 1):
            A[m + 1][r] += r * A[m][r]
            if r:
                A[m + 1][r] += m * A[m - 1][r - 1]
    return A


A = astirling_table(30)


def rising(av, r):
    v = 1
    for k in range(r):
        v *= av + k
    return v


def Q(m, av, z):
    return sp.expand(sum(rising(av, r) * A[m][r] * z**(r - 1)
                         for r in range(1, len(A[0])) if A[m][r]))


def R(m, av):
    return Q(m, av, -x)


def Wt(f, g, z=t):
    return sp.expand(sp.diff(f, z) * g - f * sp.diff(g, z))


def Bx(f, g):
    return sp.expand(f * sp.diff(g, x) - sp.diff(f, x) * g)
