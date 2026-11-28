"""Shared exact utilities for proof-j2-h-hurwitz.

Q_{n,a}(t) = sum_r (a)_r A[n][r] t^(r-1), A[n][r] = 2-associated Stirling,
A[n][r] = r A[n-1][r] + (n-1) A[n-2][r-1].  a=1 gives Q = sum k! A[n][k+1] t^k.
No finite enumeration is used as a proof in the notes.
"""
import sympy as sp

t = sp.symbols('t')


def triangle(nmax):
    width = nmax // 2 + 3
    A = [[0] * width for _ in range(nmax + 1)]
    A[0][0] = 1
    for m in range(nmax):
        for r in range(width - 1):
            A[m + 1][r] += r * A[m][r]
            if r:
                A[m + 1][r] += m * A[m - 1][r - 1]
    return A


def rising(a, r):
    v = 1
    for i in range(r):
        v *= a + i
    return v


class Family:
    def __init__(self, nmax=80):
        self.A = triangle(nmax)
        self.nmax = nmax

    def Qpoly(self, n, a=1):
        A = self.A[n]
        return sp.Poly(sum(rising(a, r) * A[r] * t**(r-1)
                           for r in range(1, len(A)) if A[r]), t)

    def Q(self, n, a=1):
        return self.Qpoly(n, a).as_expr()

    @staticmethod
    def E(p):
        """E=1+tD on a Poly or expr; returns Poly."""
        q = sp.Poly(p, t)
        return sp.Poly(q.as_expr() + t * q.diff().as_expr(), t)

    @staticmethod
    def T(p):
        """T=(tP')' on a Poly or expr; returns Poly."""
        q = sp.Poly(p, t)
        return sp.Poly(sp.expand(t * q.diff().as_expr()).diff(t), t)

    @staticmethod
    def W(p, q):
        """W_t(p,q)=p' q - p q' on Poly inputs; returns Poly."""
        p = sp.Poly(p, t); q = sp.Poly(q, t)
        return sp.Poly(sp.expand(p.diff().as_expr()*q.as_expr()
                                 - p.as_expr()*q.diff().as_expr()), t)

    def H(self, n, a=1):
        P = self.Qpoly(n, a)
        return self.W(P, self.T(P))
