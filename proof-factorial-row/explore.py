"""Exact exploration of E_r(x)=sum_k k! A_{r,k} x^k, A_{r,k}={r\brace k}_{>=2}.

No installs, no builds. Uses sympy/numpy already present.
"""
from fractions import Fraction
import sympy as sp

x = sp.symbols('x')

MAXR = 40


def table(N):
    """A[r][k] for 0<=r<=N, 0<=k<=floor(r/2)."""
    A = [[0] * (N // 2 + 2) for _ in range(N + 1)]
    A[0][0] = 1
    for r in range(1, N + 1):
        for k in range(1, r // 2 + 1):
            A[r][k] = k * A[r - 1][k] + (r - 1) * A[r - 2][k - 1]
    return A


def Epoly(A, r):
    """E_r as sympy Poly, E_0=1."""
    coeffs = [sp.Integer(0)] + [sp.factorial(k) * A[r][k] for k in range(1, r // 2 + 1)]
    return sp.Poly(sum(c * x ** i for i, c in enumerate(coeffs)), x)


A = table(MAXR)
E = [Epoly(A, r) for r in range(MAXR + 1)]

# verify recurrence exactly
for r in range(3, MAXR + 1):
    lhs = E[r]
    rhs = x * E[r - 1].diff(x) + (r - 1) * x * E[r - 2] + (r - 1) * x ** 2 * E[r - 2].diff(x)
    assert (lhs - rhs).is_zero, r
print('recurrence verified for r=3..%d' % MAXR)

for r in range(2, 11):
    print('E_%d =' % r, sp.factor(E[r].as_expr()))

# T_r = E_r/x, T_1=0
T = {1: sp.Poly(0, x)}
for r in range(2, MAXR + 1):
    T[r] = sp.Poly(sp.cancel(E[r].as_expr() / x), x)

for r in range(4, MAXR + 1):
    lhs = T[r]
    rhs = sp.Poly(sp.diff(x * T[r - 1].as_expr() + (r - 1) * x ** 2 * T[r - 2].as_expr(), x), x)
    assert (lhs - rhs).is_zero, r
print('T-recurrence T_r=(x T_{r-1}+(r-1)x^2 T_{r-2}) verified for r=4..%d' % MAXR)

# U_r = T_{r-1} + (r-1) x T_{r-2},  T_r = (x U_r)'
U = {}
for r in range(3, MAXR + 1):
    U[r] = sp.Poly(sp.expand(T[r - 1].as_expr() + (r - 1) * x * T[r - 2].as_expr()), x)
    assert sp.Poly(sp.diff(x * U[r].as_expr(), x), x) == T[r], r
print('U-recurrence verified for r=3..%d' % MAXR)
for r in range(3, 9):
    print('U_%d =' % r, U[r].as_expr())
