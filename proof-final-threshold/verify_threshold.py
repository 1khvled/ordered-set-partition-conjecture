"""Exact symbolic verification of the threshold identities used in NOTES.md."""
import sympy as sp
from math import comb

z, t = sp.symbols('z t')
a, b, c, d = sp.symbols('a b c d')


def E(r):
    return sum(z**k / sp.factorial(k) for k in range(r + 1))


def B(j):
    return sp.Rational(1) if j <= 0 else sum(z**r / sp.factorial(r) for r in range(1, j)) + t * (
        sp.exp(z) - E(j - 1))


def F_series(j, nmax):
    D = sp.series(1 - B(j), z, 0, nmax + 1).removeO()
    return sp.series(1 / D, z, 0, nmax + 1).removeO()


print('T1: B_j\' = 1 + B_{j-1}')
for j in range(2, 8):
    assert sp.simplify(sp.diff(B(j), z) - (1 + B(j - 1))) == 0, j
print('  ok j=2..7')

print('T2: pole separation tau_j - tau_{j-1}')
for j in range(2, 7):
    Pp, P = E(j - 2), E(j - 1)
    Qp, Q = 2 - Pp, 2 - P
    tau_p, tau = Qp / (sp.exp(z) - Pp), Q / (sp.exp(z) - P)
    p = z**(j - 1) / sp.factorial(j - 1)
    claimed = -p * (sp.exp(z) - 2) / ((sp.exp(z) - P) * (sp.exp(z) - Pp))
    assert sp.simplify(tau - tau_p - claimed) == 0, j
print('  ok j=2..6')

print('T3: F_j = F_{j-1}/(1-(1-t)p F_{j-1})')
for j in range(2, 7):
    p = z**(j - 1) / sp.factorial(j - 1)
    nmax = 9
    d = 1 - B(j - 1)
    lhs = sp.series(1 / (d - (1 - t) * p), z, 0, nmax + 1).removeO()
    rhs = sp.series((1 / d) / (1 - (1 - t) * p / d), z, 0, nmax + 1).removeO()
    assert sp.simplify(sp.expand(lhs - rhs)) == 0, j
print('  ok j=2..6')

print('T4: coefficient identity Q_n = P_n + (1-t)binom(n,j-1)sum_m binom(n-j+1,m) Q_m P_{n-j+1-m}')


def O(nmax, j):
    arr = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        v = sum(comb(n, r) * arr[n - r] for r in range(1, min(j - 1, n) + 1))
        v += sum(comb(n, r) * t * arr[n - r] for r in range(j, n + 1))
        arr.append(sp.expand(v))
    return arr


for j in range(2, 7):
    P, Q = O(14, j - 1), O(14, j)
    for n in range(0, 15):
        rhs = P[n]
        if n >= j - 1:
            rhs += (1 - t) * comb(n, j - 1) * sum(
                comb(n - j + 1, m) * Q[m] * P[n - j + 1 - m] for m in range(n - j + 2))
        assert sp.expand(rhs - Q[n]) == 0, (j, n)
print('  ok j=2..6, n<=14')

print('T5: j=2 outer convolution O_n = n! sum_k H_k/k!')


def H(nmax, j):
    h = [sp.Integer(1)]
    for n in range(1, nmax + 1):
        if n < j:
            h.append(sp.Integer(0)); continue
        v = sum(comb(n, k) * h[n - k] for k in range(1, min(j - 1, n) + 1))
        v += sum(comb(n, q) * t * h[n - q] for q in range(j, n + 1))
        h.append(sp.expand(v))
    return h


Q = O(14, 2)
h = H(14, 2)
for n in range(15):
    assert sp.expand(sum(sp.factorial(n) / sp.factorial(k) * h[k] for k in range(n + 1)) - Q[n]) == 0, n
print('  ok n<=14')

print('ALL EXACT CHECKS PASSED')
