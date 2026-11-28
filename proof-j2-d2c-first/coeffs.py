"""Test coefficient-level positivity of W(Q_{m,a}, Q_{m-2,a+1})."""
import sympy as sp
from explore import A, Q, W, t

t = sp.Symbol("t", real=True)


def sep(s):
    print("\n" + "=" * 68)
    print(s)
    print("=" * 68)

sep("1. coefficient pair condition (a+s)A[m][r]A[m-2][s] - (a+r)A[m][s]A[m-2][r] >= 0 for r>s")
bad = 0
total = 0
for m in range(6, 20):
    for a0 in range(1, 6):
        for r in range(2, m // 2 + 1):
            for s in range(1, r):
                if A[m][r] == 0 or A[m - 2][s] == 0:
                    continue
                lhs = (a0 + s) * A[m][r] * A[m - 2][s]
                rhs = (a0 + r) * A[m][s] * A[m - 2][r]
                total += 1
                if lhs < rhs:
                    bad += 1
                    if bad < 15:
                        print("BAD", m, a0, r, s, lhs, rhs)
print(f"total={total} bad={bad}")

sep("2. Is A[m][r]/((a+r) A[m-2][r]) increasing in r?")
for m in [8, 10, 12]:
    for a0 in [1, 2]:
        vals = []
        for r in range(1, m // 2 + 1):
            if A[m][r] and A[m - 2][r]:
                vals.append((r, sp.Rational(A[m][r], (a0 + r) * A[m - 2][r])))
        print(f"m={m} a={a0}: {[(r, float(v)) for r,v in vals]}")

sep("3. Direct small W coefficient lists")
for (m, a0) in [(8, 1), (10, 1), (12, 1)]:
    p = Q(A, m, a0, t)
    w = Q(A, m - 2, a0 + 1, t)
    pol = sp.Poly(W(p, w), t)
    print(f"m={m} a={a0}: coeffs(high->low)={pol.all_coeffs()}")
