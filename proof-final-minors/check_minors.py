"""Exact small checks for proof-final-minors/NOTES.md.

Verifies identity (4) (lag-j coefficient-minor recurrence) on actual
O_{n,j} rows and reports the generalized tail T_{n,d,k} sign.
Pure Python integer arithmetic.
"""
from math import comb


def rows(j, nmax):
    O = [[1]]
    for n in range(1, nmax + 1):
        res = [0]
        for r in range(1, min(j - 1, n) + 1):
            c = comb(n, r)
            src = O[n - r]
            for k, v in enumerate(src):
                if k >= len(res):
                    res.extend([0] * (k - len(res) + 1))
                res[k] += c * v
        for r in range(j, n + 1):
            c = comb(n, r)
            src = O[n - r]
            for k, v in enumerate(src):
                k1 = k + 1
                if k1 >= len(res):
                    res.extend([0] * (k1 - len(res) + 1))
                res[k1] += c * v
        while len(res) > 1 and res[-1] == 0:
            res.pop()
        O.append(res)
    return O


def cf(a, k):
    return a[k] if 0 <= k < len(a) else 0


def delta(O, n, d, k):
    return cf(O[n], k + 1) * cf(O[n - d], k) - cf(O[n], k) * cf(O[n - d], k + 1)


def check_identity(jmax=6, nmax=24):
    for j in range(2, jmax + 1):
        O = rows(j, nmax)
        for n in range(1, nmax + 1):
            ccoef = comb(n - 1, j - 1) if n >= j else 0
            q = O[n - j] if n >= j else [0]
            u = O[n - 1]
            deg = len(O[n])
            for k in range(-1, deg + 1):
                lhs = delta(O, n, j, k) if n >= j else None
                if lhs is None:
                    continue
                rhs = (n + 2 * k) * delta(O, n - 1, j - 1, k) if n - 1 >= j - 1 and j - 1 >= 1 else 0
                rhs += 2 * cf(u, k + 1) * cf(q, k)
                for r in range(2, min(j - 1, n) + 1):
                    if j - r >= 1:
                        rhs += comb(n, r) * delta(O, n - r, j - r, k)
                if ccoef:
                    Ttail = (k + 1) * cf(q, k) * cf(q, k) - cf(q, k) * cf(q, k + 1) - k * cf(q, k - 1) * cf(q, k + 1)
                    rhs += ccoef * Ttail
                if lhs != rhs:
                    raise AssertionError((j, n, k, lhs, rhs))
    print(f"identity OK for j<= {jmax}, n<= {nmax}")


def check_general_tail(jmax=4, nmax=20):
    bad = []
    for j in range(2, jmax + 1):
        O = rows(j, nmax)
        for n in range(1, nmax + 1):
            ccoef = comb(n - 1, j - 1) if n >= j else 0
            if not ccoef:
                continue
            q = O[n - j]
            for d in range(1, n + 1):
                if d == j:
                    continue
                b = O[n - d]
                deg = max(len(q), len(b))
                for k in range(0, deg + 1):
                    T = (k + 1) * cf(b, k) * (cf(q, k) - cf(q, k + 1)) - k * cf(b, k + 1) * (cf(q, k - 1) - cf(q, k))
                    if T < 0:
                        bad.append((j, n, d, k, T))
    print(f"general tail negative count={len(bad)}")
    for row in bad[:10]:
        print(" ", row)


if __name__ == "__main__":
    check_identity()
    check_general_tail()
