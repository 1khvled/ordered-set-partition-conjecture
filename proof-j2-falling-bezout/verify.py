"""Exact finite audits for proof-j2-falling-bezout.

Audits only; no finite check is used as a proof.
"""
from fractions import Fraction
import math
from functools import lru_cache


@lru_cache(None)
def A(m, k):
    # 2-associated Stirling numbers of the second kind
    if k == 0:
        return 1 if m == 0 else 0
    if k < 0 or m < 2 * k:
        return 0
    return k * A(m - 1, k) + (m - 1) * A(m - 2, k - 1)


@lru_cache(None)
def rising(a, k):
    r = 1
    for i in range(k):
        r *= a + i
    return r


@lru_cache(None)
def d(n, k):
    # d_{n,k} = [t^k] O_{n,2}(1,t)
    if k == 0:
        return Fraction(math.factorial(n), 1)
    s = Fraction(0, 1)
    for m in range(2 * k, n + 1):
        s += Fraction(rising(n - m + 1, k), math.factorial(m)) * A(m, k)
    return Fraction(math.factorial(n), 1) * s


@lru_cache(None)
def S2(n, k):
    if k < 0 or k > n:
        return 0
    if n == 0 and k == 0:
        return 1
    res = 0
    for j in range(k + 1):
        res += ((-1) ** (k - j)) * math.comb(k, j) * (j ** n)
    return Fraction(res, math.factorial(k))


@lru_cache(None)
def c_coeff(n, i):
    s = Fraction(0, 1)
    for k in range(i, n // 2 + 1):
        s += S2(k, i) * d(n, k)
    return s


def minor_2x2(mat, r1, r2, c1, c2):
    return mat(r1, c1) * mat(r2, c2) - mat(r1, c2) * mat(r2, c1)


def check_cauchy_binet(N=10):
    # verify identity (3): det C = sum det D det S
    for n1 in range(N):
        for n2 in range(n1 + 1, N + 1):
            maxk = n2 // 2
            maxi = maxk
            for i1 in range(maxi + 1):
                for i2 in range(i1 + 1, maxi + 1):
                    lhs = minor_2x2(c_coeff, n1, n2, i1, i2)
                    rhs = Fraction(0, 1)
                    for k1 in range(maxk + 1):
                        for k2 in range(k1 + 1, maxk + 1):
                            rhs += minor_2x2(d, n1, n2, k1, k2) * minor_2x2(S2, k1, k2, i1, i2)
                    if lhs != rhs:
                        raise AssertionError(("CB", n1, n2, i1, i2, lhs, rhs))
    print("Cauchy-Binet identity (3): exact residual 0 on checked range")


def check_tp2(N=12):
    for n1 in range(N):
        for n2 in range(n1 + 1, N + 1):
            maxk = max(n1 // 2, n2 // 2)
            for i1 in range(maxk + 1):
                for i2 in range(i1 + 1, maxk + 1):
                    det = minor_2x2(c_coeff, n1, n2, i1, i2)
                    if det < 0:
                        raise AssertionError(("TP2", n1, n2, i1, i2, det))
    print("c TP2 minors: all >=0 on checked range")


def check_samples():
    # Direction audit samples from NOTES
    n, i = 4, 1
    det = c_coeff(n, i) * c_coeff(n + 1, i + 1) - c_coeff(n, i + 1) * c_coeff(n + 1, i)
    assert (c_coeff(4, 1), c_coeff(4, 2), c_coeff(5, 1), c_coeff(5, 2)) == (51, 6, 421, 110)
    assert det == 3084
    det2 = c_coeff(3, 0) * c_coeff(4, 1) - c_coeff(3, 1) * c_coeff(4, 0)
    assert det2 == 138
    print("Direction samples: OK", det, det2)


def check_obstruction():
    # j2-base 5.4 witness
    # P=10+11t+t^2, Q=10+22t+100t^2
    # falling factorial: (t)_0=1,(t)_1=t,(t)_2=t(t-1)=t^2-t
    # P = 10 + 12*t + 1*(t)_2? 10+11t+t^2 = a0+a1*t+a2*(t^2-t)=a0 + (a1-a2)*t + a2*t^2
    # a2=1, a1-a2=11 => a1=12, a0=10
    P = [Fraction(10), Fraction(12), Fraction(1)]
    Q = [Fraction(10), Fraction(122), Fraction(100)]
    r1 = Q[1] / P[1]
    r2 = Q[2] / P[2]
    assert r1 < r2
    discQ = Fraction(22) ** 2 - 4 * Fraction(10) * Fraction(100)
    assert discQ == -3516
    print("Obstruction witness: ratios increasing", r1, r2, "Disc(Q)", discQ)


if __name__ == "__main__":
    check_samples()
    check_obstruction()
    check_cauchy_binet(N=9)
    check_tp2(N=11)
    print("ALL CHECKS PASSED")
