"""Literature-support exploration.

Tests (exactly, with rational arithmetic) whether
    E_n(x) = sum_k k! * D(n,k) * x^k
is real-rooted for the 2-associated Stirling numbers D(n,k) = {n\brace k}_{>=2}.
If E_n is real-rooted, Newton's inequalities applied to its coefficients give

    (k! D_k)^2 >= (k+1)/k * (d-k+1)/(d-k) * ((k-1)! D_{k-1}) ((k+1)! D_{k+1})

which is STRONGER than the target SLC  k^2 D_k^2 >= (k+1)^2 D_{k-1}D_{k+1}.
"""
from fractions import Fraction
from math import factorial


def table(N):
    D = [[0] * (N + 2) for _ in range(N + 2)]
    D[0][0] = 1
    for n in range(1, N + 1):
        for k in range(1, n // 2 + 1):
            D[n][k] = k * D[n - 1][k] + (n - 1) * D[n - 2][k - 1]
    return D


# ---- polynomial helpers, ascending coefficient lists over Fractions ----
def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def deriv(p):
    return trim([p[i] * i for i in range(1, len(p))])


def divmod_poly(a, b):
    a = list(a)
    b = trim(list(b))
    assert b
    q = [Fraction(0)] * max(0, len(a) - len(b) + 1)
    while len(a) >= len(b) and trim(a):
        shift = len(a) - len(b)
        c = a[-1] / b[-1]
        q[shift] = c
        for i, bi in enumerate(b):
            a[i + shift] -= c * bi
        trim(a)
    return trim(q), trim(a)


def sturm(p):
    p = trim([Fraction(x) for x in p])
    s0 = p
    s1 = deriv(p)
    seq = [s0, s1]
    while len(seq[-1]) > 1:
        q, r = divmod_poly(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-x for x in r])
    return seq


def signs_at(seq, x):
    out = []
    for p in seq:
        v = sum(c * x**i for i, c in enumerate(p))
        out.append(0 if v == 0 else (1 if v > 0 else -1))
    return out


def variations(signs):
    s = [v for v in signs if v != 0]
    return sum(1 for i in range(1, len(s)) if s[i] != s[i - 1])


def count_negative_roots(p):
    seq = sturm(p)
    # Sturm: V(-infinity) - V(0-) counts roots in (-inf,0); handle p(0)=0 separately
    v_minus = variations(signs_at(seq, Fraction(-10**60)))
    v_zero = variations(signs_at(seq, Fraction(0)))
    return v_minus - v_zero


def count_positive_roots(p):
    seq = sturm(p)
    v_zero = variations(signs_at(seq, Fraction(0)))
    v_plus = variations(signs_at(seq, Fraction(10**60)))
    return v_zero - v_plus


def main():
    N = 40
    D = table(N)
    print("n  deg  degree_expected  neg_roots  pos_roots  real_rooted?  newton_for_k!D?")
    for n in range(2, N + 1):
        d = n // 2
        coeffs = [Fraction(factorial(k) * D[n][k]) for k in range(d + 1)]
        neg = count_negative_roots(coeffs)
        pos = count_positive_roots(coeffs)
        c = coeffs
        newton = True
        for k in range(1, d):
            if k >= d:
                continue
            lhs = c[k] ** 2
            rhs = Fraction(k + 1, k) * Fraction(d - k + 1, d - k) * c[k - 1] * c[k + 1]
            if lhs < rhs:
                newton = False
        print(f"{n:2d} {len(coeffs)-1:3d} {d:3d} {neg:9d} {pos:9d} "
              f"{str(neg + pos == d):>13s} {str(newton):>17s}")


if __name__ == "__main__":
    main()
