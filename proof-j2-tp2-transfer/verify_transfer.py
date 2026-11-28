"""Exact audits for the all-n TP2 proof in NOTES.md.

Pure Python integer/Fraction arithmetic.  Finite loops are audits and finite
evidence only; the all-index proof is in NOTES.md.
"""

from fractions import Fraction
from math import factorial


def rising(a, k):
    out = 1
    for j in range(k):
        out *= a + j
    return out


def associated_stirling(N):
    A = [[0] * (N + 1) for _ in range(N + 1)]
    A[0][0] = 1
    for m in range(1, N + 1):
        for k in range(1, N // 2 + 1):
            A[m][k] = k * A[m - 1][k]
            if m >= 2:
                A[m][k] += (m - 1) * A[m - 2][k - 1]
    return A


def d_matrix(N, A):
    d = [[0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        d[n][0] = factorial(n)
        for k in range(1, N + 1):
            value = factorial(n) * sum(
                Fraction(rising(n - m + 1, k), factorial(m)) * A[m][k]
                for m in range(2, n + 1)
            )
            assert value.denominator == 1
            d[n][k] = int(value)
    return d


def direct_f_coefficients(N):
    # f[n][k] = [z^n t^k] 1/(1-z-tC(z)).
    f = [[Fraction(0) for _ in range(N + 1)] for _ in range(N + 1)]
    f[0][0] = Fraction(1)
    for n in range(1, N + 1):
        for k in range(N + 1):
            value = f[n - 1][k]
            if k:
                value += sum(
                    f[n - m][k - 1] / factorial(m)
                    for m in range(2, n + 1)
                )
            f[n][k] = value
    return f


def h(r):
    # [z^r] (1/(1-z)) * (C(z)/z), r >= 1; h(0)=0.
    if r < 1:
        return Fraction(0)
    return sum((Fraction(1, factorial(q + 1)) for q in range(1, r + 1)), Fraction(0))


def production_Q(N):
    Q = [[Fraction(0) for _ in range(N + 1)] for _ in range(N + 1)]
    for n in range(N + 1):
        Q[n][0] = Fraction(1)
        for m in range(1, n + 1):
            Q[n][m] = h(n - m)
    return Q


def matmul(A, B):
    size = len(A)
    return [
        [sum((A[i][k] * B[k][j] for k in range(size)), Fraction(0)) for j in range(size)]
        for i in range(size)
    ]


def block_identity_Q(Q, m):
    # diag(I_m, Q), with the trailing Q reindexed from zero.
    size = len(Q)
    D = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if i < m or j < m:
                D[i][j] = Fraction(int(i == j))
            else:
                D[i][j] = Q[i - m][j - m]
    return D


def left_product(Q):
    size = len(Q)
    P = [row[:] for row in Q]
    for m in range(1, size):
        P = matmul(P, block_identity_Q(Q, m))
    return P


def det(A):
    n = len(A)
    if n == 0:
        return Fraction(1)
    if n == 1:
        return A[0][0]
    total = Fraction(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        total += (-1 if j % 2 else 1) * A[0][j] * det(minor)
    return total


def minor2(A, i, j, p, q):
    return A[i][p] * A[j][q] - A[i][q] * A[j][p]


def main():
    N = 14
    A = associated_stirling(N)
    d = d_matrix(N, A)
    f = direct_f_coefficients(N)
    b = [[Fraction(d[n][k], factorial(n)) for k in range(N + 1)] for n in range(N + 1)]
    Q = production_Q(N)

    # Independent coefficient identity and exact left-production factorization.
    for n in range(N + 1):
        for k in range(N + 1):
            assert factorial(n) * f[n][k] == d[n][k]
    P = left_product(Q)
    assert P == b
    print("coefficient identity (13): zero residuals through n=14 (finite audit)")
    print("left-production factorization: zero residuals through n=14 (finite audit)")

    print("d[3] =", d[3][:3])
    print("d[4] =", d[4][:3])
    print("d[5] =", d[5][:3])

    # Direction samples.
    delta_d = d[4][1] * d[5][2] - d[4][2] * d[5][1]
    delta_b = minor2(b, 4, 5, 1, 2)
    zero_delta_d = d[3][0] * d[4][1] - d[3][1] * d[4][0]
    zero_delta_b = minor2(b, 3, 4, 0, 1)
    print("sample active Delta_d(n=4,k=1):", delta_d)
    print("sample active Delta_b(n=4,k=1):", delta_b)
    print("sample k=0 Delta_d(n=3):", zero_delta_d)
    print("sample k=0 Delta_b(n=3):", zero_delta_b)
    assert delta_d == 3084 and delta_b == Fraction(257, 240)
    assert zero_delta_d == 102 and zero_delta_b == Fraction(17, 24)

    q_boundary = minor2(Q, 2, 4, 1, 3)
    q_positive = minor2(Q, 3, 4, 1, 2)
    print("Q boundary minor rows(2,4), cols(1,3):", q_boundary)
    print("Q positive minor rows(3,4), cols(1,2):", q_positive)
    assert q_boundary == Fraction(1, 4)
    assert q_positive == Fraction(13, 144)

    # Exact monotonicity audit for h_{r+1}/h_r.
    for r in range(1, N - 1):
        difference = h(r + 2) / h(r + 1) - h(r + 1) / h(r)
        assert difference < 0
    print("h-ratio finite differences: strictly negative (finite audit)")

    # Infeasible witness: arbitrary increasing h could violate TP2, but these
    # values violate the proved decreasing-ratio constraint.
    bad_h = (Fraction(1), Fraction(2), Fraction(5))
    bad_minor = bad_h[1] * bad_h[1] - bad_h[0] * bad_h[2]
    bad_ratio_difference = bad_h[2] / bad_h[1] - bad_h[1] / bad_h[0]
    print("infeasible h=(1,2,5) would give minor:", bad_minor)
    print("its ratio difference (must be < 0):", bad_ratio_difference)
    assert bad_minor < 0 and bad_ratio_difference > 0

    # Every 2x2 minor of Q and b on the truncation.
    q_min = None
    b_min = None
    for i in range(N + 1):
        for j in range(i + 1, N + 1):
            for p in range(N + 1):
                for qcol in range(p + 1, N + 1):
                    value_q = minor2(Q, i, j, p, qcol)
                    value_b = minor2(b, i, j, p, qcol)
                    assert value_q >= 0
                    assert value_b >= 0
                    q_min = value_q if q_min is None or value_q < q_min else q_min
                    b_min = value_b if b_min is None or value_b < b_min else b_min
    print("minimum sampled Q 2x2 minor:", q_min)
    print("minimum sampled b 2x2 minor:", b_min)

    # Exact obstruction to upgrading Q from TP2 to TN.
    q4 = [row[:4] for row in Q[2:6]]
    q4_det = det(q4)
    print("Q rows(2,3,4,5), cols(0,1,2,3) determinant:", q4_det)
    assert q4_det == Fraction(-1, 4320)

    print("All exact audits passed; finite loops are not used as proof.")


if __name__ == "__main__":
    main()
