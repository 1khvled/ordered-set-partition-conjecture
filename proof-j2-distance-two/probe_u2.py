"""Probe the identity D_{j+2} = [U^2 D_j + (2j-1) U D_j + j(j-1) D_j]/(j(j+1))
and positivity of W(D_j, U^2 D_j).
"""
import sympy as sp
from explore import astirling_table, Q_coeffs, poly, W, positive_all_real

t = sp.Symbol("t", real=True)
U = lambda f: sp.expand(2 * f + t * sp.diff(f, t))

def check_identity(A, m, a):
    """D_j = Q_{m,a}; D_{j+2} = Q_{m-2,a+2}. Identity with j=a."""
    p = poly(Q_coeffs(A, m, a))
    d2 = poly(Q_coeffs(A, m - 2, a + 2))
    expr = (U(U(p)) + (2 * a - 1) * U(p) + a * (a - 1) * p) / (a * (a + 1))
    return sp.simplify(d2 - expr)

def main():
    A = astirling_table(22)
    print("=== identity residual D_{j+2} vs U^2 combo ===")
    worst = 0
    for m in range(4, 16):
        for a in range(1, 6):
            res = check_identity(A, m, a)
            if res != 0:
                print("BAD", m, a, res)
                worst += 1
    print("bad:", worst)
    print()

    print("=== W(p, U p) and W(p, U^2 p) positivity ===")
    bad1 = bad2 = 0
    for m in range(4, 18):
        for a in range(1, 6):
            p = poly(Q_coeffs(A, m, a))
            if p == 0:
                continue
            w1 = sp.expand(W(p, U(p)))
            w2 = sp.expand(W(p, U(U(p))))
            ok1 = positive_all_real(w1)[0]
            ok2 = positive_all_real(w2)[0]
            if not ok1:
                bad1 += 1
                print("W(p,Up) FAIL", m, a)
            if not ok2:
                bad2 += 1
                print("W(p,U^2p) FAIL", m, a)
    print("bad W(p,Up):", bad1, " bad W(p,U^2p):", bad2)
    print()

    print("=== sample factorizations of W(p, U^2 p) ===")
    for (m, a) in [(6, 1), (7, 1), (8, 1), (8, 2), (10, 1), (10, 3)]:
        p = poly(Q_coeffs(A, m, a))
        print(f"m={m} a={a}: W(p,U^2p) = {sp.factor(W(p, U(U(p))))}")

if __name__ == "__main__":
    main()
