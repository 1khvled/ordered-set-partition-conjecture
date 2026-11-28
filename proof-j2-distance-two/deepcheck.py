"""Deep finite check of (I2) and the two premises of its exact reduction.
Uses exact integer coefficients; Sturm root counts via sympy.
Finite evidence only -- NOT a proof.
"""
import sys
import sympy as sp
from explore import astirling_table, Q_coeffs, poly, W

t = sp.Symbol("t", real=True)
U = lambda f: sp.expand(2 * f + t * sp.diff(f, t))

def has_no_real_roots(p):
    P = sp.Poly(sp.expand(p), t)
    if P.degree() == 0:
        return True
    return P.count_roots(-sp.oo, sp.oo) == 0

def main():
    Mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 34
    A = astirling_table(Mmax)
    print(f"Mmax={Mmax}")
    fails = []
    counts = {"v": 0, "Uv": 0, "r": 0}
    for m in range(6, Mmax):
        for a in range(1, max(1, m - 1)):
            p = poly(Q_coeffs(A, m, a))
            v = poly(Q_coeffs(A, m - 2, a + 1))
            r = poly(Q_coeffs(A, m - 2, a + 2))
            if p == 0 or v == 0 or r == 0:
                continue
            wv = sp.expand(W(p, v))
            wUv = sp.expand(W(p, U(v)))
            wr = sp.expand(W(p, r))
            o1 = has_no_real_roots(wv)
            o2 = has_no_real_roots(wUv)
            o3 = has_no_real_roots(wr)
            counts["v"] += o1
            counts["Uv"] += o2
            counts["r"] += o3
            if not (o1 and o2 and o3):
                fails.append((m, a, o1, o2, o3))
                print("FAIL", m, a, o1, o2, o3)
    print("counts (no-real-root):", counts)
    print("total fails:", len(fails))
    return fails

if __name__ == "__main__":
    main()
