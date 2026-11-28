"""Exact exploration of Q_{m,a} = B_{m,a}/t for associated Stirling numbers.

Uses only stdlib for A_{m,r}; sympy (already installed) for Sturm checks.
This is exploration, NOT a proof.
"""
from fractions import Fraction
import sympy as sp

t = sp.Symbol("t", real=True)

# --- associated Stirling numbers A_{m,r} = {m \brace r}_{>=2} ---
def astirling_table(Mmax):
    # A[m][r], 0 <= m <= Mmax, 0 <= r <= Mmax//2
    R = Mmax // 2 + 2
    A = [[0] * R for _ in range(Mmax + 1)]
    A[0][0] = 1
    for m in range(0, Mmax):
        for r in range(0, R - 1):
            term = r * A[m][r]
            if r >= 1:
                term += m * A[m - 1][r - 1]
            A[m + 1][r] += term
    return A

def rising(a, r):
    v = 1
    for k in range(r):
        v *= (a + k)
    return v

def Q_coeffs(A, m, a):
    """Q_{m,a}(t) = sum_{r>=1} (a)_r A_{m,r} t^{r-1}; list c[0..d]."""
    cs = []
    for r in range(1, len(A[0])):
        if A[m][r]:
            cs.append(rising(a, r) * A[m][r])
        else:
            cs.append(0)
    # strip trailing zeros
    while cs and cs[-1] == 0:
        cs.pop()
    return cs

def poly(cs):
    return sum(c * t ** i for i, c in enumerate(cs))

def W(f, g):
    return sp.expand(sp.diff(f, t) * g - f * sp.diff(g, t))

def positive_all_real(p):
    """Strict positivity of polynomial p on all of R (exact Sturm)."""
    p = sp.Poly(p, t)
    if p.degree() == 0:
        return p.LC() > 0, []
    nreal = p.count_roots(-sp.oo, sp.oo)
    if nreal == 0:
        # odd degree always has a real root, so degree is even and
        # the constant sign equals the sign of the leading coefficient.
        return p.LC() > 0, []
    return False, p.real_roots()

def main():
    Mmax = 22
    A = astirling_table(Mmax)
    # sanity: A_{4,2}=1, A_{5,2}=10, A_{6,3}=15 etc
    print("A[4][2]", A[4][2], "A[5][2]", A[5][2], "A[6][3]", A[6][3])
    print()

    # 1) verify recurrence (M): Q_{m+1,a} = (t Q_{m,a})' + m a t Q_{m-2+1,a+1}? 
    #    exact: Q_{m+1,a} = (tQ_{m,a})' + m*a*t*Q_{m-1,a+1}
    print("=== recurrence (M) residual check ===")
    bad = 0
    for m in range(2, 12):
        for a in range(1, 5):
            lhs = poly(Q_coeffs(A, m + 1, a))
            rhs = sp.diff(t * poly(Q_coeffs(A, m, a)), t) + m * a * t * poly(Q_coeffs(A, m - 1, a + 1))
            res = sp.expand(lhs - rhs)
            if res != 0:
                bad += 1
                print("BAD", m, a, res)
    print("bad count:", bad)
    print()

    # 2) distance-two coupled interlacing tests:
    #    p = Q_{m,a}, v = Q_{m-2,a+1}, r = Q_{m-2,a+2}
    #    check W(p,v)>0, W(p, U v)>0, W(p,r)>0
    print("=== distance-two coupled interlacing ===")
    fails = []
    U = lambda f: 2 * f + t * sp.diff(f, t)
    for m in range(5, 17):
        for a in range(1, 6):
            p = poly(Q_coeffs(A, m, a))
            v = poly(Q_coeffs(A, m - 2, a + 1))
            r = poly(Q_coeffs(A, m - 2, a + 2))
            if p == 0 or v == 0 or r == 0:
                continue
            Wp_v = sp.expand(W(p, v))
            Wp_Uv = sp.expand(W(p, U(v)))
            Wp_r = sp.expand(W(p, r))
            okv = positive_all_real(Wp_v)[0]
            okUv = positive_all_real(Wp_Uv)[0]
            okr = positive_all_real(Wp_r)[0]
            tag = []
            if not okv: tag.append("v")
            if not okUv: tag.append("Uv")
            if not okr: tag.append("r")
            if tag:
                fails.append((m, a, tuple(tag)))
                print("FAIL", m, "a", a, "which:", tag)
    print("total fails:", len(fails))
    print()

    # 3) show sample W(p,r) for a few cases
    print("=== sample W(p,r) ===")
    for (m, a) in [(7, 1), (8, 1), (9, 1), (10, 1), (8, 2), (9, 2), (10, 3)]:
        p = poly(Q_coeffs(A, m, a))
        r = poly(Q_coeffs(A, m - 2, a + 2))
        print(f"m={m} a={a}: W(p,r) = {sp.factor(W(p,r))}")

if __name__ == "__main__":
    main()
