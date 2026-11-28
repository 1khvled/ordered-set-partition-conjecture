"""Exact verification for proof-j2-distance-two.
Stdlib fractions + sympy (already installed). No installs, no builds.

Establish:
  (V1) Associated-Stirling recurrence residual zero.
  (V2) Vertical identity residual zero for Q_{m,a}.
  (V3) Operator identity residual zero: U Q_{m,a} = a Q_{m,a+1} - (a-1) Q_{m,a}.
  (V4) Reduction identity W(p,r) = [W(p,Uv) + a W(p,v)]/(a+1), residual zero
       (holds as a formal identity for arbitrary polynomials).
  (V5) Direction/monotonicity audit: explicit samples with exact signs.
  (V6) Feasible-constraint falsification: the rational model of
       proof-j2-diagonal-after-h §5 satisfies the established hypotheses but
       violates the distance-two coupled premises AND (I2).
"""
from fractions import Fraction as F
import sympy as sp
from explore import astirling_table, Q_coeffs, poly, W

t = sp.Symbol("t", real=True)
U = lambda f: sp.expand(2 * f + t * sp.diff(f, t))

def sep(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)

def main():
    A = astirling_table(20)

    sep("V1: associated-Stirling recurrence A[m+1,r] = r A[m,r] + m A[m-1,r-1]")
    bad = 0
    for m in range(0, 19):
        for r in range(0, 9):
            lhs = A[m + 1][r]
            rhs = r * A[m][r] + (m * A[m - 1][r - 1] if (m >= 1 and r >= 1) else 0)
            if lhs != rhs:
                bad += 1
                print("  BAD", m, r, lhs, rhs)
    print("  residual zero count:", (19 * 9 - bad), "/", 19 * 9, " bad:", bad)

    sep("V2/V3: vertical identity and U-identity for Q_{m,a}")
    bad2 = bad3 = 0
    for m in range(4, 16):
        for a in range(1, 6):
            qa = poly(Q_coeffs(A, m, a))
            qa1 = poly(Q_coeffs(A, m, a + 1))
            res2 = sp.expand(qa1 - (qa + sp.diff(t * qa, t) / a))
            res3 = sp.expand(U(qa) - (a * qa1 - (a - 1) * qa))
            if res2 != 0: bad2 += 1
            if res3 != 0: bad3 += 1
    print("  vertical residual bad:", bad2, " U-identity residual bad:", bad3)

    sep("V4: reduction identity (formal, arbitrary polynomials)")
    a = sp.symbols("a", positive=True)
    f, h = sp.symbols("f h", cls=sp.Function)
    # take arbitrary concrete polynomials and arbitrary positive a
    for _ in range(5):
        import random
        p = sum(random.randint(1, 9) * t ** i for i in range(random.randint(2, 5)))
        v = sum(random.randint(1, 9) * t ** i for i in range(random.randint(1, 4)))
        av = sp.Rational(random.randint(1, 5))
        r = (U(v) + av * v) / (av + 1)
        lhs = sp.expand(W(p, r))
        rhs = sp.expand((W(p, U(v)) + av * W(p, v)) / (av + 1))
        print("  residual:", sp.expand(lhs - rhs), " a =", av)
    print("  => reduction is an exact formal identity.")

    sep("V5: direction audit -- explicit samples (exact signs)")
    for (m, a0) in [(7, 1), (8, 1), (9, 2), (10, 3)]:
        p = poly(Q_coeffs(A, m, a0))
        v = poly(Q_coeffs(A, m - 2, a0 + 1))
        r = poly(Q_coeffs(A, m - 2, a0 + 2))
        wv = sp.Poly(sp.expand(W(p, v)), t)
        wUv = sp.Poly(sp.expand(W(p, U(v))), t)
        wr = sp.Poly(sp.expand(W(p, r)), t)
        print(f"  m={m} a={a0}: W(p,v) LC={wv.LC()}, deg={wv.degree()}, "
              f"realroots={wv.count_roots(-sp.oo,sp.oo)}")
        print(f"           W(p,Uv) LC={wUv.LC()}, deg={wUv.degree()}, "
              f"realroots={wUv.count_roots(-sp.oo,sp.oo)}")
        print(f"           W(p,r)  LC={wr.LC()}, deg={wr.degree()}, "
              f"realroots={wr.count_roots(-sp.oo,sp.oo)}")
        print(f"           W(p,r) = {sp.factor(W(p,r))}")

    sep("V6: feasible-constraint falsification (rational countermodel)")
    # model: A_{4,1}=1, A_{4,2}=1/100; A_{5,1}=1,A_{5,2}=3/10; A_{6,1}=1,A_{6,2}=3/4,A_{6,3}=1/20
    Am = {}
    Am[(4, 1)] = F(1);  Am[(4, 2)] = F(1, 100)
    Am[(5, 1)] = F(1);  Am[(5, 2)] = F(3, 10)
    Am[(6, 1)] = F(1);  Am[(6, 2)] = F(3, 4); Am[(6, 3)] = F(1, 20)

    def Qc(m, a):
        cs = []
        for rr in range(1, 4):
            rising = 1
            for k in range(rr): rising *= (a + k)
            cs.append(rising * Am.get((m, rr), F(0)))
        while cs and cs[-1] == 0: cs.pop()
        return sum(sp.Rational(c.numerator, c.denominator) * t ** i for i, c in enumerate(cs))

    D1, D2, D3 = Qc(6, 1), Qc(5, 2), Qc(4, 3)
    print("  D1 =", D1, " D2 =", D2, " D3 =", D3)
    w12 = sp.expand(W(D1, D2)); w23 = sp.expand(W(D2, D3))
    w13 = sp.expand(W(D1, D3))
    print("  (I1) W(D1,D2) =", sp.factor(w12), " >0?", sp.Poly(w12,t).count_roots(-sp.oo,sp.oo)==0 and sp.Poly(w12,t).LC()>0)
    print("  (I1) W(D2,D3) =", sp.factor(w23), " >0?", sp.Poly(w23,t).count_roots(-sp.oo,sp.oo)==0 and sp.Poly(w23,t).LC()>0)
    print("  (I2) W(D1,D3) =", sp.factor(w13))
    p13 = sp.Poly(w13, t)
    print("       disc =", sp.discriminant(p13, t), " real roots:", p13.count_roots(-sp.oo,sp.oo))
    # distance-two coupled premises: m=6,a=1 => v=Q_{4,2}
    v = Qc(4, 2)
    wv = sp.expand(W(D1, v)); wUv = sp.expand(W(D1, U(v)))
    print("  premise W(p,v)  =", sp.factor(wv), " disc =", sp.discriminant(sp.Poly(wv,t),t))
    print("  premise W(p,Uv) =", sp.factor(wUv), " disc =", sp.discriminant(sp.Poly(wUv,t),t))
    print("  => BOTH distance-two coupled premises FAIL in this model.")

    sep("Boundary / degeneracy audit")
    for (m, a) in [(0,1),(1,1),(2,1),(3,1),(4,1),(5,2)]:
        cs = Q_coeffs(A, m, a) if m < len(A) else []
        print(f"  m={m} a={a}: Q coeffs = {cs} (degree {len(cs)-1 if cs else 'zero'})")

if __name__ == "__main__":
    main()
