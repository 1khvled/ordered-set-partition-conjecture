"""Explore the a-shift ladder: A_a = W(Q_{m,a},Q_{m-2,a}),
D_a = W(Q_{m,a},Q_{m-2,a+1}), C_a = W(Q_{m,a},UQ_{m-2,a+1}) (proved >0)."""
import sympy as sp
from explore import A, Q, W, t

t = sp.Symbol("t", real=True)
U = lambda f: sp.expand(2 * f + t * sp.diff(f, t))


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)


def posinfo(pol):
    pol = sp.Poly(pol, t)
    if pol.degree() == 0:
        return pol.LC() > 0, 0
    nr = pol.count_roots(-sp.oo, sp.oo)
    return pol.LC() > 0 and nr == 0, nr


sep("Verify D_a = (B_a + (a-1) A_a)/a exactly")
for m in range(6, 12):
    for a0 in range(1, 5):
        p = Q(A, m, a0, t)
        q0 = Q(A, m - 2, a0, t)
        q1 = Q(A, m - 2, a0 + 1, t)
        Aa = W(p, q0)
        Ba = W(p, U(q0))
        Da = W(p, q1)
        res = sp.simplify(Da - (Ba + (a0 - 1) * Aa) / a0)
        if res != 0:
            print("BAD", m, a0, res)
print("done")

sep("Positivity of A_a, B_a, C_a")
for m in range(6, 13):
    for a0 in range(1, 5):
        p = Q(A, m, a0, t)
        q0 = Q(A, m - 2, a0, t)
        q1 = Q(A, m - 2, a0 + 1, t)
        Aa = W(p, q0)
        Ba = W(p, U(q0))
        Ca = W(p, U(q1))
        Da = W(p, q1)
        okA, nrA = posinfo(Aa)
        okB, nrB = posinfo(Ba)
        okC, nrC = posinfo(Ca)
        okD, nrD = posinfo(Da)
        print(f"m={m} a={a0}: A={'OK' if okA else 'BAD'+str(nrA)} "
              f"B={'OK' if okB else 'BAD'+str(nrB)} "
              f"C={'OK' if okC else 'BAD'+str(nrC)} "
              f"D={'OK' if okD else 'BAD'+str(nrD)}")
