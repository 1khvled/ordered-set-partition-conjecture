"""Identify the square in F = x(V')^2 - a B(V,R) = c x G^2."""
import sympy as sp
from explore import A, Q, B, t

x = sp.Symbol("x", positive=True)


def sep(s):
    print("\n" + "=" * 70)
    print(s)
    print("=" * 70)

sep("Identify G in F = x(V')^2 - aB(V,R) = c x G^2")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(10,1),(11,1),(6,2),(7,2),(8,2),(9,2),(10,3)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    F = sp.expand(x * sp.diff(V, x)**2 - a0 * B(V, R))
    # square-free decomposition
    sf = sp.sqf_list(F)
    print(f"\nm={m} a={a0}: F = {sp.factor(F)}")
    print("   sqf_list =", [(sp.factor(f), e) for f, e in sf[1]])
    # try candidate polynomials
    cands = {}
    for (mm, aa, label) in [
        (n, a0+2, f"Q_{n},{a0+2}"),
        (n-1, a0+2, f"Q_{n-1},{a0+2}"),
        (n+1, a0, f"Q_{n+1},{a0}"),
        (n-1, a0+1, f"Q_{n-1},{a0+1}=H"),
        (n+1, a0+1, f"Q_{n+1},{a0+1}"),
    ]:
        cands[label] = Q(A, mm, aa, -x)
    for label, cp in cands.items():
        if cp == 0:
            continue
        # check if F/(x cp^2) is a positive constant
        rat = sp.cancel(F / (x * cp**2))
        if rat.is_number and rat > 0:
            print(f"   *** F = ({rat}) * x * [{label}(-x)]^2")
