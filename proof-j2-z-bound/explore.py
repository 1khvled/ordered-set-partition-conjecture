"""Exact symbolic exploration of the Z-bound remainder mod V."""
import sympy as sp
from ops import a, x, A, R, Q, Bx

def show(n, av):
    Rn = R(n, av)
    Vn = R(n, av + 1)
    Hn = R(n - 1, av + 1)
    D = sp.expand((3 * av + 2) * Rn - x**2 * sp.diff(Rn, x, 2))
    rem = sp.rem(D, Vn, x)
    print(f"n={n}, a={av}: deg R={sp.degree(Rn,x)}, deg V={sp.degree(Vn,x)}")
    print("  rem mod V =", sp.factor(rem))
    # remainder of x B(V,R) - a R^2 mod V
    E = sp.expand(x * Bx(Vn, Rn) - av * Rn**2)
    rem2 = sp.rem(E, Vn, x)
    print("  xB(V,R)-aR^2 mod V =", sp.factor(rem2))
    print()

for n in range(4, 9):
    for av in [1, 2, 3]:
        show(n, av)
