"""CORRECTED: remainder of F and Delta modulo V."""
import sympy as sp
from dom2 import A, Q, B, x, sep

sep("remainders mod V")
for (m, a0) in [(6,1),(7,1),(8,1),(9,1),(10,1),(11,1),(7,2),(8,2),(9,2),(10,3)]:
    n = m - 2
    R = Q(A, n, a0, -x)
    V = Q(A, n, a0 + 1, -x)
    H = Q(A, n - 1, a0 + 1, -x)
    Y = sp.expand(2*H + x*sp.diff(H, x))
    F = sp.expand(x*sp.diff(V,x)**2 - a0*B(V,R))
    Delta = sp.expand(x*sp.diff(V,x)**2 + n*x*B(V,Y) - a0*B(V,R))
    rf = sp.rem(F, V, x)
    rd = sp.rem(Delta, V, x)
    print(f"\nm={m} a={a0}")
    print("  remF  =", sp.factor(rf))
    print("  remD  =", sp.factor(rd))
