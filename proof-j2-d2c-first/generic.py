"""Generic symbolic expansion of the first-premise Wronskian B(P,V)."""
import sympy as sp

a, n, x = sp.symbols("a n x", positive=True)
R0, R1, R2, R3, R4 = sp.symbols("R0 R1 R2 R3 R4")
H0, H1, H2, H3 = sp.symbols("H0 H1 H2 H3")

# V = ((a+1) R + x R')/a
V0 = ((a + 1) * R0 + x * R1) / a
V1 = ((a + 2) * R1 + x * R2) / a
V2 = ((a + 3) * R2 + x * R3) / a
V3 = ((a + 4) * R3 + x * R4) / a

# S = a(V - R - n x H)
S0 = sp.expand(a * (V0 - R0 - n * x * H0))
S1 = sp.expand(a * (V1 - R1 - n * (H0 + x * H1)))
S2 = sp.expand(a * (V2 - R2 - n * (2 * H1 + x * H2)))
S3 = sp.expand(a * (V3 - R3 - n * (3 * H2 + x * H3)))

# P = S + x S' - (n+1) a x V
P0 = sp.expand(S0 + x * S1 - (n + 1) * a * x * V0)
P1 = sp.expand(2 * S1 + x * S2 - (n + 1) * a * (V0 + x * V1))
P2 = sp.expand(3 * S2 + x * S3 - (n + 1) * a * (2 * V1 + x * V2))

BPV = sp.expand(P0 * V1 - P1 * V0)
print("B(P,V) expanded (collected in H):")
print(sp.collect(sp.expand(BPV), [H0, H1, H2, H3]))
print()
print("number of terms:", len(sp.Add.make_args(BPV)))
