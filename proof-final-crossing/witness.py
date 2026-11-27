# Exact witness for the crossing obstruction in proof-final-crossing/NOTES.md.
# No dependencies; pure integer/rational arithmetic.
from fractions import Fraction

# j=2 large edges:
# E1: (2,0) -> (12,1), size 10; E2: (4,0) -> (8,1), size 4.
# m-coordinate linear: E1 m=(n-2)/10, E2 m=(n-4)/4.
# Solve equality.
n = Fraction(4*2 - 10*4, 4 - 10)  # solve (n-2)/10=(n-4)/4 -> 4(n-2)=10(n-4)
# Actually solve properly:
# 4(n-2)=10(n-4) -> 4n-8=10n-40 -> 6n=32 -> n=16/3
n = Fraction(16, 3)
m = (n - 2) / 10
print("crossing n =", n, "m =", m)
assert m == Fraction(1, 3)
print("m is non-integer:", m.denominator != 1)
print("Obstruction witness verified.")
