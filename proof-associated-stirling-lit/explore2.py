"""Test real-rootedness of factorial-weighted 2-associated Stirling row polynomials
    F_n^{(c)}(x) = sum_k (k!)^c D(n,k) x^k   for c = 1, 2, 3.
Only c=1 is needed: Newton on c=1 gives SLC (even stronger).
"""
from fractions import Fraction
from math import factorial
from explore import table, count_negative_roots, count_positive_roots

if __name__ == "__main__":
    N = 30
    D = table(N)
    print("n  c1  c2  c3   (True = all d roots real & negative)")
    for n in range(2, N + 1):
        d = n // 2
        res = []
        for c in (1, 2, 3):
            coeffs = [Fraction(factorial(k) ** c * D[n][k]) for k in range(d + 1)]
            neg = count_negative_roots(coeffs)
            pos = count_positive_roots(coeffs)
            res.append(neg + pos == d)
        print(f"{n:2d}  {str(res[0]):5s} {str(res[1]):5s} {str(res[2]):5s}")
