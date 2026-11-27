import numpy as np
from math import comb, factorial

def compute_O_coeffs(j, nmax):
    Os = [None] * (nmax + 1)
    Os[0] = [1]
    for n in range(1, nmax + 1):
        deg = n // j
        coeff = [0] * (deg + 1)
        for r in range(1, min(j, n + 1)):
            Onr = Os[n - r]
            c = comb(n, r)
            for k, val in enumerate(Onr):
                if k <= deg:
                    coeff[k] += c * val
        for r in range(j, n + 1):
            Onr = Os[n - r]
            c = comb(n, r)
            for k, val in enumerate(Onr):
                if k + 1 <= deg:
                    coeff[k + 1] += c * val
        Os[n] = coeff
    return Os

def Dval(j):
    return factorial(j) * ((j + 2) / (j + 1)) ** j

if __name__ == "__main__":
    j = 3
    for m in [10, 20, 40, 60]:
        n = j * m  # s=0
        Os = compute_O_coeffs(j, n)
        On = Os[n]
        cnm = On[m]
        Dm = Dval(j) * (m ** j)
        # Q_{m,s}(x) = sum_{ell=0}^m (On[m-ell]/cnm) * (-x/Dm)^ell
        # as polynomial in x: coeff of x^ell is (On[m-ell]/cnm) * (-1/Dm)^ell
        coeffs = []  # highest degree first for np.roots? np.roots takes [a_m,...,a_0]
        for ell in range(m, -1, -1):
            coeffs.append((On[m - ell] / cnm) * ((-1.0 / Dm) ** ell))
        roots = np.roots(coeffs)
        nonreal = [r for r in roots if abs(r.imag) > 1e-6 * max(1, abs(r.real))]
        print(f"j={j} m={m} n={n} nonreal={len(nonreal)} max|im|={max(abs(r.imag) for r in roots):.2e}")
        real_parts = sorted([r.real for r in roots])
        # show quantiles
        qs = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
        print("  quantiles:", {q: real_parts[int(q*(len(real_parts)-1))] for q in qs})
