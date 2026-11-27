from math import comb, factorial
from flint import acb, acb_poly, ctx

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
    ctx.dps = 80
    j = 3
    for m in [10, 20, 40, 60]:
        n = j * m
        Os = compute_O_coeffs(j, n)
        On = Os[n]
        cnm = On[m]
        Dm = Dval(j) * (m ** j)
        # Q in x: coeff of x^ell is (On[m-ell]/cnm) * (-1/Dm)^ell
        coeffs = []
        for ell in range(m + 1):
            coeffs.append(acb(On[m - ell] / cnm) * ((-acb(1) / Dm) ** ell))
        poly = acb_poly(coeffs)
        roots = poly.roots()
        nonreal = 0
        max_im = 0.0
        xs = []
        for r in roots:
            re = float(r.real.mid())
            im = float(r.imag.mid())
            rad = float(abs(r.real.rad()) + abs(r.imag.rad()))
            if abs(im) > 1e-20 or rad > 1e-20:
                # check enclosure
                if abs(r.imag.mid()) > float(r.imag.rad()) + 1e-30:
                    nonreal += 1
            max_im = max(max_im, abs(im))
            xs.append(re)
        xs.sort()
        print(f"j={j} m={m} n={n} nonreal={nonreal} max|im|={max_im:.2e}")
        qs = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
        print("  quantiles:", {q: xs[int(q * (len(xs) - 1))] for q in qs})
