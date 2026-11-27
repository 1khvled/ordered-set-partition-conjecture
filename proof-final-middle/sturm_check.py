from math import comb
from flint import fmpz_poly, fmpq_poly

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

def sturm_count(poly):
    # Standard Sturm sequence for squarefree polynomial over Q
    if isinstance(poly, fmpz_poly):
        poly = fmpq_poly(list(poly))
    seq = [poly, poly.derivative()]
    while True:
        a, b = seq[-2], seq[-1]
        if b.degree() < 0:
            break
        q, r = divmod(a, b)
        r = -r
        if r.degree() < 0 and r[0] == 0:
            break
        seq.append(r)
        if len(seq) > poly.degree() + 5:
            raise RuntimeError('too long')
    def var_at(sign_x):
        s = []
        for p in seq:
            d = p.degree()
            if d < 0:
                val = p[0]
            else:
                lc = p[d]
                val = lc if sign_x > 0 else lc * ((-1) ** d)
            if val > 0:
                s.append(1)
            elif val < 0:
                s.append(-1)
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
    return var_at(-1) - var_at(1)

if __name__ == '__main__':
    j = 3
    Os = compute_O_coeffs(j, 120)
    for n in [30, 60, 90, 120]:
        On = Os[n]
        p = fmpz_poly(On)
        g = p.gcd(p.derivative())
        sqf = p // g if g.degree() > 0 else p
        rc = sturm_count(sqf)
        print(f'j={j} n={n} deg={p.degree()} sqf_deg={sqf.degree()} real_roots={rc} all_real={rc==sqf.degree()}')
