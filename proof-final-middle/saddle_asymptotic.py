import mpmath as mp
mp.mp.dps = 60
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

def A_j(j, z):
    s = mp.mpf('0')
    for r in range(1, j):
        s += z**r / mp.factorial(r)
    return s

def C_j(j, z):
    return mp.e**z - sum(z**r / mp.factorial(r) for r in range(j))

def G_j(j, z):
    return 1 / (1 - A_j(j, z))

def H_j(j, z):
    return mp.factorial(j) * z**(-j) * C_j(j, z) * G_j(j, z)

def logH(j, z):
    return mp.log(H_j(j, z))

def psi(j, z):
    h = abs(z) * mp.mpf('1e-7')
    if h < mp.mpf('1e-12'):
        h = mp.mpf('1e-12')
    Hp = (H_j(j, z + h) - H_j(j, z - h)) / (2 * h)
    return z * Hp / H_j(j, z)

def psi_prime(j, z):
    h = abs(z) * mp.mpf('1e-6')
    if h < mp.mpf('1e-10'):
        h = mp.mpf('1e-10')
    return (psi(j, z + h) - psi(j, z - h)) / (2 * h)

def alpha_pole(j):
    f = lambda z: A_j(j, z) - 1
    return mp.findroot(f, (mp.mpf('0.1'), mp.mpf('0.99')))

def bisect_root(f, lo, hi, tol=None, maxiter=300):
    if tol is None:
        tol = mp.mpf('1e-40')
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        return None
    for _ in range(maxiter):
        mid = (lo + hi) / 2
        fm = f(mid)
        if abs(fm) < tol or abs(hi - lo) < tol:
            return mid
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2

def saddle(j, alpha):
    ap = alpha_pole(j)
    mu = j * alpha / (1 - alpha)
    f = lambda z: psi(j, z) - mu
    lo, hi = mp.mpf('1e-6'), ap - mp.mpf('1e-6')
    z0 = bisect_root(f, lo, hi)
    return z0

if __name__ == '__main__':
    j = 3
    ap = alpha_pole(j)
    print('alpha_pole', ap)
    for m in [20, 40]:
        n = j * m
        Os = compute_O_coeffs(j, n)
        On = Os[n]
        cnm = On[m]
        for ell in [1, m // 4, m // 2, 3 * m // 4, m - 1]:
            alpha = mp.mpf(ell) / m
            z = saddle(j, alpha)
            if z is None:
                print('no saddle', m, ell)
                continue
            k = m - ell
            d = j * ell  # s=0
            # saddle approx: [z^d] H^k G ~ H(z)^k G(z) / (z^d * sqrt(2*pi*k*nu(z)))
            # where nu(z)=z*psi'(z)
            nu = z * psi_prime(j, z)
            Hz = H_j(j, z)
            Gz = G_j(j, z)
            approx = (Hz**k * Gz) / (z**d * mp.sqrt(2 * mp.pi * k * nu))
            # actual coefficient
            actual_coeff = mp.mpf(On[m - ell]) * (mp.factorial(j) ** k) / mp.factorial(n)
            ratio = actual_coeff / approx
            print(f'm={m} ell={ell} alpha={float(alpha):.3f} z*={float(z):.6f} nu={float(nu):.6f} ratio={float(ratio):.6f}')
