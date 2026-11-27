import mpmath as mp
mp.mp.dps = 50

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

def psi(j, z):
    # z*H'/H
    h = abs(z) * mp.mpf('1e-7')
    if h < mp.mpf('1e-12'):
        h = mp.mpf('1e-12')
    Hp = (H_j(j, z + h) - H_j(j, z - h)) / (2 * h)
    return z * Hp / H_j(j, z)


def bisect_root(f, lo, hi, tol=None, maxiter=200):
    if tol is None:
        tol = mp.mpf('1e-30')
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        # try to expand/shrink
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

def alpha_pole(j):
    f = lambda z: A_j(j, z) - 1
    lo, hi = mp.mpf('0.1'), mp.mpf('2')
    for _ in range(100):
        if f(hi) > 0:
            break
        hi *= 2
    return mp.findroot(f, (lo, hi))

if __name__ == '__main__':
    for j in [3, 4, 5]:
        ap = alpha_pole(j)
        print(f"j={j} alpha_pole={ap}")
        for a in [0.1, 0.25, 0.5, 0.75, 0.9]:
            mu = j * a / (1 - a)
            f = lambda z: psi(j, z) - mu
            lo = mp.mpf('1e-6')
            hi = ap - mp.mpf('1e-6')
            # ensure bracket
            flo, fhi = f(lo), f(hi)
            if flo * fhi > 0:
                # try shrinking hi
                hi2 = ap - mp.mpf('1e-4')
                if f(lo) * f(hi2) <= 0:
                    hi = hi2
                else:
                    print('  fail', a, 'no bracket', flo, fhi)
                    continue
            z0 = bisect_root(f, lo, hi)
            if z0 is None:
                print('  fail', a, 'no bracket2')
                continue
            print(f"  alpha={a:.2f} mu={mu:.6f} z*={z0}")
