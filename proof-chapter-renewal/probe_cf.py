"""Exploration only: residue / Jacobi-fraction structure of H_n / O_{n-d}."""
import sympy as sp
from explore import rows, W

t = sp.symbols('t')


def pick_residues(p, q):
    """Assume q has simple real roots; return (linear coefficient, [(pole, weight)]) for p/q
    written as a*t+b+sum w/(mu-t), with w=-p(mu)/q'(mu)."""
    Q = sp.Poly(q, t)
    roots = sp.nroots(Q, n=25, maxsteps=200)
    P = sp.Poly(p, t)
    # polynomial part
    a, r = sp.div(P, Q)
    res = []
    for mu in roots:
        if abs(sp.im(mu)) > 1e-12:
            return None  # not real-rooted q
        m = sp.re(mu)
        weight = -sp.N(p.subs(t, m) / sp.diff(q, t).subs(t, m), 20)
        res.append((float(m), complex(weight)))
    return sp.expand(a.as_expr()), res


def jacobi_params(p, q):
    """Euclidean-type expansion p/q = alpha0*t+beta0 + 1/(alpha1*t+beta1 + 1/(...)).
    Returns list of (alpha,beta)."""
    P, Q = sp.Poly(p, t), sp.Poly(q, t)
    out = []
    while Q.degree() >= 1 and P.degree() == Q.degree() + 1:
        a, r = sp.div(P, Q)  # P = a*Q + r
        out.append((a.all_coeffs()[0], a.all_coeffs()[-1] if a.degree() >= 0 else 0))
        if r.is_zero:
            break
        P, Q = Q, r
    return out


if __name__ == '__main__':
    for j in [2, 3, 4]:
        O, H = rows(j, 4 * j)
        print('=== j =', j)
        for d in [1, j]:
            n = 3 * j
            p, q = H[n], O[n - d]
            out = pick_residues(p, q)
            if out is None:
                print(' d=%d  q not real-rooted' % d)
                continue
            a, res = out
            print(' d=%d n=%d  a=%s  deg diff=%d' % (d, n, a, sp.degree(p,t)-sp.degree(q,t)))
            print('   poles/weights:', [(round(mu,5), round(w.real,6), round(w.imag,6)) for mu,w in res])
            print('   all weights > 0?', all(w.real > 0 and abs(w.imag)<1e-8 for mu,w in res))
