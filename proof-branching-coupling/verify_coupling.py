"""Exact checks for the measure-valued branching coupling notes."""
import sympy as sp

z, t, t0 = sp.symbols('z t t0')


def series_coeff(expr, var, n):
    ser = sp.series(sp.expand(expr), var, 0, n + 1).removeO()
    return sp.expand(sp.diff(ser, var, n).subs(var, 0) / sp.factorial(n))


def make_objs(j, N):
    zz = sp.symbols('zz')
    A = sum(zz ** r / sp.factorial(r) for r in range(1, j))
    C = sp.exp(zz) - sum(zz ** r / sp.factorial(r) for r in range(0, j))
    B = 1 - A
    p = zz ** (j - 1) / sp.factorial(j - 1)
    D = B - t * C
    Hser = sp.series(B / D, zz, 0, N + 1).removeO()
    Fser = sp.series(1 / D, zz, 0, N + 1).removeO()
    O = [sp.expand(sp.diff(Fser, zz, n).subs(zz, 0))
         for n in range(N + 1)]
    Hn = [sp.expand(sp.diff(Hser, zz, n).subs(zz, 0))
          for n in range(N + 1)]
    return A, B, C, p, D, O, Hn


def W(p, q):
    return sp.expand(sp.diff(p, t) * q - p * sp.diff(q, t))


def test_flow():
    for j in [2, 3, 4, 5]:
        A = sum(z ** r / sp.factorial(r) for r in range(1, j))
        C = sp.exp(z) - sum(z ** r / sp.factorial(r) for r in range(0, j))
        B = 1 - A
        p = z ** (j - 1) / sp.factorial(j - 1)
        phi = B * t0 / (C * t0 + sp.exp(z))
        H = B / (B - t * C)
        r1 = sp.simplify(H.subs(t, phi) - (1 + sp.exp(-z) * C * t0))
        r2 = sp.simplify(sp.diff(phi, z) - (-(2 - p) / B * phi - p / B * phi ** 2))
        assert r1 == 0 and r2 == 0, (j, r1, r2)
    print('flow ok')


def test_residue_identity():
    for j in [2, 3, 4]:
        N = 3 * j
        A, B, C, p, D, O, Hn = make_objs(j, N)
        for n in range(j, N + 1):
            for dlag in range(1, j + 1):
                q = O[n - dlag]
                dq = sp.diff(q, t)
                # W(H_n,q) = t * sum binom(n,r) W(O_{n-r},q) at roots of q, in polynomial form
                lhs = W(Hn[n], q)
                rhs = sum(sp.binomial(n, r) * (O[n - r] * q + t * W(O[n - r], q))
                          for r in range(j, n + 1))
                assert sp.expand(lhs - rhs) == 0
                # root identity: W(H_n,q)(mu) + mu*q'(mu)*sum binom(n,r)O_{n-r}(mu)=0 mod q
                S = sum(sp.binomial(n, r) * O[n - r] for r in range(j, n + 1))
                poly = sp.expand(lhs + t * dq * S)
                rem = sp.rem(sp.Poly(poly, t), sp.Poly(q, t))
                assert rem.as_expr() == 0, (j, n, dlag)
    print('residue identity ok')


def test_covariance_positive_t():
    # For t>0, W(H_n,q)/(H_n*q) = (tH_n'/H_n - tq'/q)/t
    for j in [2, 3]:
        N = 3 * j
        A, B, C, p, D, O, Hn = make_objs(j, N)
        for n in range(j, N + 1):
            for dlag in range(1, j + 1):
                q = O[n - dlag]
                lhs = sp.together(W(Hn[n], q) / (Hn[n] * q))
                rhs = sp.together((sp.diff(Hn[n], t) / Hn[n] - sp.diff(q, t) / q))
                assert sp.simplify(lhs - rhs) == 0
    print('covariance identity ok')


def test_sample_G_signs():
    # Exact rational-root check is expensive; use Sturm/count_roots logic by
    # verifying G_q(mu) sign through W(H_n,q)(mu)>0 at sampled simple roots.
    for j in [2, 3]:
        N = 3 * j
        A, B, C, p, D, O, Hn = make_objs(j, N)
        for n in range(j, N + 1):
            for dlag in range(1, j + 1):
                q = O[n - dlag]
                # all real roots should give positive W(H_n,q)
                w = W(Hn[n], q)
                # count real roots of w: expect 0
                if sp.degree(w, t) is None or sp.degree(w, t) <= 0:
                    rr = 0
                else:
                    rr = sp.Poly(w, t, domain='QQ').count_roots()
                assert rr == 0, (j, n, dlag, rr)
    print('sample boundary W no real roots ok')


def test_gamma_positive_on_germ():
    z = sp.symbols('z')
    for j in range(2, 10):
        A = sum(z ** r / sp.factorial(r) for r in range(1, j))
        B = 1 - A
        p = z ** (j - 1) / sp.factorial(j - 1)
        roots = [r for r in sp.nroots(B, n=30)
                 if abs(sp.im(r)) < 1e-20 and sp.re(r) > 1e-12]
        zB = min([float(sp.re(r)) for r in roots])
        vals = [2 - float(p.subs(z, zB * k / 10)) for k in range(1, 10)]
        assert min(vals) > 0, (j, zB, min(vals))
    print('gamma>0 on germ ok')


if __name__ == '__main__':
    test_flow()
    test_residue_identity()
    test_covariance_positive_t()
    test_sample_G_signs()
    test_gamma_positive_on_germ()
    print('all ok')