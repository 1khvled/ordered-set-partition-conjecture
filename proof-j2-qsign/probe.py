from fractions import Fraction
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'proof-j2-distance-two'))
from explore import astirling_table, Q_coeffs, poly, W, t
import sympy as sp

U = lambda f: sp.expand(2*f + t*sp.diff(f,t))


def check(Mmax=24):
    A = astirling_table(Mmax)
    bad_uv = []
    bad_v = []
    for m in range(5, Mmax+1):
        for a in range(1, 8):
            if m-2 < 0:
                continue
            p = poly(Q_coeffs(A, m, a))
            s = poly(Q_coeffs(A, m-1, a))
            v = poly(Q_coeffs(A, m-2, a+1))
            q1 = poly(Q_coeffs(A, m-1, a+1))  # Q_{m-1,a+1}
            if p == 0 or v == 0 or s == 0 or q1 == 0:
                continue
            Uv = sp.expand(U(v))
            WpUv = sp.expand(W(p, Uv))
            okUv = sp.Poly(WpUv, t).count_roots(-sp.oo, sp.oo) == 0 and sp.Poly(WpUv, t).LC() > 0
            Wpv = sp.expand(W(p, v))
            okv = sp.Poly(Wpv, t).count_roots(-sp.oo, sp.oo) == 0 and sp.Poly(Wpv, t).LC() > 0
            if not okUv:
                bad_uv.append((m, a))
            if not okv:
                bad_v.append((m, a))
            for lam in sp.Poly(Uv, t).real_roots():
                lam_n = sp.N(lam)
                p_l = sp.N(p.subs(t, lam))
                s_l = sp.N(s.subs(t, lam))
                q1_l = sp.N(q1.subs(t, lam))
                v_l = sp.N(v.subs(t, lam))
                Uvp = sp.N(sp.diff(Uv, t).subs(t, lam))
                wanted = -p_l * Uvp
                coupled = -q1_l * Uvp
                if wanted <= 0:
                    print('BAD Uv', m, a, 'lam', lam_n, 'wanted', wanted, 'coupled', coupled)
                    print('  p', p_l, 'q1', q1_l, 's', s_l, 'v', v_l, 'Uvprime', Uvp)
    print('bad_uv', len(bad_uv), 'bad_v', len(bad_v))
    print('bad_uv list', bad_uv[:20])
    print('bad_v list', bad_v[:20])


if __name__ == '__main__':
    check(int(sys.argv[1]) if len(sys.argv) > 1 else 24)
