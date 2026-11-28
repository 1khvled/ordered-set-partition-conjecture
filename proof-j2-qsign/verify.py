"""Exact verification for proof-j2-qsign.

Checks the corrected sign identity at roots of U Q_{m-2,a+1} and
audits the gap/obstruction. Uses sympy (already installed).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'proof-j2-distance-two'))
from explore import astirling_table, Q_coeffs, poly, W, t
import sympy as sp

U = lambda f: sp.expand(2*f + t*sp.diff(f, t))


def sep(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main():
    A = astirling_table(26)

    sep("V1: Corrected identity at roots of Uv (NOT the NOTES §6 formula)")
    # NOTES §6 claimed p(lambda) = (t*s)'(lambda) at roots of Uv. False.
    # Correct: p(lambda) = (t*s)'(lambda) + (m-1)*a*lambda*v(lambda)
    # where p=Q_{m,a}, s=Q_{m-1,a}, v=Q_{m-2,a+1}.
    bad_old = 0
    bad_new = 0
    checked = 0
    for m in range(5, 16):
        for a in range(1, 6):
            p_ = poly(Q_coeffs(A, m, a))
            s = poly(Q_coeffs(A, m-1, a))
            v = poly(Q_coeffs(A, m-2, a+1))
            Uv = sp.expand(U(v))
            if Uv == 0 or sp.Poly(Uv, t).degree() < 1:
                continue
            ident = sp.expand(p_ - (sp.diff(t*s, t) + (m-1)*a*t*v))
            if ident != 0:
                bad_new += 1
            for lam in sp.Poly(Uv, t).nroots(n=80):
                if abs(sp.im(lam)) > 1e-30:
                    continue
                lam = sp.re(lam)
                checked += 1
                p_l = complex(sp.N(p_.subs(t, lam), 80))
                ts_l = complex(sp.N(sp.diff(t*s, t).subs(t, lam), 80))
                if abs(p_l - ts_l) > 1e-30:
                    bad_old += 1
    print(f"  checked={checked}  old_formula_counterexamples={bad_old}  symbolic_identity_failures={bad_new}")

    sep("V2: Correct sign decomposition of -p(lambda)*Uv'(lambda)")
    # -p*Uv' = a(-Q1*Uv' + s*Uv') + (m-1)*a*lambda^2*W(v,v')
    # where Q1=Q_{m-1,a+1}
    bad = 0
    checked = 0
    for m in range(5, 16):
        for a in range(1, 6):
            p = poly(Q_coeffs(A, m, a))
            s = poly(Q_coeffs(A, m-1, a))
            v = poly(Q_coeffs(A, m-2, a+1))
            q1 = poly(Q_coeffs(A, m-1, a+1))
            Uv = sp.expand(U(v))
            if Uv == 0 or sp.Poly(Uv, t).degree() < 1:
                continue
            Uvp_poly = sp.diff(Uv, t)
            for lam in sp.Poly(Uv, t).nroots(n=60):
                if abs(sp.im(lam)) > 1e-30:
                    continue
                lam = sp.re(lam)
                checked += 1
                Uvp = complex(Uvp_poly.subs(t, lam))
                p_l = complex(p.subs(t, lam))
                q1_l = complex(q1.subs(t, lam))
                s_l = complex(s.subs(t, lam))
                v_l = complex(v.subs(t, lam))
                vp_l = complex(sp.diff(v, t).subs(t, lam))
                vpp_l = complex(sp.diff(v, t, 2).subs(t, lam))
                Wvv = vp_l*vp_l - v_l*vpp_l
                lhs = -p_l*Uvp
                rhs = a*(-q1_l*Uvp + s_l*Uvp) + (m-1)*a*(2*v_l*v_l + (complex(lam)**2)*Wvv)
                if abs(lhs-rhs) > 1e-15 * (1+abs(lhs)):
                    bad += 1
    print(f"  checked={checked}  decomposition_failures={bad}")

    sep("V3: Coupled control on Q1 and Laguerre term")
    bad_q1 = 0
    bad_lag = 0
    checked = 0
    for m in range(5, 16):
        for a in range(1, 6):
            v = poly(Q_coeffs(A, m-2, a+1))
            q1 = poly(Q_coeffs(A, m-1, a+1))
            Uv = sp.expand(U(v))
            if Uv == 0 or sp.Poly(Uv, t).degree() < 1:
                continue
            Uvp_poly = sp.diff(Uv, t)
            for lam in sp.Poly(Uv, t).nroots(n=60):
                if abs(sp.im(lam)) > 1e-30:
                    continue
                lam = sp.re(lam)
                checked += 1
                Uvp = complex(Uvp_poly.subs(t, lam))
                q1_l = complex(q1.subs(t, lam))
                term1 = (-q1_l*Uvp).real
                if term1 <= 1e-15:
                    bad_q1 += 1
                v_l = complex(v.subs(t, lam))
                vp_l = complex(sp.diff(v, t).subs(t, lam))
                vpp_l = complex(sp.diff(v, t, 2).subs(t, lam))
                Wvv = (vp_l*vp_l - v_l*vpp_l).real
                if Wvv < -1e-15:
                    bad_lag += 1
    print(f"  checked={checked}  coupled_Q1_bad={bad_q1}  laguerre_bad={bad_lag}")

    sep("V4: Q_{m-1,a} sign at roots of Uv is NOT controlled")
    good = 0
    bad_sign = 0
    checked = 0
    for m in range(5, 16):
        for a in range(1, 6):
            s = poly(Q_coeffs(A, m-1, a))
            v = poly(Q_coeffs(A, m-2, a+1))
            Uv = sp.expand(U(v))
            if Uv == 0 or sp.Poly(Uv, t).degree() < 1:
                continue
            Uvp_poly = sp.diff(Uv, t)
            for lam in sp.Poly(Uv, t).nroots(n=60):
                if abs(sp.im(lam)) > 1e-30:
                    continue
                lam = sp.re(lam)
                checked += 1
                Uvp = complex(Uvp_poly.subs(t, lam))
                s_l = complex(s.subs(t, lam))
                val = (-s_l*Uvp).real
                if val > 1e-15:
                    good += 1
                elif val < -1e-15:
                    bad_sign += 1
    print(f"  checked={checked}  -s*Uv'>0 (helps): {good}  -s*Uv'<0 (hurts): {bad_sign}")

    sep("V5: Countermodel - obstruction survives in §5 model")
    # Use rational model from proof-j2-distance-two §5
    from fractions import Fraction as F
    Am = {(4,1):F(1),(4,2):F(1,100),(5,1):F(1),(5,2):F(3,10),(6,1):F(1),(6,2):F(3,4),(6,3):F(1,20)}
    def Qc(m,a):
        cs=[]
        for rr in range(1,4):
            rising=1
            for k in range(rr): rising*=(a+k)
            cs.append(rising*Am.get((m,rr),F(0)))
        while cs and cs[-1]==0: cs.pop()
        return sum(sp.Rational(c.numerator,c.denominator)*t**i for i,c in enumerate(cs))
    p=Qc(6,1); s=Qc(5,1); v=Qc(4,2); q1=Qc(5,2)
    Uv=sp.expand(U(v))
    print("  p=",p)
    print("  W(p,Uv)=",sp.factor(W(p,Uv)))
    print("  W(q1,Uv)=",sp.factor(W(q1,Uv)), "(coupled control holds)")
    print("  W(s,Uv)=",sp.factor(W(s,Uv)), "(s sign term)")
    # evaluate decomposition
    for lam in sp.Poly(Uv,t).real_roots():
        Uvp=sp.diff(Uv,t).subs(t,lam)
        q1_l=q1.subs(t,lam); s_l=s.subs(t,lam)
        v_l=v.subs(t,lam); vp_l=sp.diff(v,t).subs(t,lam); vpp_l=sp.diff(v,t,2).subs(t,lam)
        Wvv=sp.expand(vp_l**2 - v_l*vpp_l)
        t1=sp.expand(1*(-q1_l*Uvp))
        t2=sp.expand(1*(s_l*Uvp))
        t3=sp.expand((6-1)*1*(2*v_l**2 + lam**2*Wvv))
        total_val=sp.expand(-p.subs(t,lam)*Uvp)
        print(f"  lam={sp.N(lam)}: T1={sp.N(t1)} T2={sp.N(t2)} T3={sp.N(t3)} total={sp.N(total_val)}")

    sep("V6: W(Q_{m-1,a}, Q_{m-2,a+1}) - a candidate invariant")
    bad_sv=0; total_sv=0
    for m in range(5,18):
        for a in range(1,7):
            s=poly(Q_coeffs(A,m-1,a)); v=poly(Q_coeffs(A,m-2,a+1))
            if s==0 or v==0: continue
            total_sv+=1
            wsv=sp.expand(W(s,v))
            poly_w = sp.Poly(wsv,t)
            if poly_w.count_roots(-sp.oo,sp.oo)!=0 or poly_w.LC()<=0:
                bad_sv+=1
    print(f"  checked={total_sv}  bad={bad_sv} (on actual triangle, finite evidence)")

if __name__ == "__main__":
    main()
