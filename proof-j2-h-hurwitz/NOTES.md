# H-Wronskian Hurwitz/Bezoutian attack: quadratic certificate and exact obstruction

Date: 2026-11-28
Scope: only ordered-big-blocks/proof-j2-h-hurwitz was created/edited. No
installs, builds, external agents or models were used. SymPy exact arithmetic
is used for identities and certificates; finite enumeration and floating root
computations are labelled evidence, never proof.

## 0. Honest outcome

**No all-index proof of the target was obtained. No counterexample to the
actual associated-Stirling target was found.** What is new here:

1. An exact **quadratic-form / discriminant certificate** for the target
   (sections 1-2): under the recurrence split the Wronskian is a quadratic in
   c=(n-1)a, and h0>0, h2>0, h1^2<4 h0 h2 is sufficient for H_{n,a}>0 on all of R.
2. Exact finite verification of that certificate for the actual rows
   6<=n<=12, a in {1,2} (section 3), plus much broader but still finite
   coefficient/Hurwitz scans (section 5).
3. An **exact rational obstruction** (section 4): polynomials S,T,K satisfying
   the coupled-type interlacings (S surrounds T, K leftmost in S), with
   h0>0 and h2>0 on R, for which the discriminant condition fails and,
   already at c=1, the corresponding Wronskian has **two negative real
   roots**. Hence the previously proved interlacings plus positivity of the
   diagonal terms do **not** imply the target. A new quantitative domination
   invariant is indispensable.

The obstruction is not a counterexample to the actual P_{n,a}=Q_{n,a}
family; it is an exact obstruction to this certificate route.

## 1. Exact quadratic decomposition

Write D=t d/dt, E=1+D, U=1+E=2+tD, T(f)=(t f')'=E(f') and
W(f,g)=f'g-fg'. Use the exact recurrence (M) and identity (V) of the prior
notes:

    P_{n,a}=A+cB,  A=E P_{n-1,a},  B=t P_{n-2,a+1},  c=(n-1)a.        (1.1)

    T_{n,a}=alpha+c beta,  alpha=U T_{n-1,a}=(1+E)T_{n-1,a},
    beta=E^2 P_{n-2,a+1}.                                            (1.2)

Bilinearity of W gives the exact identity

    H_{n,a}=W(P_{n,a},T_{n,a})=h0+c h1+c^2 h2.                       (1.3)

    h0=W(A,alpha)=W(E P_{n-1,a}, U T_{n-1,a}),                       (1.4)
    h2=W(B,beta)=W(t P_{n-2,a+1}, E^2 P_{n-2,a+1}),                  (1.5)
    h1=W(A,beta)+W(B,alpha)
      =W(E P_{n-1,a}, E^2 P_{n-2,a+1})
       +W(t P_{n-2,a+1}, U T_{n-1,a}).                               (1.6)

verify.py V1-V3 print symbolic residual 0 for (1.3), (1.4), (1.5) with
generic functions S=P_{n-1,a}, K=P_{n-2,a+1}; V2/V3 identify h0,h2 as the
displayed operator Wronskians. These are algebraic identities, valid for
*any* S,K, and hold for every n>=6, a>=1.

**Sufficient certificate.** If, pointwise on R,

    h0>0,   h2>0,   h1^2 < 4 h0 h2,                                  (1.7)

then the quadratic (1.3) is strictly positive for every real c, hence in
particular for c=(n-1)a>0:

    h0+c h1+c^2 h2
      =(c sqrt(h2)+h1/(2 sqrt(h2)))^2 + h0 - h1^2/(4 h2) > 0.        (1.8)

Evaluation of the square is legitimate because h2>0 on R. The condition
h1^2<4 h0 h2 is a Routh/Bezout-type discriminant certificate for the
polynomial family in c; equivalently it says the 2x2 Bezoutian-type matrix
[[h0,h1/2],[h1/2,h2]] is positive definite pointwise.

## 2. What is known about the three pieces

* h0 is an operator transform of the induction hypothesis. If T_{n-1,a}
  surrounds P_{n-1,a} (the target at index n-1), numerical evidence on
  random interlacing pairs strongly supports W(ES,U T)>0 on R whenever S
  surrounds T (test_ops2.py, stress_op.py: 420 random exact/Sturm tests, 0
  failures). **This operator lemma is not proved here** and is not used as
  a theorem.
* h2=W(tK,E^2K) is *not* positive for a general negative-rooted K:
  test_ops.py found 101 exact failures among 125 representative
  K=prod(1+r_k t) with r_k in {1,2,3,5,8}. So positivity of h2 is a
  property of the actual row K=P_{n-2,a+1}, not of real-rootedness alone.
* Similarly h1^2<4 h0 h2 is not a formal consequence of interlacing:
  test_cauchy_generic.py found 148 exact failures among 148 sampled
  triplets satisfying the basic interlacings. Section 4 gives one explicit
  witness.

Therefore (1.7) is genuinely a **new invariant**, not a restatement of the
coupled theorem or of multiplier/coefficient positivity.

## 3. Exact finite evidence for the actual rows (not a proof)

verify.py V4 uses exact Sturm counts. For 6<=n<=12, a in {1,2}:

* h0 has no real root, h0(0)>0, positive leading coefficient;
* h2 has no real root, h2(0)>0, positive leading coefficient;
* D:=4 h0 h2-h1^2 has no real root, D(0)>0, positive leading coefficient.

Hence (1.7) holds in this range and H_{n,a}>0 on R there. This is **finite
evidence only**. The same certificate was checked to fail for generic data
(section 2), so this is a property to be proved, not a formal identity.

## 4. Exact obstruction to the certificate route

Take the rational data

    S=(1+t/8)(1+t),   T=1+t/7,   K=(1+t/9)(1+t/3),   c=1.         (4.1)

Exact root order: K roots -9<-3, S roots -8<-1, T root -7; hence

    -9 < -8 < -7 < -3 < -1 < 0,

i.e. T lies strictly between the roots of S (S surrounds T), and the roots
of K and S strictly alternate with K leftmost / S rightmost. Sturm counts
confirm W(S,T) and W(K,S) have **no real roots**, with W(S,T)(0)=55/56>0
and W(K,S)(0)=-49/72<0: exactly the coupled-type interlacing orientation
appearing in the recurrence.

For these data (1.4)-(1.6) give h0(0)=57/14>0 and h2(0)=1>0, both with no
real roots and positive leading coefficients; verify.py V5 prints the exact
forms. Nevertheless

    D=4 h0 h2-h1^2   has exactly 2 real roots,                      (4.2)

and for c=1,

    H=h0+h1+h2
     =(168 t^4+2224 t^3+14397 t^2+57834 t+102654)/13608,           (4.3)

which has roots -4.8164748397..., -4.6036334696... and one complex pair.
In particular H<0 on a subinterval of the negative axis (the quartic has
positive leading coefficient and two negative real roots, so it is negative
between them). verify.py V5 prints every count exactly.

**Consequence.** No argument using only

1. the coupled interlacings T_{n-1} prec P_{n-1} and P_{n-2,a+1}
   interlacing P_{n-1} in the realised orientation, and
2. positivity of the two diagonal Wronskians h0,h2

can prove H_{n,a}>0. The missing quantitative inequality (1.7) must be
proved from the actual associated-Stirling recurrence; it does not follow
from the previously established interlacing theorems. This obstruction is
of the same logical type as the J-reduction obstruction of
proof-j2-rootchain2 section 5 but targets a different (Bezoutian/
discriminant) inference and is not a repetition of it.

## 5. Further finite evidence and failed shortcuts (all labelled evidence)

* Exact coefficient positivity: 4<=n<=140, 1<=a<=40, every coefficient of
  H_{n,a} strictly positive (coefscan.py; 0 failures). Positive coefficients
  alone do not settle t<0.
* Exact Routh first-column positivity for 6<=n<=18; floating Hurwitz scans
  for 4<=n<=130, 1<=a<=30 found no root with nonnegative real part
  (routh.py, bigsearch.py; the 2519 "suspects" were numerical artifacts of
  underflowing imaginary parts, all with max Re(z)<0).
* No simple recurrence for H_n alone was found: H_n-t^2 H_{n-2} has
  positive coefficients and no real roots for 8<=n<=15, but its
  factorisations are irreducible and no closure was found
  (hrecur.py, factorD.py).
* Hermite-Biehler even/odd parts of H_n (evenodd.py) and Routh minors
  (routh.py) showed no factorisation into known positive families.
* An SOS ansatz H=A(t)^2+t B(t)^2 with deg A=deg H/2, deg B=deg A-1 has
  **no solution** for n=8,10 (sos.py).
* Coefficient ratios of H_n are log-concave but fail Newton's inequality
  (coeffcheck.py), and log-concavity with positive coefficients does
  **not** imply Hurwitz stability: test_lc.py found thousands of exact
  counterexamples. So this common shortcut is unavailable.

None of these items is used as an all-index proof.

## 6. Proof-policy audit

### 6.1 Direction audit

* (1.8) uses h2>0 and h1^2<4 h0 h2 with the **strict** direction H>0; at
  c=(n-1)a>0 the square is >=0 and the residual h0-h1^2/(4h2) is >0.
  Sample: c=1 in section 4 satisfies h0,h2>0 but violates h1^2<4h0h2, and
  the conclusion fails -- the direction is not reversible.
* Target direction: H_{n,a}(t)>0 for t<0; (1.7) proves the stronger statement
  on all of R. Sample n=6,a=1: H=20(1620 t^2+450 t+107), H(0)=2140>0
  (verify.py V6).

### 6.2 Monotonicity audit

No monotonicity claim is used. The only "increasing" statement is that
u=A/B has derivative -W(A,B)/B^2; this was explored but is **not** used in
any conclusion (we observed W(A,B)<0 for actual rows, but no proof and no
deduction from it). No endpoint optimisation is used.

### 6.3 Feasible-constraint falsification

Section 4 is exactly this audit: the constraints actually proved elsewhere
(interlacing orientations) and the additionally verified h0,h2>0 permit a
witness with two negative real roots of H. The witness is recorded exactly;
the argument therefore cannot be completed with those constraints alone.

### 6.4 Boundary and degeneracy

* n=4,5: deg P=1, H_{4,1}=36, H_{5,1}=400, chain vacuous (verify.py V6).
* n>=6: active range; c=(n-1)a>0 for a>=1.
* a=1 reduction to higher a via common multipliers is prior context and is
  not re-proved here.
* h0,h2 are never divided by in the proof; (1.8) divides by h2 only under
  the hypothesis h2>0. Zeros of P or K are not division points in the
  displayed identities; Wronskians are evaluated by polynomial identities,
  not by ratios at roots.

### 6.5 Symbolic residuals

verify.py prints residual 0 for (1.3), (1.4), (1.5) with generic S,K,c
(V1-V3). Exact Sturm counts are printed for V4/V5. Finite scans are labelled
finite evidence (sections 3, 5). No numerical root value is used as a proof.

### 6.6 Retraction rule

No theorem was claimed. An earlier candidate operator lemma
"S surrounds T implies W(ES,U T)>0" was stress-tested but **not** promoted
to a theorem; it remains a conjecture (section 2). A first obstruction
model with S=(1+20t)(1+12t), T=1+16t, K=1+21t had the wrong K/S
orientation and was replaced by the exact correct-orientation witness
(4.1); the earlier model is not used. No invalid all-index claim is left in
the record.

## 7. Precise remaining gap

Prove (1.7) for S=P_{n-1,a} and K=P_{n-2,a+1} for all n>=6, a>=1.
The two diagonal pieces are finite-evidence positive; the discriminant piece
is the new invariant. Equivalently, prove that the 2x2 Bezoutian-type matrix
of (A,B) against (alpha,beta) is positive definite at every real t. Section
4 shows this cannot be inferred from the coupled interlacings alone. A
recurrence-specific domination inequality (quantitative, not merely
interlacing) is still missing.

## 8. Files

* NOTES.md -- this report.
* common.py -- exact integer/SymPy construction of Q_{n,a} and operators.
* verify.py -- all exact identities, finite certificate check, obstruction.
* explore1.py through evenodd.py, test_*.py, find_obstruction*.py --
  evidence/failed-route scripts; none is a proof.

No packages were installed and no build was run.
