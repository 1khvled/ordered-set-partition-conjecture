# The remaining root chain: exact reduction, multiplier reduction, and a new obstruction

Date: 2026-11-28  
Scope: only ordered-big-blocks/proof-j2-rootchain2 was created/edited. The
pre-existing proof-j2-rootchain directory was read but not modified. No
installs, builds, external agents, or models were used. Exact SymPy arithmetic
is used for identities and rational certificates; numerical root computations
are never used as proof.

## 0. Honest outcome

**No all-index proof of the requested chain was obtained, and no actual
associated-Stirling counterexample was found.**

The exact remaining assertion is reduced to one polynomial inequality:

\[
\boxed{\quad H_{n,a}(x):=W_x\bigl(R,(xR')'\bigr)>0\qquad(x>0),\quad}
\tag{0.1}
\]

where \(R=Q_{n,a}(-x)\). Equivalently, \(R\) strictly surrounds the roots of
\((xR')'\). This is exactly the missing link \(\rho_j<\zeta_j\); all other links
in the requested chain are proved below.

Two new rigorous reductions are obtained:

1. **Parameter reduction.** The \(a=1\) case implies every \(a\ge1\), because
   both \(Q_{n,a}\) and \((tQ_{n,a}')'\) are obtained from their \(a=1\)
   analogues by common multiplier sequences.
2. **Integration reduction.** If \(G_n=Q_{n,1}\), then
   \(G_n=(tJ_n)'\), where \(J_n=t^{-1}\int_0^tG_n(u)\,du\). The coupled theorem
   makes \(G_n\) real-rooted, and Rolle then makes \(J_n\) real-rooted.
   Moreover the coefficient sequence of \(J_n\) is itself a multiplier sequence.
   Nevertheless these facts do **not** imply (0.1): §5 gives a new exact rational
   obstruction satisfying all of them.

Thus the coupled theorem in proof-u-operator/NOTES.md §7, even after these two
natural reductions, does not close the gap. A new recurrence-specific invariant
is still required. This note does not claim D2C-1 for general \(a\).

## 1. Exact reduction of the chain

Put

\[
d=\deg R=\lfloor n/2\rfloor-1,\qquad C=(xR')'.
\]

Write the positive roots as

\[
0<\rho_1<\cdots<\rho_d,\quad
0<\sigma_1<\cdots<\sigma_d,\quad
0<\eta_1<\cdots<\eta_{d-1},\quad
0<\zeta_1<\cdots<\zeta_{d-1}.
\]

For \(n=4,5\), \(d=1\), so the requested chain has no active index. Henceforth
\(n\ge6\) and \(1\le j\le d-1\).

### 1.1 The two \(V\)-links

Identity (V) in the \(x\)-variable is

\[
aV=(a+1)R+xR'.
\tag{1.1}
\]

At \(R(x)\ne0\),

\[
\frac{xR'(x)}{R(x)}=\sum_i\frac{x}{x-\rho_i},\qquad
\left(\frac{xR'}R\right)'=-\sum_i\frac{x}{(x-\rho_i)^2}<0.
\tag{1.2}
\]

The degree count and endpoint signs in (1.1) therefore give

\[
\sigma_j<\rho_j<\sigma_{j+1}.
\tag{1.3}
\]

Rolle gives

\[
\rho_j<\eta_j<\rho_{j+1}.
\tag{1.4}
\]

The last link is not a consequence of (1.4) alone. At \(x=\eta_j\), (1.1)
gives \(V(\eta_j)=(a+1)R(\eta_j)/a\). The sign of \(R(\eta_j)\) is the sign of
\(V\) on \((\sigma_j,\sigma_{j+1})\), so

\[
\eta_j<\sigma_{j+1}<\rho_{j+1}.
\tag{1.5}
\]

### 1.2 The derivative links

The polynomial \(xR'\) has roots \(0,\eta_1,\ldots,\eta_{d-1}\). Rolle once
more gives

\[
0<\zeta_1<\eta_1<\zeta_2<\eta_2<\cdots<
\zeta_{d-1}<\eta_{d-1}.
\tag{1.6}
\]

Consequently,

\[
\sigma_j<\rho_j<\zeta_j<\eta_j<\sigma_{j+1}
\quad\Longleftrightarrow\quad
\rho_j<\zeta_j.
\tag{1.7}
\]

Since \(\zeta_j<\eta_j<\rho_{j+1}\), (1.7) is precisely strict interlacing of
\(R\) and \(C\), with \(R\) surrounding \(C\). The Wronskian criterion gives

\[
R\prec C\quad\Longleftrightarrow\quad H_{n,a}>0\text{ on }x>0.
\tag{1.8}
\]

At \(\rho_j\), this is also \(K(\rho_j)>0\), where
\(K=C/R'\). The relevant root-flow intervals are exactly the intervals on
which the desired chain forces \(K>0\); no global positivity of \(K\) between
consecutive roots of \(R\) is asserted.

**Correction to the stale reduction.** proof-j2-rootchain/NOTES.md §1 wrote
\(\eta_j<\rho_{j+1}<\sigma_{j+1}\). The correct order is (1.5):
\(\eta_j<\sigma_{j+1}<\rho_{j+1}\). The stale file was not modified.

## 2. Reduction from all \(a\) to \(a=1\)

Use the \(t\)-variable and write

\[
P_{n,a}(t)=Q_{n,a}(t)=\sum_{k\ge0}p_{n,a,k}t^k,
\qquad p_{n,a,k}=(a)_{k+1}A_{n,k+1}.
\]

Identity (V) says that the step \(a\to a+1\) multiplies the coefficient of
degree \(k\) by

\[
\frac{a+1+k}{a}.
\tag{2.1}
\]

An arithmetic sequence \(k+b\) is a multiplier sequence: its exponential
generating function is \(e^z(z+b)\), which is in the Laguerre–Pólya class.
Products of multiplier sequences are multiplier sequences. Hence there is a
multiplier sequence \(M_a\) such that

\[
P_{n,a}=M_aP_{n,1}.
\tag{2.2}
\]

Define

\[
T_{n,a}=(tP_{n,a}')'.
\]

On the coefficient of \(t^{k-1}\), \(T_{n,a}\) is obtained from
\(p_{n,a,k}\) by multiplication by \(k^2\). Therefore the same parameter
multiplier applies:

\[
T_{n,a}=M_aT_{n,1},
\tag{2.3}
\]

where on degree \(k-1\) the sequence is the shifted multiplier
\(m_{j+1}\). A common multiplier preserves strict interlacing. Since

\[
C_{n,a}(x)=-T_{n,a}(-x),
\]

it follows that

\[
\boxed{\text{the }a=1\text{ case of (0.1) implies every }a\ge1.}
\tag{2.4}
\]

This reduction uses (V), not a root-flow argument.

## 3. The exact recurrence for the missing Wronskian

Let \(D=t\,d/dt\) and \(E=1+D\). For any polynomial \(f\),

\[
T(f):=(tf')'=Ef'.
\]

Applying this to (M),

\[
P_{n,a}=(tP_{n-1,a})'+(n-1)a\,tP_{n-2,a+1},
\]

gives the exact identity

\[
\boxed{
T_{n,a}=(E+1)T_{n-1,a}
 +(n-1)a\,E^2P_{n-2,a+1}.
}
\tag{3.1}
\]

verify.py V3 prints the generic symbolic residual \(0\). Equation (3.1) is
the recurrence that any all-index proof must control; the two coupled
interlacings alone do not supply the required sign of
\(W(P_{n,a},T_{n,a})\).

## 4. The integration reduction and what the coupled theorem does give

For \(a=1\), put \(G_n=P_{n,1}\) and

\[
J_n(t)=\frac1t\int_0^tG_n(u)\,du
      =\sum_{k\ge0}k!A_{n,k+1}t^k.
\]

Then

\[
(tJ_n)'=G_n,\qquad G_n=EJ_n.
\tag{4.1}
\]

The coupled theorem makes \(G_n\) real-rooted with simple negative roots.
Since \(tJ_n\) has degree \(d+1\) and the additional root \(0\), Rolle applied
to \((tJ_n)'=G_n\) shows that \(J_n\) has \(d\) simple negative roots.

Furthermore, if \(B_n(z)=\sum_rA_{n,r}z^r\), then

\[
\sum_{k\ge0}\frac{j_{n,k}}{k!}z^k=\frac{B_n(z)}z.
\tag{4.2}
\]

Bóna–Mező real-rootedness of \(B_n\), with its root at \(0\), shows that the
coefficient sequence of \(J_n\) is a multiplier sequence. Thus the actual
\(J_n\) has both natural properties used in the next obstruction. Neither is
sufficient.

## 5. New exact obstruction to the integration reduction

This is not the stale local recurrence model. It is an exact polynomial witness
showing that even **real-rootedness of \(J\), multiplier coefficients of \(J\),
and real-rootedness of \(P=EJ\) and \(T=EP'\)** do not imply \(W(P,T)>0\) on
\(t\le0\).

Take

\[
J(t)=(1+t)(1+2t)(1+3t)=1+6t+11t^2+6t^3.
\]

Then

\[
P=EJ=1+12t+33t^2+24t^3,
\qquad
T=EP'=12+132t+216t^2.
\]

Both \(P\) and \(T\) are real-rooted because \(E\) and the coefficient
multiplier \(k\) are multiplier operators. The coefficient EGF of \(J\) is

\[
1+6z+\frac{11}{2}z^2+z^3.
\]

Twice this polynomial is \(2z^3+11z^2+12z+2\), whose discriminant is
\(2024>0\); its coefficients are positive, so all three roots are negative.
Thus \(J\)'s coefficient sequence is also a multiplier sequence.

Nevertheless

\[
W(P,T)=12(432t^4+528t^3+219t^2+30t+1),
\]

and exactly

\[
W(P,T)\big|_{t=-1/10}=-\frac{2211}{625}<0.
\tag{5.1}
\]

So the proposed reduction stops exactly at (0.1). The witness is **not** a
\(Q_{n,a}\) polynomial and is not claimed to satisfy (M); it is an obstruction
to the inference from the coupled theorem plus (4.1)–(4.2), not a counterexample
to the actual associated-Stirling theorem.

## 6. Proof-policy audit

### 6.1 Direction audit

Each link is audited separately:

- (1.2) is strictly decreasing, hence \(\sigma_j<\rho_j\), not the reverse.
- Rolle gives \(\rho_j<\eta_j\).
- \(V(\eta_j)=(a+1)R(\eta_j)/a\) has the sign of \(V\) on
  \((\sigma_j,\sigma_{j+1})\), hence
  \(\eta_j<\sigma_{j+1}<\rho_{j+1}\).
- Rolle gives \(\zeta_j<\eta_j\).
- Therefore the only missing direction is \(\rho_j<\zeta_j\).

An exact feasible direction sample is \(n=6,a=1\):

\[
\sigma_1\in(1/75,1/70),\quad
\rho_1\in(1/50,1/48),\quad
\zeta_1=5/36,\quad \eta_1=5/18,
\]
\[
\sigma_2\in(2/5,41/100),\quad
\rho_2\in(53/100,54/100).
\]

Thus
\(\sigma_1<\rho_1<\zeta_1<\eta_1<\sigma_2<\rho_2\).
verify.py V5 certifies every interval by exact rational endpoint signs.

### 6.2 Monotonicity audit

The only monotonicity used in the root reduction is

\[
\left(\frac{xR'}R\right)'=-\sum_i\frac{x}{(x-\rho_i)^2}<0
\quad(x>0).
\]

No parameter monotonicity, endpoint optimization, or root-flow velocity bound
is used. Multiplier sequences are used only for their root-preservation and
interlacing theorems.

### 6.3 Feasible-constraint falsification

Section 5 is an exact rational/algebraic witness satisfying all the natural
constraints extracted from the coupled theorem and the integration reduction,
but violating the target Wronskian sign. It proves that those constraints are
insufficient. It does not falsify the actual triangle.

### 6.4 Boundary and degeneracy cases

- \(n=4,5\): \(d=1\); no \(\eta_j,\zeta_j,\sigma_{j+1}\), so the chain is
  vacuous.
- \(n\ge6\): active indices are \(1\le j\le d-1\).
- \(a=1\): the base case to which §2 reduces.
- \(a\ge2\): covered by the common-multiplier reduction, with no zero divisor.
- \(R\) and \(V\) have no common root by strict vertical interlacing, so all
  displayed divisions by \(R\) are legitimate.
- The obstruction has simple real roots and positive coefficients throughout.

### 6.5 Symbolic residuals

verify.py prints:

- V1: zero residuals for (V), generic \(a\), \(4\le n\le12\);
- V2: zero residuals for (M), generic \(a\), \(4\le n\le12\);
- V3: zero residual for (3.1);
- V4: zero residuals for the parameter-multiplier product;
- V5: exact rational direction sample;
- V6: exact obstruction value \(-2211/625\) and discriminant \(2024\);
- V7: exact retraction certificate for the numerical candidate;
- V8: zero residual for \(2R+xR'-V\) at \(a=1,n=8\).

All finite enumerations are labelled finite evidence, never an all-index proof.

### 6.6 Retraction rule

Two exploratory claims are explicitly retracted:

1. A high-precision scan at \(n=129,a=1\) appeared to compare
   \(\zeta_1=4.81025\ldots\times10^{-21}\) with
   \(\rho_1=1.92410\ldots\times10^{-20}\). Exact rational counts at
   \(q=10^{-20}\) show **two** roots of \(R\) and only **one** root of \(C\) in
   \((0,q)\). The numerical roots had different indices; this is not a
   counterexample. verify.py V7 gives the exact Descartes/sign certificate.
2. An exploratory Möbius-coefficient calculation used the indices in the
   binomial sum in reverse. It is discarded and is not used in this note.

No invalid all-index theorem is left in the record. The stale directory's
order error is corrected in §1.1 without modifying that directory.

## 7. Finite exact evidence (not a proof)

A corrected exact Sturm search found no positive zero of
\(H_{n,a}=W(R,(xR')')\) for

\[
6\le n\le60,\qquad a\in\{1,2\}.
\]

The search was stopped after those two parameter values for resource reasons.
This is finite evidence only. The script is exact_search_correct.py.

## 8. Exact remaining gap

The missing all-index statement is (0.1), or equivalently

\[
Q_{n,1}(t)\prec (tQ_{n,1}(t)')'
\quad\text{for every }n\ge6,
\tag{8.1}
\]

after which §2 gives every \(a\ge1\). The coupled theorem proves real-rootedness
and the two usual interlacings, and (M) gives (3.1), but the exact witness in
§5 shows that the natural integration/multiplier consequences do not imply
(8.1). A stronger invariant specific to the actual associated-Stirling initial
data remains to be proved.

## 9. Files

- NOTES.md — this report.
- verify.py — exact identities, rational direction sample, obstruction, and
  retraction certificate. Run py -3.14 verify.py.
- exact_search_correct.py — finite exact Sturm search; not a proof.
- verify-out.txt — saved output of verify.py.

No packages were installed and no build was run.
