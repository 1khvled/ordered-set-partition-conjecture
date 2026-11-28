# Exact obstruction to the remaining root-chain step

Date: 2026-11-28  
Scope: only ordered-big-blocks/proof-j2-rootchain was created/edited.  
No installs, builds, external agents, models, or browser use.

## 0. Honest outcome

**No all-index proof of the requested chain was obtained.** An exact rational
obstruction shows that the two coupled interlacings of
proof-u-operator/NOTES.md §7, even together with the exact recurrences (M)
and (V), do **not** imply the remaining link
\(\rho_j<\zeta_j\). Thus merely iterating the coupled theorem in its present
form cannot close this gap.

This is **not** a counterexample to the theorem for the actual
associated-Stirling triangle. The witness is a rational model for the local
recurrence step. It satisfies every local input named in the task but violates
the desired conclusion. A new family-specific invariant is therefore required.

The rigorous reduction is: three links of the requested chain are automatic,
and the entire problem is exactly \(R\prec (xR')'\), equivalently
\(K(\rho_j)>0\), where \(K=(xR')'/R'\).

## 1. Exact reduction of the chain

Let
\[
 R(x)=Q_{n,a}(-x),\qquad V(x)=Q_{n,a+1}(-x),\qquad
 C(x)=(xR')'(x),
\]
and let \(k=\deg R=\lfloor n/2\rfloor-1\). Write
\[
 0<\rho_1<\cdots<\rho_k,\qquad
 0<\sigma_1<\cdots<\sigma_k
\]
for the roots of \(R,V\). Let \(\eta_j\) be the roots of \(R'\), and
\(\zeta_j\) the roots of \(C\).

From (V),
\[
 \frac Va=\frac{a+1}{a}+\frac{xR'}{aR},
 \qquad
 \frac{xR'}R=\sum_i\frac{x}{x-\rho_i}.
\]
For \(x>0\),
\[
 \left(\frac{xR'}R\right)'
 =-\sum_i\frac{x}{(x-\rho_i)^2}<0.
\]
The degree count and endpoint signs therefore give
\[
 \sigma_j<\rho_j<\sigma_{j+1}\quad(1\le j<k).
 \tag{1.1}
\]

Rolle's theorem gives
\[
 \rho_j<\eta_j<\rho_{j+1}.
 \tag{1.2}
\]
The polynomial \(xR'\) has roots \(0,\eta_1,\ldots,\eta_{k-1}\). Differentiating
and applying Rolle once more gives
\[
 0<\zeta_1<\eta_1<\zeta_2<\eta_2<\cdots<
 \zeta_{k-1}<\eta_{k-1}.
 \tag{1.3}
\]
In particular, \(\zeta_j<\eta_j\). Combining (1.1)–(1.3),
\[
 \eta_j<\rho_{j+1}<\sigma_{j+1}.
\]
Hence
\[
 \sigma_j<\rho_j<\zeta_j<\eta_j<\sigma_{j+1}
 \quad\Longleftrightarrow\quad
 \rho_j<\zeta_j\quad(1\le j<k).
 \tag{1.4}
\]

At \(\rho_j\), equation (1.2) places \(\rho_j\) in
\((\eta_{j-1},\eta_j)\), with \(\eta_0:=0\). On this interval \(R'\) has a
fixed sign, and \(C/R'=K\) changes sign at \(\zeta_j\), from positive to
negative. Therefore
\[
 \rho_j<\zeta_j
 \quad\Longleftrightarrow\quad
 K(\rho_j)>0,
 \qquad K=\frac{(xR')'}{R'}.
 \tag{1.5}
\]
Thus the exact missing assertion is \(R\prec C\).

For \(n=4,5\), \(k=1\), so there is no \(\zeta_j,\eta_j\), or
\(\sigma_{j+1}\); the chain is vacuous. For \(n\ge6\), the active range is
\(1\le j\le k-1\).

## 2. Exact rational obstruction

Work first in the \(t\)-variable. Put
\[
 q(t)=1+\frac{26}{5}t,\qquad
 p(t)=1+\frac{103}{20}t,\qquad
 A(t)=1+\frac{203}{10}t+78t^2,
\]
and \(\mathcal U f=2f+tf'\). Define, for \(a\ge1\),
\[
 \begin{aligned}
 Q_{4,a}(t)&=a+\frac{13}{5}a(a+1)t,\\
 Q_{5,a}(t)&=a+\frac{103}{40}a(a+1)t,\\
 Q_{6,a}(t)&=a+\frac{203}{20}a(a+1)t
             +13a(a+1)(a+2)t^2.
 \end{aligned}
 \tag{2.1}
\]
Thus \(Q_{4,1}=q\), \(Q_{5,1}=p\), and \(Q_{6,1}=A\).

### 2.1 (V) holds identically in \(a\)

For each \(n=4,5,6\),
\[
 aQ_{n,a+1}=(a+1)Q_{n,a}+tQ_{n,a}'.
 \tag{2.2}
\]
This follows directly from (2.1), and the symbolic residual is zero.

### 2.2 (M) holds identically in \(a\)

At \(n=6\),
\[
 Q_{6,a}=(tQ_{5,a})'+5a\,tQ_{4,a+1}.
 \tag{2.3}
\]
Indeed,
\[
 (tQ_{5,a})'=a+\frac{103}{20}a(a+1)t,
\]
while
\[
 5a\,tQ_{4,a+1}
 =5a(a+1)t+13a(a+1)(a+2)t^2.
\]
Adding gives exactly (2.1). The symbolic residual is zero for generic \(a\).

### 2.3 Both section-7 coupled interlacings hold

In the positive \(x\)-variable set
\[
 R_0(x)=A(-x)=1-\frac{203}{10}x+78x^2.
\]
Its discriminant is
\[
 \Delta=\frac{10009}{100}>0.
\]
Its two roots \(\rho_1<\rho_2\) are positive because the leading coefficient,
constant term, and sum of roots are positive.

The roots of \(p(-x)\) and \((\mathcal Uq)(-x)\) are respectively
\[
 s_p=\frac{20}{103},\qquad s_u=\frac5{39},
 \qquad 0<s_u<s_p.
\]
Exact evaluation gives
\[
 R_0(s_u)=-\frac{25}{78}<0,\qquad
 R_0(s_p)=-\frac9{10609}<0.
\]
Since \(R_0\) is negative between its two roots,
\[
 \rho_1<s_u<s_p<\rho_2.
 \tag{2.4}
\]
Thus \(A\) strictly surrounds both \(p\) and \(\mathcal Uq\). In the
negative-root \(t\)-convention this is precisely the conclusion of the
section-7 coupled theorem at \(r=6\):
\[
 G_5\preceq G_6,\qquad \mathcal UG_4\preceq G_6,
 \quad G_4=q,\ G_5=p,\ G_6=A.
 \tag{2.5}
\]

### 2.4 Nevertheless \(\rho_1<\zeta_1\) fails

Here
\[
 C_0(x)=(xR_0')'(x)
 =-\frac{203}{10}+312x,
 \qquad
 \zeta_1=\frac{203}{3120}.
\]
Exact comparisons give
\[
 0<\zeta_1<s_u,
 \qquad
 R_0(\zeta_1)=\frac{391}{41600}>0.
\]
By (2.4), \(\zeta_1<\rho_2\). A quadratic with positive leading coefficient
is positive to the left of its smaller root and to the right of its larger
root. Since \(\zeta_1<\rho_2\) and \(R_0(\zeta_1)>0\), necessarily
\[
 \boxed{\ \zeta_1<\rho_1\ }.
 \tag{2.6}
\]
So the witness violates the exact missing link, despite (2.2), (2.3), and
(2.5).

## 3. Why this is an exact obstruction, not a theorem counterexample

The actual associated-Stirling triangle has fixed initial data. The rational
rows (2.1) are not claimed to be those actual rows. What is proved is the
logical implication failure
\[
 \text{(M),(V), and the two section-7 coupled interlacings}
 \ \not\Rightarrow\ R\prec(xR')'.
 \tag{3.1}
\]
The witness is feasible for every local hypothesis in (3.1), and all
arithmetic is exact. Therefore an all-index proof must use an additional
property of the actual triangle; it cannot be obtained by treating the
coupled theorem as a black box and iterating it unchanged.

The smallest-root result in proof-j2-z-bound is neither repeated nor
weakened. No numerical root check is used as proof. No multiplier-flow,
\(A+tB\), or WRS argument is used.

## 4. Proof-policy audit

### 4.1 Direction audit

* (1.1): \(xR'/R\) is strictly decreasing on \(x>0\), so each \(V\)-root
  lies in \((\rho_{j-1},\rho_j)\), with \(\rho_0:=0\). Direction is
  \(\sigma_j<\rho_j\), not the reverse.
* (1.2)–(1.3): Rolle gives \(\rho_j<\eta_j<\rho_{j+1}\) and
  \(\zeta_j<\eta_j\). The last chain link is
  \(\eta_j<\rho_{j+1}<\sigma_{j+1}\).
* (2.4): both exact values are negative, and the quadratic is negative
  precisely between its roots. Hence
  \(\rho_1<s_u<s_p<\rho_2\).
* (2.6): \(R_0(\zeta_1)>0\), \(\zeta_1<s_u<\rho_2\), and the leading
  coefficient is positive. Therefore \(\zeta_1\) is left of \(\rho_1\).
  The obstruction direction is \(\zeta_1<\rho_1\), opposite to the target.

### 4.2 Monotonicity audit

The only monotonicity used in the reduction is
\[
 \left(\frac{xR'}R\right)'
 =-\sum_i\frac{x}{(x-\rho_i)^2}<0.
\]
There is no parameter monotonicity or endpoint optimization. The witness
uses exact sign evaluations, not finite differences.

### 4.3 Feasible-constraint falsification

Section 2 is the required exact falsification. It satisfies (M), (V), and
both coupled interlacings, but fails \(R\prec C\). This proves that those
constraints alone are insufficient. It does not claim infeasibility for the
actual initial data.

### 4.4 Boundary and degeneracy cases

* \(n=4,5\): \(\deg R=1\); the requested chain has no active \(j\).
* \(n\ge6\): active indices are \(1\le j\le k-1\).
* The obstruction uses \(k=2\), so both roots and \(\zeta_1\) exist.
* All divisions in (2.1)–(2.3) are by \(a>0\); no zero-parameter case is
  included.
* No root is multiple in the obstruction: \(\Delta=10009/100>0\).

### 4.5 Symbolic residuals

verify.py prints zero residuals for (V) and (M) with generic \(a\), the
exact discriminant, the three exact evaluations above, and the exact strict
comparisons. It performs no finite enumeration and no numerical root solve.

### 4.6 Retraction rule

No all-index theorem is claimed. The failed implication is recorded in
(3.1), and the witness is explicitly distinguished from the actual
associated-Stirling triangle. No invalid proof remains in this directory.

## 5. Files

* NOTES.md — this report.
* verify.py — exact symbolic residuals and rational sign certificate.

Run from this directory: py -3.14 verify.py. No packages installed, no
builds, and all files are small.
