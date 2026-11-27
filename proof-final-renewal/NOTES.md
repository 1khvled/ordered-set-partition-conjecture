# Nonlinear threshold renewal: exact non-preservation theorem and the weighted gap

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-renewal` was created/edited.

---

## 0. Outcome

I attacked the **only** exact threshold relation, the bilinear renewal

$$
F_j=\frac{F_{j-1}}{1-(1-t)p_jF_{j-1}},\qquad p_j=\frac{z^{j-1}}{(j-1)!}.
\tag{R}
$$

Two rigorous results are obtained:

1. **Exact non-preservation theorem (Theorems A and B below).** The renewal is
   **not** a real-rootedness preserving map on the natural cone of input
   families. Even if the entire threshold-$(j-1)$ family
   $(P_n)=(O_{n,j-1})$ consists of real-rooted polynomials with only negative
   roots, **nonnegative (PF) coefficients, the correct degrees**
   $\deg P_n=\lfloor n/(j-1)\rfloor$, and **strictly interlacing consecutive
   rows**, the renewal can produce a row with complex roots. Exact, finite
   certificates are given at $j=3$, $n=2$ (discriminant $-3$) and at $j=3$,
   $n=8$ (Sturm count $2<4$ together with a negative coefficient).
   Consequently **no preservation theorem of the form**

   > "real-rootedness + interlacing + PF + correct degree pattern of the
   > inputs $\Rightarrow$ real-rootedness of the renewal output"

   can exist, in any generality that ignores the *specific* origin of the
   rows $O_{n,j-1}$. This is a genuine obstruction to threshold induction,
   strictly stronger than the linear/Möbius no-go of
   `proof-final-threshold/NOTES.md` because it kills the nonlinear route too.

2. **Exact reduction of the remaining gap (Section 5).** The obstruction is
   localized to a *weighted bilinear inequality* at the level of the
   convolution term: writing $C_N=(Q\star P)_N$, one has the residue identity
   $Q_n(\mu)=(1-\mu)\binom{n}{j-1}C_N(\mu)$ at every root $\mu$ of $P_n$, so the
   sign pattern of $C_N$ at the roots of $P_n$ must satisfy a non-alternating,
   global condition. The sign pattern is computed exactly for the actual
   families and is **not** a simple alternation; $C_N$ itself is not
   real-rooted. Hence neither the interlacing-family (MSS) criterion nor a
   termwise sign argument can run.

The conjecture itself is **not** disproved: the actual rows $O_{n,j-1}$
satisfy additional relations (they are themselves renewal outputs from
threshold $j-2$, and obey the universal differential recurrence) that the
counterexample families do not.

---

## 1. The exact identities (all proved here)

Conventions: $A_j=\sum_{r=1}^{j-1}z^r/r!$, $C_j=e^z-\sum_{r=0}^{j-1}z^r/r!$,
$p_j=z^{j-1}/(j-1)!$, $F_j=1/(1-A_j-tC_j)=\sum_n O_{n,j}z^n/n!$.
Write $P_n=O_{n,j-1}$, $Q_n=O_{n,j}$, $T=t\,d/dt$, and
$\binom{n}{k}=0$ for $k<0$ or $k>n$.

### Lemma R0 (reciprocal identity).

$$
\frac1{F_j}=\frac1{F_{j-1}}-(1-t)p_j .
$$

*Proof.* $A_j=A_{j-1}+p_j$ and $C_{j-1}=C_j+p_j$, hence
$1-A_j-tC_j=(1-A_{j-1}-tC_{j-1})-p_j+tp_j$. $\blacksquare$

### Lemma R1 (renewal, exact).

$$
F_j=\frac{F_{j-1}}{1-(1-t)p_jF_{j-1}}
=F_{j-1}+(1-t)p_jF_jF_{j-1}.
\tag{R1}
$$

*Proof.* Solve R0 for $F_j$; the second form is the first multiplied out. $\blacksquare$

### Lemma R2 (coefficient form).

For $n\ge j-1$, with $N:=n-j+1$,

$$
Q_n=P_n+(1-t)\binom{n}{j-1}\sum_{m=0}^{N}\binom{N}{m}Q_mP_{N-m};
\qquad Q_n=P_n\ \ (n<j-1).
\tag{R2}
$$

*Proof.* Expand (R1) in EGFs. Since $p_jF_jF_{j-1}=\frac{z^{j-1}}{(j-1)!}
\sum_{a,b}Q_aP_b z^{a+b}/(a!b!)$, the coefficient of $z^n/n!$ is
$\sum_{a+b=N}\frac{n!}{(j-1)!\,a!\,b!}Q_aP_b
=\binom{n}{j-1}\sum_{a=0}^{N}\binom{N}{a}Q_aP_{N-a}$. $\blacksquare$

Define the **binomial convolution** $(A\star B)_N=\sum_{m=0}^{N}\binom{N}{m}A_mB_{N-m}$;
it is the coefficient of $z^N/N!$ in the EGF product. Then (R2) reads
$Q_n=P_n+(1-t)\binom{n}{j-1}(Q\star P)_{N}$. Verified exactly for $j=2,\dots,6$,
$n\le14$ in `verify_all.py` (0 failures).

### Lemma R3 (triangularity).

For $j\ge2$ and $n\ge j-1$ we have $N=n-j+1\le n-1$. Hence (R2) determines
$Q_n$ recursively from $P_0,\dots,P_n$ and $Q_0,\dots,Q_{n-1}$. The map
$\mathcal R_j:(P_n)_{n\ge0}\mapsto(Q_n)_{n\ge0}$ is therefore a well-defined
triangular (nonlinear) operator, and its value at $n$ does not depend on
$Q_m$ for $m>n-1$. $\blacksquare$

### Lemma R4 ($u$-form; positivity of the multiplier).

Put $u=1-t$ and $S_n(u):=O_{n,j-1}(1-u)$, $R_n(u):=O_{n,j}(1-u)$. Then

$$
R_n=S_n+u\binom{n}{j-1}(R\star S)_{N}.
\tag{R4}
$$

Moreover $2-e^z+uC_j(z)=1-A_j(z)-tC_j(z)$, so
$R_n(u)=n![z^n]\bigl(2-e^z+uC_j(z)\bigr)^{-1}$; the roots of $O_{n,j}$ in
$(-\infty,0]$ correspond exactly to roots of $R_n$ in $[1,\infty)$. Thus the
factor that was $(1-t)$ becomes a **positive** multiplier $u\ge1$ at all
putative roots, so the difficulty is not a sign of the multiplier but the
bilinear convolution itself. $\blacksquare$

### Lemma R5 (residue identity).

Let $\mu$ be any root of $P_n$ (real or complex). Then

$$
Q_n(\mu)=(1-\mu)\binom{n}{j-1}C_N(\mu),\qquad C_N:=(Q\star P)_N .
\tag{R5}
$$

In particular for a negative real root $\mu<0$, $1-\mu>0$, so
$\operatorname{sgn}Q_n(\mu)=\operatorname{sgn}C_N(\mu)$. $\blacksquare$

### Lemma R6 (degrees).

$\deg_tP_n=\lfloor n/(j-1)\rfloor$ and $\deg_tQ_n=\lfloor n/j\rfloor$.
For $j\ge3$ there are $n$ with $\deg Q_n<\deg P_n$, and then the top
coefficients of $P_n$ in degrees $\lfloor n/j\rfloor+1,\dots,
\lfloor n/(j-1)\rfloor$ are exactly cancelled by $-(1-t)\binom{n}{j-1}C_N$.
*Proof.* Order of $C_j$ at $0$ is $j$ and order of $A_j$ is $1$, so the
$t$-degree of $n![z^n](1-A_j-tC_j)^{-1}$ is $\lfloor n/j\rfloor$; the
cancellation is the equality $Q_n=P_n+(1-t)\binom{n}{j-1}C_N$. $\blacksquare$

---

## 2. Theorem A — minimal exact obstruction ($j=3$, $n=2$)

**Theorem A.** Let $j=3$ and take
$P_0=1$, $P_1=1+t$, $P_2=1+2t+\tfrac12t^2$ (these are $L_0(-t),L_1(-t),L_2(-t)$;
they have positive coefficients, all roots real and negative, and $P_1$
strictly interlaces $P_2$: the roots are $-1$ and $-2\pm\sqrt2$, ordered
$-3.414<-1<-0.586$). Let $(Q_n)$ be generated by (R2). Then

$$
Q_2=\tfrac12 t^2+t+2,\qquad \operatorname{disc}_t(Q_2)=1-4\cdot\tfrac12\cdot2=-3<0,
$$

so $Q_2$ has nonreal roots.

*Proof.* $N=2-3+1=0$, so (R2) gives $Q_2=P_2+(1-t)\binom22 Q_0P_0
=1+2t+\tfrac12t^2+1-t=2+t+\tfrac12t^2$. The discriminant is $-3$. The
interlacing of $P_1,P_2$ is certified exactly by the vanishing of the number
of real roots of the Wronskian $P_1'P_2-P_1P_2'$ (see `verify_all.py`). $\blacksquare$

This is the smallest possible instance: at $n=2$ the convolution has only the
term $Q_0P_0=1$, so the obstruction is purely the mismatch between the
quadratic input $P_2$ and the linear output forced by $\deg Q_2=\lfloor2/3\rfloor=1$.

---

## 3. Theorem B — degree-matched, strictly interlacing, PF obstruction

Theorem A could be dismissed because $\deg P_1=1\neq\lfloor1/2\rfloor=0$ and
$\deg P_2=2\neq\lfloor2/2\rfloor=1$. The next theorem removes that loophole.

**Theorem B.** Fix $j=3$. For $n=0,\dots,8$ set

$$
P_n(t)=L_{\lfloor n/2\rfloor}^{(\alpha_n)}(-t),\qquad
\alpha_n=\begin{cases}0,&n\text{ even},\\ \tfrac12,&n\text{ odd},\end{cases}
$$

where $L_d^{(\alpha)}$ is the generalized Laguerre polynomial. Then:

1. $\deg P_n=\lfloor n/2\rfloor=\lfloor n/(j-1)\rfloor$;
2. every coefficient of $P_n$ is positive;
3. $P_n$ has only real negative roots (exactly $\deg P_n$ real roots by Sturm);
4. $P_{n-1}$ **strictly** interlaces $P_n$ for $n=1,\dots,8$ (exact certificate:
   the Wronskian $P_{n-1}'P_n-P_{n-1}P_n'$ has **no** real zero);

and yet the renewal output at $n=8$ is

$$
Q_8(t)=\frac{4817t^4+91632t^3+1562808t^2-3691536t+2032488}{24},
$$

whose coefficient of $t$ is $-153814<0$; by exact Sturm counting it has only
**2** real roots out of degree 4. Hence $Q_8$ is not real-rooted and has a
negative coefficient, so it is not in the PF cone.

*Proof sketch.* Properties 1–3 are standard for $L_d^{(\alpha)}(-t)$ with
$\alpha>-1$. Property 4 is checked exactly: for real-rooted polynomials with
simple roots and degrees differing by at most one, interlacing is equivalent
to the Wronskian being sign-definite, i.e. having no real zero; Sturm counts
give 0 for all eight consecutive pairs. The values of $Q_m$ for $m\le8$ are
generated by the triangular recursion (R2) with exact rational arithmetic;
$Q_8$ is displayed above. The negative coefficient alone already shows
$Q_8$ cannot have all roots in $(-\infty,0]$ (a real-rooted polynomial with
nonpositive roots and positive leading coefficient has all coefficients
nonnegative). Exact Sturm counting gives 2 real roots, confirming that the
failure is genuine complex roots rather than a sign pattern. All checks are in
`verify_all.py`; no floating point enters any certificate. $\blacksquare$

---

## 4. Corollary C — what is definitively ruled out

**Corollary C.** Let $\mathcal C$ be any class of polynomial sequences
$(P_n)$ specified by any **subset** of the following four properties:

* every $P_n$ is real-rooted with all roots $\le0$;
* consecutive rows interlace;
* every coefficient of $P_n$ is nonnegative (PF cone);
* $\deg P_n=\lfloor n/(j-1)\rfloor$.

Then for $j=3$ the renewal $\mathcal R_3$ does **not** map $\mathcal C$ into
the real-rooted families: Theorem B exhibits a family whose rows $P_0,\dots,P_8$
satisfy all four properties and whose renewal output $Q_8$ has nonreal roots.
Consequently:

1. There is **no** general "renewal-PF" theorem of the shape
   $P_n$ PF/interlacing $\Rightarrow Q_n$ PF/interlacing.
2. In particular, the interlacing-family (MSS) criterion, the common-interlacer
   cone closure (`proof-final-common-interlacer/NOTES.md`), PF-multiplier and
   CZDS arguments, and binomial-convolution preservation cannot prove the
   induction step from generic properties of the previous row.
3. Combined with Theorems A–B of `proof-final-threshold/NOTES.md`, **every**
   transformation class previously proposed — derivatives, integrals, shifts,
   Möbius maps in $t$, fixed-kernel convolutions — and now the nonlinear
   renewal acting on a generic input cone, is excluded as a route that uses
   only those generic properties.

Any successful threshold induction must therefore exploit a property of the
*actual* rows $O_{n,j-1}$ beyond real-rootedness, interlacing, PF-ness and the
degree pattern — for instance that they are themselves renewal outputs from
threshold $j-2$, or that they satisfy the universal differential recurrence
of `README.md`. These are exactly the routes the narrower tracks
(`proof-final-wrs`, `proof-final-inner`, `proof-final-transform`) failed to
close, so the obstruction is consistent with the conjecture being true.

---

## 5. The precise remaining weighted bilinear inequality

Here is the exact local form of what must be proved, and why it is subtle.

Fix $n$, $N=n-j+1$, $C_N:=(Q\star P)_N$, $c_n=\binom{n}{j-1}$. By (R5), at each
real root $\mu$ of $P_n$,

$$
Q_n(\mu)=(1-\mu)c_n\,C_N(\mu),\qquad 1-\mu>0 .
\tag{WBI-1}
$$

Since $O_{n,j}$ has nonnegative coefficients, $Q_n(t)>0$ for $t\ge0$; and
$\deg Q_n<\deg P_n$ when $j\ge3$. A necessary condition for $Q_n$ to have
$\deg Q_n$ negative real roots is therefore that the sign vector
$(\operatorname{sgn}C_N(\mu_i))_i$, indexed by the ordered roots
$\mu_1<\dots<\mu_{\deg P_n}<0$ of $P_n$, has exactly the right number of sign
changes and endpoint signs; sufficiency additionally requires $Q_n$ to have
no nonreal zeros. Equivalently (Aissen–Schoenberg–Whitney) $Q_n$ is
real-rooted with nonpositive roots iff the finite coefficient sequence
$(q_{n,0},\dots,q_{n,\lfloor n/j\rfloor})$ is PF$_\infty$, i.e. all minors of
its lower-triangular Toeplitz matrix are nonnegative. Writing
$q_{n,k}=p_{n,k}+c_n\bigl(c_{N,k}-c_{N,k-1}\bigr)$ turns this into an
*explicit weighted bilinear inequality on the two coefficient arrays*.

Two facts show that no cheap version of (WBI-1) suffices:

* **The convolution is not real-rooted.** For the actual families, $C_N$ has
  nonreal roots already at $j=2$, $n=3$ and $j=3$, $n=8$ (checked exactly by
  Sturm count in `probe1.py`). Hence $C_N$ is not an interlacing-family member
  and the MSS/Obreschkoff positive-combination closure does not apply.
* **The sign vector is not alternating.** For $j=3$, $n=16$ the actual
  signs of $C_N$ at the eight roots of $P_{16}=O_{16,2}$ are
  $(+,-,-,+,+,+,+,+)$; for $j=2$, $n=14$ they are
  $(-,-,-,-,-,+,-,+,+,+,+,+,+)$ (`probe_cross.py`, numerical for display only).
  Thus the required inequality is genuinely weighted/global, not a pairwise
  root-separation statement.

So the exact remaining task is: **prove, for the actual rows only, the
PF$_\infty$ (or equivalent Bezoutian-PSD) inequality for the coefficient array
generated by $q_{n,k}=p_{n,k}+c_n(c_{N,k}-c_{N,k-1})$, using the full history
of the actual rows — not merely their real-rootedness and interlacing.**
The counterexample families of Sections 2–3 satisfy all the generic
hypotheses and violate this inequality, which pinpoints why the history is
indispensable.

---

## 6. Files and verification

* `NOTES.md` — this report.
* `renewal_lib.py` — exact EGF/O-row generation, binomial convolution,
  triangular renewal, exact Sturm root counting.
* `verify_all.py` — **exact certificates** for Theorem A (discriminant $-3$),
  Theorem B (degree/PF/Wronskian interlacing/Sturm $2<4$, negative
  coefficient), and an independent check of (R2) for $j=2..6$, $n\le14$.
  Run: `py -3.14 verify_all.py`. All assertions pass.
* `verify_obstruction.py` — the detailed Theorem-B certificate with printed
  roots and coefficients.
* `probe1.py`, `probe_cross.py`, `probe_minors.py`, `probe_cancel.py`,
  `probe_cone.py`, `probe_counter.py`, `probe_sharp.py` — exploratory scripts
  (numerical displays are labelled as evidence; no float is used as a
  certificate).

No large caches, installs, model weights, or builds were created. Total size
is well under 1 MB; `__pycache__` is gitignored by the repository.
