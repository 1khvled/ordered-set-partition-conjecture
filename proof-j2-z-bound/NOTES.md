# The `Z < 3a+2` bound: a rigorous root-gap theorem at the smallest root, and the exact reduction at every root

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-j2-z-bound` was created/edited.  
No installs, builds, external agents, models, or browser use. SymPy (already
installed) is used only for exact symbolic residuals and explicitly labelled
finite checks.

---

## 0. Outcome (honest status)

**Rigorous theorem proved here (Theorem A).** Let $n\ge4$, $a\ge1$, and put
$R=Q_{n,a}(-x)$, $V=Q_{n,a+1}(-x)$. Let $\rho_1$ be the smallest positive root
of $R$ and $\sigma_1$ the smallest positive root of $V$. Then

$$
\boxed{\ \frac{\sigma_1}{\rho_1}\ \ge\ e^{-1/(a+1)}\ >\ \frac{a}{a+1}\ }
\tag{G1}
$$

and consequently (Theorem B below)

$$
Z(\sigma_1)=\frac{\sigma_1^2R''(\sigma_1)}{R(\sigma_1)}\ <\ 2a+1\ <\ 3a+2 .
\tag{Z1}
$$

This is the first rigorous root-gap estimate of the kind the preceding note
identified as the missing ingredient: the note's "natural next step"
(`proof-j2-d2c-first/NOTES.md` §5.3) is now a theorem **at the smallest root**,
which is the root at which $Z$ attains its maximum in all finite checks.

**Rigorous reduction (Theorem B).** At *every* positive root $\sigma_j$ of $V$,
if $\sigma_j\ge\frac{a}{a+1}\rho_j$ (where $\rho_j$ is the corresponding root of
$R$) then $Z(\sigma_j)<2a+1<3a+2$. So the whole remaining gap is exactly the
statement (G) $\sigma_j/\rho_j\ge a/(a+1)$ for all $j$.

**What is *not* proved here.** The estimate (G) at indices $j\ge2$. The
coefficient/continued-fraction method that proves (G1) needs the positivity of
$K=(xR')'/R'$ along the whole root-flow path; at the first root this is an
exact coefficient-ratio consequence, but at later roots it is an interlacing
statement that is **not** implied by the proved structural constraints alone
(exact rational witnesses are given in §6). Finite checks through $n=26$,
$a\le 200$ find $\sigma_j/\rho_j\ge a/(a+1)$ at every root, but that is labelled
finite evidence, not a proof. I therefore do **not** claim D2C-1 for general
$a\ge2$ here.

No retraction is needed: every displayed identity and inequality is checked in
§7, including all directions (policy item 1).

---

## 1. Setup, conventions, and exact inputs

Use the conventions of `proof-j2-d2c-first/NOTES.md` §1. For $m\ge6$, $n=m-2\ge4$, $a\ge1$, $x=-t>0$:

$$
R(x)=Q_{n,a}(-x),\quad V(x)=Q_{n,a+1}(-x),\quad Y=\mathcal UQ_{n-1,a+1}(-x),
$$

with $\mathcal Uf=2f+xf'$, $B(F,G)=FG'-F'G$, and $Z=x^2R''/R$ evaluated at a
root of $V$.

Established inputs (zero residual, `verify.py` V5):

$$
\text{(M)}\quad Q_{m,a}=(tQ_{m-1,a})'+(m-1)a\,t\,Q_{m-2,a+1},
$$
$$
\text{(V)}\quad a\,Q_{m,a+1}=(a+1)Q_{m,a}+t\,Q_{m,a}' .
$$

Also used: the associated-Stirling recurrence $A_{m+1,r}=rA_{m,r}+mA_{m-1,r-1}$
(`verify.py` V3), real-rootedness of every $Q_{m,a}$, and **Newton's inequality**
for polynomials with positive coefficients and real negative roots.

The exact reduction already established in `proof-j2-d2c-first/NOTES.md` §4 is

$$
\text{at a root of }V:\qquad xR'=-(a+1)R,\qquad xV'=-\frac{R}{a}\bigl[(a+1)(a+2)-Z\bigr],
\tag{1.1}
$$
$$
Z=(a+1)^2-\sum_i\frac{x^2}{(x-\rho_i)^2},
\tag{1.2}
$$

where $\rho_i>0$ are the roots of $R$ and $x$ is any point with $R(x)\ne0$.
(1.2) is the partial-fraction identity
$R''/R=(R'/R)^2-\sum_i(x-\rho_i)^{-2}$ together with (1.1).

---

## 2. Coefficient structure and the ratio lemma

Write

$$
\frac{R(x)}{a}=\sum_{j\ge0}(-1)^j\alpha_j x^j,\qquad
\alpha_j=(a+1)_j\,A_{n,\,j+1},\quad \alpha_0=1 .
\tag{2.1}
$$

Indeed $Q_{n,a}(t)=\sum_{r\ge1}(a)_rA_{n,r}t^{r-1}$ and $(a)_{j+1}=a(a+1)_j$, so
$Q_{n,a}(-x)=a\sum_j(-1)^j(a+1)_jA_{n,j+1}x^j$. Let

$$
r_j:=\frac{\alpha_j}{\alpha_{j+1}}=\frac{A_{n,j+1}}{(a+j+1)A_{n,j+2}} .
\tag{2.2}
$$

**Lemma 2.1.** $r_0<r_1<\cdots<r_{k-1}$ on the active range, where
$k=\deg R=\lfloor n/2\rfloor-1$.

*Proof.* Apply Newton's inequality to the coefficient sequence
$c_j=(a)_{j+1}A_{n,j+1}=a\alpha_j$ of the real-rooted $Q_{n,a}(t)$:
$c_{j+1}^2\ge c_jc_{j+2}\cdot\frac{j+2}{j+1}\cdot\frac{k-j}{k-j-1}$ for the
interior indices. Hence
$$
\frac{r_{j+1}}{r_j}=\frac{\alpha_{j+1}^2}{\alpha_j\alpha_{j+2}}
=\frac{c_{j+1}^2}{c_jc_{j+2}}
\ge\frac{(j+2)(k-j)}{(j+1)(k-j-1)}\ >\ 1 .
$$
$\blacksquare$ (`verify.py` V1 checks 500 exact instances, 0 violations.)

**Lemma 2.2.** With $\mu:=\dfrac{\alpha_1^2}{\alpha_2}
=\dfrac{(a+1)A_{n,2}^2}{(a+2)A_{n,3}}$, one has $\mu>6$ for every $n\ge6$, $a\ge1$.

*Proof.* Since $\frac{a+1}{a+2}\ge\frac23$, it suffices to show
$A_{n,2}^2>9A_{n,3}$. We use the exact recurrences (verified in `verify.py` V3)

$$
A_{m+1,2}=2A_{m,2}+m,\qquad A_{m+1,3}=3A_{m,3}+mA_{m-1,2}.
$$

Base: $m=6$: $A_{6,2}=25$, $A_{6,3}=15$, $625>135$; $m=7$:
$A_{7,2}=56$, $A_{7,3}=105$, $3136>945$. Induction step for $m\ge7$,
assuming the claim at $m$ and $m-1$:

$$
A_{m+1,2}^2=(2A_{m,2}+m)^2=4A_{m,2}^2+4mA_{m,2}+m^2,
$$
$$
9A_{m+1,3}=27A_{m,3}+9mA_{m-1,2}<3A_{m,2}^2+9mA_{m-1,2}.
$$

Since $A_{m,2}=2A_{m-1,2}+(m-1)>2A_{m-1,2}$, it is enough that
$4A_{m-1,2}^2+8mA_{m-1,2}\ge 9mA_{m-1,2}$, i.e. $4A_{m-1,2}\ge m$.
As $A_{m-1,2}=2^{m-2}-m$, this is $2^m\ge5m$, true for $m\ge5$.
$\blacksquare$ (`verify.py` V2: 2500 exact instances, 0 violations; minimum
$\mu=250/9$ at $n=6,a=1$.)

---

## 3. The first-root gap: proof of Theorem A

### 3.1 A coefficient estimate for $\rho_1$

**Lemma 3.1.** $\rho_1\le r_1/4$.

*Proof.* Put $t:=r_1/4=\alpha_1/(4\alpha_2)$, so
$\alpha_1t=\alpha_1^2/(4\alpha_2)=\mu/4$. By Lemma 2.1, $r_j\ge r_1$ for
$j\ge1$, hence $\alpha_j\le\alpha_1/r_1^{\,j-1}$ and therefore

$$
\alpha_j t^j\le \alpha_1\frac{r_1^{\,j}}{r_1^{\,j-1}4^j}
=\frac{\alpha_1r_1}{4^j}=\frac{\mu}{4^j}\qquad(j\ge2).
$$

Thus, with $p(x)=R(x)/a=1-\alpha_1x+\alpha_2x^2-\cdots$,

$$
p(t)\ \le\ 1-\frac{\mu}{4}+\mu\sum_{j\ge2}4^{-j}
=1-\frac{\mu}{4}+\frac{\mu}{12}=1-\frac{\mu}{6}<0
$$

by Lemma 2.2. Since $p(0)=1>0$ and $\rho_1$ is the first positive zero of $p$,
$p(t)<0$ forces $\rho_1<t=r_1/4$. $\blacksquare$ (`verify.py` V4.)

### 3.2 Concavity of $xR'$ on $(0,\rho_1)$

**Lemma 3.2.** $(xR')'(x)=R'(x)+xR''(x)<0$ for every $x\in(0,\rho_1)$.

*Proof.* From (2.1),
$\frac1a(xR')'=-\alpha_1+4\alpha_2x-9\alpha_3x^2+\cdots$, an alternating
series whose $j$-th magnitude is $T_j=j^2\alpha_jx^{j-1}$. For $j\ge1$,

$$
\frac{T_{j+1}}{T_j}=\Bigl(\frac{j+1}{j}\Bigr)^2\frac{x}{r_j}
\le 4\frac{\rho_1}{r_1}<1
$$

by Lemma 2.1 (so $r_j\ge r_1$) and Lemma 3.1 (so $\rho_1/r_1<1/4$).
The magnitudes are strictly decreasing, so the alternating series has the
sign of its first term, $-\alpha_1<0$. $\blacksquare$

Also $R'<0$ on $(0,\rho_1)$ (because $R'/R=\sum_i(x-\rho_i)^{-1}<0$ there and
$R>0$). Hence

$$
K(x):=\frac{(xR')'(x)}{R'(x)}>0\qquad (0<x<\rho_1).
\tag{3.1}
$$

### 3.3 The multiplier flow and the gap

For $\theta\ge0$ define

$$
P_\theta(x):=\sum_{j\ge0}(-1)^j(1+\theta j)\alpha_j x^j .
\tag{3.2}
$$

By (V), $\;V=(a+1)P_{1/(a+1)}$: indeed
$\frac{1}{a+1}V=\sum_j(-1)^j\frac{a+j+1}{a+1}\alpha_jx^j
=\sum_j(-1)^j(1+\frac{j}{a+1})\alpha_jx^j$.

The sequence $m_j=1+\theta j$ is a multiplier sequence: its exponential
generating function $\sum_j m_j z^j/j!=(1+\theta+\theta z)e^z$ lies in the
Laguerre–Pólya class. Hence every $P_\theta$ has only real zeros. Let
$z(\theta)$ be its smallest positive zero; $z(0)=\rho_1$ and
$z(1/(a+1))=\sigma_1$. On any interval where the zeros are simple (they move
continuously and cannot collide under a strictly monotone flow), implicit
differentiation gives

$$
\frac{d z}{d\theta}=-\frac{\partial_\theta P_\theta(z)}{\partial_xP_\theta(z)}
=-\frac{z}{1+\theta\,\frac{(xR')'(z)}{R'(z)}}
=-\frac{z}{1+\theta K(z)} .
$$

For $z\in(0,\rho_1)$ we have $K(z)>0$ by (3.1), hence

$$
\frac{d\log z}{d\theta}=-\frac{1}{1+\theta K(z)}\ \ge\ -1 .
$$

Integrating $\theta$ from $0$ to $1/(a+1)$:

$$
\log\frac{\sigma_1}{\rho_1}\ \ge\ -\frac{1}{a+1},
\qquad\text{i.e.}\qquad
\frac{\sigma_1}{\rho_1}\ \ge\ e^{-1/(a+1)}>\frac{a}{a+1},
\tag{G1}
$$

using $e^{-u}>1-u$ for $u>0$. This proves Theorem A.

**Direction audit.** $\partial_\theta P_\theta=xR'$; at $z\in(0,\rho_1)$, $R'<0$,
so $\partial_\theta P_\theta<0$ and $z'=-\partial_\theta P/\partial_xP<0$: roots
move **left** as $\theta$ increases, so $\sigma_1<\rho_1$, consistent with
(1.1). The inequality $z'/z\ge-1$ was checked link by link: denominator
$1+\theta K\ge1$ because $K>0$ and $\theta\ge0$; inverting preserves direction
because all quantities are positive. No reversed step. (Policy item 1.)

---

## 4. From the gap to the bound: proof of Theorem B

Let $\sigma_j$ be any positive root of $V$, and let $\rho_j$ be the $j$-th
positive root of $R$. The roots of $V$ interlace those of $R$ (from
$V/R=\frac1a[(a+1)+xR'/R]$ and the monotonic decrease of $xR'/R$ on each gap),
so $\sigma_j\in(\rho_{j-1},\rho_j)$ with $\rho_0:=0$.

Set $u_i=\frac{\sigma_j}{\sigma_j-\rho_i}$. Then $u_i>0$ for $i<j$ and
$u_i<0$ for $i\ge j$; write $v_i=-u_i=\frac{\sigma_j}{\rho_i-\sigma_j}>0$ for
$i\ge j$. By (1.1), $\sum_i u_i=-(a+1)$, hence

$$
\sum_{i\ge j}v_i=a+1+\sum_{i<j}u_i\ >\ a+1 .
\tag{4.1}
$$

If $\sigma_j\ge\frac{a}{a+1}\rho_j$, then

$$
v_j=\frac{\sigma_j}{\rho_j-\sigma_j}
=\frac{\sigma_j/\rho_j}{1-\sigma_j/\rho_j}\ \ge\ \frac{a/(a+1)}{1/(a+1)}=a .
\tag{4.2}
$$

Therefore, by (1.2),

$$
Z(\sigma_j)=(a+1)^2-\sum_i u_i^2
\le (a+1)^2-v_j^2\ \le\ (a+1)^2-a^2=2a+1<3a+2 .
\tag{4.3}
$$

This proves Theorem B. Note that (4.1) is not even needed for (4.3); it is
recorded because it shows the bound is far from tight.

Combining Theorems A and B gives (Z) at the smallest root for every $n\ge4$,
$a\ge1$, hence the first D2C premise at the root where its margin is smallest
in all finite checks.

---

## 5. Why the method does not immediately give (G) at $j\ge2$

For the $j$-th root path of the flow (3.2), implicit differentiation gives the
same formula

$$
\frac{d\log z_j}{d\theta}=-\frac{1}{1+\theta K(z_j)},\qquad
K=\frac{(xR')'}{R'} .
$$

Here $z_j(\theta)\in(\sigma_j,\rho_j)$, not $(0,\rho_1)$. Lemma 3.2 was proved
on $(0,\rho_1)$ by a *first-term-dominance* alternating-series argument, which
fails for $x>\rho_1$. At a root of $P_\theta$ one has the exact identity

$$
K=1-\theta\,\frac{x^2R''(x)}{R(x)}\Big|_{R_0\text{-quantities}},
$$

so $K>0$ at every root of $P_\theta$ along the path is equivalent to a bound on
the second log-derivative of $R$ at those intermediate points — a statement of
the same strength as (Z) itself. Numerically $K>0$ holds throughout the path
(`K_probe.py`: minima $0.89$–$0.995$ for $n\le20$), and equivalently the exact
interlacing chain

$$
\sigma_j<\rho_j<\zeta_j<\eta_j<\sigma_{j+1}
$$

holds, where $\eta_j$ are the roots of $R'$ and $\zeta_j$ those of $(xR')'$
(`eta_probe.py`, `zeta_probe.py`). But this chain is **not** a consequence of
real-rootedness, positive coefficients, and ratio monotonicity: §6 gives exact
witnesses satisfying all of those constraints but violating the chain (hence
violating $K>0$ and (G)). Consequently the recurrence (M) must enter again at
every index, and the first-root argument does not transfer automatically.

This is a genuine gap, recorded rather than papered over.

---

## 6. Feasible-constraint falsification (policy item 3)

Exact rational / algebraic witnesses (all computations finite and exact;
`test_chain.py`, `test_chain2.py`, `test_multShift.py`, `test_general.py`,
`test_general2.py`, `test_concavity.py`):

1. For a PF polynomial with positive coefficients and strictly increasing
   coefficient ratios, the multiplier flow $1+\theta j$ can move the $j$-th
   root left **faster** than $e^{-\theta}$. Example (roots
   $0.0114,0.0247,0.0535$, $\theta=0.2867$, $a=2$): first root moves to
   $0.0080$ while $e^{-\theta}\rho_1=0.0075$, i.e. $\sigma_1<e^{-\theta}\rho_1$.
   So real-rootedness + increasing ratios alone do **not** give (G1); the
   coefficient estimate of §3 (which uses the large ratio $r_1\gg\rho_1$ from
   the associated-Stirling structure) was essential.

2. For strongly separated roots (ratio $\ge4$) the chain of §5 held in 0 of
   thousands of failures, but weakly separated roots break it; hence any
   all-index proof must use the specific root separation forced by (M), not
   just structural positivity.

3. For arbitrary real-rooted $f$ with positive coefficients,
   $-xv'/f$ at roots of $v=\frac1a[(a+1)f+xf']$ can be $<a$, and $Z$ can be
   positive at non-first roots. So (Z) is not a formal consequence of (V)
   alone; (M) is essential.

These witnesses are reported as failures of the *structural* route, not as
counterexamples to the theorem (they are not feasible for (M)).

---

## 7. Proof-policy audit

### 7.1 Direction audit (item 1)

- Lemma 2.1: $\frac{r_{j+1}}{r_j}\ge\frac{(j+2)(k-j)}{(j+1)(k-j-1)}>1$ —
  verified sample-wise in `verify.py` V1 with the stated direction.
- Lemma 2.2: strict inequality $\mu>6$; the chain $A_{n,2}^2>9A_{n,3}$ is
  proved by induction with each link's direction checked; endpoints $n=6,7$ done.
- Lemma 3.1: $p(t)\le1-\mu/6<0$ with $t=r_1/4$; since $p(0)=1>0$ and $p$
  crosses zero first at $\rho_1$, $\rho_1<t$. Sample: $p(t)\approx-4.86$ vs
  bound $-4.21$ (`verify.py` V4).
- Lemma 3.2: term ratio $\le4\rho_1/r_1<1$ — direction checked; endpoint
  $j=1$ gives the largest bound factor $4$.
- §3.3 flow: $z'/z=-1/(1+\theta K)\ge-1$; integration direction $\theta\ge0$
  gives $\sigma_1/\rho_1\ge e^{-1/(a+1)}$. Sample $(m,a)=(6,1)$: $\sigma_1=1/9$,
  $\rho_1=1/(3\cdot2)=1/6$, ratio $0.6667\ge2/3$ ✓; $Z=0<3a+2=5$ ✓.
- §4: $v_j\ge a\Rightarrow Z\le(a+1)^2-a^2=2a+1<3a+2$; check $a=2$:
  $5<8$ ✓ (V6).

### 7.2 Monotonicity audit (item 2)

- $r_j$ strictly increasing: Newton's inequality supplies the finite difference
  $\frac{r_{j+1}}{r_j}-1>0$; no derivative used. Endpoint $r_0$ is the minimum.
- Alternating magnitudes $T_j$ decreasing: exact ratio bound; endpoint $j=1$.
- Multiplier $m_j=1+\theta j$ increasing in $j$; endpoint $j=0$ gives factor 1.
- $\Phi(x)=e^{-x}$ monotonicity only in the elementary inequality
  $e^{-u}>1-u$ for $u>0$ (endpoint $u=0$ equality).

### 7.3 Boundary and degeneracy cases (item 4)

- $n=4,5$: $k=1$ (or $0$). For $k=1$, $\alpha_1$ exists and the single root of
  $V$ satisfies $v_1=a+1$ exactly, so $Z=0$; the coefficient argument still
  applies (check $n=4$: $\mu=3(a+1)^2/(2(a+2))>6$? at $a=1$,
  $\mu=12/6=2$ — **not** $>6$). This case is covered directly: with one root,
  $Z=0<3a+2$ by (1.2). Lemma 2.2 is only claimed/used for $n\ge6$; for
  $n=4,5$ the bound is automatic.
- $a=1$: proved already in `proof-j2-d2c-first/NOTES.md`; the present proof
  also covers it once $n\ge6$.
- $R$(root of $V$)$\ne0$ by interlacing, so $Z$ is well defined; no division by
  zero is hidden.
- $\deg R=\deg V=k$; the multiplier family is nonsingular at $\theta=0$ and
  $\theta=1/(a+1)$ (simple roots), and the flow argument is used only on
  intervals of simplicity.

### 7.4 Symbolic residuals (item 5)

`verify.py` prints: V1 500/500 ratio inequalities with 0 violations; V2 2500
exact checks $\mu>6$ with 0 violations; V3 zero residual for both
associated-Stirling recurrences; V4 the exact certificate $p(r_1/4)\le1-\mu/6$;
V5 zero residuals for (M) and (V); V6 direction checks for $2a+1<3a+2$.
All root computations (`gap_probe.py`, `allroots_ratio.py`, `probe.py`,
`K_probe.py`, etc.) are explicitly labelled **finite evidence**.

### 7.5 Retraction rule (item 6)

No invalid theorem is left in the record: §0 states explicitly that the
all-index case $j\ge2$ is open here. §5 records the exact point where the
first-root method stops being sufficient.

---

## 8. Files

- `NOTES.md` (this report).
- `ops.py` (exact Stirling table and polynomial operations).
- `verify.py` (exact residuals, ratio checks, direction checks).
- `explore.py`, `explore2.py`, `probe.py`, `check71.py`, `gap_probe.py`,
  `allroots_ratio.py`, `zeta_probe.py`, `eta_probe.py`, `test_general.py`,
  `test_general2.py`, `test_concavity.py`, `test_terms.py`, `test_terms2.py`,
  `test_zeta_general.py`, `test_chain.py`, `test_chain2.py`,
  `test_multShift.py`, `diag_probe.py`, `dbg.py` (explorations / finite
  evidence; not proofs).

From this directory: `py -3.14 verify.py`. No packages installed, no build.

---

## 9. Addendum: non-flow smallest-root proof and an exact later-root obstruction

This addendum was prepared after a server restart. It preserves the preceding
root-flow proof, but gives an independent recurrence/interlacing proof of the
smallest-root conclusion that does not use root flow. It also replaces the
qualitative structural falsification in §6.3 with a single exact rational
witness.

### 9.1 Direct non-flow proof at the smallest root

Let `sigma` be the smallest positive root of `V`, and write the positive roots
of `R` as `0<rho_1<...<rho_d`. Define

\[
p_j=\frac{\sigma}{\rho_j-\sigma}>0.
\]

The vertical identity gives `sum p_j=a+1`. Let `h_1<...<h_e` be the positive
roots of `H=Q_{n-1,a+1}(-x)`. The coupled interlacings imply `sigma<h_k` and

\[
\rho_1<h_1<\rho_2<\cdots<h_{d-1}<\rho_d
\]

when `e=d-1`, with one extra final `H`-root when `e=d`.

At the first root of `V`, `V'(sigma)<0`. Hence `B(V,H)>0` and `B(V,Y)>0`
give `H(sigma)>0` and `Y(sigma)>0`. Applying (M) to `V` at `sigma` gives

\[
H+\sigma H'=(n-1)(a+1)\sigma K=\frac{Y}{2}>0,
\qquad K=Q_{n-2,a+2}(-\sigma).
\]

Dividing by `H(sigma)>0` and taking the logarithmic derivative of the factored
polynomial `H` yields

\[
\sum_{k=1}^e\frac{\sigma}{h_k-\sigma}<1.
\tag{9.1}
\]

Since `h_{j-1}<rho_j` for `j>=2`,

\[
\sum_{j=2}^d p_j
=\sum_{j=2}^d\frac{\sigma}{\rho_j-\sigma}
<\sum_{j=2}^d\frac{\sigma}{h_{j-1}-\sigma}
\le \sum_{k=1}^e\frac{\sigma}{h_k-\sigma}<1.
\tag{9.2}
\]

Thus `p_1=a+1-sum_{j>=2}p_j>a`. With
`u_j=sigma/(sigma-rho_j)=-p_j`, the exact identity (1.2) gives

\[
Z(\sigma)=(a+1)^2-\sum_j u_j^2
<(a+1)^2-a^2=2a+1<3a+2.
\tag{9.3}
\]

For `d=1` the tail is empty, `p_1=a+1`, and `Z=0`; for `n=4`, `H` is constant
and (9.1) is the empty sum `0<1`. Therefore (9.3) covers every `n>=4,a>=1`.

This proof uses only (M), (V), and the coupled interlacings. It is the direct
recurrence-based root-gap estimate requested in the task, at the smallest root.

### 9.2 Exact obstruction at a non-first root

The following example is not a `Q`-polynomial and does not satisfy (M). It does
satisfy the vertical identity, has positive coefficients in `t`, and has simple
positive roots with strict vertical interlacing. Nevertheless it violates the
requested bound at a non-first root.

Let `a=20` and

\[
R(x)=(1-2x)\left(1-\frac{13x}{14}\right)
     \left(1-\frac{10x}{11}\right),
\qquad
V(x)=\frac{21R(x)+xR'(x)}{20}.
\]

The roots of `R` are `1/2,14/13,11/10`. At `x=1`,

\[
(u_1,u_2,u_3)
=\left(\frac1{1-1/2},\frac1{1-14/13},\frac1{1-11/10}\right)
=(2,-13,-10),
\]

so `sum u_j=-21=-(a+1)` and therefore `V(1)=0`. Exactly,

\[
Z(1)=21^2-(2^2+13^2+10^2)=441-273=168>62=3a+2.
\tag{9.4}
\]

The signs of `V` at `0,1/2,1,14/13,11/10` are `+,-,0,+,-`. Since `V` is
cubic, these sign changes place one root in each of
`(0,1/2),(1/2,14/13),(14/13,11/10)`; `x=1` is the middle root. Thus the
failure is not caused by non-real roots or broken vertical interlacing.

This is an exact obstruction to extending §9.1 using only (V), real-rootedness,
positive coefficients, and vertical interlacing. It does not obstruct the
actual associated-Stirling theorem, because (M) is absent. An all-root proof
must therefore use an additional (M)-specific estimate at `j>=2`; no such
estimate is proved in this directory.

### 9.3 Audit of the addendum

- **Direction:** (9.1) follows from `H+sigma H'>0`; (9.2) is strict because
  `h_{j-1}<rho_j`; (9.3) follows from `p_1>a`. Exact sample `n=4,a=2`:
  `sigma=1/12`, `rho=1/9`, `p_1=3>2`, `Z=0`.
- **Monotonicity:** no parameter monotonicity is used in §9.1; only the proved
  root order and a termwise tail comparison.
- **Falsification:** (9.4) is an exact rational/algebraic witness, not numerical
  evidence.
- **Boundary:** `d=1`, constant `H`, and `a=1` are covered explicitly.
- **Residuals:** `verify.py` V7 prints `V(1)=0`, `sum u_j=-21`, `Z(1)=168`,
  and the exact sign pattern; V8 prints the `n=4,a=2` direction sample.
- **Retraction:** none. The addendum does not retract any preceding theorem; it
  strengthens the recorded status by making the structural obstruction exact.
