# First D2C premise: rigorous `a=1` case, all-index reduction, and the exact remaining gap

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-j2-d2c-first` was created/edited.
No installs, builds, external agents, models, or browser use. SymPy (already
installed) is used only for exact symbolic residuals and explicitly labelled
finite checks.

## 0. Outcome

**Rigorous theorem proved here.** Let $n\ge4$ and $a=1$. In the $x=-t$ variables
of `proof-j2-sign-inequality/NOTES.md` §7, at every root $x>0$ of
$V=Q_{n,2}(-x)$,

\[
x\,(V')^2 + n\,x\,B(V,Y) \;>\; B(V,R),
\tag{7.3-a1}
\]

with $R=Q_{n,1}(-x)$, $Y=\mathcal U Q_{n-1,2}(-x)$, $B(F,G)=FG'-F'G$.
Equivalently $W_t(Q_{n+2,1},Q_{n,2})>0$ on $\mathbb R$, i.e. the **first D2C
premise holds for $a=1$, every $m=n+2\ge6$**.

**All-index reduction (rigorous, every $a\ge1$).** For general $a$, inequality
(7.3) is *equivalent* to the single bound

\[
\boxed{\,Z:=x^2 R''/R < 3a+2\quad\text{at every root }x>0\text{ of }V\,}
\tag{Z}
\]

(derived in §4). For $a=1$ this follows from Laguerre. For $a\ge2$ the
Laguerre bound $Z<(a+1)^2$ is **strictly weaker** than (Z), and no proved
constraint strengthens it. Accordingly D2C-1 for $a\ge2$ remains open here.

**New rigorous lemma (all $a$).** Every coefficient of
$W(Q_{m,a},Q_{m-2,a+1})$ is strictly positive (§3). This is *not* sufficient
to prove (7.3) (§5 gives an exact witness).

No retraction is needed: every displayed inequality is checked in §6.

---

## 1. Setup and conventions

Use exactly the conventions of `proof-j2-sign-inequality/NOTES.md` §1 and §7.
With $m\ge6$, $n=m-2\ge4$, $a\ge1$, $x=-t>0$:

\[
R=Q_{n,a}(-x),\quad V=Q_{n,a+1}(-x),\quad H=Q_{n-1,a+1}(-x),
\]
\[
Y=\mathcal U H=2H+xH',\qquad S=Q_{n+1,a}(-x),\qquad P=Q_{n+2,a}(-x),
\]
where $\mathcal Uf=2f+xf'$ and primes are $d/dx$. The Wronskians are
$B(F,G)=FG'-F'G$, so $B(F,G)=W_t(f,g)$ when $F=f(-x)$, $G=g(-x)$.

Established inputs (zero residual, `verify.py` V1):

\[
\text{(M)}\quad Q_{m,a}=(tQ_{m-1,a})'+(m-1)a\,t\,Q_{m-2,a+1},
\]
\[
\text{(V)}\quad Q_{m,a+1}=Q_{m,a}+\tfrac1a(tQ_{m,a})',
\quad\text{i.e.}\quad aV=(a+1)R+xR',
\]
\[
\text{(H)}\quad kA_{r,k+1}A_{r-1,k}\ge(k+1)A_{r,k}A_{r-1,k+1},
\]
proved unconditionally in `proof-u-operator/NOTES.md` §7;

and the coupled interlacing (`proof-u-operator/NOTES.md` §7, after the PF
multiplier), in the $B$-orientation:

\[
B(V,H)=W_t(Q_{n,a+1},Q_{n-1,a+1})>0,
\qquad
B(V,Y)=W_t(Q_{n,a+1},\mathcal UQ_{n-1,a+1})>0.
\tag{C}
\]

For $n=4$, $H=a+1$ and $V=(a+1)-(a+1)(a+2)x$ are explicit, giving
$B(V,H)=B(V,Y)=(a+1)^2(a+2)>0$; so (C) holds for every $n\ge4$.

---

## 2. Exact reduction of the first premise (re-verified)

At a root $x>0$ of $V$, (M) applied to $P$ and $S$ gives the exact identity
(`verify.py` V4, generic residual $0$, 50/50 actual residuals $0$)

\[
\boxed{\frac Pa=aR+xV'-nxY\qquad(V=0).}
\tag{7.1}
\]

Multiplying by $V'$ and using $B(V,R)=-V'R$, $B(V,Y)=-V'Y$ at $V=0$ gives

\[
\frac{PV'}a=x(V')^2+nxB(V,Y)-aB(V,R).
\tag{7.3}
\]

Since $W_t(p,v)=B(P,V)=PV'$ at $V=0$ and $\deg p=\deg v+1$ in the
nonautomatic range, $W_t(p,v)>0$ on $\mathbb R$ is equivalent to the
positivity of the right-hand side of (7.3) at each root of $V$. This is the
target inequality.

---

## 3. New lemma: all coefficients of $W(Q_{m,a},Q_{m-2,a+1})$ are positive

**Theorem.** For $m\ge6$, $a\ge1$, every coefficient of
$W(Q_{m,a},Q_{m-2,a+1})$ is strictly positive (hence the Wronskian is positive
on $[0,\infty)$).

*Proof.* Write $A_{m,r}={m\brace r}_{\ge2}$ and $(a)_r=a(a+1)\cdots(a+r-1)$.
Since $(a+1)_s=\frac{a+s}{a}(a)_s$,

\[
W(Q_{m,a},Q_{m-2,a+1})
=\frac1a\sum_{r,s}(r-s)(a)_r(a)_s\,A_{m,r}A_{m-2,s}(a+s)\,t^{r+s-3},
\]

because $W(t^{r-1},t^{s-1})=(r-s)t^{r+s-3}$. For $r>s$ the paired
coefficient of $t^{r+s-3}$ is

\[
\frac{(a)_r(a)_s(r-s)}{a}\Big[(a+s)A_{m,r}A_{m-2,s}-(a+r)A_{m,s}A_{m-2,r}\Big].
\tag{3.1}
\]

It suffices to show $g(r):=\dfrac{A_{m,r}}{(a+r)A_{m-2,r}}$ is strictly
increasing on the active range $1\le r\le\lfloor m/2\rfloor$. By (H) applied
at rows $m$ and $m-1$,

\[
\frac{A_{m,r+1}}{A_{m,r}}\ge\frac{r+1}{r}\frac{A_{m-1,r+1}}{A_{m-1,r}},
\qquad
\frac{A_{m-1,r+1}}{A_{m-1,r}}\ge\frac{r+1}{r}\frac{A_{m-2,r+1}}{A_{m-2,r}},
\]

hence $\dfrac{A_{m,r+1}A_{m-2,r}}{A_{m,r}A_{m-2,r+1}}\ge\left(\dfrac{r+1}{r}\right)^2$.
Therefore

\[
\frac{g(r+1)}{g(r)}
\ge\Big(\frac{r+1}{r}\Big)^2\frac{a+r}{a+r+1}
\ge\Big(1+\frac2r\Big)\frac{a+r}{a+r+1}>1,
\]

since $\frac2r\ge\frac1{a+r}$ for $a\ge1,r\ge1$. Thus the bracket in (3.1)
is positive, and with $(a)_r(a)_s(r-s)>0$ every paired coefficient is
positive. For each exponent occurring in $W$ there is a pair with both
entries active, so every coefficient is strictly positive. $\blacksquare$

`verify.py` V2 checks the raw ratio inequality on 1698 exact pairs (0
violations) and V3 checks coefficient positivity in 84 exact cases.

---

## 4. Exact dominance criterion and the rigorous `a=1` proof

At a root $x>0$ of $V$, differentiate (V): $aV'=(a+2)R'+xR''$. Two exact
identities follow (no assumption beyond $V=0$):

\[
aB(V,R)=\frac{(a+1)(a+2)R^2}{x}-xRR'',
\qquad
xV'=-\frac{R}{a}\big[(a+1)(a+2)-Z\big],
\tag{4.1}
\]

where $Z=x^2R''/R$ is evaluated at the root. Put

\[
w:=(a+1)(a+2)-Z.
\]

Then

\[
x(V')^2=\frac{R^2}{a^2x}w^2,\qquad
aB(V,R)=\frac{R^2}{x}w,
\qquad
\boxed{x(V')^2-aB(V,R)=\frac{R^2}{a^2x}\,w\,(w-a^2).}
\tag{4.2}
\]

Laguerre's inequality for the real-rooted $R$,
$L[R]=(R')^2-RR''>0$ at a non-root of $R$ (roots of $R$ and $V$ are distinct),
gives $x^2L[R]/R^2=(a+1)^2-Z>0$, i.e.

\[
Z<(a+1)^2\qquad\Longrightarrow\qquad w>(a+1).
\tag{4.3}
\]

**Case $a=1$.** Then $(a+1)^2=4$ and $a^2=1$, so (4.3) gives $w>a+1=2>1=a^2$.
By (4.2), $x(V')^2-B(V,R)>0$. Adding $nxB(V,Y)>0$ from (C) proves (7.3-a1).
This is rigorous for every $n\ge4$ (the $n=4$ constants are in §1).

**General $a$.** From (4.2), $x(V')^2>aB(V,R)$ is equivalent to
$w>a^2$, i.e.

\[
Z<(a+1)(a+2)-a^2=3a+2.
\tag{Z}
\]

Laguerre only gives $Z<(a+1)^2$, and for $a\ge2$
$(a+1)^2-(3a+2)=a^2-a-1>0$; so (Z) is strictly stronger than the proved
bound. By (C) the term $nxB(V,Y)$ is positive, so (Z) is sufficient (and, by
the exact identity (4.2) together with $nxB(V,Y)>0$, the only missing
magnitude estimate).

**Integral form of $Z$.** (V) is the first-order system
$(x^{a+1}R)'=a\,x^aV$; since $R(0)=a$ is finite,

\[
R(x)=a\,x^{-(a+1)}\int_0^x s^a V(s)\,ds,\qquad x>0.
\tag{4.4}
\]

Then $Z=(a+1)(a+2)+\dfrac{axV'}{R}$ at a root of $V$, equivalently (Z) reads

\[
\left|\int_0^{x}s^aV(s)\,ds\right|
<\frac{x^{a+2}|V'(x)|}{a^2}\quad\text{at a root }x\text{ of }V.
\tag{4.5}
\]

No proved invariant supplies (4.5) for $a\ge2$.

---

## 5. Why the structural route cannot close `a>=2`

### 5.1 Coefficient positivity is insufficient

An exact witness (rational coefficients, `ratio_test.py`) shows that two
polynomials $f,g$ with only negative zeros, positive coefficients, and
strictly increasing coefficient ratios $a_k/b_k$ can still have
$W(f,g)$ with positive coefficients **and two negative real zeros**, e.g.

\[
f=t^3+8t^2+\tfrac{309}{25}t+\tfrac{54}{25},\qquad
g=t^2+\tfrac{64}{3}t+112,
\]
\[
W(f,g)=\tfrac1{75}\big(75t^4+3200t^3+37073t^2+134076t+100368\big),
\quad \operatorname{disc}\text{-type real zeros present}.
\]

So §3 does not imply (7.3). The associated-Stirling recurrence (M) is
essential; §4 is exactly where it enters.

### 5.2 Laguerre is insufficient for `a>=2`

The only proved pointwise bound on $Z$ is (4.3). For $a\ge2$ it does not
reach $3a+2$ (§4). Higher Turán inequalities introduce $V''$ through
$aV''=(a+3)R''+xR'''$ and are not controlled by the established coupled
theorems. The artificial rational countermodel of
`proof-j2-distance-two/NOTES.md` §5 satisfies all structural hypotheses but
fails the first premise; it is not feasible for (M), so it is an obstruction
to the *structural method*, not a counterexample to the theorem.

### 5.3 What remains

An all-index proof for $a\ge2$ needs a new invariant bounding the second
logarithmic derivative of $R$ at roots of $V$: precisely (Z)/(4.5). Equivalently
it needs a distance-two coupled induction stronger than the two separate
coupled interlacings in (C). No such invariant was found here.

**Natural next step (not a proof).** Writing $u_i=x/(x-\rho_i)$ for the roots
$\rho_i$ of $R$ and $x$ a root of $V$, one has $\sum u_i=-(a+1)$ and
$Z=(a+1)^2-\sum u_i^2$. At the smallest root of $V$ all $u_i$ are negative and
the nearest root of $R$ dominates; the recurrence forces the sharp gap
$\rho_1-\sigma_1\approx\sigma_1/(a+1)$, which is why numerically
$Z_{\max}\approx0.05a$ while the threshold is $3a+2$ (verified for
$a\le200$, $m\le26$). A rigorous root-gap estimate would close (Z) with a
large margin; none is established here.

---

## 6. Proof-policy audit

### 6.1 Direction audit

- (7.3) orientation: at a root of $V$, $B(P,V)=PV'$, and
  $W_t(p,v)(\lambda)=-p(\lambda)v'(\lambda)>0$ means $p,v'$ have the same
  sign. Exact sample $(m,a)=(6,1)$, $\lambda=-1/9$: `verify.py` V6 gives
  $W_t(p,v)=62>0$. At $x=1/9$: $V'= -18$, $R=1/3$, $Y=4$, hence
  $x(V')^2=36$, $B(V,R)=-V'R=6$, $B(V,Y)=-V'Y=72$, and with $n=4$,
  $a=1$ the right-hand side of (7.3) is
  $36+4\cdot(1/9)\cdot72-6=36+32-6=62=PV'/a$. Direction $>0$ confirmed.
- (4.2): $w>(a+1)>0$, so the sign of $x(V')^2-aB(V,R)$ is the sign of
  $w-a^2$. For $a=1$ this is $>0$; no reversed chain is used.
- (3.1): bracket positivity proved link-by-link from (H) at rows $m$ and
  $m-1$; each link direction checked in V2.

### 6.2 Monotonicity audit

The only monotonicity claim is that $g(r)=\frac{A_{m,r}}{(a+r)A_{m-2,r}}$ is
strictly increasing in $r$; the derivative is replaced by the exact ratio
$g(r+1)/g(r)$ and printed $\ge(1+2/r)\frac{a+r}{a+r+1}>1$ in §3. The endpoint
that gives the smallest value is $r=1$, confirming the upper/lower direction.
No other monotonicity claim is made; no endpoint optimization is used.

### 6.3 Feasible-constraint falsification

Using only the proved constraints (real-rootedness, positive coefficients,
ratio monotonicity, coupled interlacings), §5.1 gives an exact witness with
$W>0$ on $[0,\infty)$ but with real zeros, so the constraints do not force
(7.3). This is reported as a failure of the structural route, not hidden.

### 6.4 Boundary and degeneracy cases

- $m\le5$: $v=Q_{m-2,a+1}$ is constant for $m=5$; (D2C-1) is automatic; no
  roots of $V$, so (7.3) is vacuous.
- $m=6$ ($n=4$): $H,Y$ are constants; (C) holds explicitly (§1); the $a=1$
  proof applies.
- $a=1$: proved in §4.
- $a\ge2$: reduced to (Z), open here.
- $R(\text{root of }V)\ne0$ (no common roots by interlacing), so $Z$ is
  well defined; no division by zero is hidden.
- Degrees: $\deg P=\deg V+1$ in the nonautomatic range, so the Wronskian
  orientation is unambiguous; equal degrees occur only within the same $m$.

### 6.5 Symbolic residuals

`verify.py` prints: V1 zero residuals for (M) and (V) on 60 exact cases;
V2 0 violations of the ratio inequality on 1698 pairs; V3 84/84 coefficient
positivity; V4 generic residual $0$ for (7.1) and 50/50 actual residuals;
V5 and V6 are explicitly labelled finite evidence. `dom2.py` re-verifies
(7.3) numerically with a single consistent $x$ symbol (an earlier script with
two distinct $x$ symbols produced a spurious SOS; that artifact is discarded).

---

## 7. Status and files

- **Proved:** (7.3) for $a=1$, all $n\ge4$; all-index reduction to (Z);
  all-index coefficient positivity of the target Wronskian.
- **Open:** (7.3) for $a\ge2$, hence D2C/I2 in full generality.
- Files: `NOTES.md` (this report), `verify.py` (exact residuals and labelled
  finite checks), `explore.py`, `generic.py`, `ops.py`, `probe.py`,
  `coeffs.py`, `shift.py`, `ztest.py`, `dom.py`, `dom2.py`, `rem.py`,
  `sos_id.py`, `sqid.py`, `ratio_test.py` (explorations; not proofs).

From this directory: `py -3.14 verify.py`. No packages installed, no build.
