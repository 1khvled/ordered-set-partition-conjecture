# Operator lemma: $\mathcal U$ preserves interlacing (all $r$)

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-u-operator`. No installs/builds.

## 1. Statement

Let $p,q\in\mathbb R[x]$ have strictly interlacing negative zeros in one of the two
orientations realised by the factorial row:

- **Case A** ($\deg p=\deg q+1$): $p$ leftmost and rightmost,
  $\alpha_1<\beta_1<\alpha_2<\beta_2<\cdots<\alpha_n<\beta_n<\alpha_{n+1}<0$.
- **Case B** ($\deg p=\deg q$): $q$ leftmost and rightmost,
  $\beta_1<\alpha_1<\beta_2<\alpha_2<\cdots<\beta_n<\alpha_n<0$.

Set $\mathcal U f:=2f+xf'=(x^2f)'/x$. Then
\[
P:=\mathcal Up,\qquad Q:=x\,\mathcal Uq
\]
have strictly interlacing negative zeros in the corresponding orientation; equivalently
\[
\boxed{W(P,Q)(x)\le 0\quad\text{for all }x\in\mathbb R.}\tag{L}
\]
(Exact checks give $W(P,Q)(0)=-4p(0)q(0)<0$ and no real zeros, so in fact $<0$;
the inequality $\le0$ is what the induction needs.)

In the factorial-row application, $p=F_{r-2}$, $q=F_{r-3}$, so (L) is exactly the
missing input of `proof-factorial-row/NOTES.md` §5.

## 2. Differentiation preserves interlacing

**Lemma (Hermite–Biehler/Gauss–Lucas).** If $F,G\in\mathbb R[x]$ have strictly
interlacing real zeros, then $F',G'$ have strictly interlacing real zeros.

*Proof.* By the Hermite–Biehler theorem, $F+iG$ has all zeros in the open upper
half-plane $\mathbb H$ iff $F,G$ have real, strictly interlacing zeros (with the
appropriate degree/leading-coefficient orientation). If all zeros of $F+iG$ lie in
$\mathbb H$, then all zeros of $(F+iG)'=F'+iG'$ lie in the convex hull of the zeros
of $F+iG$ by Gauss–Lucas; since $\mathbb H$ is convex, they lie in $\mathbb H$.
Applying Hermite–Biehler to $F'+iG'$ gives that $F',G'$ have strictly interlacing
zeros. $\square$

A limiting version is used below: if $F_\varepsilon,G_\varepsilon$ strictly interlace
for all small $\varepsilon>0$ and converge coefficientwise to $F,G$, then
$W(F,G)\le0$ on $\mathbb R$ and $F,G$ have (possibly non-strict) interlacing zeros.

## 3. From $\mathcal U$ to differentiation

Write
\[
A(x)=x^2p(x),\qquad B(x)=x^2q(x).
\]
The roots of $A$ are the roots of $p$ together with $0$ (double); those of $B$ are
the roots of $q$ together with $0$ (double). In both cases $A,B$ have generalized
interlacing roots ending in the common double root $0$, with the same tail
orientation: the last non-zero root before $0$ is a $p$-root, so the four roots at
$0$ must alternate $q,p,q,p$ when made distinct.

For small $\varepsilon>0$ set
\[
a_1=-3\varepsilon,\ a_2=-\varepsilon,\qquad b_1=-4\varepsilon,\ b_2=-2\varepsilon,
\]
and define
\[
A_\varepsilon(x)=p(x)\prod_{j=1}^2(x-a_j),\qquad
B_\varepsilon(x)=q(x)\prod_{j=1}^2(x-b_j).
\]
Because $b_1<a_1<b_2<a_2<0$ and the original non-zero roots are strictly interlacing,
the root sets of $A_\varepsilon,B_\varepsilon$ are strictly interlacing for all
sufficiently small $\varepsilon>0$. By Lemma §2, $A'_\varepsilon$ and $B'_\varepsilon$
have strictly interlacing zeros, hence
\[
W(A'_\varepsilon,B'_\varepsilon)(x)<0\quad\text{for all }x\in\mathbb R.
\]

Letting $\varepsilon\to0$, $A'_\varepsilon\to A'$ and $B'_\varepsilon\to B'$
coefficientwise. The limiting version of §2 gives
\[
W(A',B')(x)\le0\quad\text{for all }x\in\mathbb R,
\]
so $A'$ and $B'$ have generalized interlacing zeros.

Now
\[
A'(x)=\bigl(x^2p\bigr)'=x\,(2p+xp')=xP(x),
\]
\[
B'(x)=\bigl(x^2q\bigr)'=x\,(2q+xq')=Q(x)=x\,\mathcal Uq(x).
\]
Thus $A'$ and $B'$ are exactly $xP$ and $Q=x\,\mathcal Uq$. Their generalized
interlacing (with the common root $0$ at the end) is precisely the desired
interlacing of $P$ and $Q$: deleting the common factor $x$ from $A'$ leaves $P$,
while $B'$ remains $Q=x\,\mathcal Uq$.

By the interlacing criterion (Lemma 1 of `proof-factorial-row/NOTES.md`, in its
non-strict form), $W(P,Q)\le0$ on $\mathbb R$. This proves (L).

*Remark.* One also has $W(A',B')=x^2W(\mathcal Up,\mathcal Uq)$, so the limit gives
$W(\mathcal Up,\mathcal Uq)\le0$ for $x\neq0$; but the conclusion needed for the
induction is the stronger statement $W(\mathcal Up,x\,\mathcal Uq)\le0$, which is
exactly the interlacing of $A'$ and $B'$ proved above.

## 4. Closing the factorial-row induction

With (L) proved for all admissible $p,q$, the induction in
`proof-factorial-row/NOTES.md` §3–§5 closes without gaps:

- Induction hypothesis gives $F_{r-1},F_{r-2}$ real-rooted with negative roots.
- §4.1 gives $B_1<0$; §4.2 gives $W(C,D_C)\le0$.
- (L) (proved above) gives $W(\mathcal UF_{r-2},x\mathcal UF_{r-3})\le0$.
- Hence $W_r=B_1+(r-1)[x^2W(C,D_C)+(r-2)x^2W(\mathcal UF_{r-2},x\mathcal UF_{r-3})]<0$.
- By Lemma 1 of `proof-factorial-row/NOTES.md`, $F_r$ (hence $E_r$) has only real
  non-positive zeros.
- Newton's inequalities then give (N') and (SLC) as already proved in
  `proof-associated-stirling-lit/NOTES.md` §3.

Thus the factorial-weighted 2-associated row is fully real-rooted for all $r$, and
SLC/(H) follows. $\blacksquare$

## 6. Explicit deletion step and audit of an external objection

An external audit raised two objections. Both concern a different operator than
the one the lemma actually needs.

First, it is correct that $F_r=(1+xD)F_{r-1}+(r-1)x\,\mathcal U F_{r-2}$ and
that $(1+xD)A\ne\mathcal UA$. But the missing input in
`proof-factorial-row/NOTES.md` is $W(\mathcal U F_{r-2},x\mathcal U F_{r-3})\le0$,
not interlacing of $(xF_{r-1})'$ with $x\mathcal U F_{r-2}$. The $B_2$
decomposition already isolates exactly $\mathcal U F_{r-2}$ and
$x\mathcal U F_{r-3}$, so the operator mismatch does not affect the lemma.

Second, the counterexample $q=x+1$, $p=x+9/10$ shows that $q\preceq p$ does
not imply $p\preceq x\mathcal Uq$. That implication is not used. What the
lemma needs is $\mathcal Up\preceq x\mathcal Uq$. For this counterexample
$\mathcal Up=3x+9/5$ has root $-3/5$, while $x\mathcal Uq=x(3x+2)$ has roots
$-2/3$ and $0$; hence $\mathcal Up$ does interlace $x\mathcal Uq$, and
$W(\mathcal Up,x\mathcal Uq)=-\frac95(5x^2+6x+2)<0$ with discriminant
$-324/25$. The script `check_audit.py` verifies this exactly, together with
the strict interlacing of the perturbed pair.

For completeness, here is the deletion step omitted in §3. Write
$R=\mathcal Uq$, so $B'=Q=xR$. The limit gives generalized interlacing of
$A'=xP$ and $B'=xR$. They share the root $0$, and every other root is negative.
Deleting the common root $0$ from both gives generalized interlacing of $P$
and $R$. Appending $0$ to $R$ as its rightmost root gives generalized
interlacing of $P$ and $Q=xR$ in every degree/orientation case: if $P$ has one
more root than $R$, the pattern is $p_1<r_1<\cdots<p_{d+1}<0(q)$; if the
degrees are equal, either $0(q)$ is rightmost or $P$ and $R$ already alternate
with $0(q)$ appended. Hence $W(P,Q)\le0$ by the interlacing criterion, which
is exactly (L).
