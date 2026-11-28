# Exact missing sign lemma at roots of $UQ_{m-2,a+1}$ (distance-two D2C gap)

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-j2-qsign` was created/edited.  
No installs, builds, external agents, models, or browser use. SymPy (already
installed) is used only for exact symbolic residuals and exact/numeric finite
checks.

## 0. Outcome

**No all-index proof of the missing sign lemma, and no all-index proof of D2C,
was obtained.** The exact formula stated in `proof-j2-distance-two/NOTES.md` §6
is incorrect at roots of $UQ_{m-2,a+1}$, and the corrected exact condition is
not a sign on $Q_{m-1,a}(\lambda)$ alone: it is a magnitude inequality involving
$Q_{m-1,a+1}(\lambda)$, $Q_{m-1,a}(\lambda)$, and a nonnegative Laguerre term.
The previously proved facts do not imply that inequality; the rational §5
countermodel of `proof-j2-distance-two/NOTES.md` satisfies the coupled control
but violates the corrected total sign. This is an exact obstruction under the
proved hypotheses.

What is proved here is the corrected identity/decomposition and the exact sign
condition. Finite evidence on the actual triangle (through the ranges in
`verify.py`) finds no violation of D2C; this is labelled finite evidence, never
a proof.

## 1. Setup

Use the notation of `proof-j2-distance-two/NOTES.md`:
$A_{m,r}={m\brace r}_{\ge2}$,
$Q_{m,a}(t)=\sum_{r\ge1}(a)_r A_{m,r}t^{r-1}$,
$Uf:=2f+tf'$, $W(f,g)=f'g-fg'$.

Fix $m\ge5$, $a\ge1$. Put
\[
p:=Q_{m,a},\qquad s:=Q_{m-1,a},\qquad q_1:=Q_{m-1,a+1},\qquad
v:=Q_{m-2,a+1},\qquad u:=Uv=UQ_{m-2,a+1}.
\]
The D2C target is
\[
W(p,v)>0,\qquad W(p,u)>0\qquad\text{on }\mathbb R. \tag{D2C}
\]

Established facts used (same status as `proof-j2-distance-two/NOTES.md` §5,
`proof-j2-diagonal-after-h/NOTES.md` §3, `proof-u-operator/NOTES.md` §7):
all active coefficients positive; (H); real-rootedness and negative roots;
PF multipliers $\Lambda_a$; coupled interlacing
$Q_{k-1,b}\preceq Q_{k,b}$, $UQ_{k-1,b}\preceq Q_{k,b}$ in the orientation
$W(Q_{k,b},Q_{k-1,b})>0$, $W(Q_{k,b},UQ_{k-1,b})>0$; identities (V),(U),(M),(R).
The associated-Stirling recurrence (M) is essential (obstruction in
`proof-j2-distance-two/NOTES.md` §5).

## 2. Correction to the §6 formula (policy item 6)

`proof-j2-distance-two/NOTES.md` §6 states: "Let $\lambda$ be a root of $Uv$.
By (M), $p(\lambda)=(ts)'(\lambda)$." This is **false**. (M) is
\[
p=(ts)'+(m-1)a\,t\,v. \tag{M}
\]
At a root $\lambda$ of $u=Uv$, $v(\lambda)$ is not zero in general, hence
\[
\boxed{p(\lambda)=(ts)'(\lambda)+(m-1)a\,\lambda\,v(\lambda).} \tag{C1}
\]
`verify.py` V1 gives 150/150 counterexamples to the old formula and 0 symbolic
failures for (C1) over $5\le m\le15,1\le a\le5$ (symbolic identity check per
$(m,a)$, plus 150 numeric root evaluations). This retraction is recorded here;
no result in this directory uses the old formula.

The old formula is correct at roots of $v$ (the other D2C premise), since the
$(m-1)a t v$ term vanishes there.

## 3. Exact sign decomposition at a root of $u$

Let $\lambda$ be a real root of $u=Uv$. Since $\deg p=\deg u+1$ in the
nonautomatic cases, $W(p,u)>0$ is equivalent, by exact sign count, to
$-p(\lambda)u'(\lambda)>0$ at every root $\lambda$ of $u$.

Using (U), $(ts)'=Us-s=a(q_1-s)$. Also
\[
W(tv,u)= (tv)'u-tv u' = 2v^2+t^2W(v,v'), \tag{Wtv}
\]
an exact polynomial identity (symbolic residual 0 by direct expansion:
$(v+tv')(2v+tv')-tv(3v'+tv'')=2v^2+t^2((v')^2-vv'')$). At $u(\lambda)=0$ this
gives $-\lambda v(\lambda)u'(\lambda)=2v(\lambda)^2+\lambda^2W(v,v')(\lambda)$.

Combining (M), $(ts)'=a(q_1-s)$, and (Wtv), at every root $\lambda$ of $u$:
\[
\boxed{
-p(\lambda)u'(\lambda)=
\underbrace{a\bigl(-q_1(\lambda)u'(\lambda)\bigr)}_{T_1(\lambda)}
+\underbrace{a\bigl(s(\lambda)u'(\lambda)\bigr)}_{T_2(\lambda)}
+\underbrace{(m-1)a\bigl(2v(\lambda)^2+\lambda^2W(v,v')(\lambda)\bigr)}_{T_3(\lambda)}.
}
\tag{D}
\]

`verify.py` V2: 150/150 numeric root evaluations, 0 decomposition failures
(80-digit evaluation).

## 4. What the coupled theorem controls

The coupled theorem at $(m-1,a+1)$ gives $UQ_{m-2,a+1}\preceq Q_{m-1,a+1}$ with
orientation $W(Q_{m-1,a+1},u)>0$. At a root $\lambda$ of $u$,
$W(q_1,u)(\lambda)=-q_1(\lambda)u'(\lambda)$, hence
\[
T_1(\lambda)=a(-q_1(\lambda)u'(\lambda))>0. \tag{Cpl}
\]

The Laguerre term satisfies $W(v,v')=(v')^2-vv''\ge0$ for real-rooted $v$;
strictly $>0$ at simple-root points relevant here. Since $\lambda^2\ge0$ and
$v(\lambda)^2\ge0$,
\[
T_3(\lambda)\ge0. \tag{Lag}
\]

`verify.py` V3: 150/150 checks, `coupled_Q1_bad=0`, `laguerre_bad=0`.

## 5. The actual missing condition (exact)

By (D),(Cpl),(Lag), $W(p,u)>0$ is equivalent to
\[
\boxed{
T_1(\lambda)+T_3(\lambda) > -T_2(\lambda)
\quad\text{for every root }\lambda\text{ of }u=UQ_{m-2,a+1},
}
\tag{Sign}
\]
i.e.
\[
-q_1(\lambda)u'(\lambda) + (m-1)\bigl(2v(\lambda)^2+\lambda^2W(v,v')(\lambda)\bigr)
> -s(\lambda)u'(\lambda)
\quad (\text{after dividing by }a>0).
\]

This is **not** a sign condition on $Q_{m-1,a}(\lambda)=s(\lambda)$ alone.
When $-s(\lambda)u'(\lambda)>0$ (the "good" sign for $s$), $T_2(\lambda)<0$
and opposes the target; when $-s(\lambda)u'(\lambda)<0$, $T_2(\lambda)>0$ and
helps. On the actual triangle, `verify.py` V4 finds 122/150 root evaluations
with $-s(\lambda)u'(\lambda)>0$ and 28/150 with $-s(\lambda)u'(\lambda)<0$;
nevertheless the total $-p(\lambda)u'(\lambda)$ is positive in all checked
actual cases. Thus the sign of $s(\lambda)$ alone cannot close D2C.

## 6. Exact obstruction under proved hypotheses

The conjunction $\Sigma$ of proved facts from `proof-j2-distance-two/NOTES.md`
§5 does not imply (Sign), hence does not imply $W(p,u)>0$, hence does not imply
D2C or (I2).

*Witness (same rational model as `proof-j2-distance-two/NOTES.md` §5).*
$m=6,a=1$, $p=Q_{6,1}=1+\frac32 t+\frac3{10}t^2$, $s=Q_{5,1}=1+\frac35 t$,
$q_1=Q_{5,2}=2+\frac95 t$, $v=Q_{4,2}=2+\frac3{50}t$, $u=Uv=4+\frac9{25}t$.
The unique root of $u$ is $\lambda=-100/9=-22.222\ldots$. Exact evaluation gives
\[
T_1=\frac{171}{25}=6.84>0,\qquad
T_2=-\frac{111}{50}=-2.22,\qquad
T_3=\frac{40}{3}=13.333\ldots\ge0,
\]
\[
T_1+T_2+T_3=-\frac{3127}{150}=-20.846\ldots < 0.
\]
This matches $W(p,u)=\frac3{500}(9t^2+400t+970)$, whose positive discriminant
$\frac{28143}{6250}$ gives a sign change. `verify.py` V5 prints the same
decomposition.

The lower-level coupled controls hold in this model ($W(q_1,u)=171/25>0$,
$W(s,u)=111/50>0$, $W(q_1,v)=87/25>0$, $W(s,v)=57/50>0$), so the failure is
exactly the magnitude condition (Sign), not a failure of the coupled theorem.

## 7. Direction and monotonicity audit (policy items 1–2)

- Orientation: with $\deg p=\deg u+1$, $W(p,u)=p'u-pu'$ positive on $\mathbb R$
  means $p$ surrounds $u$; at a root of $u$ this is $-p(\lambda)u'(\lambda)>0$.
  Checked on the sample table in `proof-j2-distance-two/NOTES.md` §4: LC and
  Sturm root count give $W(p,u)>0$.
- (D) is an algebraic identity; no monotonicity claim is made. The Laguerre
  inequality $W(v,v')=(v')^2-vv''\ge0$ is the standard real-rootedness
  inequality; it is used only with direction $\ge0$.
- No endpoint monotonicity argument is used. The §5 countermodel is an explicit
  witness that the inequality can be reversed under $\Sigma$.

## 8. Boundary/degeneracy and finite evidence (policy items 4–5)

- $m\le3$: $Q$ constant/zero, D2C vacuous or automatic as in
  `proof-j2-distance-two/NOTES.md` §7.
- $m=4$: $Q_{4,a}$ linear, $Q_{2,a+1}$ constant, automatic.
- $a=1$: (U) loses the $(a-1)$ term but (M),(D),(Cpl),(Lag) remain valid.
- $\deg p=\deg u+1$ in nonautomatic cases; equal-degree obstruction does not
  occur for this pair.
- Symbolic residuals: (M) checked by `explore.py`/`verify.py`; (Wtv) checked by
  direct expansion above and used in V2; (D) checked by V2 (0 failures).
- Finite evidence (not a proof): `proof-j2-distance-two/deepcheck.py` checked
  D2C on the actual triangle through $m\le29$ with 0 failures. V6 here checks
  $W(Q_{m-1,a},Q_{m-2,a+1})>0$ on the actual triangle for $5\le m\le17,1\le a\le6$:
  78/78, 0 bad. This candidate invariant is not proved.

## 9. Conclusion

The exact missing lemma at a root $\lambda$ of $UQ_{m-2,a+1}$ is not
"$Q_{m-1,a}(\lambda)$ has a fixed sign"; the NOTES §6 formula that suggested that
form is false. The correct exact condition is the magnitude inequality (Sign).
It is implied by no combination of the proved facts in $\Sigma$ (exact rational
witness above). A rigorous all-index proof of D2C would therefore have to prove
(Sign) using the actual associated-Stirling recurrence/triangle in an essential
new way (e.g. a TP closure strong enough to control the $T_2$ term, or a
distance-two coupled induction that propagates a stronger invariant). That
proof was not found here.
