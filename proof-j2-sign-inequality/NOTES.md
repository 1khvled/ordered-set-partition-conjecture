# All-index proof of the corrected magnitude inequality `(Sign)`

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-j2-sign-inequality` was created/edited.  
No installs, builds, external agents, models, or browser use.

## 0. Outcome

**Theorem (all-index magnitude lemma).** Let $m\ge6$ and $a\ge1$, and put

\[
p=Q_{m,a},\qquad s=Q_{m-1,a},\qquad v=Q_{m-2,a+1},\qquad
u=\mathcal Uv,\qquad \mathcal U f=2f+tf'.
\]

At every real root $\lambda$ of $u$,

\[
\boxed{T_1(\lambda)+T_3(\lambda)>-T_2(\lambda),}
\tag{Sign}
\]

where

\[
T_1=a(-Q_{m-1,a+1}u'),\qquad
T_2=a(su'),\qquad
T_3=(m-1)a\bigl(2v^2+t^2W(v,v')\bigr).
\]

The proof uses the full associated-Stirling recurrence, not the retracted
identity $p=(ts)'$ at roots of $u$. Its new ingredient is a sum-of-squares
certificate for a refined residue residual. It proves the **second** D2C
premise only:

\[
W(Q_{m,a},\mathcal UQ_{m-2,a+1})>0\qquad\text{on }\mathbb R.
\]

**Retraction (policy item 6).** An earlier version of this note described
this as closing D2C/I2. That was unjustified: D2C also requires the first
premise $W(Q_{m,a},Q_{m-2,a+1})>0$, which is not proved below. No D2C or I2
closure is claimed. The first premise is investigated separately in §7.

## 1. Inputs and conventions

Write $n=m-2\ge4$. The established inputs used are:

1. the associated-Stirling polynomial recurrence
   \[
   Q_{m,a}=(tQ_{m-1,a})'+(m-1)a\,tQ_{m-2,a+1}; \tag{M}
   \]
2. the vertical identity
   \[
   Q_{m,a+1}=Q_{m,a}+\frac1a(tQ_{m,a})'; \tag{V}
   \]
3. real-rootedness of every $Q_{m,a}$, with simple negative roots in the
   nonconstant range;
4. the coupled theorem of `proof-u-operator/NOTES.md` §7, after applying the
   PF multiplier $\Lambda_{a+1}$:
   \[
   W(Q_{n,a+1},Q_{n-1,a+1})>0\qquad(n\ge5). \tag{C}
   \]

The case $n=4$ of (C) is automatic and is checked explicitly below.

For the proof set $x=-t>0$ and define

\[
R(x)=Q_{n,a}(-x),\quad V(x)=Q_{n,a+1}(-x),\quad
H(x)=Q_{n-1,a+1}(-x),
\]
\[
S(x)=Q_{n+1,a}(-x),\quad P(x)=Q_{n+2,a}(-x),\quad
U(x)=2V+xV'.
\]

Primes in the rest of the proof mean $d/dx$. Define the opposite-oriented
Wronskian

\[
B(F,G):=FG'-F'G.
\]

Thus, if $F=f(-x)$ and $G=g(-x)$, then

\[
B(F,G)=W_t(f,g). \tag{B-translation}
\]

## 2. Correct root calculus and a retraction

At a root $x>0$ of $U$,

\[
xV'=-2V. \tag{U-root}
\]

Also $V(x)\ne0$: at a simple root of $V$, $U=xV'\ne0$.

A tempting but false shortcut is $u'(\lambda)=-3v(\lambda)/\lambda$.
It is **not used**. The correct identity is

\[
-tv\,u'_t=2v^2+t^2W_t(v,v'). \tag{L-root}
\]

In the $x$-variable, at a root of $U$ this becomes

\[
-xVU'=2V^2+x^2(V'^2-VV'').
\]

Since $V$ has only real roots, Laguerre gives $V'^2-VV''\ge0$. Hence

\[
-xVU'>0,\qquad\text{so } VU'<0. \tag{2.1}
\]

The original Wronskian value at $\lambda=-x$ is

\[
W_t(p,u)(\lambda)=-p(\lambda)u'_t(\lambda)=P(x)U'(x).
\]

By (2.1), it is therefore enough to prove

\[
P(x)V(x)<0. \tag{2.2}
\]

**Retraction audit.** In the explicit sample $(m,a)=(6,1)$ used in §6,
$\lambda=-2/27$, $v(\lambda)=2/3$, and $u'(\lambda)=54$, whereas the false
shortcut would give $27$. The error was differentiating the root equation
$tv'/v=-2$ as though it held on an interval. No step below uses that error.

The earlier false formula $p(\lambda)=(ts)'(\lambda)$ is likewise not used:
the full extra term $(m-1)a\lambda v(\lambda)$ in (M) is retained.

## 3. Refined residue identity

Translating (M) for $s=Q_{n+1,a}$ gives

\[
S=R+xR'-naxH.
\]

Using (V), $a(V-R)=R+xR'$, so

\[
\boxed{S=a(V-R-nxH).} \tag{3.1}
\]

Similarly, translating (M) for $p=Q_{n+2,a}$ gives

\[
P=S+xS'-(n+1)axV. \tag{3.2}
\]

At a root of $U$, (U-root) yields the exact identity

\[
\begin{aligned}
PV
&=-xB(S,V)-SV-(n+1)axV^2. \tag{3.3}
\end{aligned}
\]

Now substitute (3.1). Since

\[
B(xH,V)=xB(H,V)-HV=-xB(V,H)-HV,
\]

one obtains

\[
B(S,V)=a\bigl(B(V,R)+nxB(V,H)+nHV\bigr).
\]

The $nHV$ terms cancel against $SV/a=V^2-RV-nxHV$, leaving

\[
\boxed{
xB(S,V)+SV+(n+1)axV^2
=a\Bigl(G+n x^2B(V,H)\Bigr),
} \tag{3.4}
\]

where

\[
\boxed{
G=xB(V,R)+(n+1)xV^2+V^2-RV.
} \tag{3.5}
\]

Combining (3.3) and (3.4),

\[
\boxed{
PV=-a\Bigl(G+n x^2B(V,H)\Bigr).
} \tag{3.6}
\]

This is the refined residue identity. It keeps the full coupled Wronskian
term instead of trying to dominate $T_2$ by a ratio bound alone.

## 4. The sum-of-squares certificate

Put

\[
M=1+(n+1)x(a+1)>0.
\]

The vertical identity in the $x$-variable is

\[
aV=(a+1)R+xR'. \tag{4.1}
\]

Direct expansion of (3.5), using (4.1), gives the polynomial identity

\[
\boxed{
\begin{aligned}
a^2G={}&ax^2(R'^2-RR'')\\
&+(a+1)M\left(R+\frac{xR'}{a+1}\right)^2
+\frac{ax^2}{a+1}(R')^2.
\end{aligned}
} \tag{SOS}
\]

The residual of (SOS) is identically zero; `verify.py` checks it for
completely generic $R,R',R''$ and symbolic $a,n,x$.

Every term on the right is nonnegative:

- $R'^2-RR''\ge0$ by Laguerre because $R$ is real-rooted;
- $M>0$ and the square is nonnegative;
- $a,x>0$.

In fact $G>0$ strictly. Otherwise every term would vanish. The last term
would force $R'=0$. The first would then force $RR''=0$. Since $R$ has only
simple roots, $R'=0$ implies $R\ne0$, so $R''=0$. The middle term would
then be $(a+1)MR^2>0$, a contradiction.

It remains to note that $B(V,H)>0$. For $n\ge5$ this is exactly (C), by
(B-translation). For $n=4$,

\[
H=a+1,\qquad V=(a+1)-3(a+1)(a+2)x,
\]

so

\[
B(V,H)=-V'H=3(a+1)^2(a+2)>0.
\]

Therefore (3.6) gives $PV<0$, proving (2.2). By (2.1), $PU'>0$, which is
precisely $W_t(p,u)(\lambda)>0$.

Finally, the corrected decomposition from `proof-j2-qsign/NOTES.md` is

\[
W_t(p,u)(\lambda)=T_1(\lambda)+T_2(\lambda)+T_3(\lambda).
\]

Hence

\[
T_1+T_3>-T_2,
\]

which proves `(Sign)`. $\blacksquare$

Because $u$ has simple negative roots and $\deg p=\deg u+1$, the same root
sign count gives $W_t(p,u)>0$ on all of $\mathbb R$. This proves only the
second D2C premise; the first premise $W_t(p,v)>0$ remains a separate
question.

## 5. Why the artificial countermodel is not feasible here

The rational model in `proof-j2-distance-two/NOTES.md` §5 deliberately
omits (M). It therefore does not satisfy the translated recurrence (3.1),
the refined identity (3.4), or the sum-of-squares identity (SOS). It remains
a valid obstruction to proofs using only the older structural hypotheses,
but it is not a counterexample to the theorem above.

In particular, the proof uses (M) essentially: (3.1) and (3.2) are exactly
where the associated-Stirling coefficients enter the residue calculation.

## 6. Proof-policy audit

### 6.1 Direction audit

Use the actual feasible sample $(m,a)=(6,1)$:

\[
p=1+50t+90t^2,\quad s=1+20t,\quad Q_{5,2}=2+60t,
\]
\[
v=2+18t,\quad u=4+54t,\quad \lambda=-\frac2{27}.
\]

Exact evaluation gives

\[
T_1=132>0,\qquad T_2=-26,\qquad T_3=\frac{40}{3}\ge0,
\]

and

\[
T_1+T_2+T_3=\frac{358}{3}>0.
\]

Thus the required direction is indeed $T_1+T_3>-T_2$, not the reverse.
At the same point $G=166/243>0$ and $B(V,H)=36>0$, so (3.6) has the
claimed strict direction.

### 6.2 Monotonicity audit

No monotonicity claim is used. The only differential inequalities are
Laguerre's $R'^2-RR''\ge0$ and the proved coupled Wronskian $B(V,H)>0$.
There is no endpoint optimization or reversed finite difference.

### 6.3 Feasible-constraint falsification

The certificate (SOS) is a sum of nonnegative terms for every real-rooted
$R$ and every $a,x>0$. A zero would force $R=R'=R''=0$, contradicting
simplicity of the roots of $R$. Thus the proved feasible constraints do
not permit a counterexample to $G>0$. The older artificial model is not
feasible for (M), as recorded in §5.

### 6.4 Boundary and degeneracy cases

- $m\le4$: the D2C pair is constant/linear and automatic or vacuous.
- $m=5$: $v=Q_{3,a+1}$ is constant, so $u$ is a nonzero constant and has
  no roots; `(Sign)` is vacuous.
- $m=6$, equivalently $n=4$: handled explicitly in §4.
- $a=1$: (V) loses its $(a-1)$ term, but every displayed identity remains
  valid; the sample above is this case.
- Inactive coefficients do not create zero active leading coefficients;
  the degree and simple-root statements are those of the established
  $Q$-family.
- $R$ and $V$ have the same degree; $H$ has one lower degree. The target
  pair $p,u$ has $\deg p=\deg u+1$, so the Wronskian orientation is
  unambiguous.
- At a root of $U$, $V\ne0$, so no division by $V$ is hidden in the proof.

### 6.5 Symbolic residuals

`verify.py` prints zero residuals for:

- (M) and (V) on an exact finite test range;
- the translations (3.1), (3.2), and (3.4);
- the generic sum-of-squares identity (SOS);
- the exact sample decomposition in §6.1.

Finite Sturm checks, when printed, are labelled finite evidence only. The
all-index proof is the algebraic argument in §§2–4.

## 7. Status of the first D2C premise

The first premise is

\[
\boxed{W_t(p,v)>0\qquad\text{on }\mathbb R,} \tag{D2C-1}
\]

where $v=Q_{n,a+1}$ and $p=Q_{n+2,a}$. **No all-index proof of (D2C-1)
was obtained.** Accordingly D2C and I2 remain open in this note.

### 7.1 Exact root reduction

At a root $x>0$ of $V$, the target is

\[
B(P,V)=PV'>0.
\]

Put

\[
Y:=2H+xH'=\mathcal U Q_{n-1,a+1}(-x).
\]

The full recurrence (M), not the retracted root-of-$u$ formula, gives the
exact identity modulo $V$

\[
\boxed{\frac Pa=aR+xV'-nxY\qquad\text{at }V=0.} \tag{7.1}
\]

For completeness, (3.1) gives $S=a(V-R-nxH)$, and (3.2) gives
$P=S+xS'-(n+1)axV$. Substituting these and cancelling modulo $V$ yields
(7.1); `verify.py` V9 checks both the generic residual and the actual
triangle reduction.

The established interlacings give

\[
B(V,R)>0,\qquad B(V,H)>0,\qquad B(R,H)>0,\qquad B(V,Y)>0. \tag{7.2}
\]

Here the first is vertical interlacing, the second and fourth are the
coupled theorem, and the third is the adjacent diagonal step (I1). At a
root of $V$, these say that $R,H,Y$ all have sign opposite to $V'$.
Multiplying (7.1) by $V'$ therefore gives the exact magnitude condition

\[
\boxed{
x(V')^2+n x B(V,Y)>a B(V,R).
} \tag{7.3}
\]

Indeed,

\[
\frac{PV'}a
=x(V')^2-aB(V,R)+nxB(V,Y).
\]

Thus the coupled theorem controls every sign in (7.3), but it does not
supply the required magnitude domination between the three terms.

### 7.2 Why the second-premise SOS certificate does not close this

The SOS proof in §4 used the root-of-$U$ relation $xV'=-2V$. At a root of
$V$ itself that relation is unavailable, and the residual in (7.3) contains
the genuinely competing negative term $-aB(V,R)$. The sum-of-squares
identity (SOS) is therefore not a certificate for (7.3).

No alternative all-index SOS, three-step induction, or residue domination
proving (7.3) was found. The exact Sturm checks in `verify.py` V11 are
**finite evidence only**, not a proof. No counterexample on the actual
associated-Stirling triangle was found in those checks.

### 7.3 Direction and boundary audit for (D2C-1)

For $(m,a)=(6,1)$,

\[
v=2+18t,\qquad p=1+50t+90t^2,\qquad \lambda=-\frac19,
\]

so

\[
W_t(p,v)(\lambda)=-p(\lambda)v'(\lambda)
=-\left(-\frac{31}{9}\right)18=62>0.
\]

This confirms the direction $>0$ at an explicit feasible point. For $m=5$,
$v=Q_{3,a+1}$ is constant and (D2C-1) is automatic. The nonautomatic
all-index range $m\ge6$ is exactly the range left open by (7.3). No
monotonicity claim is used in this reduction.

## 8. Files and reproduction

- `NOTES.md` — proof of `(Sign)` and the open status of (D2C-1).
- `verify.py` — exact symbolic residuals, direction samples, and explicitly
  labelled finite evidence.

From this directory:

```powershell
py -3.14 verify.py
```

No packages are installed and no build is performed.
