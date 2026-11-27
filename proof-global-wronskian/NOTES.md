# Global negative-axis Wronskian positivity: complex-analytic / Pick / resolvent attack

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-global-wronskian` was created or edited.  
Target: inequality (15) of `proof-final-wrs-corrected/NOTES.md`,
\[
\mathcal W_{n,j}(t)=W(O_{n,j},O_{n-j,j})(t)>0\quad
\text{for every }j\ge2,\ n\ge j,\ t\le0 .
\tag{15}
\]

---

## 0. Outcome

**No unconditional all-index proof of (15) was obtained.** The known open
conjecture is confirmed open by the live problem list
(<https://symmetricfunctions.com/realRootedOpenProblems.htm>, entry
"Conjecture (Alexandersson, 2026)" for ordered set partitions; status: open).

What *is* proved here is a set of exact structural facts that sharpen the gap
and, in particular, show that the representations requested by this track
cannot be *shortcuts*: they are equivalent reformulations of the same missing
real-rootedness/Pick statement, not weaker conditions.

1. **Proposition A (Pick equivalence).** Under the standing simple-real-root
   induction hypothesis, positivity of $W(O_n,O_{n-j})$ on $t\le0$ is
   equivalent to: (i) positivity on **all** of $\mathbb R$; (ii) the ratio
   $R_{n,j}=O_n/O_{n-j}$ being a Pick/Nevanlinna function; (iii) the Bezoutian
   matrix $M(O_{n-j},O_n)$ being positive definite; (iv) the WRS residue
   inequality. So a "Stieltjes/Pick/Bezoutian proof" is *the same* problem in
   another language (§2, §3).
2. **Proposition B (no free lunch).** If (15) holds for all indices, then
   **every** $O_{n,j}$ is real-rooted (simple negative roots) and
   $O_{n-j,j}\interl O_{n,j}$. Hence any proof of (15) must itself contain a
   real-rootedness argument for the actual rows (§4).
3. **Proposition C (exact cross-family reduction).** (15) follows from a
   single *common-interlacer* condition on the family
   $\{H_{n,j}\}\cup\{\binom nr O_{n-r,j}:1\le r\le j-1\}$ relative to
   $O_{n-j,j}$. This is a clean sufficient condition, but it is conditional
   on (H) and on lag interlacing, both open (§5).
4. **Proposition D (double-resolvent kernel).** There is an exact Cauchy
   representation
   \[
   W(O_n,O_{n-j})(t)=\frac{n!(n-j)!}{(2\pi i)^2}\oint\oint
   \frac{\Phi(z,w)}{D(z)^2D(w)^2}\,\frac{dz\,dw}{z^{n+1}w^{n-j+1}},
   \quad \Phi=C(z)B(w)-C(w)B(z),
   \]
   with $B=1-A$, $D=1-A-tC$. The kernel $\Phi=(z-w)\Psi(z,w)$ is **not**
   positive: for $j=2$, $\Psi(z,z)=(2-z)e^z-2$ is positive at $z=1$ and
   negative at $z=3$. So no termwise-positive contour-kernel expansion of
   this type exists (§6).
5. **Proposition E (why the "resolvent is Pick" idea stalls).** For each fixed
   real $z>0$, $F_j(z,t)=1/(B(z)-tC(z))$ **is** a Pick function of $t$, being
   a resolvent at the real number $B(z)/C(z)$ with weight $1/C(z)>0$.
   But (a) coefficient extraction in $z$ uses a complex contour measure and is
   not a positive operation, (b) $B/C$ is not monotone on $(0,\infty)$
   (exact derivative computed in §7), so the two-variable kernel $F_j$ is
   not totally positive, and (c) $O_{n,j}$ itself is not Pick once
   $\lfloor n/j\rfloor\ge2$, so no positive combination of these resolvents
   can equal it (§7).
6. **Finite evidence for a stronger statement.** Exact Sturm counts verify
   that $W(O_n,O_{n-d})$ has **no real zero** for all $j=2,\dots,8$,
   $1\le d\le j$, $n\le\min(6j,40)$ (972 pairs), and likewise
   $W(H_{n,j},O_{n-j,j})$ has no real zero (172 pairs). This supports the
   stronger **all-lag Pick Conjecture P** (§8), which is also open and which
   implies the full real-rootedness/interlacing conjecture.
7. **Precise remaining gap** is stated in §9: it is exactly the WRS margin
   (equivalently Bezoutian domination), now identified as Pick positivity of
   $R_{n,j}$. Nothing in the requested list of techniques supplies it without
   already proving real-rootedness.

---

## 1. Setting and standing conventions

Fix $j\ge2$ and write (as in `proof-final-wrs-corrected`, §1)
\[
A(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\quad
E(z)=\sum_{r=0}^{j-1}\frac{z^r}{r!},\quad
C(z)=e^z-E(z),\quad
B(z)=1-A(z)=2-E(z),
\]
\[
D(z,t)=1-A(z)-tC(z)=B(z)-tC(z),\quad
F(z,t)=\frac1{D(z,t)},\quad
O_n=n!\,[z^n]F,\quad O_0=1 .
\]
Set $q=O_{n-j}$, $p=O_n$, $W(p,q)=p'q-pq'$ (derivative in $t$), and
\[
H_{n,j}=n![z^n]\frac{B(z)}{D(z,t)}
      =O_n-\sum_{r=1}^{j-1}\binom nr O_{n-r}
      =t\sum_{r=j}^n\binom nr O_{n-r}.
\tag{1}
\]
The identities used below are proved in `proof-final-inner/NOTES.md`,
Lemma 1, and `proof-final-wrs-corrected/NOTES.md`, §4. Degree facts:
$\deg_t O_n=\lfloor n/j\rfloor$ and $\deg_t H_{n,j}=\lfloor n/j\rfloor$,
with positive leading coefficients; hence $\deg p=\deg q+1$, and
$a_0=\operatorname{lc}(p)/\operatorname{lc}(q)>0$.

Throughout §2–§5 we impose the **standing induction hypothesis**: $q$ has
$e=\lfloor n/j\rfloor-1$ distinct negative roots
$\mu_1<\dots<\mu_e<0$, and $p$ has positive leading coefficient. The root
$\mu=-\infty$ case $e=0$ ($n=j$) is trivial: $O_j$ is linear with a negative
root.

---

## 2. Proposition A — (15) on $t\le0$ is Pick positivity

Let $p,q$ satisfy the standing hypothesis and write the partial fraction
decomposition of $R=p/q$:
\[
R(t)=at+b+\sum_{i=1}^{e}\frac{r_i}{t-\mu_i},\qquad
a=\frac{\operatorname{lc}p}{\operatorname{lc}q}>0,\qquad
r_i=\frac{p(\mu_i)}{q'(\mu_i)} .
\tag{2}
\]
Then, identically as polynomials,
\[
W(p,q)=q^2R'(t)=q^2\Bigl(a-\sum_{i=1}^{e}\frac{r_i}{(t-\mu_i)^2}\Bigr),
\tag{3}
\]
and evaluating the right side at $t=\mu_i$ (the singularities cancel) gives
\[
\boxed{\,W(p,q)(\mu_i)=-r_i\,q'(\mu_i)^2\,}.
\tag{4}
\]

**Proposition A.** With $p,q$ as above, the following are equivalent:

* (A1) $W(p,q)(t)>0$ for every $t\in\mathbb R$;
* (A2) $W(p,q)(t)>0$ for every $t\le0$;
* (A3) $r_i<0$ for every $i$;
* (A4) $R=p/q$ is a Pick (Nevanlinna) function: $\Im R(z)>0$ for $\Im z>0$,
  i.e. $R(z)=\alpha+\beta z+\sum_{i} w_i/(\mu_i-z)$ with $\beta\ge0$,
  $w_i>0$, $\mu_i\in\mathbb R$;
* (A5) the WRS inequality
  $-p(\mu_i)/q'(\mu_i)>0$ holds at every root $\mu_i$ of $q$.

*Proof.*
(A1)$\Rightarrow$(A2) is trivial. (A2)$\Rightarrow$(A3): each $\mu_i<0$ lies
in the tested half-line, so (4) and (A2) give $-r_iq'(\mu_i)^2>0$, hence
$r_i<0$ (roots are simple, $q'(\mu_i)\ne0$).
(A3)$\Rightarrow$(A1): by (3), $a-\sum_i r_i/(t-\mu_i)^2\ge a>0$ for every
real $t\notin\{\mu_i\}$, and $W$ is a polynomial, so $W>0$ everywhere.
(A3)$\Leftrightarrow$(A4): if $r_i<0$, put $w_i=-r_i>0$ and rewrite (2) as
$R(t)=at+b+\sum_i w_i/(\mu_i-t)$; for $z\in\mathbb C^+$,
$\Im(\mu_i-z)^{-1}=\Im z/|\mu_i-z|^2>0$, so $\Im R(z)=a\Im z+\sum_i
w_i\Im z/|\mu_i-z|^2>0$. Conversely, a rational Pick function has the
Nevanlinna representation $R(z)=\alpha+\beta z+\int(\frac1{s-z}-\frac{s}{1+s^2})d\rho(s)$
with $\beta\ge0$ and $\rho$ a positive measure; for a rational function
$\rho=\sum w_j\delta_{s_j}$, $w_j>0$, so near a pole $s_j$,
$R(z)\sim w_j/(s_j-z)$ has residue $-w_j<0$, i.e. $r_j<0$.
(A3)$\Leftrightarrow$(A5) is the definition of $r_i$. $\qquad\square$

**Consequence.** The requested "global negative-axis Wronskian positivity"
is *exactly* the statement that the ratio $O_n/O_{n-j}$ is a Pick function.
It is not a weakening of WRS that avoids real-rootedness: (A2) forces
$r_i<0$ at the already-existing negative poles, and (A3) then upgrades it
to positivity on all of $\mathbb R$.

---

## 3. Proposition A′ — Pick = Bezoutian positivity

For the same $p,q$ define the Bezoutian kernel
$\mathcal B(t,s)=\frac{p(t)q(s)-p(s)q(t)}{t-s}$ and let $M(q,p)$ be its
coefficient matrix in the monomial basis (`proof-bezout/NOTES.md`, §2,
uses the same convention). Then for simple real roots of $q$,
\[
\mathcal B(\mu_i,\mu_i)=-q'(\mu_i)p(\mu_i)=q'(\mu_i)^2\,w_i,
\qquad w_i=-r_i ,
\tag{5}
\]
so $M(q,p)\succ0$ iff all $w_i>0$ (the extra leading direction contributes
one more positive square). Thus (A3)–(A5) are the same as the Bezoutian
domination target (16) of `proof-bezout/NOTES.md`:
\[
M(q,A_n)+c\,v(q)v(q)^T\;\succ\;c\,M(q,h_q),\qquad h_q=t(t-1)q'.
\tag{6}
\]
**The complex-analytic/Pick route and the Bezoutian route are one and the
same inequality.** No new positivity is gained by changing languages.

---

## 4. Proposition B — (15) implies real-rootedness (so cannot be free)

**Proposition B.** Assume (15) for every $j\ge2$, $n\ge j$. Then for each
$j\ge2$ and every $n\ge0$, $O_{n,j}$ is real-rooted with simple negative
roots, and for $n\ge j$, $O_{n-j,j}\interl O_{n,j}$.

*Proof.* Fix $j$ and induct on $n$. For $n<j$ the $t$-sum in the first-block
recurrence is empty, so $O_{n,j}$ is a positive constant (no roots). For
$n=j$, $O_j=\text{const}+t$ has one negative root. Let $n>j$ and assume
$O_{n-j}$ has $e=\lfloor n/j\rfloor-1$ simple negative roots. By
Proposition A, (15) gives (A3), hence (A1): $W(O_n,q)>0$ on $\mathbb R$.
By (4), $\operatorname{sgn}O_n(\mu_i)=-\operatorname{sgn}q'(\mu_i)$, so the
signs of $O_n$ at the ordered roots $\mu_1<\dots<\mu_e$ alternate. Since
$\deg O_n=e+1$ with positive leading coefficient, $O_n(-\infty)$ has sign
$(-1)^{e+1}$, while $\operatorname{sgn}O_n(\mu_1)=-\operatorname{sgn}q'(\mu_1)
=(-1)^e$ (since $q'(\mu_1)$ has sign $(-1)^{e-1}$); these are opposite, and
$O_n(0)>0$ has sign opposite to $O_n(\mu_e)=-q'(\mu_e)<0$. The intermediate
value theorem gives one root in each of
\[
(-\infty,\mu_1),(\mu_1,\mu_2),\dots,(\mu_{e-1},\mu_e),(\mu_e,0),
\]
i.e. $e+1=\deg O_n$ distinct simple negative roots, interlacing $q$.
This completes the induction. $\qquad\square$

(This is the Wronskian-language form of Theorem 2 of
`proof-weighted/NOTES.md`; the point here is the *logical* consequence:
proving (15) already proves real-rootedness of every actual row. Hence any
"resolvent representation" proof of (15) that merely assumes the rows are
real-rooted is circular.)

---

## 5. Proposition C — exact cross-family reduction (conditional)

Let $\mathcal F_{n,j}=\{\,H_{n,j}\,\}\cup\{\binom nr O_{n-r,j}:1\le r\le j-1\}$
(the term $\binom nj\cdot O_{n-j}=\binom nj$ is a constant and contributes
nothing to a Wronskian). Say that a real-rooted polynomial $f$ of degree
$\deg q$ or $\deg q+1$, with positive leading coefficient, is **aligned with
$q$** if its roots $\nu_1<\dots<\nu_d$ satisfy the "newer row" interlacing
orientation
\[
\nu_1<\mu_1<\nu_2<\mu_2<\dots<\mu_e<\nu_{e+1}
\quad(d=e+1),\qquad
\mu_1<\nu_1<\mu_2<\nu_2<\dots<\mu_e<\nu_e
\quad(d=e).
\tag{7}
\]

**Lemma C1.** If $f$ is aligned with $q$, then $r_f=f(\mu_i)/q'(\mu_i)<0$
for all $i$, hence $W(f,q)>0$ on $\mathbb R$ by Proposition A.

*Proof.* For $d=e+1$, $f(\mu_i)$ has sign $(-1)^{e+1-i}$ and $q'(\mu_i)$
has sign $(-1)^{e-i}$. For $d=e$, $f(\mu_i)$ has sign $(-1)^{e-i+1}$ and
$q'(\mu_i)$ has sign $(-1)^{e-i}$. In both cases the ratio is negative.
$\qquad\square$

**Proposition C.** If every member of $\mathcal F_{n,j}$ is aligned with
$q=O_{n-j,j}$, then $\mathcal W_{n,j}=W(O_n,q)>0$ on $\mathbb R$, hence (15).

*Proof.* By bilinearity of $(\cdot,q)\mapsto W(\cdot,q)$ and identity (1),
\[
\mathcal W_{n,j}=\sum_{r=1}^{j-1}\binom nr W(O_{n-r},q)+W(H_{n,j},q),
\tag{8}
\]
which is a **positive** linear combination of the Wronskians in Lemma C1.
$\qquad\square$

This is the exact cross-family content requested by
`proof-final-wrs-corrected/NOTES.md`, §8: the whole older tail is the single
cross-family Wronskian $W(H_{n,j},q)$. Proposition C isolates precisely what
is needed — that $H_{n,j}$ and all lag members $O_{n-r}$, $1\le r\le j-1$,
share one aligned interlacer $q$. The hypothesis is **not** known: it
includes (H) (real-rootedness of $H_{n,j}$) and lag interlacing of the $O$
rows, both open. In particular for $j=2$ it reduces to the open boundary
statement "$H_{n,2}$ interlaces $O_{n-2,2}$" (cf.
`proof-final-common-interlacer/NOTES.md`, §6.2).

---

## 6. Proposition D — the double-resolvent representation and its fatal kernel

For $p=O_n$, $q=O_{n-j}$ we have $\partial_tF=C F^2$ (`proof-weighted`,
eq. (1)), so
\[
W(p,q)=n!(n-j)![z^nw^{n-j}]
\bigl(C(z)F(z)^2F(w)-F(z)C(w)F(w)^2\bigr).
\]
Since $F=1/D$, factoring gives the **exact double-resolvent identity**
\[
\boxed{
W(O_n,O_{n-j})(t)=\frac{n!(n-j)!}{(2\pi i)^2}\oint\oint
\frac{\Phi(z,w)}{D(z)^2D(w)^2}\,\frac{dz\,dw}{z^{n+1}w^{n-j+1}},
}
\tag{9}
\]
where (contours are small circles around $z=w=0$)
\[
\Phi(z,w)=C(z)D(w)-C(w)D(z)
         =C(z)-C(w)+C(w)A(z)-C(z)A(w)
         =C(z)B(w)-C(w)B(z).
\tag{10}
\]
The last two equalities were verified symbolically for $j=2,\dots,5$
(`verify_structures.py`). The $t$-dependence cancels from $\Phi$, which is
an entire function of $(z,w)$; $\Phi(z,z)=0$, hence $\Phi=(z-w)\Psi(z,w)$
with $\Psi$ symmetric.

**Proposition D.** The kernel $\Psi$ is not sign-definite on $(0,\infty)^2$.

*Proof.* For $j=2$, $A(z)=z$, $C(z)=e^z-1-z$ and a direct expansion gives
\[
\Phi(z,w)=(1-w)e^z-(1-z)e^w+2(w-z),
\qquad
\Psi(z,z)=\partial_z\Phi(z,z)=(2-z)e^z-2 .
\]
At $z=1$, $\Psi(1,1)=e-2>0$; at $z=3$, $\Psi(3,3)=-e^3-2<0$. By symmetry
$\Psi(z,w)$ takes both signs. $\qquad\square$

So (9) cannot be turned into a positive-kernel/positive-measure proof by any
choice of real positive contours: the integrand's numerator changes sign
already on the diagonal. (The factor $D(z)^2D(w)^2$ does not repair this;
on the relevant contours it is complex.) This rules out the most direct
"positive kernel expansion" interpretation of the requested resolvent route.

---

## 7. Proposition E — pointwise Pick of the EGF, and why it does not transfer

**(a) Pointwise Pick.** For any fixed real $z>0$ with $B(z)\ne0$,
$B(z),C(z)$ are real and $C(z)=\sum_{r\ge j}z^r/r!>0$. Hence
\[
F(z,t)=\frac{1/C(z)}{B(z)/C(z)-t}
      =\frac{w(z)}{\alpha(z)-t},\qquad
w(z)=\frac1{C(z)}>0,\quad \alpha(z)=\frac{B(z)}{C(z)}\in\mathbb R,
\tag{11}
\]
and for $\Im t>0$, $\Im F(z,t)=w(z)\Im t/|\alpha(z)-t|^2>0$: $F(z,\cdot)$
is a Pick function (a rank-one resolvent). $\alpha(z)$ has simple poles at
the zeros of $C(z)$ but is real wherever defined.

**(b) Coefficient extraction is not positive.** Formally
$O_n(t)=n![z^n]F(z,t)=\frac{n!}{2\pi i}\oint F(z,t)z^{-n-1}dz$; the weight
$z^{-n-1}dz$ is a complex measure on a circle around $0$, not a positive
measure on $(0,\infty)$. Pick functions are closed under **positive**
combinations only; this is not one.

**(c) No total positivity of the two-variable kernel.** For
$z_1,z_2>0$ and real $t_1,t_2$, the $2\times2$ minor of the kernel (11) is
\[
\det\begin{pmatrix}F(z_1,t_1)&F(z_1,t_2)\\F(z_2,t_1)&F(z_2,t_2)\end{pmatrix}
=\frac{w(z_1)w(z_2)\,(\alpha_1-\alpha_2)(t_2-t_1)}
 {(\alpha_1-t_1)(\alpha_1-t_2)(\alpha_2-t_1)(\alpha_2-t_2)},
\qquad \alpha_i=\alpha(z_i),
\]
as an exact expansion (verified from (11)). Whenever the four denominator
factors have a fixed sign, the minor's sign is that of
$(\alpha_1-\alpha_2)(t_2-t_1)$, so a fixed kernel sign needs $\alpha$
monotone. It is **not**:
\[
\alpha'(z)=\frac{e^z\bigl(p(z)-2\bigr)+2E_{j-2}(z)}{C(z)^2},
\qquad p(z)=\frac{z^{j-1}}{(j-1)!},\quad E_{j-2}=\sum_{r=0}^{j-2}\frac{z^r}{r!},
\tag{12}
\]
verified symbolically for $j=2,\dots,6$ (`verify_structures.py`). For $j=2$
the numerator $e^z(z-2)+2$ equals $2-e<0$ at $z=1$ and $2>0$ at $z=2$. So
$\alpha$ decreases then increases on $(0,\infty)$. For $j=2$, with
$t_1=-10<t_2=-1$, the minor is positive for $(z_1,z_2)=(1,2)$ (where
$\alpha_1=0>\alpha_2$) and negative for $(z_1,z_2)=(2,3)$ (where
$\alpha_1<\alpha_2<0$), all four denominator factors being positive in both
cases: the kernel $F_j(z,t)$ is not totally positive on $(0,\infty)^2$.
(Also $D(z,-1/10)=1-z+\frac1{10}(e^z-1-z)$ changes sign between $z=2$ and
$z=4$, so even the positivity region of $D$ is not stable under $t\le0$.)

**(d) $O_n$ itself is not Pick.** A rational Pick function grows at most
linearly at infinity, while $O_{n,j}(t)$ has degree $\lfloor n/j\rfloor\ge2$
for $n\ge2j$. Hence $O_n$ is not a positive combination of the rank-one
resolvents (11). The Pick property can only live at the level of the
**ratio** $O_n/O_{n-j}$ (degree difference one) — which is Proposition A,
i.e. the open WRS.

---

## 8. Finite evidence: the stronger all-lag Pick conjecture

`probe_pick.py` generates the actual rows over $\mathbb Q[t]$ and uses exact
Sturm root counting. Result: for all $j=2,\dots,8$, all
$1\le d\le j$ and all $j\le n\le\min(6j,40)$, the polynomial
$W(O_{n,j},O_{n-d,j})$ has **zero real roots** (972 pairs, 0 failures).
Under the standing simple-root hypothesis, Proposition A says each ratio
$O_n/O_{n-d}$ is Pick. `probe_cross.py` checks the cross-family Wronskian
$W(H_{n,j},O_{n-j,j})$ for $j=2,\dots,8$, $n\le\min(6j,40)$: 172 pairs, 0
real roots, 0 failures.

This suggests the stronger statement

> **Conjecture P (all-lag Pick).** For every $j\ge2$, $n\ge j$ and
> $1\le d\le j$, the rational function $O_{n,j}(t)/O_{n-d,j}(t)$ is Pick;
> equivalently $W(O_{n,j},O_{n-d,j})(t)>0$ for all real $t$.
> (The restriction $d\le j$ is necessary: only then is
> $\deg_tO_n-\deg_tO_{n-d}=1$, as required for a Pick function.)

Conjecture P immediately implies the full real-rootedness and consecutive
interlacing conjecture: the same IVT induction as in Proposition B, run with
$d=1$ (compare $O_n$ directly to $O_{n-1}$ at every root of $O_{n-1}$), gives
simple negative roots of $O_n$ interlacing those of $O_{n-1}$. It is at least
as hard as the target, so it cannot serve as an induction hypothesis without
an independent proof. These counts are evidence only; no float enters them.

---

## 9. Precise remaining gap

Let $q=O_{n-j,j}$, $c=\binom{n-1}{j-1}$, $A_n=nO_{n-1}+2TO_{n-1}
+\sum_{r=2}^{j-1}\binom nr O_{n-r}$. The exact all-index statement needed is

\[
\boxed{\;\mathcal W_{n,j}(t)>0
\quad\text{for every }j\ge2,\ n\ge j,\ t\le0.\;}
\tag{15}
\]

By §2–§4 the following are equivalent for the actual family, under the
simple-root induction:

1. (15);
2. $W(O_n,q)(t)>0$ for **all** real $t$;
3. $R_{n,j}=O_n/q$ is a Pick/Nevanlinna function;
4. $M(q,O_n)\succ0$ (Bezoutian positive definite), i.e. (6);
5. at each root $\mu$ of $q$,
   \[
   -\frac{A_n(\mu)}{O_{n-j}'(\mu)}
   >\binom{n-1}{j-1}\mu(\mu-1)
   \quad\Longleftrightarrow\quad \frac{O_n(\mu)}{q'(\mu)}<0,
   \tag{13}
   \]
   the WRS margin of `proof-final-wrs-corrected`, §8.

The requested techniques do not close it, for the following exact reasons:

* **Pick/Stieltjes/resolvent = Bezoutian = WRS** (Props. A, A′): they are
  reformulations of the gap, not hypotheses weaker than it.
* **Pointwise resolvent Pickness of $F$ (11) does not transfer** to the
  coefficients: extraction is complex, $B/C$ is not monotone (12), $F$ is
  not TP, and $O_n$ is not Pick (Prop. E).
* **The exact Cauchy kernel (9) is indefinite** (Prop. D), so there is no
  termwise-positive kernel expansion of that form.
* **Common-interlacer suffices but is open** (Prop. C): it requires (H) and
  aligned lag interlacing, which are themselves the conjectured conclusion.
* The local/coefficient/PF/TN/edge obstructions recorded in
  `proof-final-wrs-corrected`, §§5–7 and `proof-final-renewal`, §4 remain in
  force; §6 adds the exact failure of the natural two-variable kernel.

Consequently an unconditional proof of (15) would already prove
real-rootedness of all actual rows (Prop. B) and is a full research-level
result. No such proof was obtained in this track; the gap is precisely (13)
/ (6) / item 3 above.

---

## 10. Files and reproduction

* `NOTES.md` — this report.
* `verify_structures.py` — exact symbolic checks of (10), (12), the $j=2$
  sign change of $\Psi$, and the explicit Pick decomposition of $O_5/O_3$.
  Run: `py -3.14 verify_structures.py` (all residuals 0).
* `probe_pick.py` — 972 exact Sturm counts of real zeros of
  $W(O_n,O_{n-d})$. Run: `py -3.14 probe_pick.py`.
* `probe_cross.py` — 172 exact Sturm counts for $W(H_{n,j},O_{n-j,j})$.
  Run: `py -3.14 probe_cross.py`.

No packages were installed, no builds were run, and no large files or caches
were created. `__pycache__` is gitignored by the repository. Total new size
well under 1 MB.
