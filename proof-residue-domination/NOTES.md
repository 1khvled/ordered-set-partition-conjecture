# Global signed residue-mass domination — exact obstruction

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-residue-domination` was created or edited.  
Outcome: **no proof of $G_q(\mu)>0$ was obtained.** What *is* proved is an exact
equivalence (Theorem 1): the residue-mass domination at all roots of
$q=O_{n-d,j}$ is **precisely** the statement that $H_{n,j}$ is real-rooted with
roots aligned-interlacing $q$. Consequently (Theorem 2) the domination for all
indices implies the main Wronskian/WRS conjecture (15) and hence full
real-rootedness of the ordered big-block rows — a known open problem. Each of
the three requested routes (uniform contour/saddle estimates, global signed
domination, renewal theorem for the reciprocal kernel) is then shown in
Theorem 3 to be either a restatement of that open interlacing theorem or
structurally inapplicable. The exact induction gap is stated in §5. All
supporting identities and sign patterns are certified exactly/high-precision
in `verify_obstruction.py` (all checks PASS).

---

## 0. Setting and target

Fix $j\ge2$. With $A=\sum_{r=1}^{j-1}z^r/r!$, $C=e^z-\sum_{r=0}^{j-1}z^r/r!$,
$B=1-A$, $D=B-tC$, $F=1/D=\sum_n O_n(t)z^n/n!$, $g=C/B$,
$H=1/(1-tg)=BF=\sum_n H_n(t)z^n/n!$, we have the exact recurrences

$$O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}+t\sum_{r=j}^n\binom nr O_{n-r},
\qquad H_n=O_n-\sum_{r=1}^{j-1}\binom nr O_{n-r}=t\sum_{r=j}^n\binom nr O_{n-r}. \tag{0.1}$$

Degrees: $\deg_t O_n=\lfloor n/j\rfloor$, $\deg_t H_n=\lfloor n/j\rfloor$, both
with positive leading coefficients, and $H_n(0)=0$ for $n\ge j$.

For $1\le d\le j$ put $q:=O_{n-d}$ and, for a root $\mu$ of $q$ with
$q'(\mu)\ne0$,

$$G_q(\mu):=\sum_{r=j}^{n}\binom nr\frac{O_{n-r}(\mu)}{q'(\mu)}. \tag{0.2}$$

The standing induction hypothesis (IH) is: for all $m<n$, $O_m$ has
$\lfloor m/j\rfloor$ simple negative roots and $W(O_m,O_{m-k})>0$ on
$\mathbb R$ for $1\le k\le j$ (all-lag Pick). Under IH every root $\mu$ of
$q$ is negative and simple. Target: $G_q(\mu)>0$ at every such $\mu$.

Write $S_n:=\sum_{r=j}^n\binom nr O_{n-r}$, so $H_n=tS_n$.

---

## 1. The exact residue identity

**Lemma 1.1.** For every root $\mu$ of $q=O_{n-d}$,

$$\boxed{\;G_q(\mu)=\frac{S_n(\mu)}{q'(\mu)}=\frac{H_n(\mu)}{\mu\,q'(\mu)},\qquad
W(H_n,q)(\mu)=-H_n(\mu)q'(\mu)=-\mu\,q'(\mu)^2\,G_q(\mu).\;} \tag{1.1}$$

*Proof.* $H_n=tS_n$ gives $H_n(\mu)=\mu S_n(\mu)$, hence
$S_n(\mu)/q'(\mu)=H_n(\mu)/(\mu q'(\mu))$; this is $G_q(\mu)$ by (0.2). At a
zero of $q$, $W(H_n,q)(\mu)=H_n'(\mu)q(\mu)-H_n(\mu)q'(\mu)=-H_n(\mu)q'(\mu)$.
Substituting $H_n(\mu)=\mu q'(\mu)G_q(\mu)$ gives the second identity. $\qquad\square$

Equivalently, the polynomial identity (exact, verified in `verify_obstruction.py`
V1) is

$$W(H_n,q)+t\,q'\,S_n\equiv0\pmod{q},\qquad H_n-tS_n\equiv0. \tag{1.2}$$

Since $\mu<0$ under IH, (1.1) gives the chain used throughout:

$$G_q(\mu)>0\iff \frac{H_n(\mu)}{q'(\mu)}<0\iff W(H_n,q)(\mu)>0. \tag{1.3}$$

---

## 2. Theorem 1 — the domination is aligned interlacing

**Theorem 1.** Let $q\in\mathbb R[t]$ have $e\ge1$ simple negative roots
$\mu_1<\cdots<\mu_e<0$ and positive leading coefficient. Let $p\in\mathbb R[t]$
have positive leading coefficient, $p(0)=0$, and $\deg p\in\{e,e+1\}$ (these
hold for $p=H_n$, $q=O_{n-d}$ by §0). Then the following are equivalent:

1. **(D)** $p(\mu_i)/q'(\mu_i)<0$ for every $i$ (equivalently, for $p=H_n$,
   $G_q(\mu_i)>0$ by (1.3));
2. **(I)** $p$ has only simple real nonpositive roots and they are **aligned**
   with the $\mu_i$: writing the roots of $p$ as
   $\nu_1<\cdots<\nu_{\deg p-1}<0$ together with $\nu_{\deg p}=0$,
   - if $\deg p=e+1$: $\nu_1<\mu_1<\nu_2<\mu_2<\cdots<\mu_e<\nu_{e+1}=0$;
   - if $\deg p=e$: $\mu_1<\nu_1<\mu_2<\nu_2<\cdots<\nu_{e-1}<\mu_e<\nu_e=0$;
3. **(P)** $W(p,q)=p'q-pq'>0$ on $\mathbb R$ (as a polynomial, i.e. also at the
   $\mu_i$);
4. **(R)** $p/q=a t^{\deg p-e}+b+\sum_i r_i/(t-\mu_i)$ has $r_i<0$ for all $i$,
   with $a=\operatorname{lc}(p)/\operatorname{lc}(q)>0$; hence $p/q$ is a
   Pick/Nevanlinna function when $\deg p=e+1$.

*Proof.*
(D)$\Rightarrow$(R) is the definition $r_i=p(\mu_i)/q'(\mu_i)$.
(R)$\Rightarrow$(P): $W(p,q)=q^2(p/q)'$ and $(p/q)'=a(\deg p-e)-\sum_i r_i/(t-\mu_i)^2
\ge a(\deg p-e)>0$ everywhere off the poles (and $W$ is a polynomial, so
$>0$ everywhere); at $\mu_i$ one has $W(\mu_i)=-p(\mu_i)q'(\mu_i)>0$.
(P)$\Rightarrow$(D): evaluate $W(\mu_i)=-p(\mu_i)q'(\mu_i)>0$.
(D)$\Rightarrow$(I): $q'(\mu_i)$ has sign $(-1)^{e-i}$, so (D) forces
$\operatorname{sgn}p(\mu_i)=(-1)^{e-i+1}$. With $p(0)=0$, $p(-\infty)$ of sign
$(-1)^{\deg p}$, and positive leading coefficient, the intermediate value
theorem gives exactly one root of $p$ in each gap displayed in (I) (details:
for $\deg p=e+1$ the sign of $p$ alternates across the ordered points
 $-\infty,\mu_1,\ldots,\mu_e,0$ starting with $(-1)^{e+1}$ and ending at $0$;
 for $\deg p=e$ it starts with $(-1)^{e}$ and the first gap has no root). This
produces $\deg p$ distinct simple roots, so $p$ is real-rooted and aligned.
(I)$\Rightarrow$(D): at $\mu_i$, the number of **negative** roots $\nu_j$
exceeding $\mu_i$ is $e-i$ in both cases (the root $\nu_{\deg p}=0$ is excluded
from the product in (2.1)); hence
$\operatorname{sgn}\prod_{j<\deg p}(\mu_i-\nu_j)=(-1)^{e-i}$ in both cases, and
dividing by $q'(\mu_i)$ of sign $(-1)^{e-i}$ gives
$p(\mu_i)/q'(\mu_i)<0$. $\qquad\square$

**Root formula (used in §4).** If $\nu_1<\cdots<\nu_{e'}<0$ are the nonzero
roots of $H_n$ ($e'=\deg H_n-1$), then

$$G_q(\mu_i)=\frac{\operatorname{lc}(H_n)\prod_{j=1}^{e'}(\mu_i-\nu_j)}{q'(\mu_i)}. \tag{2.1}$$

Its sign at $\mu_i$ depends on the **whole** interlacing pattern (how many
$\nu_j$ lie to the right of $\mu_i$), not on local data at $\mu_i$ alone.

---

## 3. Theorem 2 — no free lunch

Let **(D)** denote the domination at index $n$, $d=j$:
$G_{O_{n-j}}(\mu)>0$ at every real root $\mu$ of $O_{n-j}$.

**Theorem 2 (no free lunch, conditional form).** Assume the all-lag
induction hypothesis of the other tracks: for all $m<n$, $O_m$ is real-rooted
with $W(O_m,O_{m-k})>0$ on $\mathbb R$ for $1\le k\le j$. Then (D) implies the
main WRS at index $n$,

$$W(O_{n,j},O_{n-j,j})(t)>0\quad\text{on }\mathbb R\ (\text{hence on }t\le0), \tag{15$_n$}$$

and therefore, by the IVT argument of Proposition B of
`proof-global-wronskian/NOTES.md`, $O_{n,j}$ is real-rooted with simple
negative roots and $O_{n-j,j}\interl O_{n,j}$. Since (D) is to be proved at
*every* index, proving it for all $n$ under the standing IH closes exactly the
induction that proves the main WRS conjecture (15) and full real-rootedness.
Hence (D) is at least as hard as the open Alexandersson conjecture.

*Proof.* Put $q=O_{n-j}$. By (D) and Theorem 1, $W(H_n,q)>0$ on $\mathbb R$.
The exact decomposition (`proof-all-lag-pick/NOTES.md` eq. (4); verified
exactly in V4) at $d=j$ is

$$W(O_n,q)=W(H_n,q)+\sum_{r=1}^{j-1}\binom nr W(O_{n-r},q). \tag{3.1}$$

Each lag $j-r\in[1,j-1]$ and index $n-r<n$, so by IH every term in the sum
is $>0$ on $\mathbb R$. Hence $W(O_n,O_{n-j})>0$ on $\mathbb R$, giving
(15$_n$); the IVT argument then gives the stated roots/interlacing of $O_n$.
$\qquad\square$

**Remark 2.2 (what (D) does *not* by itself give).** For $d<j$ the all-lag
decomposition (eq. (5) of `proof-all-lag-pick/NOTES.md`) subtracts IH-positive
terms, so $W(H_n,O_{n-d})>0$ (equivalently (D) at lag $d$) is only one
summand of the quantitative domination (6) of that note; it does not alone
imply all-lag Pick at index $n$. The $d=j$ step above is the one lag with no
subtracted terms, which is why (D) closes the main-WRS induction.

**Corollary 3.1.** Any proof of (D) under the standing IH necessarily contains
a real-rootedness argument for $O_n$ (Theorem 2). A "shortcut" proof by
contour estimates, a positive kernel, or a renewal theorem that merely
*assumes* the rows are real-rooted is circular: the assumption is the
conclusion of Theorem 2.

---

## 4. Theorem 3 — the three requested techniques

### 4.1 Uniform contour / saddle-point estimates cannot work

For fixed $\mu<0$, Cauchy's formula and (1.1) give

$$G_q(\mu)=\frac{n!}{2\pi i\,\mu\,q'(\mu)}\oint \frac{B(z)}{B(z)-\mu C(z)}\,
\frac{dz}{z^{n+1}},\qquad \text{poles at } g(z)=1/\mu<0. \tag{4.1}$$

Two independent obstructions:

1. **Global vs local.** By (2.1), the sign of $G_q(\mu_i)$ at a root $\mu_i$ of
   $q$ is determined by the position of $\mu_i$ relative to *all* roots
   $\nu_j$ of $H_n$. A saddle-point/steepest-descent estimate of the integral
   (4.1) at a single $\mu_i$ uses only the local geometry of $B-\mu_i C$ near a
   saddle and cannot determine that global count. Certifying the sign at all
   roots simultaneously is, by Theorem 1, exactly proving the aligned
   interlacing of $H_n$ with $q$ — a real-rootedness statement (Corollary 3.1).

2. **No uniform contour.** As $\mu$ runs over the $e\sim n/j$ distinct roots of
   $q$, the values $1/\mu$ are distinct negative numbers, and the pole set
   $\{z:g(z)=1/\mu\}$ of the integrand in (4.1) changes (in number and
   position) with $\mu$; $g=C/B$ is not globally monotone (equivalently
   $\alpha=B/C$ is not, §4.3). There is no fixed saddle/contour whose residue
   contribution has a uniform sign across the root set. Any "uniform over all
   roots" estimate would have to track the argument-principle count of zeros
   of $1-\mu g(z)$ inside the contour for each $\mu$, i.e. again the root
   configuration of $H_n$.

### 4.2 Global signed domination is the conjecture itself

By Theorem 1 the domination at all roots of $q$ is equivalent to
$W(H_n,q)>0$ on $\mathbb R$; by Theorem 2 this inequality at $d=j$, supplied
at every index, is exactly the missing input that closes the induction
proving (15) and full real-rootedness. So a "global signed domination
theorem" producing $G_q>0$ *is* a real-rootedness/interlacing theorem, not a
weaker input. In particular, under IH the induction step is **not** closed:
IH controls only $O_m$, $m<n$, whereas aligned interlacing of $H_n$ with $O_{n-d}$
is the separate conjecture (H) of `proof-final-inner/NOTES.md` (see §5).

### 4.3 Renewal theorem for the reciprocal kernel is structurally inapplicable

At EGF level $H(z,\cdot)$ is a Pick function and the renewal $H=1/(1-tg)$ is a
positive Möbius flow (`proof-chapter-renewal/NOTES.md`, Theorem A). But:

* **Coefficient extraction is not a positive operation.** By Theorem C1 of
  that note, for $n\ge2j$ there is no positive measure $\mu$ with
  $H_n(t)=\int H(z,t)\,d\mu(z)$, since a positive integral of Pick slices is
  Pick of at most linear growth while $\deg H_n\ge2$. Hence flow positivity
  cannot be pushed to (D).

* **The reciprocal kernel is not Stieltjes/Pick.** With $\alpha=B/C$ one has
  exactly (verified in V5)

  $$\alpha'(z)=\frac{e^z(p(z)-2)+2E_{j-2}(z)}{C(z)^2},
  \qquad E_{j-2}=\sum_{r=0}^{j-2}\frac{z^r}{r!}, \tag{4.2}$$

  whose numerator is negative near $z=0$ (e.g. $-0.0998$ at $z=0.1$) and
  positive for large $z$ (e.g. $8.7\times10^9$ at $z=20$, $j=2$). So $\alpha$
  decreases then increases: $1/g=B/C$ is not a globally monotone/Stieltjes
  reciprocal kernel, and no standard renewal domination theorem (which needs
  a fixed total-positivity/Stieltjes orientation) applies to it.

* **Termwise positivity fails.** The residue mass is genuinely signed for far
  lags: for $j=2$, $n=8$, $q=O_6$, the residues $O_k(\mu)/q'(\mu)$ take **both**
  signs for lags $3,4,5,6$ (V3), while only lags $1,2$ are IH-positive. Any
  renewal theorem for (D) must therefore be a *global* domination, which by
  §4.2 is the conjecture.

---

## 5. The exact induction gap

Assume IH (all-lag Pick for $O_m$, $m<n$). To close the induction at $(n,d)$
one must prove $G_{O_{n-d}}(\mu)>0$ at roots $\mu$ of $q=O_{n-d}$. By Theorem 1
this is equivalent to $H_n$ being real-rooted with roots aligned-interlacing
$q$. But:

* $H_n$ is **not** among the polynomials covered by IH (it involves $O_n$
  itself: $H_n=O_n-\sum_{r<j}\binom nr O_{n-r}$);
* the close-lag part of $G_q$ is IH-positive (older-by-$\le j$ residues $>0$),
  while the far-lag part has no fixed sign (§4.3, V3);
* hence IH supplies no sign information about $H_n$'s roots, and the required
  domination cannot be derived from IH by termwise sign reasoning.

For $d=j$ the missing statement is exactly the $H$-boundary
$W(H_n,O_{n-j})>0$ (`proof-cross-boundary/NOTES.md`, §9.1 / (9.1)) and, by
Theorem 2, closing it for all $n$ proves (15). For $d<j$ the all-lag
decomposition (`proof-all-lag-pick/NOTES.md`, eq. (5)) contains IH-positive
**subtracted** terms, so even granted the $d=j$ boundary one still needs the
quantitative domination (6) of that note.

Thus the gap is irreducible to termwise positivity, positive-kernel transfer,
global Stieltjes monotonicity of $1/g$, or IH alone. It is exactly the
aligned-interlacing/real-rootedness conjecture (H)+all-lag Pick for the family
$\{H_{n,j}\}\cup\{O_{n,j}\}$.

---

## 6. Exact verification (`verify_obstruction.py`, all PASS)

* **V1** — 355 exact triples $(j,n,d)$ ($2\le j\le6$, $n\le4j$):
  $H_n=tS_n$ and $W(H_n,q)+t q'S_n\equiv0\pmod q$. Establishes (1.1)/(1.2).
* **V2** — same 355 triples: at every negative root $\mu$ of $q=O_{n-d}$,
  $G_q(\mu)>0$ (40-digit root approximations), and the roots of $H_n$
  aligned-interlace those of $q$ in the sense of Theorem 1(2). Establishes
  that the domination *is* the interlacing statement empirically.
* **V3** — $j=2,n=8,q=O_6$: 4 far lags have residues of both signs;
  termwise positivity is false (exact witness).
* **V4** — exact decomposition (3.1), 56 cases.
* **V5** — symbolic numerator of $\alpha'$ changes sign (Stieltjes obstruction).

`explore_residue.py` contains the underlying row generator and sample residue
tables. No floats are used as certificates in any theorem above; root-based
checks use high-precision approximations and are labelled as evidence only.

---

## 7. Conclusion

The global signed residue-mass domination $G_q(\mu)>0$ was **not proved**. The
exact result of this track is:

1. **Theorem 1** — $G_q(\mu)>0$ at all roots $\mu$ of $q=O_{n-d}$ $\iff$
   $H_n$ is real-rooted with roots aligned-interlacing $q$ (equivalently
   $W(H_n,q)>0$ on $\mathbb R$).
2. **Theorem 2** — the domination at $d=j$, supplied at every index under the
   standing IH, closes exactly the induction proving (15) and hence full
   real-rootedness of $\{O_{n,j}\}$, a known open conjecture; so any proof of
   the domination is already a real-rootedness proof (Corollary 3.1).
3. **Theorem 3** — the three requested techniques fail for exact reasons:
   contour/saddle estimates are local while the sign is globally determined
   by interlacing (§4.1); global signed domination *is* the conjecture (§4.2);
   and the reciprocal kernel $1/g=B/C$ is not Stieltjes while coefficient
   extraction is not positive (§4.3).

No counterexample was found. All exact identities and sign patterns are
certified in `verify_obstruction.py`. No installs, builds, or large caches;
new data is well under 1 MB and `__pycache__` is gitignored.
