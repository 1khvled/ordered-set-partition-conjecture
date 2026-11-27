# Chapter / renewal attack on the cross-family Pick boundary

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-chapter-renewal` was created or edited.  
Outcome: **no proof of $W(H_{n,j},O_{n-d,j})>0$ was obtained.** What *is* proved is a new exact
**positive-matrix (Markov/branching) flow representation** of the chapter EGF (Theorem A),
three exact renewal identities for the boundary Wronskian (Theorem B), and a sharp
obstruction showing why this representation cannot be transferred to the coefficient ratio
by any positive-kernel/coupling argument (Theorem C). The residue-level gap is stated
exactly in §5; it is the same open domination as in `proof-cross-boundary/NOTES.md` §9,
now in "renewal-weighted residue" form.

---

## 0. Setting (chapter dictionary)

Fix $j\ge2$ and write

```
A(z) = sum_{r=1}^{j-1} z^r/r!,   C(z) = e^z - sum_{r=0}^{j-1} z^r/r!,
B(z) = 1 - A(z),   p(z) = z^{j-1}/(j-1)!,   g(z) = C(z)/B(z),
D(z,t) = B(z) - t C(z),   F(z,t) = 1/D(z,t) = sum_n O_n(t) z^n/n!,
H(z,t) = 1/(1-t g(z)) = B(z)F(z,t) = sum_n H_n(t) z^n/n!.
```

Chapter interpretation (`proof-final-inner/NOTES.md` Lemma 5):

* $H_n(t)/n!$ is the weighted count of ordered set partitions of $[n]$ into **chapters**,
  each chapter = one block of size $\ge j$ followed by a string of blocks of size $<j$;
  $h_{n,m}=[t^m]H_n$ counts those into exactly $m$ chapters.
* $H=\sum_{m\ge0}t^m g^m$ is the renewal equation $H=1+tgH$; equivalently
  $H_n=\delta_{n,0}+t\sum_k\binom nk a_k H_{n-k}$ with $a_k=k![z^k]g$, $a_k=0$ for $k<j$.
* $F=H/(1-A)$, i.e. $O_n=\sum_k\binom nk d_{n-k}H_k$ with
  $d_r=r![z^r](1-A)^{-1}\ge0$ (positive outer renewal).
* First-block form: $H_n=t\sum_{q=j}^n\binom nq O_{n-q}$, $O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}+H_n$.

Degrees: $\deg_tO_n=\deg_tH_n=\lfloor n/j\rfloor$ (positive leading coefficients); for
$1\le d\le j$, $\deg O_{n-d}\in\{\deg H_n-1,\deg H_n\}$.

Target: $W(H_n,O_{n-d})>0$ on $\mathbb R$, $W(p,q)=p'q-pq'$ (derivative in $t$).

Elementary facts used throughout (all proved in `proof-final-inner/NOTES.md` Lemma 2/3):

$$A'=1+A-p,\qquad C'=C+p,\qquad B'=B-2+p,\qquad B\,g'=(2-p)g+p. \tag{0.1}$$

---

## 1. Theorem A — the chapter renewal is a positive-matrix Möbius flow

This is the requested Markov/branching representation. It is exact and all-index.

### Theorem A

Let $z_B\in(0,\infty]$ be the least positive zero of $B(z)$ (so $B>0$ on $(0,z_B)$; e.g.
$z_B=1$ for $j=2$). For $z\in(0,z_B)$ define

$$\boxed{\;\phi_z(t_0)=\frac{B(z)\,t_0}{C(z)\,t_0+e^z}\;},\qquad
M(z)=\begin{pmatrix}B(z)&0\\ C(z)&e^z\end{pmatrix}. \tag{A1}$$

Then:

1. **(Positive-matrix / Markov form.)** $M(z)$ has nonnegative entries and
   $\det M(z)=B(z)e^z>0$. Hence $\phi_z$ is an automorphism of the upper half-plane:
   for $\Im t_0>0$,
   $$\Im\phi_z(t_0)=\frac{\det M(z)\cdot\Im t_0}{|C(z)t_0+e^z|^2}>0. \tag{A2}$$
   Equivalently $\phi_z$ is the likelihood ratio of the two positive linear forms given by
   the rows of $M(z)$: a two-state positive Markov (branching) transfer.
2. **(Flow equation.)** $\phi_z$ is the flow of the Riccati generator
   $$\frac{d}{dz}\phi_z(t)=-\frac{2-p(z)}{B(z)}\phi_z(t)-\frac{p(z)}{B(z)}\phi_z(t)^2,
   \qquad \phi_0(t)=t. \tag{A3}$$
3. **(Exact renewal identity.)** For every $t_0$,
   $$\boxed{\;H\!\left(z,\phi_z(t_0)\right)=1+e^{-z}C(z)\,t_0,\qquad
   F\!\left(z,\phi_z(t_0)\right)=\frac{1+e^{-z}C(z)t_0}{B(z)}.\;} \tag{A4}
   $$
4. **(Pick semigroup at EGF level.)** Writing $\alpha(z)=B(z)/C(z)>0$ on $(0,z_B)$,
   $$H(z,t)=\frac{\alpha(z)}{\alpha(z)-t},\qquad F(z,t)=\frac{1/C(z)}{\alpha(z)-t},
   \tag{A5}$$
   so $H(z,\cdot)$ and $F(z,\cdot)$ are rank-one Pick (Nevanlinna) functions of $t$ for
   every $z\in(0,z_B)$. By (A2) the flow $t_0\mapsto\phi_z(t_0)$ and the positive affine
   map $t_0\mapsto1+e^{-z}Ct_0$ both preserve Pick functions; hence the whole chapter
   renewal transports Pick functions to Pick functions in the EGF variable.

*Proof of 1.* $B,C,e^z>0$ on $(0,z_B)$, $\det M=Be^z>0$, and for a real Möbius map
$(at+b)/(ct+d)$ one has $\Im\frac{at+b}{ct+d}=\frac{(ad-bc)\Im t}{|ct+d|^2}>0$. Here
$a=B,b=0,c=C,d=e^z$. $\square$

*Proof of 2.* Put $y=1/t$; then (A3) is equivalent to the linear inhomogeneous equation
$y'=\frac{2-p}{B}y+\frac pB$, which $g=C/B$ solves by (0.1). Hence $y(z)=g(z)+K\,e^z/B(z)$
for a constant $K$; with $y(0)=1/t_0$ and $g(0)=0$, $K=1/t_0$. Therefore
$y(z)=\frac{Ct_0+e^z}{Bt_0}$, i.e. $t(z)=Bt_0/(Ct_0+e^z)=\phi_z(t_0)$. $\square$

*Proof of 3.* Direct: $B-\phi_z(t_0)C=\dfrac{B\,e^z}{Ct_0+e^z}$, hence
$H(z,\phi_z(t_0))=\dfrac{B}{B-\phi_z C}=\dfrac{Ct_0+e^z}{e^z}=1+e^{-z}Ct_0$, and
$F=H/B$. $\square$

*Proof of 4.* $H=B/(B-tC)=\bigl(B/C\bigr)/\bigl(B/C-t\bigr)$; $F=1/(B-tC)=(1/C)/(\alpha-t)$.
For $\Im t>0$, $\Im(\alpha-t)^{-1}=\Im t/|\alpha-t|^2>0$ when $\alpha$ is real. $\square$

Exact symbolic verification of (A3)–(A5), of the Riccati equation
$F'=(2+(t-1)p)F^2-F$, and of the PDE of `proof-final-inner` Lemma 3, for
$j=2,3,4,5$: `verify_flow.py`, all six residuals identically $0$.

### Remark A.6 (what Theorem A is and is not)

(A4) says the entire chapter renewal is a **single positive-matrix Möbius flow**: the
"state" $t_0$ is transported by the Markov transfer $M(z)$, and the value of $H$ along
the transported state is affine with **positive** slope $e^{-z}C(z)>0$. This is the
strongest Markov/branching structure the EGF has. It proves Pickness of $H(z,\cdot)$
for fixed small $z$, *not* of the coefficient polynomials $H_n(t)$, to which we now turn.

---

## 2. Theorem B — exact renewal identities for the boundary Wronskian

### Theorem B1 (product-rule renewal)

For $q=O_{n-d}$ and $1\le d\le j$,

$$W(H_n,q)=\sum_{r=j}^{n}\binom nr\Bigl[\,O_{n-r}\,q+t\,W(O_{n-r},q)\Bigr]. \tag{B1}$$

*Proof.* $W$ is $\mathbb R$-bilinear, so expand $H_n=t\sum_{r\ge j}\binom nr O_{n-r}$ and use
$W(tf,q)=f\,q+t\,W(f,q)$. $\square$

### Theorem B2 (renewal-weighted residue)

Assume $q=O_{n-d}$ has simple roots. At each root $\mu$ of $q$,

$$-\frac{H_n(\mu)}{q'(\mu)}=\frac{-\mu}{q'(\mu)^2}\sum_{r=j}^{n}\binom nr
W(O_{n-r},q)(\mu)=-\mu\sum_{r=j}^{n}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}. \tag{B2}$$

Consequently, under the standing simple-negative-root hypothesis, the following are
equivalent (this is Proposition 2.1 of `proof-cross-boundary/NOTES.md` in renewal form):

$$W(H_n,q)>0\text{ on }\mathbb R
\iff \sum_{r=j}^{n}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}>0
\quad\text{for every root }\mu\text{ of }q. \tag{B3}$$

*Proof.* (B1) at a zero of $q$ kills the $O_{n-r}q$ terms, giving
$W(H_n,q)(\mu)=\mu\sum_{r\ge j}\binom nr W(O_{n-r},q)(\mu)$; and
$W(O_{n-r},q)(\mu)=-O_{n-r}(\mu)q'(\mu)$. Divide by $q'(\mu)^2$ and note
$W(H_n,q)(\mu)=-H_n(\mu)q'(\mu)$; then use the Pick equivalence. $\square$

### Theorem B3 (truncated inner renewal)

Let $d_r=r![z^r](1-A)^{-1}\ge0$. Then for $n\ge j$,

$$\frac{H_n}{t}=\sum_{s=0}^{n-j}\binom ns H_s\sum_{l=0}^{n-j-s}\binom{n-s}{l}d_l. \tag{B4}$$

*Proof.* Insert $O_k=\sum_{s\le k}\binom ks d_{k-s}H_s$ (from $F=H(1-A)^{-1}$) into
$S_n=\sum_{k\le n-j}\binom nk O_k$ and use
$\binom nk\binom ks=\binom ns\binom{n-s}{k-s}$. $\square$

### The summands are not termwise sign-definite

(B1) and (B4) are positive-coefficient renewals, but their summands are not positive:

* Exact witness (`probe_lagcoeff.py`): $j=2$, $n=8$, $q=O_6$, $r=5$ (lag $3$),
  $W(O_3,O_6)$ has a coefficient $-233660$; for $r=5$, even
  $O_3O_6+tW(O_3,O_6)$ has a coefficient $-1080$ (`probe_prodrule.py`). Thus neither
  the Wronskian nor the product-rule bracket is coefficientwise positive for far lags.
* The same phenomenon occurs for every tested $j\ge3$ (scripts in this directory).

This reproduces, in the renewal language, the far-lag obstruction of
`proof-cross-boundary/NOTES.md` §4: only close lags $1\le n-d-k\le j$ are IH-positive.

---

## 3. Theorem C — precise obstruction to the requested transfer

The Markov flow of Theorem A is a genuine positive-kernel representation at the EGF
level. Theorem C shows that **no positive-kernel, coupling or stochastic-monotonicity
argument based on it can prove Pickness of the coefficient ratio**.

### Theorem C1 (no positive-measure representation of the coefficients)

Let $I\subset(0,z_B)$ be any interval and $\mu$ a finite positive measure on $I$. Then
for $n\ge2j$ there is **no** representation

$$H_n(t)=\int_I H(z,t)\,d\mu(z)\quad\text{or}\quad O_n(t)=\int_I F(z,t)\,d\mu(z)
\qquad(\mu\ge0). \tag{C1}$$

*Proof.* By Theorem A(4), each $H(z,\cdot)$ (resp. $F(z,\cdot)$) is a Pick function of $t$;
positive integrals of Pick functions are Pick (with the Nevanlinna representation
obtained by integrating the representing measures). A Pick function rational in $t$
has at most linear growth at infinity, whereas $\deg H_n=\deg O_n=\lfloor n/j\rfloor\ge2$
for $n\ge2j$. Contradiction. $\square$

Thus coefficient extraction — the Cauchy functional $f\mapsto n![z^n]f$ — is an
**irreducibly non-positive** operation on the Pick cone generated by the Markov flow.
This is a structural statement about the requested route, not a finite check.

### Theorem C2 (the kernel has no global orientation)

For $z_1,z_2\in(0,\infty)$ with $C(z_i)\ne0$, set $\alpha_i=B(z_i)/C(z_i)$. The $2\times2$
minor of the rank-one resolvent kernel $F(z,t)$ is

$$\det\begin{pmatrix}F(z_1,t_1)&F(z_1,t_2)\\F(z_2,t_1)&F(z_2,t_2)\end{pmatrix}
=\frac{(\alpha_1-\alpha_2)(t_2-t_1)}
{C(z_1)C(z_2)(\alpha_1-t_1)(\alpha_1-t_2)(\alpha_2-t_1)(\alpha_2-t_2)}. \tag{C2}$$

Whenever the denominators have a fixed sign, the sign of the minor is
$\operatorname{sgn}\bigl[(\alpha_1-\alpha_2)(t_2-t_1)\bigr]$. Now:

* on $(0,z_B)$ we have $B>0$, $B'=B-2+p<0$ (since $B<1$ there), $C>0$, $C'=C+p>0$,
  hence $\alpha'=(B'C-BC')/C^2<0$: $\alpha$ is **strictly decreasing**;
* $\alpha$ crosses $0$ at $z_B$ and, for large $z$,
  $\alpha(z)\sim -z^{j-1}e^{-z}/(j-1)!\to0^-$, so $\alpha$ is eventually **increasing**
  on $(z_B,\infty)$ (a minimum exists; compare `proof-global-wronskian/NOTES.md`
  Prop. E, eq. (12)).

Hence the orientation required by (C2) reverses: there is no subinterval of $(0,\infty)$
containing both a decreasing and an increasing piece on which $F$ is a totally positive
kernel, and in particular no positive-kernel proof on the natural domain. (On the small-$z$
germ $(0,z_B)$ the kernel *is* TP for $t\le0$; but by C1 that germ cannot produce the
coefficients by a positive measure.)

### Theorem C3 (coupling cannot replace the residue inequality)

By (B3), $W(H_n,O_{n-d})>0$ is **equivalent** to the signed renewal inequality
$\sum_{r\ge j}\binom nr O_{n-r}(\mu)/O'_{n-d}(\mu)>0$ at every root $\mu$. Under the
all-lag Pick induction hypothesis, the terms with $1\le d-r\le j$ (close lags, i.e.
$O_{n-r}$ older than $q=O_{n-d}$ by at most $j$) have **negative** sign, while far-lag
terms have no fixed sign (Theorem B's witness above). A coupling/monotone-kernel
argument would have to show termwise dominance of the positive mass over the negative
mass; termwise dominance is false. Hence any successful proof must establish a
**global weighted domination** among residues, which is exactly the open gap.

---

## 4. What the chapter/renewal route achieves

1. **New positive structure (Theorem A).** The EGF $H$ is not merely a renewal series:
   it is a **one-parameter positive-matrix Möbius flow**, with $H(z,\phi_z(t_0))$ affine
   of positive slope in $t_0$. This is a rigorous Markov-chain/branching representation
   of the chapters, and it gives Pickness of every EGF slice $H(z,\cdot)$, $z\in(0,z_B)$,
   by a semigroup argument.
2. **Exact renewal form of the boundary (Theorem B).** The boundary Wronskian is a
   positive binomial renewal of product-rule brackets, and its Pick property is exactly
   a renewal-weighted residue sign condition (B3). The truncated inner renewal (B4) gives
   a second positive representation.
3. **Sharp impossibility (Theorem C).** The flow's positivity cannot be transferred to
   $H_n,O_{n-d}$ by a positive measure/coupling: coefficient extraction is not a positive
   operation on the Pick cone (C1), the kernel has no global orientation (C2), and the
   residue condition is not termwise (C3).

These are consistent with `proof-cross-boundary` §9 and `proof-global-wronskian` §9: the
boundary Pick statement remains open and is at least as hard as full real-rootedness.

---

## 5. The exact remaining gap in renewal form

Under the simple-negative-root and all-lag Pick induction hypothesis for $q=O_{n-d}$:

$$\boxed{\;
W(H_n,q)>0\ \text{on }\mathbb R
\iff
\sum_{r=j}^{n}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}>0
\quad\text{for every root }\mu<0\text{ of }q.
} \tag{G}
$$

For $d=j$ this is the clean domination
$\sum_{k=0}^{n-j-1}\binom nk W(O_{n-j},O_k)>0$ of
`proof-cross-boundary/NOTES.md` (9.1), now read as: the **renewal-weighted residue mass**
of the chapters strictly dominates the far-lag mass at every pole of $O_{n-j}$.
For $d<j$ it is the mixed-sign condition (5.2) of that note.

What a successful chapter/renewal proof would need, beyond Theorem A: a
**measure-valued martingale or branching coupling** whose cross-variation at the poles of
$O_{n-d}$ is exactly the left side of (G), with strictly positive value. Equivalently a
renewal theorem for the *reciprocal* kernel $g(z)$ that does not require $g$ to be PF/CZDS
(the routes excluded in `proof-final-inner/NOTES.md` §2). No such coupling was found.

---

## 6. Files, verification, storage

* `NOTES.md` — this report.
* `verify_flow.py` — exact symbolic residuals for (A3)–(A5), the Riccati equation and the
  inner PDE, $j=2,3,4,5$; all $0$. Run `py -3 verify_flow.py`.
* `verify_renewal.py` — exact verification of (B1) and (B4) for $j=2,\dots,5$,
  $n\le4j$; all pass. Run `py -3 verify_renewal.py`.
* `explore.py`, `probe_cf.py`, `probe_bez.py`, `probe_lagcoeff.py`, `probe_prodrule.py` —
  exploration only: exact row generation, Pick residues/Bezoutians, and far-lag sign
  witnesses. No float is used as a certificate in any claim above.
* No installs, builds, model weights or large caches. Total new size well under 1 MB;
  `__pycache__` is gitignored by the repository.
