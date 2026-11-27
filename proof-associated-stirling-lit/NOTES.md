# Associated-Stirling SLC: literature search, applicable theorems, and the exact gap

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-associated-stirling-lit` was created/edited.
No installs, builds, external agents, or models were used. No browser was opened.

## 0. Target, normalized correctly

Let $A_{r,k}:=\left\{{r\atop k}\right\}_{\ge 2}$ (partition of $[r]$ into $k$ blocks, all of size
$\ge 2$), with recurrence
$A_{r,k}=kA_{r-1,k}+(r-1)A_{r-2,k-1}$, active range $1\le k\le d:=\lfloor r/2\rfloor$.
The target is the **strong log-concavity (SLC)** inequality

$$
\boxed{\;k^2 A_{r,k}^2\ \ge\ (k+1)^2 A_{r,k-1}A_{r,k+1}\;}\tag{SLC}
$$

Equivalently $b_k:=(k!)^2A_{r,k}$ is log-concave, since

$$
\frac{b_k^2}{b_{k-1}b_{k+1}}
=\frac{(k!)^4 A_{r,k}^2}{\bigl((k-1)!(k+1)!\bigr)^2 A_{r,k-1}A_{r,k+1}}
=\frac{k^2}{(k+1)^2}\,\frac{A_{r,k}^2}{A_{r,k-1}A_{r,k+1}},
$$

so $b_k^2\ge b_{k-1}b_{k+1}\iff$ (SLC). (Check: $(k!)^4/((k-1)!^2(k+1)!^2)=k^2/(k+1)^2$.)

This is the `SLC(r,k)` of `proof-associated-stirling/NOTES.md` §1.2, eq. (5).

## 1. Outcome

**No published theorem was found that implies (SLC) for all $r,k$ for the 2-associated
Stirling row.** Every theorem found is either (a) too weak (ordinary Newton / PF / LC
only), (b) about a different polynomial (ordinary Stirling, unrestricted large-block
counts, or a polynomial with a singleton term), or (c) a direct-induction statement for
the *ordinary* Stirling recurrence (which the present task excludes duplicating).

The closest sufficient condition is identified exactly in §3: real-rootedness of the
**factorial-weighted ordered-partition polynomial** $E_r(x)=\sum_k k!A_{r,k}x^k$.
Exact Sturm-sequence computations (§4) show $E_r$ is real-rooted for all $r\le 40$,
but no published proof of this fact was located; the Symmetric Functions Catalog
records related ordered-lift results that stop short of this polynomial.

## 2. Theorems found, with exact statements, sources, and applicability

### 2.1 Bóna–Mező (2016): real-rootedness of the *unweighted* row

> **Theorem 1 (Bóna–Mező).** Let $n\ge2$. The polynomial
> $D_n(x)=\sum_{k\ge1}D(n,k)x^k$ of degree $\lfloor n/2\rfloor$, where $D(n,k)$ counts
> partitions of $[n]$ into $k$ blocks each of size $\ge2$, has real roots only, all
> simple and non-positive. Moreover roots of $D_n$ and $D_{n-1}$ interlace in the
> prescribed sense.

Source: M. Bóna, I. Mező, *Real zeros and partitions without singleton blocks*,
European J. Combin. **51** (2016) 500–510; arXiv:0705.2734v3,
<https://arxiv.org/abs/0705.2734>, DOI <https://doi.org/10.1016/j.ejc.2015.07.021>.

Same paper, §3.1: for a real-rooted polynomial $\sum_{k=1}^n b_kx^k$ with positive
coefficients, **strong log-concavity**
$b_j^2\ge b_{j-1}b_{j+1}\frac{j+1}{j}\frac{n-j+1}{n-j}$ holds.

**Why insufficient.** With $b_k=A_{r,k}$ and $d=\lfloor r/2\rfloor$, Newton's factor is
$\frac{k+1}{k}\frac{d-k+1}{d-k}$, which dominates $(k+1)^2/k^2$ **iff** $2k\ge d$.
(SLC) requires $(k+1)^2/k^2$; in the lower half $2k<d$ Newton is strictly weaker.
This is already noted in `proof-associated-stirling/NOTES.md` §2.1.
Also SLC is not a formal consequence of real-rootedness: for $(1+x)^d$ at $k=1$,
$c_1^2/(c_0c_2)=2d/(d-1)<4$ for $d\ge3$, whereas SLC needs $\ge4$.

### 2.2 Wang–Yeh / Ahuja / Brenti: PF∞ of associated Stirling rows

> Wang–Yeh, *Polynomials with real zeros and Pólya frequency sequences*, JCTA **109**
> (2005) 63–74, Example 3: the associated Stirling numbers of both kinds satisfy
> recurrences and "Kurtz pointed out both of them are LC, and Ahuja showed that they
> are PF."

Source: <https://arxiv.org/abs/math/0611825>, <https://arxiv.org/html/math/0611825v1>.
Brenti, *Unimodal, log-concave and Pólya frequency sequences in combinatorics*,
Mem. AMS **413** (1989), also lists associated Stirling numbers of the second kind.

**Why insufficient.** PF∞ of a finite row is equivalent to real-rootedness of the
unweighted generating polynomial, so it gives exactly the same Newton bound as §2.1.
PF∞ alone cannot imply SLC: the binomial row $\binom{n}{k}$ is PF∞, yet
$k^2\binom nk^2\ge(k+1)^2\binom n{k-1}\binom n{k+1}$ is equivalent to
$k(n-k+1)\ge(k+1)(n-k)$, i.e. $2k\ge n$, which fails in the lower half.

### 2.3 Canfield–Pomerance (2002): a *stronger-than-SLC* inequality for ordinary Stirling only

> Using $S(n+1,k)=kS(n,k)+S(n,k-1)$ one shows by induction on $n$ that
> $$S(n,k)^2\ \ge\ \Bigl(1+\frac{3}{k}\Bigr)S(n,k-1)S(n,k+1),\qquad 1\le k\le n. \tag{2}$$

Source: E. R. Canfield, C. Pomerance, *On the problem of uniqueness for the maximum
Stirling number(s) of the second kind*, INTEGERS **2** (2002) #A01,
<https://emis.muni.cz/journals/INTEGERS/papers/c1/c1.pdf>.

Since $1+3/k\ge (k+1)^2/k^2$ for $k\ge1$, (2) implies SLC for **ordinary** Stirling
numbers. But the proof uses the ordinary two-term recurrence and is not a theorem
about $A_{r,k}={r\brace k}_{\ge2}$ (three-term recurrence with coefficients $k$ and
$r-1$). Reproducing it for the associated case would be a direct induction proof,
which is explicitly out of scope here.

### 2.4 Ma–Wang (2007): factorial-weighted real-rootedness for *ordinary* Stirling

> If $\mathbf n=(1,1,\dots,1)$, then $\mathcal O(\mathbf n,k)=k!S(n,k)$, and
> $F_n(x)=\sum_{k=1}^n k!S(n,k)x^k$ has only real *simple* zeros in $[-1,0]$
> (indeed $F_{n+1}=xF_n+x(1+x)F_n'$).

Source: S.-M. Ma, Y. Wang, *q-Eulerian polynomials and polynomials with only real
zeros*, arXiv:math/0611826v2, §3.1–3.2 (props. 5–6),
<https://arxiv.org/html/math/0611826v2>. Simion's multiset-composition Sturm
sequence is the underlying result: R. Simion, JCTA **36** (1984) 15–22.

This is exactly the *analogue* of the sufficient condition in §3, but for ordinary
Stirling. No corresponding published statement for $k!A_{r,k}$ was found.

### 2.5 Alexandersson (2026): ordered lifts of "singleton-or-large-block" polynomials

> Let $S_{n,k}(t)$ count set partitions of $[n]$ whose blocks have size $1$ or at
> least $k$, with $S_{n,k}(t)=\sum_r s_{n,k,r}t^r$, and let the **ordered lift** be
> $\widehat S_{n,k}(t)=\sum_r r!\,s_{n,k,r}t^r$. Then for every $k\ge2,n\ge0$,
> both $\widehat S_{n,k}$ and $S_{n,k}$ have only real non-positive zeros, with
> consecutive interlacing. Recurrence:
> $\widehat S_{n,k}=\mathcal T\widehat S_{n-1,k}+\binom{n-1}{k-1}t\,\mathcal T\widehat S_{n-k,k}$,
> $\mathcal T f=(tf)'$.
> The catalog also records that $k$-associated *block-counting* polynomials
> $\sum_k \{n\brace k\}_{\ge k'}$ are real-rooted (the $k'=2$ case is Bóna–Mező).

Source: P. Alexandersson, *Real-rooted polynomials from set partitions*, The
Symmetric Functions Catalog, 18 Nov 2026, "Big blocks" example,
<https://www.symmetricfunctions.com/realRootedSetPartitions.htm>.

**Why insufficient.** Our polynomial $E_r(x)=\sum_k k!A_{r,k}x^k$ counts ordered
partitions **with no singleton blocks**. Alexandersson's ordered lift allows an
arbitrary set of singletons; its EGF is
$e^z\cdot \frac{t(e^z-1-z)}{1-t(e^z-1-z)}$, whereas ours is
$\frac{t(e^z-1-z)}{1-t(e^z-1-z)}$. Indeed
$\widehat S_{n,2}=\sum_{j}\binom{n}{j}E_j$ is the *binomial transform* of $E$,
equivalently $E_n=\sum_j(-1)^{n-j}\binom nj \widehat S_{j,2}$. The inverse binomial
transform is not a Pólya–Schur multiplier (it does not preserve real-rootedness in
general), so this theorem cannot be invoked as a black box to get $E_r$ real-rooted.
The catalog explicitly leaves the analogous unrestricted ordered-family
real-rootedness as a conjecture ("open problems" page,
<https://www.symmetricfunctions.com/realRootedOpenProblems.htm>).

### 2.6 Fang (2026): $k!a_k$ log-concave for unrestricted large-block counts

> For $a_k=[t^k]P_n^{(b)}(t)$ (all partitions of $[n]$, counted by number of blocks
> of size $\ge b$), the sequence $k!a_k$ is log-concave. The paper also observes
> that the factorial-weighted row polynomial itself can have non-real zeros and
> gives the first threshold-four counterexample at $n=11$.

Source: C. Fang, *Eventual interlacing of large-block set-partition polynomials*,
arXiv:2611.01129v1, §8, <https://arxiv.org/abs/2611.01129>.

**Why insufficient.** Log-concavity of $k!a_k$ is Newton-strength:
$a_k^2\ge\frac{k+1}{k}a_{k-1}a_{k+1}$, not the factor $((k+1)/k)^2$ required by
(SLC). It is about a different triangle anyway. It is useful as a caution: the
factorial-weighted analogue is not automatic and must be proved for our triangle.

### 2.7 Other checked sources (no applicable theorem)

* R. Simion (1984) — real-rootedness for multiset compositions; gives ordinary
  Stirling factorial-weighted result but not the associated one.
* Wang–Yeh (2005) general recurrence criterion $A(n,k)=(rn+sk+t)A(n-1,k-1)+
  (an+bk+c)A(n-1,k)$ — associated Stirling rows are covered, but only up to PF/LC.
* Belbachir–Tebtoub (2017), *t-successive associated Stirling numbers*,
  Turkish J. Math 41 (2017) 1279–1291 — Newton + strict LC only; a different
  ("successive", consecutive-elements) family.
* Moll–Ramírez–Villamizar (2018), *Combinatorial and arithmetical properties of the
  restricted and associated Bell and factorial numbers*, arXiv:1707.08138 —
  identities and LC via real-rootedness; no factorial-weighted SLC.
* Qi (2016), *Diagonal recurrence relations, inequalities and monotonicity related
  to the Stirling numbers of the second kind*, Math. Inequal. Appl. 19 — ordinary
  Stirling diagonal inequalities, not our row inequality.
* Mező (2010), *The r-Bell numbers*, arXiv:0909.4417 — real negative zeros of
  $r$-Bell polynomials; different "r" (distinguished elements), not associated.

## 3. The exact sufficient condition and a rigorous derivation of (SLC) from it

Let $d=\lfloor r/2\rfloor$ and define the **ordered 2-associated partition polynomial**

$$
E_r(x)\ :=\ \sum_{k=1}^{d} k!\,A_{r,k}\,x^k .
$$

It counts ordered set partitions of $[r]$ into $k$ blocks, all of size $\ge2$,
weighted by $x^k$; its EGF is $t(e^z-1-z)/(1-t(e^z-1-z))$. From
$A_{r,k}=kA_{r-1,k}+(r-1)A_{r-2,k-1}$ it satisfies the exact recurrence

$$
E_r(x)=x\,E_{r-1}'(x)+(r-1)x\bigl(E_{r-2}(x)+xE_{r-2}'(x)\bigr). \tag{R}
$$

**Proposition.** If $E_r(x)$ has only real zeros for a given $r$, then (SLC) holds
for that $r$ (in fact with an extra factor $\ge1$).

*Proof.* Write $c_k:=k!A_{r,k}\ge0$, so $E_r(x)=\sum_{k=0}^d c_kx^k$ with $c_0=0$.
Newton's inequalities for real-rooted polynomials with non-negative coefficients say

$$
c_k^2\ \ge\ c_{k-1}c_{k+1}\,\frac{k+1}{k}\,\frac{d-k+1}{d-k},
\qquad 1\le k\le d-1 . \tag{N}
$$

(Standard; see Bóna–Mező §3.1 or any Newton-inequality reference.) Now

$$
\frac{c_k^2}{c_{k-1}c_{k+1}}
=\frac{(k!)^2 A_{r,k}^2}{(k-1)!\,(k+1)!\,A_{r,k-1}A_{r,k+1}}
=\frac{k}{k+1}\,\frac{A_{r,k}^2}{A_{r,k-1}A_{r,k+1}} .
$$

Substituting into (N),

$$
\frac{k}{k+1}\frac{A_{r,k}^2}{A_{r,k-1}A_{r,k+1}}
\ge \frac{k+1}{k}\frac{d-k+1}{d-k}
\ \Longrightarrow\
\frac{A_{r,k}^2}{A_{r,k-1}A_{r,k+1}}
\ge \Bigl(\frac{k+1}{k}\Bigr)^2\frac{d-k+1}{d-k}
\ge\Bigl(\frac{k+1}{k}\Bigr)^2 ,
$$

since $1\le k\le d-1$ gives $(d-k+1)/(d-k)>1$. Multiplying by
$k^2/(k+1)^2$ gives (SLC). $\blacksquare$

Thus a proof of $E_r\in\mathrm{RZ}$ for all $r$ closes the target; it is a strictly
stronger statement than (SLC) (it also implies the strict version). Alternatively,
real-rootedness of $\sum_k (k!)^2A_{r,k}x^k$ would suffice (with an even larger
Newton factor), but no published theorem covers that polynomial either.

A weaker published theorem that would still suffice is a *direct* Newton-strength
inequality for the coefficients of $E_r$, namely

$$
(k!A_{r,k})^2\ \ge\ \frac{k+1}{k}\frac{d-k+1}{d-k}
\,(k-1)!A_{r,k-1}\,(k+1)!A_{r,k+1}. \tag{N'}
$$

Note (N') is equivalently $A_{r,k}^2\ge((k+1)/k)^2\frac{d-k+1}{d-k}
A_{r,k-1}A_{r,k+1}$; it is stronger than (SLC) but weaker than full real-rootedness.

## 4. Exact computational evidence (support, not proof)

`explore.py` builds the triangle exactly and uses rational Sturm sequences to count
real roots of $E_r$ ($c=1$ below) and of $\sum_k(k!)^c A_{r,k}x^k$ for $c=2,3$:

* $c=1$: every row $2\le r\le40$ has **all $d=\lfloor r/2\rfloor$ roots real and
  non-positive** (one root at $0$, the remaining $d-1$ strictly negative, no
  positive roots). Hence Newton applies and (SLC) holds there. (The Sturm counter
  in `explore.py` includes the root at $0$ in the "negative" count.)
* $c=2$: same for $2\le r\le30$.
* $c=3$: already fails at $r=16$ (non-real roots) — so the factorial weighting
  cannot be pushed arbitrarily far; $c=1$ is the sharp useful case here.

Run: `py -3.14 explore.py` and `py -3.14 explore2.py` from this directory.

No counterexample to (SLC), to (N'), or to real-rootedness of $E_r$ ($c=1$) was found.

## 5. Exact gap to record

1. The published real-rootedness theorem for this triangle
   (Bóna–Mező 2016, Thm. 1) applies to $\sum_k A_{r,k}x^k$, and Newton there gives
   SLC only for $2k\ge d$; for $2k<d$ the factor $(d-k+1)/(d-k)<(k+1)/k$ is too small.
2. PF∞ (Ahuja / Wang–Yeh / Brenti) and plain log-concavity results for the same row
   give no improvement; binomial coefficients show PF∞ alone cannot imply SLC.
3. Stronger factorial-normalized results exist for **ordinary** Stirling
   (Canfield–Pomerance 2002 inequality (2); Ma–Wang 2007 real-rootedness of
   $\sum k!S(n,k)x^k$), but their proofs use the ordinary two-term recurrence and
   do not transfer to $\ge2$-block partitions with the given three-term recurrence.
   Re-proving them for the associated triangle is a direct induction proof and is
   out of scope for this agent.
4. The 2026 ordered-lift theorem of Alexandersson real-rootedness covers
   singleton-or-large-block partitions; our polynomial is its inverse binomial
   transform, and the inverse binomial transform is not a multiplier sequence, so
   the theorem does not apply. The catalog records no theorem for the no-singleton
   ordered case.
5. Therefore **the missing ingredient is exactly**: a proof that
   $E_r(x)=\sum_k k!A_{r,k}x^k$ is real-rooted (or just the Newton-strength
   coefficient inequality (N')). With that, (SLC) follows rigorously by the
   Proposition in §3. No published source for this ingredient was found in the
   searches performed; the finite-$r$ verification in §4 is evidence only.
