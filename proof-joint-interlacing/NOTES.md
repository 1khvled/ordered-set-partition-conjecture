# Joint inner/outer interlacing attack

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-joint-interlacing` was created or edited.
Outcome: **no proof of the joint interlacing target was obtained.**  The work
isolates an exact two-family structure (same full-history recurrence, different
initial data), proves a coupled-induction reduction, and records an exact
obstruction: closing the joint step is equivalent to the open boundary
Wronskian, and the far-lag domination required by that step is not supplied by
the induction hypotheses.  No counterexample was found in the exact checks
reported below.

## 1. Target and notation

Fix $j\ge2$.  Let
$$
A=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad C=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!},\qquad
B=1-A,\qquad F=\frac{1}{B-tC}=\sum_{n\ge0}O_n(t)\frac{z^n}{n!},
$$
$$
H=\frac{B}{B-tC}=\frac{1}{1-t g},\qquad g=C/B,\qquad H=\sum_{n\ge0}H_n(t)\frac{z^n}{n!}.
$$
Then
$$
H_n=O_n-\sum_{r=1}^{j-1}\binom nr O_{n-r}=tS_n,\qquad S_n:=\sum_{k=0}^{n-j}\binom nk O_k,
\qquad n\ge j. \tag{1.1}
$$
The joint target is: for relevant $n,d$ ($1\le d\le j$),
$$
H_{n,j}\ \text{and}\ O_{n-d,j}\ \text{are aligned interlacing with negative residues,}
\tag{1.2}
$$
i.e. $W(H_n,O_{n-d})=H_n'O_{n-d}-H_nO_{n-d}'>0$ on $\mathbb R$, equivalently
$H_n(\mu)/O_{n-d}'(\mu)<0$ at every real root $\mu$ of $O_{n-d}$.  As recorded
in the existing notes, (1.2) implies the $d=j$ boundary Wronskian and, with the
all-lag decomposition, the main WRS.  This note does not rederive those
equivalences; it attacks the joint structure.

## 2. Same recurrence, two initial conditions

**Lemma 2.1 (common recurrence).** Both families satisfy the same full-history
recurrence
$$
X_n=\sum_{k=1}^{j-1}\binom nk X_{n-k}+t\sum_{r=j}^{n}\binom nr X_{n-r},\qquad n\ge j,
\tag{2.1}
$$
with initial data
$$
O_0=1,\quad O_n=d_n=n![z^n](1-A)^{-1}\ (1\le n<j),\qquad
H_0=1,\quad H_1=\cdots=H_{j-1}=0. \tag{2.2}
$$

*Proof.* For $O_n$ this is the first-block recurrence (1.1).  For $H_n$ it is
the positive recurrence (R4) of `proof-final-inner/NOTES.md`, obtained by
extracting coefficients from $B H=t C H$: the identity $B g^m=C g^{m-1}$
($m\ge1$) gives the same binomial convolution after summing over $m$.  The
initial conditions follow from the orders of $C$ and $1-A$.  The identity was
also checked exactly by symbolic coefficient extraction for the ranges below.
$\square$

**Remark 2.2.** The two families are therefore two solutions of the same linear
$t$-dependent renewal operator, differing only in the length-$j$ initial window.
The joint interlacing question asks whether this operator preserves an aligned
interlacing relation between the $H$-solution and lagged $O$-solutions.

## 3. A coupled induction package

Define, for fixed $j$, the joint hypothesis $J(n)$:
- $J_1(n)$: $O_m$ has $\lfloor m/j\rfloor$ simple negative roots for all $m\le n$.
- $J_2(n)$: $W(O_m,O_{m-d})>0$ on $\mathbb R$ for all $m\le n$, $1\le d\le j$, $m-d\ge0$.
- $J_3(n)$: $W(H_m,O_{m-d})>0$ on $\mathbb R$ for all $m\le n$, $1\le d\le j$, $m-d\ge j$.

The lower bound $m-d\ge j$ in $J_3$ excludes the trivial constant-denominator
cases.  By Proposition B of `proof-global-wronskian/NOTES.md`, $J_2(n)$ implies
$J_1(n)$; by Proposition C of the same note, $J_3(n)$ and $J_2(n-1)$ imply the
$d=j$ case of $J_2(n)$.

**Proposition 3.1 (coupled reduction).** Assume $J(n-1)$.  Then:
1. For $d=j$, $J_2(n,d=j)$ follows from $J_3(n,d=j)$ and $J_2(n-1)$.
2. For $d<j$, $J_2(n,d)$ follows from $J_3(n,d)$, $J_2(n-1)$, and the
   quantitative domination
   $$\sum_{r=d+1}^{j-1}\binom nr W(O_{n-d},O_{n-r})(t)\;<\;
     W(H_n,O_{n-d})(t)+\sum_{r=1}^{d-1}\binom nr W(O_{n-r},O_{n-d})(t),
     \tag{3.1}
   $$
   which is the all-lag decomposition identity rearranged.
3. The only missing input for the whole package is $J_3(n)$ (for all relevant
   $d$), together with (3.1) when $d<j$.

*Proof.* (1) is the positive decomposition
$W(O_n,O_{n-j})=W(H_n,O_{n-j})+\sum_{r=1}^{j-1}\binom nr W(O_{n-r},O_{n-j})$.
(2) is the all-lag decomposition (4) of `proof-all-lag-pick/NOTES.md` with the
IH-positive terms moved to the right.  (3) is immediate. $\square$

Thus the coupled induction is sign-closed only after supplying $J_3(n)$; it
does not manufacture $J_3(n)$ from $J(n-1)$.

## 4. Why $J_3(n)$ does not follow from the same recurrence

The common recurrence (2.1) suggests trying to transfer interlacing from the
initial window to all $n$.  Two natural mechanisms fail.

### 4.1 Positive-cone transfer fails

Write $S_n=\sum_{k=0}^{n-j}\binom nk O_k$, $q_d=O_{n-d}$.  At a root $\mu$ of
$q_d$, the residue identity gives
$$
W(H_n,q_d)(\mu)>0 \iff \sum_{k=0}^{n-j}\binom nk \frac{O_k(\mu)}{q_d'(\mu)}>0.
\tag{4.1}
$$
For terms with $k<n-d$ and close lag $(n-d-k)\le j$, $J_2(n-1)$ implies
$W(q_d,O_k)>0$, hence $O_k(\mu)/q_d'(\mu)=W(q_d,O_k)(\mu)/(q_d')^2>0$.  For far
lags $(n-d-k)>j$, the sign of $W(q_d,O_k)$ is not fixed by $J_2(n-1)$ and is
indefinite in exact tests (see `proof-cross-boundary/NOTES.md` §4).  Therefore
(4.1) is a signed domination, not a positive cone.

### 4.2 Common-interlacer / Jacobi-matrix route is not available

A shared Jacobi matrix would give three-term recurrences and automatic
interlacing for orthogonal polynomial sequences.  The degree of $O_n,H_n$ in $t$
is $\lfloor n/j\rfloor$, so the natural index is not the orthogonal degree; no
Jacobi matrix with $O_n,H_n$ as consecutive characteristic polynomials was found.
The HKO/common-interlacing condition for the pair $(H_n,O_{n-d})$ is equivalent
to (1.2) itself, so it cannot be used as an independent hypothesis.

### 4.3 EGF-level Pick transport does not descend

The EGF slices $H(z,\cdot),F(z,\cdot)$ are Pick functions and the chapter
renewal has a positive-matrix Möbius flow (Theorem A of
`proof-chapter-renewal/NOTES.md`), but coefficient extraction is not a positive
operation on that Pick cone (Theorem C1 there).  Hence no EGF-level positive
kernel can certify $J_3(n)$.

## 5. Exact obstruction

The joint interlacing target (1.2) is equivalent, by the residue identity, to
the weighted domination (4.1) at all roots of $q_d$.  For $d=j$ this is exactly
the boundary domination
$$
\sum_{k=0}^{n-j-1}\binom nk W(O_{n-j},O_k)>0 \quad\text{on }\mathbb R, \tag{5.1}
$$
whose close-lag part is IH-positive and far-lag part is sign-indefinite.  For
$d<j$ it is the mixed-sign domination (3.1)/(5.2 of `proof-cross-boundary`).
Thus the joint structure reduces the problem to a single explicit signed-domination
inequality, but does not close it.

A proof of $J_3(n)$ for all $n$ under $J(n-1)$ would, by Proposition 3.1 and the
known Wronskian induction, prove the main WRS conjecture.  Consequently the
joint interlacing target is at least as hard as the open Alexandersson
conjecture; it is not a weaker intermediate lemma.

## 6. Finite exact checks

The script `probe_joint.py` (SymPy exact polynomials, `count_roots`/`real_roots`)
checked, for $j=2,\ldots,5$, $n\le \min(5j,24)$, $1\le d\le j$ (241 pairs):
- $W(H_n,O_{n-d})$ has no real roots;
- $H_n$ and $O_{n-d}$ are real-rooted and interlace with the aligned pattern
  $H,O,\ldots,H$ when $\deg H_n=\deg O_{n-d}+1$ and $O,H,\ldots,H$ when the
  degrees are equal;
- sampled HKO linear combinations $aH_n+bO_{n-d}$ are real-rooted.

The script `verify_joint.py` (SymPy exact arithmetic) checks the common
recurrence (2.1) for both families, the residue identity, and the decomposition
(3.1)/(4.1) over $j=2,\ldots,6$, $n\le 4j$, all $d$.

These are finite evidence only; no float is used as a certificate.

## 7. Files

- `NOTES.md` — this report.
- `probe_joint.py` — interlacing/HKO/Wronskian probe.
- `verify_joint.py` — exact recurrence and residue/decomposition checks.

No packages were installed, no builds were run, and no large caches were
created.  Any `__pycache__` is excluded by the repository `.gitignore`.
