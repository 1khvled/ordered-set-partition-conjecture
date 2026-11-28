# j=2 coefficient TP2: an all-index production-matrix proof

Date: 2026-11-28
Scope: only `proof-j2-tp2-transfer` was created/edited.
No installs/builds. No other agents/models.

## Outcome

**Theorem.** For every \(n\ge0\) and every \(k\ge0\),

\[
\boxed{\;
 d_{n,k}d_{n+1,k+1}-d_{n,k+1}d_{n+1,k}\ge0,
\qquad d_{n,k}=[t^k]O_{n,2}(1,t).
\;}
\]

Thus the full monomial coefficient matrix \(d\), including column \(k=0\),
is TP2. Structural zero minors are allowed.

The proof normalizes by \(n!\), writes the resulting array in the standard
degenerate Riordan/left-production form, proves that its production matrix
\(Q\) is TP2, and transfers TP2 through the finite left-production product by
Cauchy--Binet. The key scalar sequence is

\[
 h_r=\sum_{q=1}^{r}\frac1{(q+1)!}\quad(r\ge1),
\]

whose successive ratios decrease exactly. Neither (H) nor (17) is needed for
this transfer.

This proves only the requested coefficient TP2 statement. It does not prove
TP2 of the falling-basis matrix \(c_{n,i}\), and it does not close the
Bezoutian/interlacing target: `proof-j2-base/NOTES.md` §5.4 already shows that
ratio monotonicity alone is insufficient for that further conclusion.

---

## 1. Normalization and the exact production matrix

Put

\[
 C(z)=e^z-1-z,\qquad F(z,t)=\frac1{1-z-tC(z)},\qquad
 b_{n,k}=\frac{d_{n,k}}{n!}.
\]

Then

\[
 \sum_{n,k\ge0}b_{n,k}z^nt^k=F(z,t).
\]

For later use, write the rising factorial explicitly:

\[
 a^{\overline{k}}=a(a+1)\cdots(a+k-1),\qquad a^{\overline0}=1.
\]

Expansion of \(F\) gives

\[
 d_{n,k}=n!\sum_{m=2k}^{n}
 \frac{(n-m+1)^{\overline{k}}}{m!}A_{m,k},
\]

with the \(k=0\) term interpreted as \(d_{n,0}=n!\). This is formula (13) of
`proof-j2-base`; the overline makes clear that its Pochhammer symbol is
rising, not falling.

Define

\[
 d(z)=\frac1{1-z},\qquad K(z)=\frac{C(z)}z
       =\sum_{q\ge1}\frac{z^q}{(q+1)!}.
\]

The exact identity

\[
 \boxed{F(z,t)=\frac{d(z)}{1-tz\,d(z)K(z)}} \tag{1}
\]

holds because \(z\,d(z)K(z)=C(z)/(1-z)\).

For an array with generating function \(d(z)/(1-tz\,d(z)K(z))\), the standard
left-production matrix is

\[
 Q_{n,0}=d_n,\qquad
 Q_{n,m}=[z^{n-m}]d(z)K(z)\quad(m\ge1).
\]

This is the degenerate-Riordan left-production formula recorded in
`literature.md` §§7.3 and 8.1 (Chen--Fu--Ruan, equation (1.3) and its
Riordan-array corollary). It applies also when the array is noninvertible.
For a finite truncation \(B_N=(b_{n,k})_{0\le n,k\le N}\), the exact
factorization is

\[
 \boxed{
 B_N=Q_N\operatorname{diag}(I_1,Q_{N-1})
       \operatorname{diag}(I_2,Q_{N-2})\cdots
       \operatorname{diag}(I_N,Q_0).
 } \tag{2}
\]

Here \(Q_r\) is the top-left \((r+1)\times(r+1)\) truncation. In each block
diagonal factor the trailing copy of \(Q\) is reindexed from zero.

For the present \(d,K\), set

\[
 h_0=0,\qquad
 h_r=[z^r]d(z)K(z)=\sum_{q=1}^{r}\frac1{(q+1)!}\quad(r\ge1).
\]

Since \(d_n=1\), the production matrix is simply

\[
 \boxed{
 Q_{n,0}=1,\qquad Q_{n,m}=h_{n-m}\quad(m\ge1),
 \quad h_r=0\text{ for }r\le0.
 } \tag{3}
\]

Its first rows are

\[
 \begin{pmatrix}
 1&0&0&0&0\\
 1&0&0&0&0\\
 1&1/2&0&0&0\\
 1&2/3&1/2&0&0\\
 1&17/24&2/3&1/2&0
 \end{pmatrix}.
\]

---

## 2. The scalar monotonicity lemma

**Lemma.** For every \(r\ge1\),

\[
 \frac{h_{r+1}}{h_r}=1+\frac1{(r+2)!h_r}
\]

is strictly decreasing in \(r\).

**Proof.** Since

\[
 h_{r+1}=h_r+\frac1{(r+2)!},
\]

the exact finite difference is

\[
 \frac{h_{r+2}}{h_{r+1}}-\frac{h_{r+1}}{h_r}
 =
 \frac1{(r+3)!h_{r+1}}-\frac1{(r+2)!h_r}<0. \tag{4}
\]

The left denominator is larger: \(h_{r+1}>h_r>0\) and
\((r+3)!>(r+2)!\). This proves the claimed direction. \(\square\)

Consequently, if \(a>b\ge1\) and \(f\ge1\), then

\[
 \frac{h_{a+f}}{h_a}
 =\prod_{r=a}^{a+f-1}\frac{h_{r+1}}{h_r}
 \le
 \prod_{r=b}^{b+f-1}\frac{h_{r+1}}{h_r}
 =\frac{h_{b+f}}{h_b}. \tag{5}
\]

The endpoint audit is important: the product starting at the larger index
\(a\) uses no larger ratios, so it is the lower product. Rearranging (5) gives
\(h_a h_{b+f}\ge h_b h_{a+f}\).

---

## 3. \(Q\) is TP2: complete minor case audit

Take row indices \(i<j\) and column indices \(p<q\). We prove

\[
 Q_{i,p}Q_{j,q}-Q_{i,q}Q_{j,p}\ge0. \tag{6}
\]

### Case 1: \(p=0\)

Here \(Q_{i,0}=Q_{j,0}=1\).

* If \(i\le q\), then \(Q_{i,q}=0\), so the minor is \(Q_{j,q}\ge0\).
* If \(i>q\), both entries in column \(q\) are positive and the minor is
  \[
  h_{j-q}-h_{i-q}>0,
  \]
  because \(j>i\) and \(h_r\) is strictly increasing for \(r\ge1\).

### Case 2: \(p\ge1\) and \(i\le p\)

Both entries in row \(i\) vanish, so the minor is zero.

### Case 3: \(p\ge1\) and \(p<i\le q\)

The minor has the form

\[
 \det\begin{pmatrix}
 h_{i-p}&0\\
 h_{j-p}&h_{j-q}
 \end{pmatrix}
 =h_{i-p}h_{j-q}\ge0.
\]

### Case 4: \(q<i<j\)

All four entries are positive. Put

\[
 a=i-p,\qquad b=i-q,\qquad f=j-i.
\]

Then \(a>b\ge1\), \(f\ge1\), and the minor is

\[
 h_a h_{b+f}-h_b h_{a+f}\ge0
\]

by (5). In fact it is strictly positive.

These cases exhaust all row and column pairs. Therefore \(Q\) is TP2.
\(\square\)

---

## 4. TP2 transfer through the production product

We use the elementary TP2 form of Cauchy--Binet. If \(P,Q\) are nonnegative
TP2 matrices, then for any two rows \(I\) and two columns \(J\),

\[
 \det(PQ)[I,J]
 =\sum_{\substack{K\\|K|=2}}
   \det P[I,K]\det Q[K,J]\ge0. \tag{7}
\]

Thus \(PQ\) is TP2.

Every truncation \(Q_r\) is TP2 by §3. A block diagonal matrix
\(\operatorname{diag}(I_m,Q_r)\) is also TP2: a \(2\times2\) minor lies either
in one identity block, in one \(Q_r\) block, or is zero across blocks.
Equation (2) and (7) therefore show that every finite \(B_N\) is TP2.

Any \(2\times2\) minor of the infinite matrix \(B\) lies in some \(B_N\), so
\(B\) is TP2. Finally,

\[
 d_{n,k}=n!\,b_{n,k},
\]

and positive row scaling multiplies the target minor by \(n!(n+1)!>0\).
Hence \(d\) is TP2 for all \(n,k\). \(\blacksquare\)

---

## 5. Proof-policy audits

### 5.1 Direction audit

The target orientation is exactly

\[
 d_{n,k}d_{n+1,k+1}-d_{n,k+1}d_{n+1,k}\ge0.
\]

The production-minor orientation in (6) is the same row/column orientation.
Exact feasible samples:

\[
 d_3=(6,7,0),\qquad d_4=(24,45,6),\qquad d_5=(120,311,110).
\]

For \(n=4,k=1\),

\[
 45\cdot110-6\cdot311=3084>0.
\]

For \(n=3,k=0\),

\[
 6\cdot45-7\cdot24=102>0.
\]

After normalization these are respectively

\[
 \frac{257}{240}>0,\qquad \frac{17}{24}>0.
\]

Production-matrix samples in the boundary and all-positive cases are

\[
 \det Q[\{2,4\},\{1,3\}]=\frac14>0,
\]

\[
 \det Q[\{3,4\},\{1,2\}]
 =\frac49-\frac{17}{48}=\frac{13}{144}>0.
\]

### 5.2 Monotonicity audit

The only scalar monotonicity used is (4). Its finite difference is printed
exactly and is strictly negative. The product endpoint direction is audited
in (5): for \(a>b\), the product beginning at \(a\) is no larger, yielding
\(h_a h_{b+f}\ge h_b h_{a+f}\), which is the required minor direction.

### 5.3 Feasible-constraint falsification

The proof uses the exact formula for \(h_r\), not merely positivity or
monotonicity of \(h_r\). An arbitrary increasing positive sequence would not
suffice. For example, \(h_1,h_2,h_3=(1,2,5)\) would give the negative minor

\[
 h_2^2-h_1h_3=4-5=-1,
\]

but it has

\[
 \frac{h_3}{h_2}-\frac{h_2}{h_1}
 =\frac52-2=\frac12>0,
\]

so it violates the proved constraint (4). It is therefore infeasible. For the
actual factorial-partial-sum sequence, (4) excludes every such violation.
The subsequent transfer uses only nonnegative \(2\times2\) minors through
(7); there is no unproved optimization step.

### 5.4 Boundary and degeneracy audit

* \(h_0=0\) is used only in zero-entry cases. The ratio argument starts at
  \(r=1\), where \(h_r>0\).
* \(Q\) has diagonal \((1,0,0,\ldots)\) and is noninvertible. Formula (2) is
  the noninvertible left-production factorization; no inverse of \(Q\) or
  \(B\) is used.
* Rows with \(i\le p\), the boundary \(p=0\), the case \(p<i\le q\), and the
  all-positive case \(q<i\) are separately covered in §3.
* Structural zero minors of \(B\) and \(d\) are harmless. The \(k=0\) column
  is included and has the positive sample above.
* Infinite products cause no convergence issue for a fixed minor: it is
  already determined by a finite truncation and finitely many factors in (2).
* \(0!=1\), so positive row scaling also covers \(n=0\).

### 5.5 Symbolic residuals

The generating-function residual is exactly zero:

\[
 \frac{d(z)}{1-tz\,d(z)K(z)}
 =\frac{1/(1-z)}{1-tC(z)/(1-z)}
 =\frac1{1-z-tC(z)}.
\]

The coefficient residual defining \(h_r\) is zero by

\[
 [z^r]\frac1{1-z}\sum_{q\ge1}\frac{z^q}{(q+1)!}
 =\sum_{q=1}^{r}\frac1{(q+1)!}.
\]

`verify_transfer.py` independently checks the coefficient formula and the
finite factorization (2) with exact `Fraction` arithmetic. Those finite loops
are labelled audits, not used as the all-index proof.

---

## 6. Exact obstruction to a stronger TN claim

The production matrix \(Q\) is TP2 but not TN. Already

\[
 \det Q[\{2,3,4,5\},\{0,1,2,3\}]
 =
 \det\begin{pmatrix}
 1&1/2&0&0\\
 1&2/3&1/2&0\\
 1&17/24&2/3&1/2\\
 1&43/60&17/24&2/3
 \end{pmatrix}
 =-\frac1{4320}<0. \tag{8}
\]

Thus the full-TN production theorem cannot be invoked, and this proof does
not claim that \(Q\), \(B\), or \(d\) is TN. The exact negative \(4\times4\)
minor is an obstruction to that stronger certificate, not to the requested
TP2 theorem.

---

## 7. Retraction log required by `AGENT_PROOF_POLICY.md`

1. **Arithmetic retraction.** An early conversational claim used the incorrect
   value \(d_{3,1}=12\) and asserted a \(k=0\) obstruction \(-18\). The exact
   value is \(d_{3,1}=7\), and the minor is \(102>0\). The claim that the
   falling-basis matrix had the same obstruction was also wrong:
   \(c_{4,1}=d_{4,1}+d_{4,2}=51\), giving \(138\), not \(-18\). Both claims
   were retracted immediately.
2. **Kernel-pairing retraction.** A proposed transfer paired \(m<\ell\)
   using \(K_{\ell,m}=-K_{m,\ell}\). This identity is false for the actual
   three-index convolution weight. An exact residual was nonzero already at
   \(n=4,k=1\), so that claimed decomposition was retracted. The valid proof
   above uses a different composition: the degenerate-Riordan left-production
   matrix and TP2 Cauchy--Binet.

No invalid theorem from either retracted attempt remains in this note.

---

## 8. Reproduction

From this directory:

```powershell
py -3.14 verify_transfer.py
```

The script uses only the Python standard library. It prints exact direction
samples, zero factorization residuals, finite TP2 audits, the infeasible
witness rejection, and the exact negative \(4\times4\) minor (8).
