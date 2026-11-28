# Falling-basis TP2 to real roots: conditional transfer and exact obstruction

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-falling-to-realroots` was created/edited.  
No installs/builds. No other agents/models.

## Outcome

Two rigorous results are obtained.

1. **Conditional Bezout transfer.** In the falling-factorial basis
   \(f_i(t)=(t)_i\), write
   \[
   P=\sum_i c_i f_i,\qquad Q=\sum_i e_i f_i,\qquad
   \Delta_{ij}=c_i e_j-c_j e_i\quad(i<j).
   \]
   The root-positive Bezoutian has the exact decomposition
   \[
   \boxed{\quad M(P,Q)=\sum_{i<j}\Delta_{ij}K_{ij},\qquad
   K_{ij}(x,y)=\frac{f_j(x)f_i(y)-f_j(y)f_i(x)}{x-y}.\quad} \tag{1}
   \]
   Thus falling-basis TP2 gives precisely nonnegative weights in (1).
   If the basis kernels satisfy \(K_{ij}\succeq0\), then TP2 implies
   \(M(P,Q)\succeq0\); if all adjacent weights are strictly positive, the
   adjacent rank-one terms already give a positive-definite matrix.

2. **The required basis-kernel hypothesis is false for falling factorials.**
   Adjacent kernels are positive semidefinite:
   \[
   K_{i,i+1}=f_i(x)f_i(y),
   \]
   but every lag at least two is indefinite. Already
   \[
   K_{0,2}=x+y-1,\qquad
   [K_{0,2}]=\begin{pmatrix}-1&1\\1&0\end{pmatrix},
   \]
   whose eigenvalues are \((-1\pm\sqrt5)/2\).

The exact missing input is therefore **Loewner domination of the lag-\(\ge2\)
Bezout kernels by the adjacent positive kernels**, not ratio monotonicity.
In degree two this is the explicit inequality
\[
\boxed{\quad
\Delta_{12}>0,\qquad
\Delta_{01}\Delta_{12}>
\Delta_{02}(\Delta_{12}+\Delta_{02}).
\quad} \tag{2}
\]

A stronger feasible obstruction than `proof-j2-base` §5.4 is proved below.
The pair
\[
P=6+5t+t^2=(t+2)(t+3),\qquad
Q=4+9t+2t^2=(t+\tfrac12)(t+4)
\]
has:

- simple real negative roots and positive coefficients;
- strictly positive monomial TP2 minors \((34,8,1)\);
- strictly positive falling-basis TP2 minors \((42,8,1)\);
- PF (real-nonpositive-rooted) falling-coefficient row polynomials;
- no common root;

but
\[
M(P,Q)=\begin{pmatrix}34&8\\8&1\end{pmatrix},
\qquad \det M(P,Q)=-30<0,
\]
and \(Q(-3)/P'(-3)=5>0\). Hence even the conjunction of the two proved
coefficient-TP2 inputs, real-rootedness of both rows, and row PF does **not**
imply the required Bezout positivity or interlacing.

As a positive finite application, the exact \(O_{4,2},O_{5,2}\) pair is
certified by (2):
\[
M(O_{4,2},O_{5,2})=
\begin{pmatrix}2064&1920\\1920&3084\end{pmatrix},
\qquad \det=2678976>0.
\]
This is a base-case certificate only. No all-\(n\) closure is claimed.

---

## 1. Bezout sign convention and the target

For real polynomials \(P,Q\), use the root-positive kernel

\[
\mathcal K_{P,Q}(x,y)
 =\frac{Q(x)P(y)-Q(y)P(x)}{x-y}. \tag{3}
\]

Let \(M(P,Q)\) be its coefficient matrix in the monomial basis, truncated to
the size relevant for \(\deg P=d\), \(\deg Q\in\{d,d+1\}\). If \(P\) has a
simple root \(\mu\), then

\[
\mathcal K_{P,Q}(\mu,\mu)=-Q(\mu)P'(\mu). \tag{4}
\]

Consequently, for positive leading coefficients and the degree pattern above,

\[
M(P,Q)\succ0
\quad\Longleftrightarrow\quad
\frac{Q(\mu)}{P'(\mu)}<0
\text{ at every root }\mu\text{ of }P. \tag{5}
\]

For completeness, necessity follows from (4). For sufficiency, the sign
condition forces one root of \(Q\) in each alternating interval determined by
the roots of \(P\), with the appropriate exterior interval when
\(\deg Q=d+1\). The degree count exhausts all roots of \(Q\). The standard
Bezout/Obreschkoff criterion then gives positive definiteness. This is the
same orientation as `proof-final-consecutive/NOTES.md` Lemma 2.

For the actual family, \(P=O_{n,2}\), \(Q=O_{n+1,2}\), all coefficients are
positive, \(P(0)>0\), and
\(\deg Q-\deg P\in\{0,1\}\). Thus (5) is exactly the desired consecutive
interlacing criterion once real-rootedness and simplicity are supplied.

---

## 2. Exact falling-basis Bezout decomposition

Put \(f_i(t)=(t)_i=t(t-1)\cdots(t-i+1)\), with \(f_0=1\). For \(i<j\),
define

\[
K_{ij}(x,y)=\mathcal K_{f_j,f_i}(x,y)
 =\frac{f_j(x)f_i(y)-f_j(y)f_i(x)}{x-y}. \tag{6}
\]

The kernel is antisymmetric in its two polynomial arguments, so
\(K_{ji}=-K_{ij}\), while \(K_{ii}=0\). Bilinearity of (3) gives

\[
\begin{aligned}
\mathcal K_{P,Q}
 &=\sum_{i,j}c_i e_j\mathcal K_{f_j,f_i}\\
 &=\sum_{i<j}\bigl(c_i e_j-c_j e_i\bigr)K_{ij}\\
 &=\sum_{i<j}\Delta_{ij}K_{ij}.
\end{aligned} \tag{7}
\]

All sums are finite because the rows have finite degree. This proves (1).

### Theorem 1 (basis-positive transfer)

Let \(D=\max(\deg P,\deg Q)\), so the relevant Bezout matrix is
\(D\times D\). Suppose \(\Delta_{ij}\ge0\) for every \(i<j\).

* If \(K_{ij}\succeq0\) for every \(i<j\), then \(M(P,Q)\succeq0\).
* If in addition \(\Delta_{i,i+1}>0\) for every adjacent pair
  \(i=0,\ldots,D-1\), then \(M(P,Q)\succ0\), even if the higher-lag
  kernels are merely positive semidefinite.

**Proof.** The first assertion is the nonnegative sum (1). For the second,
the adjacent kernels satisfy

\[
K_{i,i+1}=f_i(x)f_i(y), \tag{8}
\]

because \(f_{i+1}=(t-i)f_i\). If \(v_i\) is the coefficient vector of
\(f_i\), their contribution is

\[
\sum_{i=0}^{D-1}\Delta_{i,i+1}v_i v_i^T.
\]

The vectors \(v_0,\ldots,v_{D-1}\) are linearly independent: \(f_i\) is monic
of degree \(i\). Positive adjacent weights therefore make this sum positive
definite. Adding the remaining positive-semidefinite kernels preserves
positive definiteness. \(\square\)

This is a rigorous additional condition under which falling-basis TP2
transfers to Bezout positivity. It is not a condition of the actual falling
basis, as the next section shows.

---

## 3. Why the falling-factorial basis is not a Bezout-positive chain

For \(r=j-i\ge1\), put \(g_r(t)=(t-i)_r\). Since \(f_j=f_i g_r\),

\[
K_{i,i+r}(x,y)
 =f_i(x)f_i(y)\frac{g_r(x)-g_r(y)}{x-y}. \tag{9}
\]

For \(r=1\), the quotient is \(1\), giving (8). For \(r=2\),

\[
K_{i,i+2}(x,y)
 =f_i(x)f_i(y)(x+y-2i-1). \tag{10}
\]

Already \(K_{0,2}=x+y-1\) has coefficient matrix

\[
\begin{pmatrix}-1&1\\1&0\end{pmatrix},
\qquad \det=-1,
\]

so it is indefinite.

More generally, every \(r\ge2\) gives an indefinite kernel. Indeed,
\(g_r\) has \(r\) distinct real roots \(i,i+1,\ldots,i+r-1\). Its derivative
takes both positive and negative values. Choose \(s\ge i\) with
\(g_r'(s)<0\); then \(f_i(s)\ne0\), and

\[
K_{i,i+r}(s,s)=f_i(s)^2g_r'(s)<0.
\]

For sufficiently large positive \(s\), \(g_r'(s)>0\), so the same kernel
takes a positive value. Therefore it cannot be positive semidefinite.

**Exact obstruction to Theorem 1's hypothesis:** the falling-factorial basis
is not a Bezout-positive chain at any lag \(\ge2\). TP2 supplies nonnegative
weights for indefinite kernels, which is not enough.

---

## 4. Degree two: the exact additional inequality

For \(d=2\), the only basis kernels are

\[
K_{01}=1,\qquad K_{02}=x+y-1,\qquad K_{12}=xy.
\]

Equation (1) therefore becomes

\[
\boxed{
M(P,Q)=
\begin{pmatrix}
\Delta_{01}-\Delta_{02}&\Delta_{02}\\
\Delta_{02}&\Delta_{12}
\end{pmatrix}.} \tag{11}
\]

Its determinant is the exact residual identity

\[
\det M(P,Q)
 =\Delta_{01}\Delta_{12}
  -\Delta_{02}(\Delta_{12}+\Delta_{02}). \tag{12}
\]

Hence, for a symmetric \(2\times2\) matrix,

\[
M(P,Q)\succ0
\quad\Longleftrightarrow\quad
\Delta_{12}>0
\ \text{and}\
\Delta_{01}\Delta_{12}>
\Delta_{02}(\Delta_{12}+\Delta_{02}). \tag{13}
\]

This proves (2). TP2 alone gives only
\(\Delta_{01},\Delta_{02},\Delta_{12}\ge0\); it does not control the
right-hand side of (12).

A useful weaker-than-full-target sufficient certificate is adjacent-kernel
domination. In degree two, let

\[
A=\begin{pmatrix}\Delta_{01}&0\\0&\Delta_{12}\end{pmatrix},
\qquad
H=\Delta_{02}\begin{pmatrix}-1&1\\1&0\end{pmatrix}.
\]

If \(\Delta_{01},\Delta_{12}>0\) and \(A\succ H\), then \(A+H=M(P,Q)\succ0\).
Condition (13) is exactly \(A\succ H\). Thus the missing input can be stated
without assuming the desired Bezout positivity: it is a Loewner bound on the
single lag-two term.

---

## 5. Exact feasible obstruction to all currently proved TP2 inputs

Consider

\[
\begin{aligned}
P(t)&=(t+2)(t+3)=6+5t+t^2,\\
Q(t)&=2(t+\tfrac12)(t+4)=4+9t+2t^2.
\end{aligned} \tag{14}
\]

Both have positive coefficients, positive leading coefficient, degree two,
and simple negative roots:

\[
\operatorname{roots}(P)=\{-3,-2\},\qquad
\operatorname{roots}(Q)=\{-4,-\tfrac12\}.
\]

Their discriminants are \(1\) and \(49\), respectively, and they have no
common root.

### 5.1 Monomial TP2 is strict

The monomial coefficient rows are

\[
d(P)=(6,5,1),\qquad d(Q)=(4,9,2).
\]

The three \(2\times2\) minors are

\[
\begin{aligned}
d_0(P)d_1(Q)-d_1(P)d_0(Q)&=54-20=34,\\
d_0(P)d_2(Q)-d_2(P)d_0(Q)&=12-4=8,\\
d_1(P)d_2(Q)-d_2(P)d_1(Q)&=10-9=1.
\end{aligned} \tag{15}
\]

All are strictly positive.

### 5.2 Falling-basis TP2 is strict

Since \(t^2=(t)_2+(t)_1\),

\[
P=6+6(t)_1+(t)_2,\qquad
Q=4+11(t)_1+2(t)_2.
\]

Thus the falling coefficient rows are

\[
c(P)=(6,6,1),\qquad c(Q)=(4,11,2),
\]

with minors

\[
\begin{aligned}
c_0(P)c_1(Q)-c_1(P)c_0(Q)&=66-24=42,\\
c_0(P)c_2(Q)-c_2(P)c_0(Q)&=12-4=8,\\
c_1(P)c_2(Q)-c_2(P)c_1(Q)&=12-11=1.
\end{aligned} \tag{16}
\]

Again all three are strictly positive.

### 5.3 Even row PF does not repair the implication

The falling-coefficient row polynomials are

\[
C_P(y)=6+6y+y^2,\qquad C_Q(y)=4+11y+2y^2.
\]

Their discriminants are \(12\) and \(89\), so both have simple real
nonpositive roots. Thus both falling coefficient rows are PF in the finite
Aissen--Schoenberg--Whitney sense.

### 5.4 The Bezoutian is nevertheless indefinite

For (14), \((\Delta_{01},\Delta_{02},\Delta_{12})=(42,8,1)\). By (11),

\[
M(P,Q)=\begin{pmatrix}34&8\\8&1\end{pmatrix},
\qquad \det M(P,Q)=34-64=-30<0. \tag{17}
\]

At the first root \(\mu=-3\) of \(P\),

\[
\frac{Q(-3)}{P'(-3)}
 =\frac{-5}{-1}=5>0, \tag{18}
\]

which violates the required negative residue. The TP2 weights are all
positive, but the adjacent matrix
\(\operatorname{diag}(42,1)\) does not dominate the indefinite lag-two term

\[
8K_{0,2}=\begin{pmatrix}-8&8\\8&0\end{pmatrix}.
\]

This is an exact feasible-constraint falsification. It shows that no
implication whose hypotheses are only

* monomial TP2,
* falling-basis TP2,
* positivity of coefficients,
* real-rootedness of the two rows,
* PF of their falling-coefficient rows,
* the correct equal-degree pattern, and
* absence of common roots

can imply consecutive Bezout positivity or interlacing.

### 5.5 Even a global TP2 array with the degree pattern can contain the pair

The obstruction is not an artifact of looking at an isolated \(2\times3\)
matrix. Pad the following rows with zeros:

\[
\begin{array}{c|c}
n&\text{falling coefficient row}\\ \hline
0&(1)\\
1&(1)\\
2&(4,1)\\
3&(4,1)\\
4&(6,6,1)\\
5&(4,11,2)\\
6&(1,4,4,4)\\
7&(1,4,4,4)
\end{array}
\]

Every \(2\times2\) minor is nonnegative, and the row degrees are
\(0,0,1,1,2,2,3,3=\lfloor n/2\rfloor\). Rows 4 and 5 are exactly the
noninterlacing pair (14). `verify.py` checks all minors of this finite block
exactly. This finite block is an audit/witness, not an all-index proof.

---

## 6. Application to the actual \(O_{4,2},O_{5,2}\) base case

The exact rows are

\[
O_{4,2}(t)=24+45t+6t^2,\qquad
O_{5,2}(t)=120+311t+110t^2.
\]

Their falling coefficient vectors are

\[
c_4=(24,51,6),\qquad c_5=(120,421,110).
\]

Therefore

\[
\begin{aligned}
\Delta_{01}&=24\cdot421-51\cdot120=3984,\\
\Delta_{02}&=24\cdot110-6\cdot120=1920,\\
\Delta_{12}&=51\cdot110-6\cdot421=3084.
\end{aligned}
\]

Equation (11) gives

\[
M(O_{4,2},O_{5,2})
 =\begin{pmatrix}2064&1920\\1920&3084\end{pmatrix}.
\]

The leading principal minors are

\[
2064>0,\qquad
\det=2064\cdot3084-1920^2=2678976>0.
\]

Equivalently, the banded residual in (12) is

\[
3984\cdot3084-1920(3084+1920)=2678976>0.
\]

Thus the actual first degree-two pair is positive definite and interlaces.
This proves a finite base case by exact algebra; it is not used as a finite
enumeration proof for larger \(n\).

---

## 7. Literature audit: which standard theorems do not close the gap

### 7.1 Falling-factorial preservation runs in the wrong direction

Fisk, *Polynomials, roots, and interlacing*, arXiv:math/0612833v2,
Proposition 7.31, states for the transform \(T(t^i)=(t)_i\):

* \(T\) maps the positive-root/alternating-coefficient class to itself;
* \(T^{-1}\) maps the negative-root/nonnegative-coefficient class to itself.

The second statement is the falling-basis-to-monomial direction. It says
that if a polynomial written in the falling basis has real nonpositive
roots, then the polynomial with the same coefficients in the monomial basis
does too. It does **not** say that TP2 falling coefficients make the
falling-basis polynomial real-rooted.

The converse direction fails even for one row. For example,

\[
6+5t+t^2=(t+2)(t+3)
\]

is real-rooted with nonnegative coefficients, but

\[
T(6+5t+t^2)=6+5(t)_1+(t)_2=6+4t+t^2
\]

has discriminant \(-8\). Thus row PF plus the falling-factorial transform is
not a sufficient real-rootedness certificate.

### 7.2 Brändén's TP2 matrix theorem needs an interlacing input sequence

Brändén, *Unimodality, log-concavity, real-rootedness and beyond*,
arXiv:1410.6601, Theorem 8.5 and Corollary 8.6, characterize constant matrices
mapping an already interlacing sequence of polynomials with nonnegative
coefficients to another such sequence. TP2 of the matrix is the criterion.

Here the natural input sequence would be \(f_i=(t)_i\). It is not
interlacing. By Obreschkoff's pencil characterization, it is enough to note

\[
f_0+\lambda f_2=1+\lambda(t^2-t),
\qquad
\operatorname{Disc}_t=\lambda(\lambda-4).
\]

At \(\lambda=1\) this is \(t^2-t+1\), which is not real-rooted. Hence
\(f_0,f_2\) do not interlace, and the domain hypothesis of Brändén's theorem
fails. This is the same failure expressed by the indefinite \(K_{0,2}\).

### 7.3 Hermite--Biehler does not create the missing part-interlacing

Hermite--Biehler converts a half-plane zero statement for \(P+iQ\) into
interlacing of its real and imaginary parts, with an orientation condition.
Falling-basis TP2 does not establish that required part-interlacing. The
pair (14) is a direct reminder that coefficient TP2 and real-rootedness of
\(P,Q\) separately do not supply the Hermite--Biehler hypothesis.

### 7.4 Pólya--Schur multiplier sequences are not falling-basis TP2

A standard multiplier sequence \(\Lambda=(\lambda_i)\) acts on **monomial**
coefficients and is tested by

\[
G_\Lambda(z)=\sum_i\frac{\lambda_i}{i!}z^i.
\]

Increasing componentwise ratios are not this condition. For example, the
strictly increasing sequence \((1,2,3)\) has normalized test polynomial

\[
1+2z+\frac32z^2,
\qquad
\operatorname{Disc}_z=4-6=-2<0,
\]

so it is not a multiplier sequence.

There is also a basis mismatch in the present application. The strict
falling-basis ratios of (14) are

\[
\lambda=\left(\frac23,\frac{11}{6},2\right),
\]

but the standard monomial multiplier would send the monomial row
\(d(P)=(6,5,1)\) to

\[
\lambda\,d(P)=\left(4,\frac{55}{6},2\right),
\]

not to \(d(Q)=(4,9,2)\). Thus even checking a Pólya--Schur condition on the
falling-basis ratios would not apply the standard theorem to the actual
change of basis. No Pólya--Schur preservation theorem can be invoked from
falling-basis TP2 alone.

### 7.5 Cauchy--Binet proves TP2 transfer, not Bezout transfer

The proved Cauchy--Binet argument in `proof-j2-falling-bezout` transfers
nonnegative \(2\times2\) minors through the Stirling matrix. Equation (1)
is a different operation: a bilinear Bezout decomposition whose basis
kernels are indefinite. Nonnegative weights cannot transfer positive
definiteness through indefinite summands.

Sources used in this audit:

* Steve Fisk, *Polynomials, roots, and interlacing*, Proposition 7.31:
  <https://arxiv.org/abs/math/0612833v2>.
* Petter Brändén, *Unimodality, log-concavity, real-rootedness and beyond*,
  Theorems 7.1, 8.5 and Corollary 8.6:
  <https://arxiv.org/abs/1410.6601>.
* Olga Holtz, *Hermite-Biehler, Routh-Hurwitz, and total positivity*:
  <https://arxiv.org/abs/math/0512591>.
* Petter Brändén and Leonardo Saud Maia Leite, *Totally nonnegative
  matrices, chain enumeration and zeros of polynomials*, Theorems 3.6--3.7:
  <https://arxiv.org/abs/2412.06595>. Their chain-polynomial theorem requires
  a lower-unitriangular TN matrix and its resolvability structure; the
  present falling coefficient TP2 matrix is not such an input.

---

## 8. Exact remaining input for \(O_{n,2}\)

For each actual pair, define

\[
\Delta^{(n)}_{ij}
 =c_{n,i}c_{n+1,j}-c_{n,j}c_{n+1,i}.
\]

The proved TP2 theorem gives \(\Delta^{(n)}_{ij}\ge0\). What is still needed
is the weighted certificate

\[
\boxed{\quad
\sum_{i<j}\Delta^{(n)}_{ij}K_{ij}\succ0.
\quad} \tag{19}
\]

This is equivalent to the target by (1), but it isolates the exact algebraic
gap. A practically weaker sufficient condition is adjacent domination, where
\(D=\max(\deg O_{n,2},\deg O_{n+1,2})\):

\[
\sum_{i=0}^{D-1}\Delta^{(n)}_{i,i+1}v_i v_i^T
\ \succ\
\sum_{0\le i<j-1}\Delta^{(n)}_{ij}K_{ij}. \tag{20}
\]

If all adjacent \(\Delta^{(n)}_{i,i+1}>0\), the left side of (20) is already
positive definite. Thus only the higher-lag terms need domination. In degree
two, (20) is exactly (13), and the strict witness (14) satisfies every TP2
and row-real-rooted hypothesis while violating it.

No proved input in `proof-j2-falling-bezout`, `proof-j2-tp2-transfer`,
`proof-j2-base`, or `proof-u-operator` establishes (19) or (20) for all
\(n\). In particular, (H) and (17) are associated-Stirling inequalities; no
positive Cauchy--Binet or Loewner identity transferring them to the
indefinite \(K_{ij}\) has been proved.

---

## 9. Proof-policy audits

### 9.1 Direction audit

* The target TP2 minor is
  \(\Delta_{ij}=c_i e_j-c_j e_i\ge0\), with the lower row first.
* In (7), the same \(\Delta_{ij}\) multiplies \(K_{ij}\); no sign reversal
  occurs.
* The Bezout kernel (3) has diagonal value \(-Q(\mu)P'(\mu)\), so (5) has the
  required direction \(Q/P'<0\).
* Exact feasible positive sample: for \(O_4,O_5\),
  \((\Delta_{01},\Delta_{02},\Delta_{12})=(3984,1920,3084)\), and
  \(\det M=2678976>0\).
* Exact feasible negative sample: for (14),
  \((\Delta_{01},\Delta_{02},\Delta_{12})=(42,8,1)\), and
  \(\det M=-30<0\).

### 9.2 Monotonicity audit

No scalar monotonicity is used to prove (1), (11), or the obstruction.
Ratio monotonicity is only the positive-denominator equivalent of TP2; it is
never used to infer Bezout positivity. The only multiplier test is an exact
quadratic discriminant, not a monotonicity claim.

### 9.3 Feasible-constraint falsification

Section 5 satisfies all the generic coefficient and row hypotheses listed
there and violates the conclusion. The witness is therefore feasible for the
constraints actually proved. It blocks any transfer argument using only
those constraints.

### 9.4 Boundary and degeneracy audit

* Degree zero: there is no root condition; the empty Bezout statement is
  vacuous.
* Degree one, equal degrees: \(M=[\Delta_{01}]\), so TP2 with
  \(\Delta_{01}>0\) is exactly the required root order.
* Degree two: all three basis kernels and all boundary cases in (11) are
  included.
* \(\deg Q=\deg P+1\): use the next larger coefficient matrix; the algebraic
  decomposition (1) is unchanged, with the appropriate truncation.
* Structural zero coefficients cause no division: (1) and (7) use only
  polynomial products and finite sums.
* \(K_{ii}=0\) and \(K_{ji}=-K_{ij}\) are explicitly used.
* Non-strict TP2 allows zero weights; Theorem 1 then concludes only positive
  semidefiniteness unless strict adjacent weights are separately supplied.
* Infinite arrays cause no convergence issue: each row and each Bezoutian is
  finite.

### 9.5 Symbolic residual audit

`verify.py` checks, with exact integer/rational arithmetic:

* the falling/monomial basis conversions;
* the full decomposition (1), not only adjacent minors;
* the degree-two determinant identity (12);
* all strict minors and discriminants in the obstruction (14)--(18);
* the actual \(O_4,O_5\) matrix and determinant;
* the indefinite basis kernel \(K_{0,2}\);
* all \(2\times2\) minors of the finite global-TP2 embedding in §5.5.

These finite checks are audits of symbolic identities and exact witnesses;
they are not used as an all-index enumeration proof.

### 9.6 Retraction log

An early version of §7.4 claimed that the normalized falling-basis ratio
vector of the strict witness had a negative Pólya--Schur test discriminant.
That arithmetic claim was wrong. After normalizing
\((2/3,11/6,2)\) by its first entry, the test polynomial is

\[
1+\frac{11}{4}z+\frac32z^2,
\qquad
\operatorname{Disc}_z=\frac{25}{16}>0.
\]

The claim is retracted. The corrected §7.4 uses the separate increasing
sequence \((1,2,3)\), whose test discriminant is \(-2\), and relies on the
basis mismatch rather than on a false discriminant for the witness. None of
Theorems 1--2, the decomposition (1), the strict obstruction (14)--(18), or
the \(O_4,O_5\) certificate is affected.

---

## 10. Status and files

**Proved:** the conditional basis-positive transfer theorem, the exact
falling-basis Bezout decomposition, the degree-two criterion, the strict
feasible obstruction, and the \(O_4,O_5\) base certificate.

**Not proved:** all-\(n\) Bezout positivity, real-rootedness, or interlacing
for \(O_{n,2}\). The missing input is domination of the indefinite
lag-\(\ge2\) falling-basis Bezout kernels, equivalently (19) or the weaker
sufficient condition (20).

Files:

* `NOTES.md` — this proof/obstruction record.
* `verify.py` — exact standard-library audit of the identities and witnesses.
