# PSD / symmetric-matrix realization attack

Date: 2026-11-28  
Scope: only `ordered-big-blocks/proof-psd-realization` was created or edited.

## Outcome

**No all-index proof was obtained, and no counterexample was found.**

The track produced three exact results:

1. A new Cramer determinant representation
   \[
   O_{n,j}(t)=\det(A_n+tB_n)
   \]
   comes directly from the first-block recurrence.  Thus
   \[
   O_{n,j}(t)=O_{n,j}(0)\det(I+tM_n),\qquad M_n=A_n^{-1}B_n,
   \]
   is an explicit coefficient-only candidate, not a construction from the unknown roots.
2. This natural candidate has an all-index rank-one update, but it is a
   **nonsymmetric** update.  It fails symmetry, diagonalizability, and a
   rational-entry test by explicit small examples.  Hence it is not the
   required PSD rank-one/bordering chain.
3. The requested `j=2` kill-test does give finite coefficient-only symmetric
   PSD rank-one updates and PSD borderings through `m=5`.  Exact rational
   Lanczos produces a canonical rational tridiagonal matrix at every step.
   All computed Lanczos coefficients are positive, but their denominators
   explode (73 digits for an odd-step diagonal entry at `m=5`), with no
   recognizable closed form.  This kills the hoped-for *simple structured
   formula*, not the conjecture or the existence of less explicit PSD matrices.

The exact verification command is:

```powershell
py -3.14 proof-psd-realization/verify.py
```

It uses only integer/`Fraction` arithmetic and performs no root finding.

---

## 1. Exact coefficient conditions for a PSD realization

Fix `j,n`, put `q=floor(n/j)`, and write

\[
O_{n,j}(t)=\sum_{k=0}^q a_k t^k,\qquad a_0>0,\qquad
c_k=\frac{a_k}{a_0}.
\]

The determinant model must be normalized:

\[
\widehat O_{n,j}(t):=\frac{O_{n,j}(t)}{O_{n,j}(0)}
=\sum_{k=0}^q c_k t^k.
\]

Trying to use `O_n(0) det(I+t M_n)` without this normalization is harmless,
but comparing determinants for different `n` requires the normalized
polynomials because `O_n(0)` varies.

### 1.1 Individual rows

A `q x q` real symmetric matrix `M` satisfies

\[
\det(I+tM)=\sum_{k=0}^q e_k(M)t^k.
\]

Consequently:

* `M` PSD is sufficient for all roots of `widehat O_n` to be real and negative.
* Conversely, if `widehat O_n` has `q` negative real roots
  `rho_1,...,rho_q`, then the diagonal matrix with entries
  `lambda_i=-1/rho_i` is positive definite and has the required determinant.
* Therefore **individual PSD realization is equivalent to real-rootedness**.
  It is not a weaker theorem.

Necessary coefficient conditions include Newton's inequalities

\[
\left(\frac{c_k}{\binom qk}\right)^2
\ge
\frac{c_{k-1}}{\binom q{k-1}}
\frac{c_{k+1}}{\binom q{k+1}},
\qquad 1\le k<q.
\]

These are not sufficient in general.

A complete coefficient-only test is the following Hankel criterion.  Let
`s_k=sum_i rho_i^k` be the power sums of the roots `rho_i` of `widehat O_n`,
computed by Newton identities, and put

\[
\mu_k=(-1)^k s_k=\sum_i\lambda_i^k,\qquad
H_q=(\mu_{r+s})_{0\le r,s\le q}.
\]

Then

\[
\widehat O_n\text{ has }q\text{ negative real roots}
\quad\Longleftrightarrow\quad H_q\succeq0.
\]

For simple roots the leading `q x q` Hankel determinant is positive and the
full `(q+1) x (q+1)` determinant is zero; in general the rank is the number
of distinct roots.  To see sufficiency, a nonreal conjugate pair `a+-ib` would allow the real
polynomial `(x-a)` times the product over the real roots other than `a` to
give a negative value of `f^T H_q f`.  This criterion is exact but still
equivalent to real-rootedness.

### 1.2 Rank-one steps

If `M^+=M+vv^T` with `M` PSD, the determinant lemma gives

\[
\det(I+tM^+)=\det(I+tM)
+t\,v^T\operatorname{adj}(I+tM)v.
\]

In normalized coefficient form,

\[
\frac{\widehat O_{n+1}}{\widehat O_n}
=1+t\,v^T(I+tM_n)^{-1}v.
\]

A necessary condition is coefficient monotonicity.  If `C_k(M)` is the
`k`-th additive compound, then

\[
e_k(M+vv^T)=e_k(M)+v^TC_k(M)v\ge e_k(M),
\]

because `C_k(M)` is PSD.  Thus `c_k(O_{n+1}) >= c_k(O_n)` is necessary.
The same conclusion holds for a PSD bordering
`[[M,v],[v^T,alpha]]`, using

\[
e_k\begin{pmatrix}M&v\\v^T&\alpha\end{pmatrix}
=e_k(M)+\alpha e_{k-1}(M)+v^TC_{k-1}(M)v.
\]

All normalized coefficients were nondecreasing in the exact checked ranges;
this is much weaker than interlacing.

### 1.3 PSD borderings and equivalence of the full chain

For a degree jump, a PSD bordering gives

\[
\det\begin{pmatrix}I+tM&tv\\tv^T&1+\alpha t\end{pmatrix}
=\det(I+tM)(1+\alpha t)
-t^2v^T\operatorname{adj}(I+tM)v.
\]

Cauchy interlacing applies to the principal block `M`.  A rank-one PSD
update gives Weyl interlacing.

Conversely, strict interlacing of two equal-degree real-rooted polynomials
is equivalent to a rank-one PSD realization: in the eigenbasis of the first
matrix, the positive residues of the rational determinant ratio determine
the squared components of the update vector.  Strict interlacing across a
degree jump is similarly equivalent to a PSD bordering, by the secular
equation and Cauchy interlacing.

Hence an all-index chain of the proposed form is **equivalent to the
conjecture**, not a reduction to a weaker positivity statement.  Any proof
that such matrices exist must prove the same root/interlacing content.

---

## 2. New exact Cramer construction from the first-block recurrence

Define the infinite unit lower-triangular matrix

\[
K_{n,m}=\delta_{n,m}
-\mathbf 1_{1\le n-m<j}\binom n{n-m}
-t\,\mathbf 1_{n-m\ge j}\binom n{n-m}.
\]

The first-block recurrence is exactly

\[
K\,(O_0,O_1,\ldots)^T=e_0.
\]

For `n>=1`, Cramer's rule gives a finite determinant.  Let `H_n(t)` be the
`n x n` matrix with rows `r=1,...,n` and columns `m=0,...,n-1`:

\[
(H_n)_{r,m}=
\begin{cases}
-1,&r=m,\\
\binom r{r-m},&1\le r-m<j,\\
t\binom r{r-m},&r-m\ge j,\\
0,&r-m<0.
\end{cases}
\]

Then

\[
\boxed{O_{n,j}(t)=\det H_n(t).}
\]

Write `H_n=A_n+tB_n`.  Since `det A_n=O_{n,j}(0)>0`,

\[
\boxed{M_n=A_n^{-1}B_n},\qquad
O_{n,j}(t)=O_{n,j}(0)\det(I+tM_n).
\]

Every entry of this candidate is explicit in binomial coefficients and
earlier constant terms.  No unknown roots are used.

### 2.1 All-index rank-one update of the Cramer candidate

The matrices nest exactly:

\[
H_{n+1}=
\begin{pmatrix}
H_n&h_n\\
g_n^T&c_n
\end{pmatrix}.
\]

Equivalently,

\[
A_{n+1}=\begin{pmatrix}A_n&a_n\\b_n^T&c_n\end{pmatrix},
\qquad
B_{n+1}=\begin{pmatrix}B_n&0\\e_n^T&0\end{pmatrix}.
\]

The zero last column of `B_{n+1}` is important.  With
`E_n=diag(M_n,0)`, one obtains

\[
A_{n+1}(M_{n+1}-E_n)
=e_{n+1}(e_n^T-b_n^TM_n,\ 0).
\]

Therefore

\[
\boxed{
M_{n+1}=
\begin{pmatrix}M_n&0\\0&0\end{pmatrix}
+u_nw_n^T,
}
\]

where

\[
u_n=A_{n+1}^{-1}e_{n+1},\qquad
w_n^T=(e_n^T-b_n^TM_n,\ 0).
\]

This is an all-index rank-one identity.  The vector `u_n` has the simple
positive formula

\[
\boxed{(u_n)_r=\frac{O_{r,j}(0)}{O_{n+1,j}(0)},\qquad 0\le r\le n.}
\]

For `j=2`, `b_n^TM_n=0` because the only nonzero entry of `b_n` multiplies
the last column of `M_n`, which is zero.  Hence

\[
(w_n)_k=\binom{n+1}{k+1}\quad(0\le k<n),\qquad (w_n)_n=0.
\]

So even in the simplest case the update is `u_n w_n^T`, not `u_n u_n^T`.

### 2.2 Exact obstruction: this rank-one update is not PSD

A PSD rank-one matrix must be symmetric.  The Cramer update is already
nonsymmetric at `j=2,n=1 -> 2`:

\[
M_{2,2}=\begin{pmatrix}1/2&0\\1/2&0\end{pmatrix}.
\]

There is a stronger diagonalizability obstruction at `j=2,n=4`:

\[
M_{4,2}=
\begin{pmatrix}
17/24&2/3&1/4&0\\
17/24&2/3&1/4&0\\
5/12&4/3&1/2&0\\
1/4&1&3/2&0
\end{pmatrix}.
\]

Its determinant is

\[
\det(I+tM_{4,2})=1+\frac{15}{8}t+\frac14t^2.
\]

Thus zero has algebraic multiplicity two in its characteristic polynomial.
Its last column is zero, while its rank is three, so zero has geometric
multiplicity one.  The zero eigenvalue has a nontrivial Jordan block.
Therefore `M_{4,2}` is not similar to any real symmetric matrix, PSD or
otherwise.  No change of basis can repair this Cramer realization into the
desired symmetric chain.

### 2.3 Exact obstruction to rational symmetric entries

Already

\[
\widehat O_{4,2}(t)=1+\frac{15}{8}t+\frac14t^2.
\]

If a rational symmetric `2 x 2` matrix `[[a,b],[b,c]]` realized this
polynomial, then

\[
(a-c)^2+4b^2
=\left(\frac{15}{8}\right)^2-4\left(\frac14\right)
=\frac{161}{64}.
\]

The left side is a sum of two rational squares.  But `161=7*23`, and both
`7` and `23` are primes congruent to `3 mod 4`, each occurring to an odd
power.  Therefore `161/64` is not a sum of two rational squares.  For
example, clearing denominators in `x^2+y^2=161 z^2` and reducing modulo `7`
forces `7|x,y`; dividing by `49` gives a valuation contradiction.

So no rational symmetric realization exists for this row.  Padding by zero
rows does not help: a rational PSD rank-two matrix has a rational Gram
factor `B`, and its two nonzero eigenvalues are those of the rational
symmetric `2 x 2` matrix `B^T B`.  Irrational symmetric realizations do
exist; the obstruction is specifically to the integer/rational-entry model.

---

## 3. Requested `j=2` odd-step kill-test

For `m>=1`, define

\[
r_m(t)=
\frac{\widehat O_{2m+1,2}(t)-\widehat O_{2m,2}(t)}
{t\widehat O_{2m,2}(t)}.
\]

The proposed odd-step rank-one update is equivalent to

\[
r_m(t)=v_m^T(I+tM_{2m})^{-1}v_m.
\]

Thus `r_m` must be a Stieltjes transform.  Its Taylor coefficients give the
moments

\[
\mu_k=(-1)^k[t^k]r_m(t)=v_m^TM_{2m}^kv_m.
\]

Exact moment Lanczos, with no roots, gives

\[
r_m(t)=
\frac{w_m}{1+\alpha_0t-
\dfrac{\beta_1t^2}{1+\alpha_1t-
\dfrac{\beta_2t^2}{1+\alpha_2t-\cdots}}}.
\]

The rational tridiagonal matrix `T_{2m}` has diagonal `alpha_i`, unit
superdiagonal, and subdiagonal `beta_i`.  It is diagonally similar to the
symmetric Jacobi matrix with off-diagonal `sqrt(beta_i)`.  Its nonzero block
has size `m=\deg\widehat O_{2m,2}`; it may be padded by zeros to size `2m`.
Exact checks give

\[
\det(I+tT_{2m})=\widehat O_{2m,2}(t).
\]

Consequently, whenever all `alpha_i,beta_i,w_m` are positive,

\[
M_{2m+1}=J_m+\sqrt{w_m}\,e_0e_0^T
\]

is an explicit symmetric PSD rank-one update, where `J_m` is the symmetric
form of `T_{2m}`.  Indeed `J_m` is symmetric and its determinant is
`widehat O_{2m,2}`, whose coefficients are positive; therefore all its real
eigenvalues are positive, so `J_m` is PSD.

### 3.1 Exact ratios for `m=1,...,5`

Up to a nonzero scalar the computed rational functions are:

\[
r_1=\frac{4}{6+3t},
\]

\[
r_2=\frac{86+80t}{120+225t+30t^2},
\]

\[
r_3=\frac{3620+8372t+2520t^2}
{5040+16681t+10430t^2+630t^3},
\]

\[
r_4=\frac{260650+976488t+816648t^2+120960t^3}
{362880+1722321t+2143134t^2+661500t^3+22680t^4},
\]

\[
r_5=\frac{28671512+148599308t+217815840t^2+94506720t^3+8316000t^4}
{39916800+246798321t+471163902t^2+309047970t^3+56548800t^4+1247400t^5}.
\]

### 3.2 Canonical Lanczos data

For `m=1`:

\[
w_1=\frac23,\qquad \alpha_0=\frac12.
\]

For `m=2`:

\[
w_2=\frac{43}{60},\qquad
\alpha=\left(\frac{325}{344},\frac{40}{43}\right),\qquad
\beta_1=\frac{4651}{7396}.
\]

For `m=3`:

\[
w_3=\frac{181}{252},
\]

\[
\alpha=\left(
\frac{129931}{130320},
\frac{169445152769}{99529848365},
\frac{67114080}{109977733}
\right),
\]

\[
\beta_1=\frac{109977733}{117939600},\qquad
\beta_2=\frac{4146123129706890}{12095101755819289}.
\]

All entries are positive.  At `m=4` and `m=5` the exact output is large;
representative entries and the growth statistics are:

| `m` | representative `alpha` | max digits in an `alpha` denominator | max digits in a `beta` denominator |
|---:|---|---:|---:|
| 4 | `alpha_0=1050833369/1050940800`; `alpha_1=16707176617557906192/8685652905705081475` | 35 | 42 |
| 5 | `alpha_0=619303119349/619304659200`; `alpha_1=425007086779366944124205951/217409566404568167596735214` | 73 | 83 |

The complete exact `alpha,beta` lists are printed by
`py -3.14 proof-psd-realization/kill_test.py`.

The `m=4` denominator

\[
8685652905705081475
=5^2\cdot13\cdot401\cdot147163\cdot452872661
\]

contains the prime `452872661`.  A denominator at `m=5` factors as

\[
728716681025945109053362151353865476905844274186973
=3\cdot31\cdot41\cdot139\cdot2346336497886019
\cdot585985418520645336323397927481,
\]

with `2346336497886019` prime.  These are not the small factorial/binomial
denominators suggested by the original recurrence.

**Kill-test conclusion.** The finite canonical matrices exist and are PSD
through `m=5`, but there is no recognizable low-complexity closed form.
The rapid appearance of large irregular factors by `m=4,5` kills the
proposed simple structured-matrix formula.  It does not kill all possible
matrix models: an all-index proof could still establish positivity of the
Lanczos coefficients by a non-closed-form argument.  That positivity is
equivalent to the odd-step interlacing statement once the denominator row
is real-rooted.

---

## 4. Analogous exact even-step bordering test

For the even step, write

\[
R_m(t)=\frac{\widehat O_{2m,2}(t)}{\widehat O_{2m-1,2}(t)}
=1+\gamma_m t+O(t^2).
\]

A PSD bordering requires

\[
R_m(t)=1+t\alpha_m-t^2s_m(t),\qquad
s_m(t)=u_m^T(I+tM_{2m-1})^{-1}u_m.
\]

Cancellation of the pole at zero forces `alpha_m=gamma_m`, and then

\[
s_m(t)=\frac{1+\alpha_mt-R_m(t)}{t^2}.
\]

Exact Lanczos on `s_m` again gives positive data through `m=5`:

| `m` | border `alpha_m` | `u_m^Tu_m` | max diagonal denominator digits |
|---:|---:|---:|---:|
| 1 | `1/2` | `0` | 0 |
| 2 | `17/24` | `83/144` | 1 |
| 3 | `517/720` | `61187/86400` | 7 |
| 4 | `28961/40320` | `48626287/67737600` | 24 |
| 5 | `2606501/3628800` | `945831411911/1316818944000` | 56 |

For example,

\[
s_2(t)=\frac{83}{144+168t},\qquad T_3=[7/6],
\]

so

\[
M_4=
\begin{pmatrix}
7/6&\sqrt{83}/12\\
\sqrt{83}/12&17/24
\end{pmatrix}
\]

has determinant `1+15/8 t+1/4 t^2`, as required.

Thus both odd and even finite steps admit coefficient-only symmetric PSD
constructions through `n=11` for `j=2`.  The even-step denominators also
explode, and no all-index positivity proof is supplied by the differential
recurrence.

---

## 5. No standard Jacobi/Favard chain in `n`

A single Jacobi/Favard chain would imply a three-term determinant recurrence

\[
D_{q+1}(t)=(1+a_qt)D_q(t)-b_qt^2D_{q-1}(t).
\]

The degree-grouped `j=2` rows already fail this.  For

\[
D_0=1,\quad D_1=\widehat O_{2,2},\quad
D_2=\widehat O_{4,2},\quad D_3=\widehat O_{6,2},
\]

the first two steps force

\[
a_1=\frac{11}{8},\quad b_1=\frac7{16},\qquad
a_2=\frac{1033}{720},\quad b_2=\frac{1003}{1152}.
\]

The resulting predicted `t^3` coefficient of `D_3` is
`-883/11520`, whereas the true coefficient is `1/8=1440/11520`.  The
residual is

\[
-\frac{2323}{11520}\ne0.
\]

`jacobi_chain_probe.py` checks every residue class `n mod j` for
`j=2,...,12` and finds the same obstruction (the first failure is at degree
three).  Therefore the differential recurrence in `n` does not provide the
orthogonality/Favard coefficients needed for a standard Jacobi matrix
generating the whole row family.  The tridiagonal matrices in Sections 3-4
are separate spectral realizations of individual rational functions, not
one Favard chain in `n`.

---

## 6. Exact computational audit

`verify.py` checks:

* 140 Cramer determinant identities for `j=2,...,8`, `1<=n<=20`;
* 133 all-index-form rank-one **nonsymmetric** Cramer updates in that range;
* 98 coefficient-derived Hankel/Jacobi identities for `j=2,...,8`,
  `1<=n<=18`;
* the requested odd-step kill-test for `m=1,...,5`;
* the analogous even-step bordering test for `m=1,...,5`;
* the Jordan, rational-square, and Favard obstructions above.

All checks pass.  The Hankel/Jacobi rows have positive leading Hankel
minors, zero full Hankel determinant, and positive Lanczos coefficients in
this finite range.  This is evidence only; no all-index positivity or
interlacing theorem is claimed.

---

## 7. Bottom line

* The Cramer route is a genuine new recurrence-based determinant identity.
* Its natural inverse matrix has a simple all-index rank-one update, but the
  update is nonsymmetric and the matrix has a nontrivial Jordan block.  It
  cannot be gauged into the proposed PSD chain.
* The `j=2` Stieltjes/Lanczos kill-test constructs exact finite PSD updates
  and borderings without roots.  The canonical entries are positive but
  become numerically and combinatorially unstructured by `m=4,5`.
* Individual PSD realizability is equivalent to real-rootedness; the full
  rank-one/bordering chain is equivalent to the conjecture.  The exact
  obstruction is therefore not a lack of a determinant identity, but the
  absence of an independent all-index positivity proof for the canonical
  Lanczos data.

## Files

* `explore.py` - coefficient power sums, complete Hankel test, exact Jacobi/Lanczos.
* `cramer_probe.py` - Cramer matrices and the nonsymmetric rank-one update.
* `kill_test.py` - requested odd-step `r_m` J-fraction kill-test.
* `even_border_probe.py` - exact even-step Stieltjes/bordering test.
* `jacobi_chain_probe.py` - exact failure of a three-term Favard chain.
* `verify.py` - consolidated exact regression checks.
