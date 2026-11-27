# Bezoutian proof track for the original ordered big-block family

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-bezout` was created or edited.

## Status

**Proved:** an exact all-index Bezoutian decomposition of the correct recurrence; positivity of the Bezoutian is equivalent to WRS/interlacing; and an exact compatible-window counterexample in which the final Bezoutian is indefinite.

**Open:** proving the required quantitative domination for the actual rows `O_{n,j}`. The recurrence and a common compatible window alone do not prove it. No all-index sum-of-previous-positive-kernels, Loewner-positivity, or sum-of-squares proof was obtained.

The computations in `exact_checks.py` are exact SymPy checks. They are evidence for the small original-family cases, not a substitute for the missing all-index inequality.

---

## 1. Conventions and the correct recurrence

Fix `j >= 2`. The original family is

\[
 B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
 +t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right),
 \qquad F_j(z,t)=\frac1{1-B_j(z,t)},
\]

so `B_j(0,t)=0`, `F_j(0,t)=1`, and

\[
 O_n=n![z^n]F_j(z,t),\qquad O_0=1.
\]

Write `T=t d/dt`, `Q(p)=p+(t-1)p'`, and

\[
 A_n=nO_{n-1}+2TO_{n-1}
       +\sum_{r=2}^{j-1}\binom nr O_{n-r},
 \qquad c=\binom{n-1}{j-1}.
\]

The correct recurrence is

\[
 \boxed{O_n=A_n+c\,tQ(O_{n-j}).} \tag{1}
\]

The binomial coefficients are `C(n,r)`, not `C(n-1,r-1)`.

For the root step put `q=O_{n-j}`, `e=deg q`. Then `deg O_n=e+1`.
At a root `mu` of `q`,

\[
 tQ(q)(\mu)=\mu(\mu-1)q'(\mu). \tag{2}
\]

Thus WRS is

\[
 -\frac{A_n(\mu)}{q'(\mu)}>c\mu(\mu-1). \tag{3}
\]

---

## 2. Root-positive Bezoutian

Use the symmetric kernel

\[
 \boxed{
 \mathcal B(q,p)(x,y)
 =\frac{q(y)p(x)-q(x)p(y)}{x-y}.} \tag{4}
\]

This is the negative of the convention
`(q(x)p(y)-q(y)p(x))/(x-y)`. The choice is made so that

\[
 \mathcal B(q,p)(\mu,\mu)=-q'(\mu)p(\mu)
 \quad\text{when }q(\mu)=0. \tag{5}
\]

Let `M(q,p)` be its coefficient matrix in the monomial basis. If
`deg q=e` and `deg p<=e+1`, this is an `(e+1) x (e+1)` symmetric matrix
(a smaller matrix if the degree drops).

The standard Bezoutian-interlacing criterion applies:

* if `q` has `e` distinct negative roots and `p` has positive leading
  coefficient and degree `e+1`, then `M(q,p)` is positive definite exactly
  when `p` has one root in each of
  `(-infinity,mu_1),...,(mu_e,0)`;
* equivalently, exactly when
  `sgn p(mu_i)=-sgn q'(mu_i)` for every root `mu_i` of `q`.

Indeed, those signs give the `e+1` distinct roots by IVT and exhaust the
degree. In a basis adapted to the roots of `q`, the Bezoutian is diagonal
with positive entries (the leading direction supplies the additional
positive entry). This is the usual Bezoutian form of the
real-rooted/interlacing criterion.

Consequently, for the actual root step,

\[
 \boxed{M(q,O_n)\succ0\quad\Longleftrightarrow\quad\text{WRS (3)}.} \tag{6}
\]

The forward implication also follows immediately from (5). Thus positive
definiteness is a matrix formulation of WRS, not a consequence of ordinary
compatibility.

---

## 3. All-index decomposition of the intrinsic term

Define

\[
 h_q(t)=t(t-1)q'(t).
\]

Since `Q(q)=q+(t-1)q'`, one has the exact polynomial identity

\[
 tQ(q)=tq+h_q. \tag{7}
\]

Bilinearity of (4) gives

\[
 \mathcal B(q,tQ(q))
 =\mathcal B(q,tq)+\mathcal B(q,h_q). \tag{8}
\]

The first term has a rank-one formula valid for every polynomial `q`:

\[
 \begin{aligned}
 \mathcal B(q,tq)(x,y)
 &=\frac{q(y)xq(x)-q(x)yq(y)}{x-y}\\
 &=q(x)q(y).
 \end{aligned} \tag{9}
\]

Therefore

\[
 \boxed{
 M(q,tQ(q))=v(q)v(q)^T+M(q,h_q),
 \qquad v(q)=(\text{coefficients of }q).} \tag{10}
\]

Combining (1), (8), and (10) gives the all-index recurrence decomposition

\[
 \boxed{
 M(q,O_n)=M(q,A_n)
       +c\,v(q)v(q)^T+c\,M(q,h_q).} \tag{11}
\]

At every root `mu` of `q`, the rank-one term vanishes and

\[
 \boxed{
 M(q,h_q)(\mu,\mu)
 =-\mu(\mu-1)q'(\mu)^2<0. \tag{12}
\]

So the intrinsic term is emphatically **not** a positive Bezout kernel.
Equation (11) exhibits the exact quantitative gap: `M(q,A_n)` must
dominate the negative root-direction contribution from `c M(q,h_q)`.

---

## 4. Exact obstruction: a compatible window is insufficient

Take the same formal data as the weighted-note obstruction:

\[
 j=2,\quad n=5,\quad c=4,\quad
 q=t+3,\quad p=(t+2)(t+4).
\]

The root `-3` of `q` lies between the roots `-4,-2` of `p`, so `(q,p)` is
an ordinary compatible pair with the required degrees. Put

\[
 A_5=5p+2tp'=9t^2+42t+40.
\]

Its roots are `-10/3` and `-4/3`, so `-3` also lies between them. Hence
`q` is compatible with `A_5`, and direct calculation gives

\[
 M(q,A_5)=
 \begin{pmatrix}86&27\\27&9\end{pmatrix},
 \qquad \det M(q,A_5)=45>0. \tag{13}
\]

Thus even the stronger statement “the main part has a positive definite
Bezoutian with `q`” holds.

Nevertheless,

\[
 h_q=t(t-1),\qquad tQ(q)=2t^2+2t,
\]

and

\[
 M(q,tQ(q))=
 \begin{pmatrix}6&6\\6&2\end{pmatrix},
 \qquad \det M(q,tQ(q))=-24<0. \tag{14}
\]

The final test row is

\[
 O_5^{\rm test}=A_5+4tQ(q)=17t^2+50t+40,
\]

whose discriminant is `-220`. Its Bezoutian is exactly

\[
 M(q,O_5^{\rm test})
 =M(q,A_5)+4M(q,tQ(q))
 =\begin{pmatrix}110&51\\51&17\end{pmatrix},
 \qquad \det=-731<0. \tag{15}
\]

At `mu=-3`, the root-direction values make the failure explicit:

\[
 M(q,A_5)(-3,-3)=5,\qquad
 4M(q,tQ(q))(-3,-3)=-48,
\]

so the final value is `-43`, whereas WRS would require a positive value.

**Obstruction.** No theorem whose hypotheses are only

1. the common compatible window,
2. the structural coefficient `c=C(n-1,j-1)`, and
3. compatibility (or even positive definiteness of the Bezoutian) of `q`
   with `A_n`

can imply `M(q,O_n) succ 0`. The displayed exact data satisfy those
hypotheses and violate the conclusion.

This also blocks the specific proposed route that treats
`M(q,tQ(q))` as another positive previous Bezout kernel: (14) is indefinite.
A universal Loewner congruence cannot repair that term by itself, because
inertia is invariant under congruence. A sum-of-squares proof would first
need an additional identity expressing the negative `M(q,h_q)` contribution
in terms of special properties of the actual preceding rows; ordinary
compatibility supplies no such identity.

---

## 5. What the exact original-family checks show

`exact_checks.py` verifies (10), (11), and the counterexample matrices. It
also constructs the original rows with the correct `C(n,r)` recurrence and
prints the leading principal minors of `M(O_{n-j},O_n)`.

For the checked ranges

\[
 j=2,\ n=2,\ldots,8;\qquad
 j=3,\ n=3,\ldots,9;\qquad
 j=4,\ n=4,\ldots,10,
\]

all those minors are positive. For example, at `j=2,n=5`,

\[
 q=6+7t,\qquad O_5=120+367t+54t^2,
\]

and

\[
 M(q,O_5)=
 \begin{pmatrix}1026&660\\660&770\end{pmatrix},
 \qquad \det=354420>0.
\]

These exact checks confirm the target and the sign convention, but they do
not prove the required domination in (11) for all `n,j`.

---

## 6. Precise remaining target

A Bezoutian induction must prove, for every actual root step,

\[
 \boxed{
 M(q,A_n)+c\,v(q)v(q)^T\succ c\,M(q,h_q),
 \qquad h_q=t(t-1)q'.} \tag{16}
\]

Equivalently it must prove WRS (3), or produce a target-specific identity
that dominates the negative intrinsic kernel using more than the common
compatible window. The recurrence alone, previous positive Bezoutians
without quantitative control, and a universal Loewner transformation of
the intrinsic term do not currently provide that missing inequality.
