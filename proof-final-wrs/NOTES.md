# Final WRS proof track for the original ordered big-block family

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-wrs` was created or edited.

## Status

No unconditional all-index proof of WRS was obtained. The strongest new exact results are:

1. a positive recurrence for the lag-`j` adjacent coefficient minor;
2. an exact formal recurrent example in which the lagged row is compatible with the main part, the final row is real-rooted, and both relevant coefficient minors are positive, but WRS still fails;
3. an exact actual-family obstruction to the tempting proof that treats all older rows as remaining in one global compatible window.

These results show that neither the first-block recurrence plus local compatibility nor the natural coefficient-ratio refinement is sufficient. The missing input is a genuine lower bound on the main-part residues (equivalently, quantitative Bezoutian domination), using the full history of the actual rows.

## 1. Family and target

Fix `j >= 2`. The original series is

\[
B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
+t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right),
\qquad F_j(z,t)=\frac1{1-B_j(z,t)},
\]
\[
O_n=n![z^n]F_j(z,t),\qquad O_0=1.
\]

Put `T=t d/dt`, `Q(p)=p+(t-1)p'`, and

\[
A_n=R_n(O_{n-1})+S_n,\qquad
R_n(p)=np+2Tp,\qquad
S_n=\sum_{r=2}^{j-1}\binom nr O_{n-r},
\]
\[
c=\binom{n-1}{j-1}.
\]

The exact differential recurrence is

\[
O_n=A_n+c\,tQ(O_{n-j}).
\tag{1}
\]

For `q=O_{n-j}` and every root `mu` of `q`, WRS is

\[
-\frac{A_n(\mu)}{q'(\mu)}>c\mu(\mu-1).
\tag{WRS}
\]

Equivalently, since `tQ(q)(mu)=mu(mu-1)q'(mu)`,

\[
\frac{O_n(\mu)}{q'(\mu)}<0.
\tag{2}
\]

Thus WRS is exactly negativity of every residue of `O_n/q` at a root of `q`.

## 2. What the first-block recurrence proves, and what it does not

The actual first-block recurrence is

\[
O_n=\sum_{r=1}^{j-1}\binom nr O_{n-r}
+t\sum_{r=j}^{n}\binom nr O_{n-r}.
\tag{3}
\]

Let `mu` be a root of `q=O_{n-j}`. Suppose the recent rows `O_{n-j},...,O_{n-1}` lie in a strict compatible window. Then:

* for `1 <= r <= j-1`, the row `O_{n-r}` is newer than `q`, so
  `sgn O_{n-r}(mu)=-sgn q'(mu)`;
* the `r=j` term vanishes;
* for `r>j`, the factor `t` contributes `sgn mu=-1`. Hence the desired sign
  `-sgn q'(mu)` follows exactly when
  `sgn O_{n-r}(mu)=sgn q'(mu)`.

Therefore (3) proves (2) if the lagged row is compatible with **all older rows occurring in (3)**. Local compatibility of the recent window does not supply those older-row signs.

This older-row compatibility cannot be replaced by a global-window invariant: it is false for the actual family. Already for `j=2`,

\[
O_2=2+t,\qquad O_2(-2)=0,
\]
but
\[
O_3(-2)=-8,\quad O_4(-2)=-42,\quad O_5(-2)=-62,\quad
O_6(-2)=1194>0.
\]

Thus the sign at an old root changes after several newer rows. A termwise proof from (3) needs a different way to control the entire older tail in the second sum.

## 3. Positive lag-`j` coefficient-minor identity

Write

\[
a_{m,k}=[t^k]O_m,\qquad
\Delta^{(d)}_{m,k}
=a_{m,k+1}a_{m-d,k}-a_{m,k}a_{m-d,k+1}.
\]

The coefficient form of (1) is

\[
a_{n,k}=(n+2k)a_{n-1,k}
+\sum_{r=2}^{j-1}\binom nr a_{n-r,k}
+c\,k(a_{n-j,k-1}-a_{n-j,k}).
\tag{4}
\]

Put `q_k=a_{n-j,k}` and `u_k=a_{n-1,k}`. Direct substitution into
`Delta^(j)_{n,k}=a_{n,k+1}q_k-a_{n,k}q_{k+1}` gives the exact identity

\[
\begin{aligned}
\Delta^{(j)}_{n,k}={}&
(n+2k)\Delta^{(j-1)}_{n-1,k}
+2u_{k+1}q_k\\
&+\sum_{r=2}^{j-1}\binom nr
   \Delta^{(j-r)}_{n-r,k}\\
&+c\bigl((k+1)q_k^2-q_kq_{k+1}
             -kq_{k-1}q_{k+1}\bigr).
\end{aligned}
\tag{5}
\]

The first-block recurrence for `u=O_{n-1}` has the term

\[
\binom{n-1}{j-1}O_{n-j}=c q,
\]

so coefficientwise

\[
u_{k+1}\ge c q_{k+1}.
\tag{6}
\]

Using (6), the last two positive-looking contributions in (5) are bounded below by

\[
\begin{aligned}
2u_{k+1}q_k
+c\bigl((k+1)q_k^2-q_kq_{k+1}-kq_{k-1}q_{k+1}\bigr)
&\ge c\bigl((k+1)q_k^2+q_kq_{k+1}
-kq_{k-1}q_{k+1}\bigr)\\
&=c\Bigl((k+1)(q_k^2-q_{k-1}q_{k+1})
+q_{k+1}(q_k+q_{k-1})\Bigr)>0,
\end{aligned}
\tag{7}
\]

where the last strict positivity uses ordinary log-concavity of the coefficients of the real-rooted polynomial `q`.

Consequently, if the lower-lag minors in (5) are positive, then
`Delta^(j)_{n,k}>0`. This is a genuine positive coefficient recurrence. It is not WRS: positivity of these minors gives no lower bound on the residues in (2).

## 4. Exact obstruction to the coefficient-minor refinement

The following is formal recurrent data, not a claimed counterexample to the actual family. It satisfies (1), the threshold lower bound (6), recent compatibility, transformed-window compatibility, and has a real-rooted final row.

Take `j=2`, `n=5`, `c=4`, and

\[
q=t+1,\qquad
p=(t+\tfrac{23}{5})(t+\tfrac{24}{25})
=t^2+\frac{139}{25}t+\frac{552}{125}.
\]

The root `-1` of `q` lies strictly between the roots `-23/5` and `-24/25` of `p`. Moreover

\[
p_0=\frac{552}{125}>4q_0,\qquad
p_1=\frac{139}{25}>4q_1,
\]

so (6) holds. The main part is

\[
A_5=R_5(p)=5p+2tp'
=9t^2+\frac{973}{25}t+\frac{552}{25}.
\]

At `mu=-1`,

\[
A_5(-1)=-\frac{196}{25}<0<A_5(0),
\]

and the leading coefficient is positive. Hence `A_5` has one root on each side of `-1`; thus `q` is also compatible with `A_5`, not merely with `p`.

The final recurrent row is

\[
r=A_5+4tQ(q)
=17t^2+\frac{973}{25}t+\frac{552}{25}.
\]

It is real-rooted, since

\[
\operatorname{disc}_t(r)
=\frac{8329}{625}>0,
\]

and its coefficients are positive. Both adjacent coefficient minors are positive:

\[
p_1q_0-p_0q_1=\frac{143}{125}>0,
\qquad
r_1q_0-r_0q_1=\frac{421}{25}>0.
\]

Nevertheless,

\[
r(-1)=\frac4{25}>0,
\]

so the required sign in (2) fails. Explicitly,

\[
-\frac{A_5(-1)}{q'(-1)}=\frac{196}{25}<8
=4(-1)(-2),
\]

which is exactly WRS failure.

Thus even all of the following are insufficient:

* compatibility of `q` with `p`;
* compatibility of `q` with `A_5`;
* the coefficient lower bound (6);
* positivity of the lower-lag and lag-`j` coefficient minors;
* real-rootedness of the final row.

The example satisfies the differential coefficient recurrence (4), but not the full first-block history of the actual `O_m`. That distinction is essential: the missing proof must use more than the structural recurrence and these coefficient minors.

## 5. Precise remaining gap

For actual rows, WRS is equivalent to either of the following root-direction statements:

\[
-\frac{A_n(\mu)}{q'(\mu)}>c\mu(\mu-1)
\quad\text{for every root }\mu\text{ of }q,
\tag{8}
\]

or

\[
M(q,A_n)+c\,v(q)v(q)^T\succ c\,M(q,h_q),
\qquad h_q=t(t-1)q',
\tag{9}
\]

where `M` is the root-positive Bezoutian and `v(q)` is the coefficient vector of `q`.

Ordinary compatibility of `q` with `A_n` gives only positivity of the left side of (8), not the strict quantitative domination by `c mu(mu-1)`. The positive minor identity (5) controls coefficient ratios, while (8) controls residues; the exact example in Section 4 shows that these are fundamentally different requirements.

An all-index proof therefore still needs a target-specific lower bound on the main-part residues, or an equivalent quantitative interlacing/Bezoutian invariant, derived from the full actual coefficient history. No such bound was obtained in this track.
