# Final minors / coefficient-ratio proof track

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-final-minors` was created or edited.

## Status

**Not proved.** No all-index proof of real-rootedness, WRS, or row PF for the
actual ordered big-block family was obtained.

What this track establishes rigorously:

1. the lag-`j` coefficient-minor identity (5) from `proof-final-wrs/NOTES.md`
   is algebraically correct;
2. it yields a conditional positivity propagation for the lag-`j` minors against
   the base row `q=O_{n-j}`;
3. the same identity does **not** close into an all-lag / all-row-pair minor
   bootstrap by itself;
4. the exact formal obstruction of `proof-final-wrs/NOTES.md` Section 4 also
   blocks promoting the resulting coefficient-minor positivity into WRS or
   Bezoutian domination.

The remaining gap is still a quantitative main-part residue / Bezoutian
domination bound that uses actual-family history beyond the local compatible
window and the coefficient recurrences used here.

## 1. Setup

Fix `j>=2`. Let

$$
O_n(t)=\sum_{k>=0} a_{n,k} t^k,\qquad a_{n,k}=[t^k]O_n,
\qquad \deg O_n=e(n)=\lfloor n/j\rfloor.
$$

Use the exact differential recurrence

$$
O_n=A_n+c\,tQ(O_{n-j}),\qquad c=\binom{n-1}{j-1},
\tag{1}
$$

where

$$
A_n=R_n(O_{n-1})+S_n,\qquad R_n(p)=np+2Tp,\qquad S_n=\sum_{r=2}^{j-1}\binom nr O_{n-r},
\qquad Q(p)=p+(t-1)p'.
$$

The coefficient form is

$$
a_{n,k}=(n+2k)a_{n-1,k}
+\sum_{r=2}^{j-1}\binom nr a_{n-r,k}
+c\,k(a_{n-j,k-1}-a_{n-j,k}).
\tag{2}
$$

For rows `m>ell`, define the adjacent-column cross-row minor

$$
\Delta^{(m-\ell)}_{m,k}
:=a_{m,k+1}a_{\ell,k}-a_{m,k}a_{\ell,k+1}.
\tag{3}
$$

The sign convention is: `Delta^(d)>0` means the coefficient ratio of the newer
row dominates the older row at index `k`.

## 2. Verification of identity (5)

Let `q=O_{n-j}`, `u=O_{n-1}`, `q_k=a_{n-j,k}`, `u_k=a_{n-1,k}`. Substituting
(2) into `Delta^(j)_{n,k}=a_{n,k+1}q_k-a_{n,k}q_{k+1}` gives

$$
\begin{aligned}
a_{n,k+1}q_k-a_{n,k}q_{k+1}
={}&(n+2k)(u_{k+1}q_k-u_kq_{k+1})+2u_{k+1}q_k\\
&+\sum_{r=2}^{j-1}\binom nr
   (a_{n-r,k+1}q_k-a_{n-r,k}q_{k+1})\\
&+c\bigl((k+1)(q_k-q_{k+1})q_k-k(q_{k-1}-q_k)q_{k+1}\bigr).
\end{aligned}
$$

Since `u_{k+1}q_k-u_kq_{k+1}=Delta^(j-1)_{n-1,k}` and
`a_{n-r,k+1}q_k-a_{n-r,k}q_{k+1}=Delta^(j-r)_{n-r,k}`, this is

$$
\begin{aligned}
\Delta^{(j)}_{n,k}={}&
(n+2k)\Delta^{(j-1)}_{n-1,k}
+2u_{k+1}q_k\\
&+\sum_{r=2}^{j-1}\binom nr \Delta^{(j-r)}_{n-r,k}\\
&+c\bigl((k+1)q_k^2-q_kq_{k+1}-kq_{k-1}q_{k+1}\bigr).
\end{aligned}
\tag{4}
$$

This is identity (5) of `proof-final-wrs/NOTES.md`. The algebra was also checked by exact integer computation against the actual
family for `j<=6, n<=24` via `proof-final-minors/check_minors.py`.

## 3. Positive tail from actual first-block lower bound

The first-block recurrence for `u=O_{n-1}` contains the term

$$
\binom{n-1}{j-1}O_{n-j}=c q,
$$

hence coefficientwise

$$
u_{k+1}\ge c q_{k+1},\qquad k>=0.
\tag{5}
$$

Therefore the two tail contributions in (4) satisfy

$$
\begin{aligned}
2u_{k+1}q_k
+c\bigl((k+1)q_k^2-q_kq_{k+1}-kq_{k-1}q_{k+1}\bigr)
&\ge c\bigl((k+1)q_k^2+q_kq_{k+1}-kq_{k-1}q_{k+1}\bigr)\\
&=c\Bigl((k+1)(q_k^2-q_{k-1}q_{k+1})+q_{k+1}(q_k+q_{k-1})\Bigr).
\end{aligned}
\tag{6}
$$

If `q` has nonnegative log-concave coefficients, then
`q_k^2-q_{k-1}q_{k+1}>=0`, and the remaining term
`q_{k+1}(q_k+q_{k-1})` is strictly positive on the active coefficient range.
Thus (6) is strictly positive on the active range.

## 4. Conditional lag-j minor propagation

### Proposition 4.1

Assume, for a fixed actual root step `(n,j)`:

- `q=O_{n-j}` is coefficientwise log-concave;
- for every `r=1,...,j-1`, `Delta^(j-r)_{n-r,k}>=0` on all relevant `k`.

Then `Delta^(j)_{n,k}>0` on all relevant `k`.

**Proof.** In (4), the first term is `(n+2k)` times a nonnegative lower minor,
the middle sum is a nonnegative combination of nonnegative lower minors, and the
tail is strictly positive by (6).  `QED`

### What this gives

The proposition is a genuine positive recurrence for the **lag matching the
block threshold**. It is stronger than mere nonnegativity of the final tail: it
propagates positivity along the `j`-window to the same base row `q`.

It does **not** by itself prove the hypotheses for all `n,j`.

## 5. Why this does not bootstrap all lags / row PF

A tempting all-index invariant is

$$
\Delta^{(d)}_{n,k}>0\qquad\text{for all }n>=1,\ 1<=d<=n,\ k>=0.
\tag{H}
$$

(H) is much stronger than lag-`j` positivity and would be a powerful
coefficient-ratio chain. However:

1. **Identity (4) is lag-specific.** It uses the differential coefficient
   recurrence (2), which has the `q=O_{n-j}` tail. Replacing the base row by an
   arbitrary `b=O_{n-d}` gives
   $$
   \Delta^{(d)}_{n,k}=(n+2k)\Delta^{(d-1)}_{n-1,k}
   +2u_{k+1}b_k+\sum_{r=2}^{j-1}\binom nr\Delta^{(d-r)}_{n-r,k}
   +c\,T_{n,d,k},
   $$
   with
   $$
   T_{n,d,k}=(k+1)(q_k-q_{k+1})b_k-k(q_{k-1}-q_k)b_{k+1}.
   $$
   This is algebraically correct, but for `d!=j` the tail `T_{n,d,k}` is not
   positive in general. Exact checks on the actual family found many negative
   values, e.g. `j=2,n=5,d=1,k=0` gives `T=b_0(q_0-q_1)<0` cases and many more
   for `j=2,3`.

2. **Lag-`j` positivity does not imply WRS/RR.** `proof-final-wrs/NOTES.md`
   Section 4 gives a formal recurrent example satisfying (2), threshold lower
   bound (5), compatibility, positive lower/lag minors, and a real-rooted final
   row, but violating WRS. Therefore coefficient minors, even including the
   lag-`j` minor, cannot be promoted to WRS without additional actual-family
   input.

3. **Row PF is not implied.** Row PF (real-rootedness of a row) is equivalent to
   the finite Pólya frequency property of its coefficient sequence. Cross-row
   minors `Delta^(d)` are different objects from the Hankel/Toeplitz minors that
   control row PF.

Consequently identity (4) cannot be turned into an all-index proof by a local
coefficient-minor induction alone.

## 6. Residue / Bezoutian route: what identity (4) does and does not give

The target WRS inequality is

$$
-\frac{A_n(\mu)}{q'(\mu)}>c\mu(\mu-1),\qquad q(\mu)=0,
\tag{7}
$$

equivalently Bezoutian domination

$$
M(q,A_n)+c\,v(q)v(q)^T \succ c\,M(q,h_q),\qquad h_q=t(t-1)q'.
\tag{8}
$$

Identity (4) controls coefficient ratios. WRS controls residues. The Section 4
formal counterexample shows these are different invariants: it has positive
coefficient minors but fails (7).

A full-history route would need a lower bound on the main-part residue
`-A_n(mu)/q'(mu)` from the actual rows
`O_{n-1},...,O_{n-j+1}` and, through the full first-block tail, from older
rows. The ordinary compatible-window hypothesis only gives positivity of this
residue, not the quantitative margin `c mu(mu-1)`. No such margin was derived
here.

## 7. Exact checks performed

Exact integer/rational scripts (run ad hoc in this session, not committed as
separate files) checked:

- identity (4) for actual `O_{n,j}`, `j<=6`, `n<=24`: no mismatches;
- all cross-row adjacent-column minors `Delta^(d)_{n,k}>=0` for `j<=5`,
  `n<=29`: no negative values found (evidence only);
- the generalized tail `T_{n,d,k}` is negative for some actual `d!=j` cases
  (`j=2,3`), confirming it is not a general positivity certificate.

## 8. Precise gap

Write the first-block recurrence at a root `mu` of `q=O_{n-j}`:

$$
\frac{O_n(\mu)}{q'(\mu)}=
\sum_{r=1}^{j-1}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}
+\mu\sum_{r=j+1}^{n}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}.
\tag{9}
$$

Under the recent compatible window, Lemma 1 gives
`O_{n-r}(mu)/q'(mu)<0` for `1<=r<=j-1`, so the first sum is strictly negative.
The dangerous tail is the second sum: because `mu<0`, a term is harmful only
when `O_{n-r}(mu)/q'(mu)<0`, i.e. when an older row has already flipped sign at
an old root. `proof-final-wrs/NOTES.md` Section 2 shows such flips occur in the
actual family (`j=2, mu=-2, O_6(-2)>0`).

Thus the missing all-index input is exactly a tail-domination bound of the form

$$
\sum_{r=1}^{j-1}\binom nr \Bigl(-\frac{O_{n-r}(\mu)}{q'(\mu)}\Bigr)
\;>\;
\mu\sum_{r=j+1}^{n}\binom nr \frac{O_{n-r}(\mu)}{q'(\mu)}
\tag{10}
$$

for every old root `mu`, proved from actual-family history rather than from
local compatibility. Equivalently one needs the Bezoutian domination (8) with a
main-part lower bound that sees the full first-block tail.

Identity (4) gives a positive lag-`j` coefficient-minor certificate, but it does
not control the signed older-tail residues in (9)-(10). This is the precise gap.
