# Consecutive interlacing: sign condition for O_{n,j} interl O_{n+1,j}

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-final-consecutive` was created/edited.

## Status

No unconditional all-index proof of consecutive interlacing was obtained.

What is proved here:

1. A necessary and sufficient sign condition for consecutive interlacing of
   two real-rooted rows with positive leading coefficients (Lemma 2).
2. An exact recurrence evaluation of that condition at roots of `O_{n,j}`
   (Theorem 4), giving a two-row/consecutive residue inequality.
3. Equivalent Wronskian / Bezout formulations of the same condition
   (Section 5).
4. An exact reduction of the open gap to a quantitative domination inequality
   at roots of the current row (Section 6).

What is not proved: the domination inequality for all `j>=2, n>=0`.

Sanity check: `verify_consecutive.py` verifies (R1), (E1), (FBroot), (Wroot),
and NSC<=>interlacing for `j<=6, n<=13` (SymPy, high-precision root evaluation;
not a proof).

## 1. Setup

Fix `j>=2`. Let

```text
B_j(z,t)=sum_{r=1}^{j-1} z^r/r! + t (e^z - sum_{r=0}^{j-1} z^r/r!),
F_j(z,t)=1/(1-B_j(z,t)),   O_n(t)=O_{n,j}(t)=n![z^n]F_j(z,t),   O_0=1.
```

Degrees and coefficients (proved from the EGF / first-block recurrence):

```text
deg O_n = d_n := floor(n/j),
O_n(0)=n![z^n]1/(1-sum_{r=1}^{j-1} z^r/r!) > 0,
lc(O_n) > 0,
all coefficients of O_n are nonnegative integers.
```

Thus every real root of a real-rooted `O_n` is <=0, and if `d_n>=1` then the
roots are <0 because `O_n(0)>0`.

The universal differential recurrence is (README / proof-weighted)

```text
O_n = R_n(O_{n-1}) + S_n + C_{n,j} t Q(O_{n-j}),                 (R)
R_n(p)=n p + 2T p,  T=t d/dt,
S_n = sum_{r=2}^{j-1} binom(n,r) O_{n-r},   (empty for j=2)
C_{n,j}=binom(n-1,j-1),
Q(p)=p+(t-1)p'.
```

The first-block recurrence is

```text
O_n = sum_{r=1}^{j-1} binom(n,r) O_{n-r} + t sum_{r=j}^{n} binom(n,r) O_{n-r}.  (FB)
```

We also use the identity, for a root `mu` of a polynomial `q`,

```text
t Q(q)(mu) = mu q(mu) + mu(mu-1)q'(mu) = mu(mu-1)q'(mu).         (Qroot)
```

## 2. Interlacing convention

Let `P,Q` be real-rooted polynomials with positive leading coefficients and
real roots

```text
P: mu_1 < ... < mu_d,      d=deg P,
Q: nu_1 < ... < nu_e,      e=deg Q,   e in {d,d+1}.
```

We say `P interlaces Q` (written `P interl Q`) if

```text
e==d+1:  nu_1 <= mu_1 <= nu_2 <= ... <= mu_d <= nu_{d+1},
e==d:    mu_1 <= nu_1 <= mu_2 <= ... <= mu_d <= nu_d.
```

This is the convention used by `interlace.py`. For strict interlacing, all
inequalities are strict and all listed roots are simple.

For the ordered big-block family, consecutive interlacing means
`O_n interl O_{n+1}` for all `n>=0`.

## 3. Necessary and sufficient sign condition

### Lemma 1 (root-sign of a real-rooted polynomial)

Let `P(t)=c prod_{i=1}^d (t-mu_i)`, `c>0`, `mu_1<...<mu_d` simple. Then

```text
sgn(P'(mu_i)) = (-1)^(d-i).                                       (S1)
```

If `x in (mu_i,mu_{i+1})`, with `mu_0:=-infinity`, `mu_{d+1}:=+infinity`,
then

```text
sgn(P(x)) = (-1)^(d-i).                                           (S2)
```

Proof. `P'(mu_i)=c prod_{k!=i}(mu_i-mu_k)`. The factors with `k<i` are
positive, `d-i` factors with `k>i` are negative. For (S2), exactly `d-i`
linear factors are negative on `(mu_i,mu_{i+1})`.  QED.

### Lemma 2 (consecutive interlacing sign condition)

Let `P,Q` be real-rooted with positive leading coefficients, `deg P=d>=0`,
`deg Q=e in {d,d+1}`. Assume `P` has simple roots `mu_1<...<mu_d`.

#### Proposition 2a (degree e=d+1)

The following are equivalent:

1. `P interl Q` strictly: `nu_1<mu_1<nu_2<...<mu_d<nu_{d+1}`.
2. `Q(mu_i)/P'(mu_i) < 0` for `i=1,...,d`.  (C)
3. `sgn Q(mu_i)=(-1)^(d+1-i)` for `i=1,...,d`.

#### Proposition 2b (same degree e=d)

The following are equivalent:

1. `P interl Q` strictly in the repo convention: `mu_1<nu_1<mu_2<...<mu_d<nu_d`.
2. `Q(mu_i)/P'(mu_i) < 0` for `i=1,...,d`.  (C)
3. `sgn Q(mu_i)=(-1)^(d+1-i)` for `i=1,...,d`.

We prove both cases.

#### Proposition 2a (degree +1)

Assume `e=d+1`. The following are equivalent:

1. `P interl Q` strictly: `nu_1<mu_1<nu_2<...<mu_d<nu_{d+1}`.
2. `Q(mu_i)/P'(mu_i) < 0` for `i=1,...,d`.
3. `sgn Q(mu_i)=(-1)^(d+1-i)` for `i=1,...,d`.

Proof. (1)=>(3): in `(nu_i,nu_{i+1})` the sign of `Q` is `(-1)^(d+1-i)` by
Lemma 1(S2), and `mu_i` lies in this interval. (3)=>(1): `Q(mu_d)` has sign
`-1` while `Q(+infinity)>0`, so a root `nu_{d+1}>mu_d` exists;
`sgn Q(mu_i)` alternates and `Q(-infinity)` has sign `(-1)^(d+1)`, giving a
root in `(-infinity,mu_1)` and one in each `(mu_i,mu_{i+1})`; total `d+1`
simple roots in the interlacing positions. (2)<=> (3) by Lemma 1(S1):
`sgn(Q(mu_i)/P'(mu_i)) = sgn Q(mu_i)*(-1)^(d-i)`, so it is negative iff
`sgn Q(mu_i)=(-1)^(d+1-i)`.  QED.

#### Proposition 2b (same degree)

Assume `e=d`. The following are equivalent:

1. `P interl Q` strictly in the repo convention: `mu_1<nu_1<mu_2<...<mu_d<nu_d`.
2. `Q(mu_i)/P'(mu_i) < 0` for `i=1,...,d`.
3. `sgn Q(mu_i)=(-1)^(d+1-i)` for `i=1,...,d`.

Proof. Identical sign bookkeeping. (1)=>(3): `mu_i in (nu_{i-1},nu_i)` in the
relevant indexing gives `d-i+1` roots of `Q` to the right, hence sign
`(-1)^(d+1-i)`. (3)=>(1): `Q(mu_d)<0` and `Q(+infinity)>0` gives
`nu_d>mu_d`; alternation gives `nu_i in (mu_i,mu_{i+1})`; count gives all
`d` roots. (2)<=> (3) is the same as 2a.  QED.

### Remark 3 (endpoint at 0 for the O-family)

For `O_n`, all real roots are <0 and `O_n(0)>0`, `lc(O_n)>0`. Therefore, once
all roots of `O_n` are simple real negative, the endpoint conditions are
automatic in the direction needed: `O_{n+1}(0)>0` is already known. The
interlacing of `O_n` by `O_{n+1}` is equivalent to (C) together with
real-rootedness of `O_{n+1}$ and the degree condition `e-d in {0,1}`.

## 4. Recurrence evaluation at roots of O_n

Let `mu` be a simple root of `P=O_n`, `d=d_n=deg P`. By (R) at index `n+1$,

```text
O_{n+1} = R_{n+1}(O_n) + S_{n+1} + C_{n+1,j} t Q(O_{n+1-j}),      (R1)
R_{n+1}(p)=(n+1)p+2Tp,  C_{n+1,j}=binom(n,j-1),
S_{n+1}=sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}.
```

Evaluating at `mu$, where $O_n(mu)=0$:

```text
O_{n+1}(mu)=2mu O_n'(mu) + S_{n+1}(mu) + C_{n+1,j} (tQ(O_{n+1-j}))(mu). (E1)
```

Dividing by `O_n'(mu)$ gives the consecutive residue

```text
Rho_{n+1}(mu) := O_{n+1}(mu)/O_n'(mu)
 = 2mu
 + sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}(mu)/O_n'(mu)
 + C_{n+1,j} (tQ(O_{n+1-j}))(mu)/O_n'(mu).                         (Rho)
```

### Theorem 4 (necessary and sufficient consecutive sign condition)

Assume $O_n$ is real-rooted with simple roots $mu_1<...<mu_{d_n}<0$, and
$O_{n+1}$ has positive leading coefficient and degree
$d_{n+1}in{d_n,d_n+1}$. Then

```text
O_n interl O_{n+1}
  <=>  Rho_{n+1}(mu_i) < 0   for all i=1,...,d_n.                  (NSC)
```

In particular, (NSC) implies that $O_{n+1}$ is real-rooted.

Equivalently,

```text
2mu_i + sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}(mu_i)/O_n'(mu_i)
      + binom(n,j-1) (tQ(O_{n+1-j}))(mu_i)/O_n'(mu_i) < 0.          (NSC')
```

Proof. This is Proposition 2a/2b applied to $P=O_n$, $Q=O_{n+1}$, combined
with (Rho). The endpoint conditions are automatic because $O_{n+1}(0)>0$ and
$O_{n+1}$ has positive leading coefficient.  QED.

### Corollary 4.1 (two-row domination form)

Since $mu_i<0$, (NSC') is equivalent to

```text
sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}(mu_i)/O_n'(mu_i)
 + binom(n,j-1) (tQ(O_{n+1-j}))(mu_i)/O_n'(mu_i)
 < 2|mu_i|.                                                        (DOM)
```

The term $2mu_i= -2|mu_i|$ is the favorable derivative contribution from the
consecutive row itself. All other terms are a domination remainder.

### Corollary 4.2 (j=2 specialization)

For $j=2$, $S_{n+1}=0$, $C_{n+1,2}=n$, and (NSC') is

```text
2mu_i + n (tQ(O_{n-1}))(mu_i)/O_n'(mu_i) < 0.                      (NSC-j2)
```

This is a direct consecutive analogue of the lag-2 WRS inequality, but the
$tQ(O_{n-1})$ term is evaluated at a root of $O_n$, not at a root of
$O_{n-1}$.

## 5. Wronskian and Bezout equivalents

Let $P=O_n, Q=O_{n+1}$, degrees $d<=e$, positive leading coefficients.

### 5.1 Wronskian

Define

```text
W(P,Q)(t) = P(t)Q'(t) - P'(t)Q(t).
```

At a simple root $mu_i$ of $P$, $W(P,Q)(mu_i) = -P'(mu_i)Q(mu_i)$. Therefore

```text
Q(mu_i)/P'(mu_i) < 0  <=>  W(P,Q)(mu_i) > 0.                     (Wroot)
```

By Proposition 2a/2b and (Wroot):

> Under the hypotheses of Theorem 4, $O_n interl O_{n+1}$ iff $W(O_n,O_{n+1})$
> is positive at every root of $O_n$.

This is a two-row sign invariant on the root set of the lower row. It is not
the same as global positivity of the Wronskian polynomial, which is stronger
and false in general for same-degree interlacing unless additional conditions
hold; the root-set condition is the exact one.

### 5.2 Bezout / residue matrix

Write the partial fraction residue

```text
Q/P = polynomial + sum_{i=1}^d rho_i/(t-mu_i),   rho_i=Q(mu_i)/P'(mu_i),
```

when $P$ has simple roots and $e<=d+1$ (the polynomial part is constant when
$e==d$, linear absent in monic-ratio form; in general for $e=d+1$ it is a
constant). The condition (C) is exactly negativity of all residues:

```text
rho_i < 0, i=1,...,d.                                             (RES)
```

Equivalently, the Bezoutian $B(P,Q)$ in the appropriate basis gives the residue
quadratic form. For adjacent-degree real-rooted polynomials with positive
leading coefficients, strict interlacing is equivalent to positive definiteness
of the corresponding Bezout/Hermite form; in the residue normalization this is
the same as all $rho_i<0$ after fixing the sign convention. We do not re-prove
the full Bezoutian equivalence here; (RES) is the exact condition used.

### 5.3 First-block decomposition form

Using (FB) at index `n+1` and evaluating at a root `mu` of `O_n`, the `r=1`
term vanishes because `(n+1)O_n(mu)=0`. The terms `r=2,...,j-1` remain in the
unweighted sum, and the terms `r=j,...,n+1` remain in the `t`-weighted sum:

```text
O_{n+1}(mu) = sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}(mu)
             + mu sum_{r=j}^{n+1} binom(n+1,r) O_{n+1-r}(mu).      (FBroot)
```

The condition (NSC) becomes

```text
sum_{r=2}^{j-1} binom(n+1,r) O_{n+1-r}(mu_i)/O_n'(mu_i)
 + mu_i sum_{r=j}^{n+1} binom(n+1,r) O_{n+1-r}(mu_i)/O_n'(mu_i) < 0. (FB-NSC)
```

This is a two-row sign invariant in the form "recent-row residue sum plus
signed older-tail sum is negative".

## 6. Gap analysis

### 6.1 What the common-window hypothesis gives

Assume the recent window $O_{n+1-j},...,O_n$ is compatible in the standard
interlacing order. At a root $mu_i$ of $O_n$:

- For $2<=r<=j-1$, $O_{n+1-r}$ is older than $O_n$, so
  $sgn O_{n+1-r}(mu_i)=sgn O_n'(mu_i)$ by Lemma 1 of proof-weighted style;
  hence $O_{n+1-r}(mu_i)/O_n'(mu_i)>0$.
- The factor $mu_i<0$ makes the $t$-sum terms involve $mu_i/q'...$ signs.

The recent compatibility alone does not control (FB-NSC), because the $r=j$
term uses $O_{n+1-j}$ at a root of $O_n$ (not a root of $O_{n+1-j}$), and the
$r>=j+1$ older tail has the same global sign problem identified in
proof-final-wrs Section 2.

### 6.2 Precise remaining gap

**Gap CONSEC-SIGN (consecutive domination).** Prove (DOM) for all $j>=2$,
$n>=0$, and all roots $mu_i$ of $O_n$, assuming the real-rootedness/interlacing
hypothesis needed for the intermediate rows.

Equivalent forms:

1. Residue form: $Rho_{n+1}(mu_i)<0$.
2. Wronskian form: $W(O_n,O_{n+1})(mu_i)>0$.
3. First-block form: (FB-NSC).
4. Bezout/residue form: all residues of $O_{n+1}/O_n$ are negative.

The condition is necessary and sufficient (Theorem 4). Any all-index proof must
establish one of these equivalent inequalities; local compatibility and
coefficient-ratio positivity are insufficient in isolation (cf.
proof-final-wrs obstruction, which is for lag j but shows the same residue-vs-
coefficient-ratio phenomenon).

### 6.3 Relation to WRS

The WRS target in proof-final-wrs evaluates $O_n/q'$ at roots of $q=O_{n-j}$.
The consecutive condition here evaluates $O_{n+1}/O_n'$ at roots of $O_n$.
They are adjacent-row analogues with different lags. Neither implies the other
directly. The retracted tail-residue counterexample (proof-final-tail-residue)
does not give a WRS counterexample; its x*q' bug was corrected and the exact
certificate shows WRS holds at that root.

## 7. Files

- `NOTES.md` - this note.
- `verify_consecutive.py` - exact SymPy checks for (R1), (E1), (NSC), (FB-NSC),
  and the Wronskian/root-sign equivalence on small cases.
