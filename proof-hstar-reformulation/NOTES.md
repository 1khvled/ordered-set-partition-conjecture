# The `v=t/(1-t)` and `h*` reformulation: exact identities and the remaining obstruction

Date: 2026-11-28.  Scope: only `proof-hstar-reformulation/` in the named worktree.
No installs, builds, external agents, or models were used.

## Status in one paragraph

The proposed denominator, `R`/`W`, and `h*` generating-function identities are
correct exactly as stated: `O_{n,j}` has all roots in `(-infinity,0)` iff
`W_{n,j}` has all roots in `(1,infinity)`.  The apparent factor
`(1-y)^{d_n}` does not create a root at `y=1`; it cancels the pole of
`O_{n,j}(1/(1-y))`, and `W_{n,j}(1)` is the positive leading coefficient of
`O_{n,j}`.  The exact transformed universal recurrence is derived below.  The
stable denominator proves nonnegativity of the binomial basis coefficients,
but it does not prove real-rootedness.  Even nonnegativity, gamma-positivity,
and real-rootedness of `g_{n,j}(m)` below `-1` are jointly insufficient:
`g(m)=(m+2)^2` satisfies all three and gives `W(y)=y^2-3y+4`, with
discriminant `-7`.  No all-index root proof is claimed.  The exact missing
statement is recorded in Section 6.

## 1. Conventions and degrees

Fix `j>=2` and put

\[
 A_j(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad
 C_j(z)=e^z-1-A_j(z),\qquad
 N_j(z)=1-A_j(z),\qquad E(z)=2-e^z.
\]

Thus `C_j=N_j-E`, `C_j` has positive coefficients and starts in degree `j`,
and

\[
 F_j(z,t)=\frac1{1-A_j(z)-tC_j(z)}
        =\sum_{n\ge0}O_{n,j}(t)\frac{z^n}{n!}.
\]

Every big block uses at least `j` labels, so
`deg O_{n,j} <= d_n=floor(n/j)`.  Conversely, for every `0<=m<=d_n`, take `m`
disjoint blocks of size `j` and make all remaining labels singleton blocks;
ordering the blocks gives a positive coefficient of `t^m`.  Hence

\[
 \deg O_{n,j}=d_n,\qquad [t^m]O_{n,j}>0\quad(0\le m\le d_n).
\]

In particular `O_{n,j}(0)>0` and `O_{n,j}(1)>0`.

## 2. Exact `v` and `W` identities; direction audit

### 2.1 Algebra

Direct expansion gives

\[
 1-A_j-tC_j=(1+t)-(1-t)A_j-te^z.
\]

With `v=t/(1-t)`, multiplication by `1+v` gives

\[
 F_j(z,t)=\frac{1+v}{N_j(z)-v(e^z-2)}
        =\frac{1+v}{N_j(z)+vE(z)}.
\]

Define

\[
 R_{n,j}(v)=(1+v)^{d_n}O_{n,j}\!\left(\frac v{1+v}\right).
\]

Then, exactly,

\[
\begin{aligned}
 W_{n,j}(y)
 &=(-1)^{d_n}y^{d_n}R_{n,j}(-1/y)\\
 &=(1-y)^{d_n}O_{n,j}\!\left(\frac1{1-y}\right).
\end{aligned}
\]

For `y!=1`, the factor `(1-y)^{d_n}` is nonzero, so the zero correspondence is

\[
 t=\frac1{1-y},\qquad y=1-\frac1t.
\]

### 2.2 Direction and endpoints

On `t<0`,

\[
 \frac{dy}{dt}=\frac1{t^2}>0,
\]

so the map is increasing.  Its endpoint directions are

\[
 t\to-\infty:\ y\to1^+,
 \qquad t\to0^-:\ y\to+\infty.
\]

Thus `t<0` corresponds to `y>1`, not to `y<1`.  Explicit samples printed by
the audit are

\[
 t=-2\mapsto y=3/2,\quad t=-1\mapsto y=2,\quad
 t=-1/2\mapsto y=3,\quad t=1/2\mapsto y=-1.
\]

Similarly, `v=t/(1-t)` has derivative `1/(1-t)^2>0` for `t<0`, with
`t->-infinity` giving `v->-1+` and `t->0-` giving `v->0-`.  Hence negative
`t`-roots correspond to `v`-roots in `(-1,0)`, not to arbitrary negative
`v`-roots.

**Boundary at `y=1`; no forced root.**  Let `a_{d_n}>0` be the leading
coefficient of `O_{n,j}`.  Since `deg O=d_n`,

\[
 W_{n,j}(y)=(1-y)^{d_n}O_{n,j}\!\left(\frac1{1-y}\right)
 \longrightarrow a_{d_n}>0
 \quad(y\to1).
\]

The apparent factor `(1-y)^{d_n}` cancels the pole of the substituted `O`;
it does not give a root at `y=1`.  The analogous cancellation at `v=-1` gives
`R_{n,j}(-1)=(-1)^{d_n}a_{d_n}`, also nonzero.  Because `deg W=d_n`, the exact equivalence is
therefore literally

\[
\boxed{
 O_{n,j}\text{ has all roots in }(-\infty,0)
 \iff
 W_{n,j}\text{ has all roots in }(1,\infty).
}
\]

For `d_n=0` both conditions are vacuous.

## 3. Exact transformed universal recurrence for `R`

The original universal recurrence is

\[
\begin{aligned}
 O_n={}&nO_{n-1}+2t\partial_tO_{n-1}
 +\sum_{r=2}^{j-1}\binom nr O_{n-r}\\
 &+\binom{n-1}{j-1}t\bigl[O_{n-j}+(t-1)\partial_tO_{n-j}\bigr],
\end{aligned}
\]

where negative-index terms and the last term for `n<j` are omitted.  Write
`d_k=floor(k/j)` and

\[
 \sigma_{n,r}(v)=(1+v)^{d_n-d_{n-r}}.
\]

For `1<=r<j`, the exponent is either `0` or `1`; in particular
`sigma_{n,1}=1+v` exactly when `j|n`, and otherwise `sigma_{n,1}=1`.

For any polynomial `p`, set `p~(v)=p(v/(1+v))`.  The required chain rules are

\[
 t\partial_tp\longmapsto v(1+v)p~'(v),
 \qquad
 t(t-1)\partial_tp\longmapsto -vp~'(v),
 \qquad
 tp\longmapsto \frac v{1+v}p~(v).
\]

The first identity is `t d_t = v(1+v)d_v`; replacing it by `v d_v` would
give an incorrect recurrence.

Applying these identities and the normalization `(1+v)^{d_n}` gives the
following exact all-index formula:

\[
\boxed{
\begin{aligned}
 R_{n,j}(v)={}&
 \sigma_{n,1}\Bigl[
 nR_{n-1,j}+2v(1+v)R'_{n-1,j}
 -2d_{n-1}vR_{n-1,j}\Bigr]\\
&+\sum_{r=2}^{j-1}\binom nr\sigma_{n,r}R_{n-r,j}\\
&+\binom{n-1}{j-1}v
 \Bigl[d_nR_{n-j,j}-(1+v)R'_{n-j,j}\Bigr].
\end{aligned}}
\]

Here the sum only includes `n>=r`, and the last line only applies for `n>=j`.

For the last line, if `e=d_n-1=d_{n-j}` and
`R=(1+v)^eO_{n-j}(v/(1+v))`, then

\[
 (1+v)^{d_n}t\bigl[O_{n-j}+(t-1)\partial_tO_{n-j}\bigr]
 =v\bigl[d_nR-(1+v)R'\bigr].
\]

This proves the formula directly from the universal recurrence, independently
of any root argument.

### 3.1 The positive-root trap survives

For `j=2,n=5`, the lag row is `O_{3,2}=6+7t`.  With
`Q(p)=p+(t-1)p'`,

\[
 tQ(O_{3,2})=14t^2-t=t(14t-1).
\]

After multiplication by `(1+v)^{d_5}=(1+v)^2` and substitution `t=v/(1+v)`,

\[
 (1+v)^2tQ(O_{3,2})=13v^2-v=v(13v-1).
\]

Thus the nonzero positive root moves from `t=1/14` to `v=1/13`; it is not
removed by the `v` reformulation.  The final row can still have negative roots
because the complete recurrence combination is different from this isolated
lag term.  No `A+tB` split is used or duplicated here.

## 4. Exact `h*` identities and what the stable denominator proves

At `t=1/(1-y)`,

\[
 (1-y)(1-A_j)-C_j=E-yN_j,
\]

and hence

\[
 F_j\!\left(z,\frac1{1-y}\right)
 =\frac{1-y}{E-yN_j}
 =\frac{1-y}{E}\sum_{m\ge0}y^m\left(\frac{N_j}{E}\right)^m.
\]

Extracting `n![z^n]` proves

\[
 \boxed{
 \sum_{m\ge0}g_{n,j}(m)y^m
 =\frac{W_{n,j}(y)}{(1-y)^{d_n+1}},
 }
\]

where

\[
 g_{n,j}(m)=n![z^n]\frac{N_j(z)^m}{E(z)^{m+1}}.
\]

Now use

\[
 \frac{N_j}{E}=1+\frac{C_j}{E}.
\]

Define

\[
 c_{n,j,r}=n![z^n]\frac{C_j(z)^r}{E(z)^{r+1}}.
\]

The binomial theorem gives the exact polynomial identity

\[
 \boxed{
 g_{n,j}(m)=\sum_{r=0}^{d_n}c_{n,j,r}\binom mr.
 }
\]

Equivalently, `c_{n,j,r}` is the `r`-th finite difference
`Delta^r g_{n,j}(0)`.

### 4.1 All-index nonnegativity

Since

\[
 \frac1E=\frac1{1-(e^z-1)}=\sum_{q\ge0}(e^z-1)^q,
\]

all coefficients of `1/E` are nonnegative.  `C_j` also has nonnegative
coefficients and starts in degree `j`.  Therefore

\[
 c_{n,j,r}\ge0,\qquad c_{n,j,r}=0\quad(r>d_n),
\]

and in fact `c_{n,j,r}>0` for every `0<=r<=d_n`.  Thus `g_{n,j}(m)>=0` for
all integers `m>=0`, and `deg g_{n,j}=d_n`.

Using `sum_{m>=0} binom(m,r)y^m=y^r/(1-y)^{r+1}` gives

\[
 \boxed{
 W_{n,j}(y)=\sum_{r=0}^{d_n}c_{n,j,r}y^r(1-y)^{d_n-r}.
 }
\]

Let

\[
 h_{n,j}(x)=\sum_{r=0}^{d_n}c_{n,j,r}x^r.
\]

Then

\[
 W_{n,j}(y)=(1-y)^{d_n}h_{n,j}\!\left(\frac y{1-y}\right).
\]

Cancelling the common factor in the two formulas for `W` yields the especially
simple identity

\[
 \boxed{h_{n,j}(x)=O_{n,j}(1+x).}
\]

For `y>1`, `x=y/(1-y)` has derivative `1/(1-y)^2>0`; its endpoint directions
are `y->1+` giving `x->-infinity` and `y->+infinity` giving `x->-1-`.
Consequently the `W`-root condition is exactly that `h_{n,j}` has all roots
in `(-infinity,-1)`.

**This is not a reduction.**  The equality `h_{n,j}(x)=O_{n,j}(1+x)` shows
that the desired root statement for `h` is the original conjecture after the
translation `x=t-1`.  The stable denominator has proved coefficient
nonnegativity, which is equivalent to the already known positivity of the
coefficients of `O`; it has not proved negative real roots.

### 4.2 Exact stable-denominator differential equation

Let

\[
 D_j(z,x)=E(z)-xC_j(z),\qquad
 \mathcal H_j(z,x)=\frac1{D_j(z,x)}
 =\sum_{n\ge0}h_{n,j}(x)\frac{z^n}{n!}.
\]

With `w_j=z^{j-1}/(j-1)!`, one has `C_j'=C_j+w_j` and therefore

\[
 D_j'=D_j-2-xw_j,
 \qquad
 \boxed{\partial_z\mathcal H_j=(2+xw_j)\mathcal H_j^2-\mathcal H_j.}
\]

This identity is exact, but positivity of its series coefficients does not
provide the missing compatible-window or root-preservation argument.  No such
implication is assumed.

## 5. Why the proposed `g`-properties are insufficient: exact obstruction

The following implication is false, even if all three hypotheses are imposed
on the same polynomial:

\[
\begin{aligned}
 &g(m)\ge0\quad(m\ge0),\\
 &g\text{ is gamma-positive},\\
 &\text{all roots of }g\text{ lie in }(-\infty,-1)
\end{aligned}
 \quad\not\Longrightarrow\quad
 \text{all roots of }W\text{ lie in }(1,\infty).
\]

Take

\[
 g(m)=(m+2)^2.
\]

It is nonnegative on `m>=0`, its only root is `-2<-1`, and

\[
 g(m)=4(1-m)^2+12m(1-m)+9m^2,
\]

so it is gamma-positive with coefficients `(4,12,9)`.  Its binomial-basis
coefficients are

\[
 c_0=4,\qquad c_1=5,\qquad c_2=2,
\]

all positive.  Nevertheless

\[
 h(x)=2x^2+5x+4,
 \qquad \operatorname{disc}_x h=25-32=-7,
\]

and

\[
 W(y)=4(1-y)^2+5y(1-y)+2y^2=y^2-3y+4,
 \qquad \operatorname{disc}_y W=9-16=-7.
\]

Thus `W` has two nonreal roots.  This is a feasible-constraint
falsification of the proposed `g`-only implication.  It is **not** a
counterexample to the ordered big-block conjecture: it is a model showing that
the proposed sufficient conditions do not control the finite-difference
polynomial `h`.

For the actual family, the finite exact box `2<=j<=8`, `0<=n<=12` has no
failures of gamma-positivity or of the property "all roots of `g` below `-1`".
That is finite evidence only.  It cannot replace the missing all-index
statement, and the model above shows why even the all-index `g` statement
alone would not suffice.

## 6. Exact remaining obstruction

The reformulation isolates the missing content as follows:

\[
 c_{n,j,r}=\Delta^r g_{n,j}(0),\qquad
 h_{n,j}(x)=\sum_{r=0}^{d_n}c_{n,j,r}x^r=O_{n,j}(1+x).
\]

Nonnegativity and gamma-positivity of `g` do not imply that this binomial
finite-difference transform `h` is real-rooted.  The unresolved all-index
lemma is precisely

\[
 \boxed{
 \text{for every }j\ge2,\ n\ge0,\quad
 h_{n,j}(x)=O_{n,j}(1+x)
 \text{ has all its roots in }(-\infty,-1).
 }
\]

By the exact identities in Sections 2 and 4, this lemma is equivalent to the
original negative-real-rootedness conjecture and to the corrected `W`
statement.  Therefore the `v`/`h*` reformulation supplies exact recurrences,
positive coefficients, and useful tests, but no strict proof reduction from
the stated `g`-properties or from stable-denominator nonnegativity alone.

## 7. Proof-policy audit

### 7.1 Direction audit

All root maps and their derivatives were checked at explicit points and at
both endpoints; see Section 2.2 and the printed audit.  There is no hidden
inequality direction in the algebraic equivalences.  The only monotonicity
claims are the displayed rational maps, whose derivatives are printed and
positive on the stated intervals.

### 7.2 Monotonicity audit

The claimed increasing maps are

\[
 y(t)=1-\frac1t,\quad \frac{dy}{dt}=\frac1{t^2}>0;
\]

\[
 v(t)=\frac{t}{1-t},\quad \frac{dv}{dt}=\frac1{(1-t)^2}>0;
\]

\[
 x(y)=\frac y{1-y},\quad \frac{dx}{dy}=\frac1{(1-y)^2}>0.
$$

For each map the audit states which endpoint gives the lower and upper end
of the target interval.  No optimization or endpoint bound is used to prove
a root statement.

### 7.3 Feasible-constraint falsification

The polynomial `(m+2)^2` in Section 5 satisfies every proposed row-level
constraint simultaneously and violates the desired conclusion.  The exact
gamma coefficients, finite differences, and both discriminants are printed.
This is the required witness that the proposed `g`-only argument is
incomplete.

### 7.4 Boundary and degeneracy enumeration

- `n=0`: `d=0`, `O=R=W=h=1`, `g=1`; no nonconstant root issue.
- `0<n<j`: `d=0`; all transformed root conditions are vacuous and `g` is a
  positive constant.
- `n=j`: `d=1`; the lag term first appears.
- `n` a multiple of `j`: the degree rises and `sigma_{n,1}=1+v`; this
  normalization jump is explicitly present in the recurrence.
- Equal-degree consecutive rows inside a `j`-block use `sigma_{n,1}=1`.
- Inactive `r>d_n` have `c_{n,j,r}=0`; all active `c_{n,j,r}` are positive.
- `v=-1` and `y=1` are regular cancellation points:
  `R_{n,j}(-1)=(-1)^{d_n}a_{d_n}` and `W_{n,j}(1)=a_{d_n}>0`; neither is a root.
- `y=0` means `t=1`; since `O_{n,j}(1)>0`, it is not a root.
- The degree proof excludes cancellation of the top coefficient, so the
  root counts are exact.

### 7.5 Retraction recorded during the boundary audit

An earlier draft incorrectly treated `(1-y)^{d_n}` as a forced root at `y=1`.
The exact boundary residual exposed the error: substitution of `t=1/(1-y)`
gives a pole of the same order, so `W_{n,j}(1)=a_{d_n}>0`.  That qualification
is retracted.  The literal equivalence stated in the task is restored, and no
invalid theorem remains in this note.

### 7.6 Printed symbolic residuals

`verify.py` prints every residual individually.  `AUDIT.txt` is the captured
run.  It checks, for `j=2,...,8` and `n=0,...,12`:

- the original EGF definition against the first-block recurrence;
- the transformed universal recurrence for `R`;
- the `R`/`W` identity and the exact cancellations
  `R(-1)=(-1)^{d_n}a_{d_n}`, `W(1)=a_{d_n}`;
- the `h*` EGF formula directly for `m=0,...,d_n+4`;
- the `c` EGF and finite-difference formulas;
- the `W`/`h` substitution and `h(x)=O(1+x)`;
- the chain rules, denominator identities, stable-denominator Riccati identity,
  positive-root trap, and obstruction model.

All 1,491 printed algebraic residuals are exactly zero.  The finite root and
gamma checks are labelled finite evidence, not proof.

Reproduce from the worktree root with:

```powershell
py -3.14 proof-hstar-reformulation/verify.py
```

## 8. Small exact rows

The following are direct exact outputs, not numerical approximations.

| `j,n` | `O_{n,j}(t)` | `R_{n,j}(v)` | `W_{n,j}(y)` | `g_{n,j}(m)` |
|---|---|---|---|---|
| `2,0` | `1` | `1` | `1` | `1` |
| `2,1` | `1` | `1` | `1` | `1` |
| `2,2` | `2+t` | `2+3v` | `3-2y` | `m+3` |
| `2,3` | `6+7t` | `6+13v` | `13-6y` | `7m+13` |
| `2,4` | `24+45t+6t^2` | `24+93v+75v^2` | `75-93y+24y^2` | `3m^2+54m+75` |
| `2,5` | `120+311t+110t^2` | `120+551v+541v^2` | `541-551y+120y^2` | `55m^2+476m+541` |
| `3,0` | `1` | `1` | `1` | `1` |
| `3,1` | `1` | `1` | `1` | `1` |
| `3,2` | `3` | `3` | `3` | `3` |
| `3,3` | `12+t` | `12+13v` | `13-12y` | `m+13` |
| `3,4` | `66+9t` | `66+75v` | `75-66y` | `9m+75` |
| `3,5` | `450+91t` | `450+541v` | `541-450y` | `91m+541` |
| `3,6` | `3690+973t+20t^2` | `3690+8353v+4683v^2` | `4683-8353y+3690y^2` | `10m^2+1003m+4683` |

## 9. Files

- `NOTES.md` — this proof/obstruction record.
- `verify.py` — exact rational/symbolic verification and policy audits.
- `AUDIT.txt` — full printed zero-residual log.
