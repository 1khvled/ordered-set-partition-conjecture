# `j=2` parametric nonvanishing: exact structure, rigorous reductions, and obstructions

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-j2-resultants` was created/edited.
No installs, no builds, no other agents/models.

## Executive summary

Target statements (from `proof-j2-base/NOTES.md`, §4):

* **(D)** `Disc_t O_n(x,t) != 0` for every `n` and every `x>0`;
* **(R)** `Res_t(O_n(x,t), O_{n+1}(x,t)) != 0` for every `n` and every `x>0`.

What this directory establishes rigorously:

1. **(Exact equivalence.)** Given the `x -> infinity` Fibonacci limit proved in
   `proof-j2-base/NOTES.md` §3, statement (D) is *equivalent* to: for every `n`,
   `O_n(x,t)` has only simple real roots for every `x>0`. Statement (D)+(R) is
   equivalent to: the rows `O_n(x,·)` strictly interlace for every `x>0`, in the
   Fibonacci orientation. In particular (D)+(R) at `x=1` *is* the `j=2`
   real-rootedness/interlacing conjecture, which is open (SymCat, issue #108).
   So the requested statements are not finite checks; any proof must resolve an
   open all-`x` problem. This is the principal exact obstruction. (§2)

2. **(Rigorous component theorem.)** The exact `x`-expansion
   `O_n(x,t) = sum_k x^k P_n^{(k)}(t)` has components
   `P_n^{(k)}(t) = (n!/m!) B_{m,a}(t)`, `m=n-k`, `a=k+1`, where
   `B_{m,a}(t) = sum_r (a)_r {m\brace r}_{>=2} t^r`. Each `B_{m,a}` has only
   simple nonpositive real roots (Theorem A), and the identity
   `B_{m,a+1} = B_{m,a} + (1/a) t B'_{m,a}` holds exactly (Theorem B), so the
   components interlace in `a` for fixed `m`. (§3)

3. **(Exact obstruction to the two most natural sufficient conditions.)**
   * The deformation is *not* a Pólya-frequency multiplier on the Fibonacci
     limit: there is an exact negative Toeplitz minor at `n=6, x=1`,
     `-1939704113/5832000000`. (§4.1)
   * The interlacing-chain components admit no low-order differential
     recurrence with constant coefficients: at `n=8` even a representation
     `D_2 = c_0 D_1 + c_1 t D_1' + c_2 t^2 D_1''` has no solution. (§4.2)
   * Strict interlacing is not transitive, so horizontal+vertical component
     interlacing cannot give the diagonal step (explicit quadratics). (§4.3)

4. **(Finite exact verification, explicitly labelled as finite.)**
   `Disc_t O_n` and `Res_t(O_n,O_{n+1})`, viewed as polynomials in `x`, become
   constant-sign polynomials after removing their lowest `x`-monomial, for all
   `2<=n<=9` (disc) and `1<=n<=8` (res), and Sturm counts show zero positive
   roots for `n<=11`/`n<=10`. This is evidence only, not an all-index proof. (§5)

The all-index gap is exactly the *diagonal* interlacing of the
`B_{M-a,a}` (§3.4): it is verified numerically through `n=12` and is equivalent
to real-rootedness of `O_n(x,t)` for all `x>0`, hence to the open conjecture.

---

## 1. Setup and exact coefficient formula

Write `C(z)=e^z-1-z`, `F(x,z,t)=1/(1-xz-tC(z))`, `O_n(x,t)=n![z^n]F`.
Expanding `1/(1-xz-tC)=sum_{r>=0} t^r C^r/(1-xz)^{r+1}` gives

```
a_{n,r}(x) := [t^r] O_n(x,t)
           = n! r! sum_{m=2r}^n {m\brace r}_{>=2} x^{n-m}/m! * binom(n-m+r, n-m).   (1)
```

Every coefficient is a polynomial in `x` with **positive integer
coefficients**, and the top coefficient is `n!/2^r`, independent of `x`.
Hence for `x>0` all `a_{n,r}(x)>0`, `O_n(x,0)>0`, the degree in `t` is
`d=floor(n/2)` for all `x>0`, and all roots lie off the positive real axis.
These are the only facts about positivity used below.

The fundamental PDE identities (derived in `proof-j2-base/NOTES.md` §1, valid
for general `x`):

```
t O_n' = U_n - n x U_{n-1} - O_n,        U_n = n![z^n]F^2,
O_{n+1}= (1+x)U_n + n(t-x)U_{n-1} - O_n,
O_{n+1}= x O_n + (1+x) t O_n' + n(x^2+t) U_{n-1}.                       (2)
```

At a hypothetical double root `O_n=O_n'=0` one gets `U_n=nx U_{n-1}`,
`O_{n+1}=n(x^2+t)U_{n-1}` and `dO_n/dx = n U_{n-1}`. These identities are exact
but do not by themselves exclude a collision (they are consistent with one).

---

## 2. Exact equivalence of (D)+(R) with the all-`x` conjecture

**Limit input (proved in `proof-j2-base/NOTES.md` §3).** With `t=x^2 u`,

```
x^{-n} O_n(x,x^2 u) -> Phi_n(u),   Phi_n(u)=n! sum_r binom(n-r,r)/2^r u^r,
```

and `Phi_n(v/(-2))` is the Fibonacci polynomial `P_n(v)=sum_j (-1)^j binom(n-j,j) v^j`.
`Phi_n` has `d=floor(n/2)` simple negative roots, and
`Phi_{n+1}(u)/Phi_n'(u)<0` at each root of `Phi_n`.

**Lemma 1 (D `<=>` simple real roots for all `x>0`).**
For fixed `n`, `O_n(x,t)` has degree `d` in `t` with `x`-independent nonzero
leading coefficient (§1). For large `x` its roots are `x^2` times the roots of
`Phi_n` (up to `o(1)`), all real and simple. A real polynomial's number of
nonreal roots changes only at a real multiple root (discriminant zero); roots
cannot escape to infinity since the leading coefficient is constant. Hence if
`Disc_t O_n(x,·) != 0` for all `x>0`, all roots are real and simple for all
`x>0`. Conversely, simple real roots give `Disc>0`. So (D) `<=>` "all roots
real and simple for all `x>0`".

**Lemma 2 (D+R `<=>` strict interlacing for all `x>0`).** Assume (D). The sign
of `Res_t(O_n,O_{n+1})` can change only at a common root, excluded by (R). At
`x->infinity` the rows interlace in the Fibonacci orientation. Continuous
deformation with no collision and no common root preserves the alternation
pattern. Hence (D)+(R) gives strict consecutive interlacing for every `x>0`,
and conversely interlacing gives (D) and (R).

**Consequence (exact obstruction).** Evaluating at `x=1`, (D) implies the
`j=2` real-rootedness conjecture and (D)+(R) implies the `j=2` consecutive
interlacing conjecture. Both are open: SymCat's
`realRootedSetPartitions.htm` records that "when all block sizes are allowed
but only the large ones are counted ... the analogous statement for ordered
set partitions" is open; the repo `literature.md` §5 and issue #108 agree.
Therefore no argument based on finitely many `n`/`x` checks can prove (D) or
(R); an all-index all-parameter idea is unavoidable. This is the exact
structural obstruction to the requested task as a routine computation.

---

## 3. Rigorous component results

### 3.1 Decomposition

From `F = sum_{k>=0} x^k (1-tC)^{-(k+1)}`,

```
O_n(x,t) = sum_{k=0}^{n} x^k P_n^{(k)}(t),
P_n^{(k)}(t) = n![z^{n-k}](1-tC(z))^{-(k+1)}.
```

Put `m=n-k`, `a=k+1`. Since `[z^m] C^r = r! {m\brace r}_{>=2}/m!`,

```
P_n^{(k)}(t) = (n!/m!) B_{m,a}(t),    B_{m,a}(t) = sum_{r>=1} (a)_r {m\brace r}_{>=2} t^r.   (3)
```

(`(a)_r=a(a+1)...(a+r-1)`.) Note `B_{m,1}=sum_r r!{m\brace r}_{>=2}t^r`
is the ordered no-singleton Bell polynomial `D_m` of Bóna–Mező.

### 3.2 Theorem A: each component is real-rooted (proved)

**Theorem A.** For every integer `m>=2` and `a>=1`, `B_{m,a}` has only simple
nonpositive real roots (and a simple root at `0`).

*Proof.* Bóna–Mező (EJC 51 (2016) 500–510; SymCat "Set partitions without
singleton blocks") prove `B_{m,1}=D_m` has only real nonpositive zeros. From
(3), `B_{m,a}` is obtained from `B_{m,1}` by multiplying the coefficient of
`t^r` by `c_r=binom(a+r-1,r)=(a)_r/r!`. The generating function is
`sum_r c_r z^r=(1-z)^{-a}=prod_{i=1}^{a}(1-z)^{-1}`, which by the
Aissen–Edrei–Schoenberg–Whitney characterization is a Pólya-frequency formal
power series (parameters `beta_1=...=beta_a=1`, `sum beta_i=a<infinity`).
A PF sequence is a multiplier sequence (Pólya–Szegő), so it preserves
real-rootedness with nonpositive roots; all `c_r>0`, so no coefficient is
annihilated and the roots stay simple (strict Toeplitz total positivity of the
finite truncation). ∎

### 3.3 Theorem B: exact interlacing in the `a`-direction (proved)

**Theorem B.** `B_{m,a+1} = B_{m,a} + (1/a) t B'_{m,a}` for all `a>=1`, `m>=2`.
Consequently `B_{m,a}` interlaces `B_{m,a+1}` (and both are real-rooted by
Theorem A).

*Proof.* `(a+1)_r=(a)_r (a+r)/a=(a)_r(1+r/a)`; the term `r(a)_r t^r` is exactly
the coefficientwise action of `t d/dt` on `B_{m,a}`. This proves the identity.
For the interlacing: the operator `1+c t d/dt` is a nonnegative combination of
the multiplier sequences `(1)_r` and `(r)_r` (both PF, the latter because
`sum_r r z^r=z/(1-z)^2` is PF), hence a multiplier sequence for all `c>=0`;
and `B_{m,a}` interlaces `t B'_{m,a}` by Rolle plus `B_{m,a}(0)=0`. By
Obreschkoff/Wagner's lemma the positive combination `B_{m,a}+(1/a)t B'_{m,a}`
interlaces `B_{m,a}`. ∎

### 3.4 The missing diagonal step

`O_n(x,t)/t = sum_k x^k (n!/m!) B_{m,k+1}/t` with `m+k=n+1=:M` is a positive
combination along the **diagonal** `m+a=M` of the grid `(B_{m,a})`. By the
interlacing-cone theorem (Obreschkoff), (D) would follow if the diagonal chain
`B_{M-a,a}/t` (`a=1,...,M-2`) strictly interlaced. It does so numerically for
every `n<=12` (`interlace_probe.py`), and both directions of grid interlacing
hold:

* **vertical**, `B_{m,a}` vs `B_{m,a+1}`: Theorem B (proved);
* **horizontal**, `B_{m,a}` vs `B_{m-1,a}`: the *same* PF multiplier
  `c_r=binom(a+r-1,r)` applied to the Bóna–Mező chain `D_m` vs `D_{m-1}`;
  PF multiplier sequences preserve interlacing of real-rooted polynomials
  (variation diminution). Verified exactly for `a<=3, m<=11`.

But:

* horizontal+vertical interlacing does **not** imply diagonal interlacing in
  general (§4.3);
* the diagonal neighbors admit no low-order differential recurrence with
  constant coefficients (§4.2), so the step cannot be obtained by the usual
  "positive combination of f and derivatives" induction.

This diagonal statement is precisely real-rootedness of `O_n(x,·)` for all
`x>0` (Lemma 1), i.e. the open conjecture. It is the exact remaining gap.

---

## 4. Exact obstructions to natural sufficient conditions

All numbers in this section are exact rationals; scripts reproduce them.

### 4.1 The PF-multiplier route on the full row fails

Writing `O_n(x,x^2 u)/n! = sum_r [binom(n-r,r)/2^r] m_{n,r}(x) u^r`, the limit
polynomial `Phi_n(u)/n!` is real-rooted, so real-rootedness would follow if
`(m_{n,r}(x))_r` were PF for each `x>0`. It is not. At `n=6, x=1`,

```
(m_{6,r}(1))_{r=0}^{3} = (1, 2383/1800, 149/108, 1),
```

and the `3x3` lower-triangular Toeplitz minor on rows `(1,2,3)`, columns
`(0,1,2)` is

```
det = -1939704113/5832000000 < 0.
```

(`witnesses.py`.) So the multiplier method cannot prove even real-rootedness,
let alone (D). This also shows why the proved component theorem (§3.2) does not
transfer to the full row.

### 4.2 No constant-coefficient differential recurrence on the diagonal

Let `D_a = B_{M-a,a}` be consecutive diagonal elements. At `M=n+1=9`
(`n=8`), the first step `D_2` vs `D_1` has no representation
`D_2=c_0 D_1 + c_1 t D_1'` with constants `c_0,c_1`, and *even*
`D_2=c_0 D_1 + c_1 t D_1' + c_2 t^2 D_1''` has no solution (`witnesses.py`).
Hence the standard Rolle/Wagner induction cannot produce the diagonal
interlacing from the component identities of §3.

### 4.3 Strict interlacing is not transitive (hence grid `!=` diagonal)

Take

```
f=(t+1)(t+4)  roots {1,4},
g=(t+2)(t+3)  roots {2,3},
h=(t+5/2)(t+7/2)  roots {2.5,3.5}.
```

Then `f` strictly interlaces `g` (`1<2<3<4`) and `g` strictly interlaces `h`
(`2<2.5<3<3.5`), but `f` does **not** interlace `h`: one needs
`1 < h_1 < 4 < h_2`, whereas `h_2=3.5<4`. Thus even the two proven component
interlacings — horizontal (`B_{m,a}` vs `B_{m-1,a}`, from Bóna–Mező plus the
PF multiplier of Theorem A) and vertical (Theorem B) — cannot be composed to
give the required diagonal interlacing `B_{M-a,a}/t` vs
`B_{M-a-1,a+1}/t`; a genuinely diagonal argument is required.

### 4.4 Earlier obstructions (from `proof-j2-base`)

The falling-basis ratio monotonicity, the simple disjoint-quadratic
certificate for the lag-two Wronskian, and the association-Bell interlacing
inference are all blocked there by explicit exact counterexamples; those
obstructions are unaffected by the present work.

---

## 5. Finite exact verification (NOT a proof)

`probe.py`, `sturm_probe.py`:

* `Disc_t O_n(x,t)` / `x^{k_n}` has all positive coefficients for `2<=n<=9`;
* `Res_t(O_n,O_{n+1})` / `x^{ell_n}` has constant-sign coefficients for
  `1<=n<=8` (sign negative for `n<=5`, positive for `n>=6`, reflecting the
  resultant normalization convention, not a root phenomenon);
* Sturm sequences give **zero** roots in `(0,infinity)` for the disc cofactors
  through `n=11` and the resultant cofactors through `n=10`. (`sympy`'s
  `count_roots(0,oo)` counts the removed `x=0` factor; the scripts strip that
  monomial first.)

These confirm the statements on a finite range and display the exact
positive-coefficient pattern, but by Lemma 1 they cannot constitute an
all-index proof of (D): doing so would settle the open `j=2` conjecture.

---

## 6. What would close the gap

Prove the diagonal interlacing of `B_{M-a,a}/t` for all `M` (equivalently,
real-rootedness of `O_n(x,·)` for all `x>0`), e.g. by:

* a planar-network / LGV proof of total nonnegativity of the coefficient array
  `(a_{n,r}(x))` (the natural left-production factorization is known to fail,
  `literature.md` §8, but another network may exist);
* a diagonal Wagner-cone argument stronger than §3.3, using the exact identity
  of Theorem B together with the Bóna–Mező recurrence for `B_{m,1}`;
* an eventual/asymptotic argument in `n` combined with an explicit finite
  threshold (as in Fang's 2026 eventual-interlacing paper for `b=4`), which
  would have to cover the middle roots.

Given Lemma 1, any such proof automatically proves the open `j=2` ordered
big-block real-rootedness conjecture; that is the exact obstruction to a
shorter answer.

---

## 7. Files

* `NOTES.md` — this report.
* `probe.py`, `sturm_probe.py` — exact sign-pattern / Sturm checks (finite).
* `cone_probe.py`, `interlace_probe.py`, `interlace2_probe.py` — component
  interlacing probes.
* `pf_probe.py`, `row_pf_probe.py` — PF tests (one gives the §4.1 witness).
* `recur_probe.py`, `relation_probe.py` — differential-relation searches
  (Theorem B identity and the §4.2 failure).
* `witnesses.py` — exact rational witnesses for §4.1–4.3.

Run from this directory, e.g. `py -3.14 witnesses.py`.
No packages installed; no builds; all files small.
