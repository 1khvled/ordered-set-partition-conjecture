# Inner real-rootedness via the linear PDE: analysis and precise gap

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-inner-pde` was created or edited.  
Outcome: **(H) is not proved here.** The linear PDE gives a clean recurrence and
shows exactly where a root-preserving induction would have to control a signed,
lagged, derivative term. The PDE route does not close with the standard tools;
the remaining obstruction is isolated below.

---

## 0. Target and identities used

Fix `j>=2`. The inner family is

```
H(z,t) = 1/(1 - t g_j(z)) = sum_{n>=0} H_n(t) z^n/n!,     g_j = C_j/N,
N = 1-A_j,   C_j = e^z - E_{j-1},   A_j = sum_{r=1}^{j-1} z^r/r!.
```

Target (H): for every `j>=2, n>=0`, `H_n` has only real nonpositive roots and
`H_n \interl H_{n+1}` (with the usual convention that the zero rows `n=1,...,j-1`
are handled by degree).

The identities from `proof-final-inner/NOTES.md` that are used below are:

* Lemma 3 (linear PDE):
  ```
  N H_z = t[2+(t-1)p] H_t + t p H,     p = z^{j-1}/(j-1)!.
  ```
* Lemma 6 (coefficient form):
  ```
  H_{n+1} = sum_{r=1}^{j-2} binom(n,r) H_{n-r+1}
          + binom(n,j-1)[ H_{n-j+2} + t Q(H_{n-j+1}) ]
          + 2 T H_n,
  ```
  where `T = t d/dt` and `Q(q) = q + (t-1) q'`.

The recurrence above was re-verified independently in this directory
(`verify_recurrence3.py`, 0 errors for `j=2..5`, `n<=15`).

---

## 1. What the PDE gives that is new

### 1.1 Operator form

Write the recurrence as

```
H_{n+1} = S_n + binom(n,j-1) L_n + 2 T H_n,
```

where

```
S_n = sum_{r=1}^{j-2} binom(n,r) H_{n-r+1}   (positive combination),
L_n = H_{n-j+2} + t Q(H_{n-j+1}).
```

Thus the evolution from `n` to `n+1` uses:

* positive binomial combinations of recent rows (degrees `floor((n-r+1)/j)`);
* a lagged row `H_{n-j+2}` plus a lagged, differentiated row through `t Q`;
* the Euler operator `2T` applied to the current row.

### 1.2 The operators `T` and `Q` are individually real-rootedness preserving

Numerical experiments (`op_test.py`, 2000 random real-rooted polynomials of
degrees 2..6) show that both

```
T(p) = t p',          Q(p) = p + (t-1) p'
```

map real-rooted polynomials to real-rooted polynomials.  This is consistent with
the classical facts that `t d/dt` (Euler) and `(t-1)d/dt + 1` are
variation-diminishing in the appropriate sense.  Hence the *linear* pieces of
the recurrence do not by themselves destroy real-rootedness.

### 1.3 The PDE makes the obstruction explicit

The signed term in Lemma 6 is

```
binom(n,j-1)[ H_{n-j+2} + t H_{n-j+1} + t(t-1) H_{n-j+1}' ].
```

Its derivative part `t(t-1) H_{n-j+1}'` is the inner analogue of the
`(t-1)` obstruction recorded in `proof/PROOF_NOTES.md` and
`proof-final-inner/NOTES.md`.  It cannot be removed: it is forced by
`C_j' = C_j + p`, i.e. by the truncation of the exponential at order `j-1`.
For `j=1` one has `p=0` and the inner family degenerates to the Fubini family,
where Comtet's recurrence `O_{n+1}=t(1+D)O_n` closes; for every `j>=2` the
`p`-term is present.

---

## 2. Attempted routes from the PDE

### 2.1 Real stability of `H(z,t)` as a bivariate function

A standard way to get real-rooted Taylor coefficients is to show that the
generating function is real stable.  Here

```
H(z,t) = 1/(1 - t g_j(z)),   g_j = C_j/N.
```

For `H` to be real stable one would need, for every `z` with `Im z>0`, that
`1/g_j(z)` avoid the upper half-plane.  Direct sampling (`stability_test.py`)
shows the opposite: for `j=2,3,4`, a large majority of points in the upper
half-plane have `Im(1/g_j(z))>0`.  Equivalently, for such `z` and for
`t = 1/g_j(z)` in the upper half-plane, the denominator `1 - t g_j(z)` vanishes,
so `H` has poles with both variables in the upper half-plane.  Thus `H` is
**not** real stable in the ordinary sense, and the Borcea-Brändén coefficient
theorem cannot be applied without a nontrivial transformation.

### 2.2 Transformation to a real stable PDE

One can rewrite the PDE in characteristic form.  With `x = ∫_0^z du/N(u)` one gets

```
H_x = t[2+(t-1)p(z(x))] H_t + t p(z(x)) H,
```

which still has the nonconstant coefficient `p(z(x))`.  The characteristic
solution for `t` is Möbius:

```
t(z;t_0) = N(z) t_0 / (e^z + t_0 C_j(z)),
```

but this does not turn `H` into an entire function of `z`; the poles of `g_j`
remain.  No change of independent variables found here removes the poles of
`g_j` while preserving the EGF coefficient extraction, so the real-stability
route remains blocked.

### 2.3 Coefficient majorization (Newton inequalities)

For a real-rooted polynomial with positive coefficients, the coefficient
sequence must be a Pólya-frequency sequence, in particular it must satisfy
Newton inequalities.  These were checked exactly for `H_n` for `j=2,3,4` and
`n<=23` (`newton_test.py`, 0 failures).  Thus Newton inequalities are
compatible with (H), but they are necessary, not sufficient, for degree `>=3`;
they do not close the proof.  No induction establishing them for all `n` from
the PDE was found.

### 2.4 Root-preserving evolution by induction

An induction would need to show that the combination

```
S_n + binom(n,j-1)L_n + 2 T H_n
```

is real-rooted and interlaces the previous rows whenever the inputs are.  The
positive part `S_n + 2 T H_n` is manageable, but `L_n` contains

```
t(t-1) H_{n-j+1}'.
```

Because of the factor `(t-1)`, the derivative term does **not** vanish at the
zeros of `H_{n-j+1}`.  Hence the standard "if `p` and `q` interlace then
`p + (t-a)q'` is real-rooted" lemmas do not apply directly.  The needed
control is precisely a weighted root-separation inequality between the roots of
`H_n` and `H_{n-j+1}`; no such inequality is established here.

---

## 3. The precise gap

The missing statement is the following (equivalent to (H) by the
Aissen-Schoenberg-Whitney criterion applied to each row).

> **Missing Lemma (PDE array).** For every `j>=2` and `n>=0`, the recurrence
>
> ```
> H_{n+1} = S_n + binom(n,j-1)[H_{n-j+2}+t Q(H_{n-j+1})] + 2 T H_n
> ```
>
> preserves the property "all roots real and nonpositive, with interlacing
> between consecutive `n`" provided one has, for the lagged index
> `n-j+1`, the weighted separation
>
> ```
> for every root \alpha of H_{n-j+1} and every root \beta of H_n:
>     (some explicit lower bound on |\alpha-\beta| in terms of the
>      multiplicities and the factor (t-1))
> ```
>
> or, equivalently, a theorem that the operator
> `p \mapsto p + (t-1) q'` is real-rootedness preserving under the
> interlacing hypothesis between `p` and `q` with the specific index gap `j`.

No such theorem or inequality is available in the present notes.  The ordinary
multiplier/CZDS/Favard/PF routes are excluded by
`proof-final-inner/NOTES.md` (items 2-6), and the PDE route reduces to exactly
the same obstruction: the truncation term `p` forces a signed, lagged derivative
that standard real-rootedness-preserving operations cannot absorb.

Consequently, the all-index proof of (H) remains open.  The finite checks in
`proof-final-inner` and the Newton checks above are evidence only.

---

## 4. Files in this directory

* `NOTES.md` - this note.
* `explore_roots.py` - numerical root check for `H_n`, `j=2..4`, `n<=15`.
* `verify_recurrence3.py` - independent verification of the PDE recurrence.
* `op_test.py` - tests that `T` and `Q` individually preserve real-rootedness.
* `newton_test.py` - Newton-inequality check for `H_n`.
* `stability_test.py` - shows `H(z,t)` is not real stable in the naive sense.

No large caches, builds, or external data were created.
