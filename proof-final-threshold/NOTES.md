# Threshold induction for ordered big blocks: clean identities and a no-go theorem

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-final-threshold` was created or edited.

**Outcome.** An induction on the threshold `j` using derivatives/integrals/shifts/
Möbius transforms/convolution **cannot** work: the pole of `F_j` in the `t`-plane
moves with `z` in a way no operator of that class can reproduce. The only exact
threshold relation is a nonlinear Möbius renewal; its coefficient form contains
a product of rows of the two families. This note proves the no-go theorem, records
the exact identities, and isolates the remaining weighted-separation obstruction.
No PF/CZDS/PDE-stability route is repeated.

---

## 0. Conventions

For `s >= 1` write

```
E_s(z)   = sum_{r=0}^s z^r/r!,
A_j(z)   = E_{j-1}(z) - 1 = sum_{r=1}^{j-1} z^r/r!,
C_j(z)   = e^z - E_{j-1}(z),
N_j(z)   = 1 - A_j(z) = 2 - E_{j-1}(z),
p_j(z)   = z^{j-1}/(j-1)!,
B_j(z,t) = A_j(z) + t C_j(z),
F_j(z,t) = 1/(1 - B_j(z,t)) = sum_{n>=0} O_{n,j}(t) z^n/n!,
tau_j(z) = N_j(z)/C_j(z)   (the pole of F_j in t).
```

Set `P_n = O_{n,j-1}`, `Q_n = O_{n,j}`. All identities are formal in `z` and
rational/exact in `t`.

---

## 1. Exact threshold identities (all proved)

### Lemma 1 (derivative identity).

```
B_j'(z,t) = 1 + B_{j-1}(z,t),                        (T1)
equivalently   (1/F_j)' = 1/F_{j-1} - 2.
```

*Proof.* `A_j' = 1 + A_{j-1}` (differentiate `sum_{r=1}^{j-1} z^r/r!`) and
`C_j' = C_j + p_j`. Hence

```
B_j' = (1 + A_{j-1}) + t(C_j + p_j) = 1 + A_{j-1} + t C_{j-1} = 1 + B_{j-1},
```

using `C_{j-1} = C_j + p_j`. The reciprocal form follows from `1/F_s = 1 - B_s`. ∎

This is the only *linear* relation between the two generating functions, and it
is a relation between **reciprocals**, whose coefficients are trivial.

### Lemma 2 (pole separation). With `x = e^z`,

```
tau_j(z) - tau_{j-1}(z)
    = - p_j(z) (e^z - 2) / [ (e^z - E_{j-1}(z)) (e^z - E_{j-2}(z)) ].   (T2)
```

*Proof.* `E_{j-1} = E_{j-2} + p_j` and `N_j = N_{j-1} - p_j`. Write
`P=E_{j-2}`, `Q=2-P`; then `tau_{j-1}=Q/(x-P)`, `tau_j=(Q-p_j)/(x-P-p_j)`.
The numerator of the difference is
`(Q-p_j)(x-P) - Q(x-P-p_j) = -p_j x + p_j(P+Q) = -p_j(x-2)`. ∎

In particular `tau_j(z) != tau_{j-1}(z)` for all `z>0` with `z != log 2`, and
`tau_j` is nonconstant: as `z -> 0+`, `N_j(z) -> 1` while `C_j(z) ~ z^j/j!`,
so `tau_j(z) ~ j!/z^j -> +infty`.

### Lemma 3 (Möbius renewal, the exact threshold relation).

```
F_j = F_{j-1} / ( 1 - (1-t) p_j F_{j-1} )
    = F_{j-1} ( 1 + (1-t) p_j F_j ).                                  (T3)
```

*Proof.* `1/F_j = 1/F_{j-1} + (t-1)p_j = 1/F_{j-1} - (1-t)p_j`; solve for `F_j`.
The second form is the identity `F_j - F_{j-1} = (1-t)p_j F_j F_{j-1}` rearranged. ∎

Coefficient form (exact, for `n >= 0`, with the second term absent for `n<j-1`):

```
Q_n = P_n + (1-t) binom(n,j-1)
          * sum_{m=0}^{n-j+1} binom(n-j+1,m) Q_m P_{n-j+1-m}.        (T4)
```

Thus `Q_n - P_n` is divisible by `(1-t)` (as a polynomial in `t`), and the
correction is **bilinear** in the two families and contains the signed factor
`(1-t)`. This is the only coefficient-level threshold recurrence; it is not a
linear operator on `(P_0,P_1,...)`.

### Lemma 4 (base `j=2`, outer convolution). For `j=2`, `A_2=z`, hence
`D_2 = 1/(1-A_2) = sum_{s>=0} z^s` and

```
O_{n,2}(t) = n! sum_{k=0}^n H_{k,2}(t)/k!,                            (T5)
```

where `H_{*,2}` is the inner family `1/(1-t C_2/(1-A_2)) = sum H_k z^k/k!`.
Note the weights are `n!/k!`, **not** binomial coefficients; the usual
binomial-transform interlacing theorem does not apply verbatim.

All five identities were checked exactly for `j<=7`, `n<=14`
(`verify_threshold.py`, 0 failures).

---

## 2. Obstruction A: no linear/Möbius threshold transformation exists

### Theorem A (linear operators in `z`).

Let `L` be any finite differential operator

```
L f = sum_{r=0}^R a_r(z,t) ∂_z^r f,
```

where each `a_r` is analytic in `z` and rational in `t` with its pole set
contained in a fixed finite set `S` independent of `z`. Then `L F_{j-1} != F_j`
for every `j >= 2`.

*Proof.* Fix `z>0` outside the (discrete) zero sets of `C_{j-1},C_j` and of the
coefficients of `L`. As a function of `t`, `F_{j-1}` has the single pole
`tau_{j-1}(z)` and `F_j` has the single pole `tau_j(z)`. Every term of
`L F_{j-1}` has its `t`-poles inside `S ∪ {tau_{j-1}(z)}`, since differentiating
in `z` does not move a singularity in the `t`-variable. Since `tau_j` is
nonconstant while `S` is finite, choose `z` with `tau_j(z) notin S` and
`z != log 2`; by Lemma 2 also `tau_j(z) != tau_{j-1}(z)`. Then the right-hand
side is analytic at `t = tau_j(z)` while the left-hand side has a pole there.
Contradiction. ∎

The same argument excludes integrals `∫_0^z b(u) f(u,t) du` and convolutions in
`z` with `t`-independent kernels **for any `z` such that**
`tau_j(z) notin {tau_{j-1}(u) : 0<=u<=z} ∪ S`. Such `z` exist: for `j=2` any
`z>1` works (`tau_1(u)=1/(e^u-1)>0` for `u>0` while `tau_2(z)<0` for `z>1`); for
`j=3`, `z=1.3` works (`tau_3(1.3) ≈ -2.184`, whereas `tau_2(u) >= -0.22` on
`[0,1.3]`). A convolution kernel with `z`-dependent `t`-poles is not a
transformation determined by `F_{j-1}` alone: its pole must already be the
target pole `tau_j(z)`, so it presupposes the answer. Shifts `z -> z+c` would
produce the pole `tau_{j-1}(z+c)`, not `tau_j(z)`, and are excluded by the same
pole-locus comparison.

### Theorem B (Möbius transformations in `t`).

There are **no** constants `a,b,c,d` with `ad-bc != 0` such that

```
tau_j(z) = (a tau_{j-1}(z) + b) / (c tau_{j-1}(z) + d)   for all z>0,  j>=2.
```

*Proof.* Put `x=e^z`, `P=E_{j-2}`, `p=p_j`, so `E_{j-1}=P+p`,
`tau_{j-1} = Q/(x-P)` with `Q=2-P` and `tau_j=(Q-p)/(x-P-p)`. Cross-multiplying
the assumed identity and multiplying by `(x-P)(x-P-p)` gives an equality of two
polynomials in `x` whose coefficients are functions of `z`:

```
(Q-p)[c Q + d(x-P)] = (x-P-p)[a Q + b(x-P)].
```

The `x^2` coefficient of the left side is `0`, of the right side is `b`; hence
`b=0`. Comparing the `x` coefficients with `b=0` gives
`d (Q-p) = a Q` identically in `z`. But `(Q-p)/Q` is nonconstant because
`p=z^{j-1}/(j-1)!` and `Q=2-E_{j-2}` are linearly independent functions; hence
`a=d=0`, which forces `tau_j=0`, impossible. ∎

Checked symbolically for `j=2,...,8` (`moebius_obstruction.py`).

**Consequence.** None of the transformations listed in the brief — derivatives,
integrals, shifts, Möbius transforms in `t`, or fixed-kernel convolutions — maps
the threshold `j-1` family to the threshold `j` family. The exact map is the
nonlinear renewal (T3)/(T4), which involves the self-convolution of `F_{j-1}` to
all orders (the geometric series in `(1-t)p_j F_{j-1}`) and therefore cannot be
replaced by a root-preserving linear operator on the coefficient vectors.

---

## 3. Obstruction B: even the known inner theorem does not close the transform

The inner real-rootedness statement (H), which the repo notes list as open, is in
fact a published theorem: the `k`-associated Bell polynomials are real-rooted and
any `k+1` consecutive rows interlace pairwise (Alexandersson, 2026; SymCat page
*Set partitions*, example “Set partitions without singleton blocks” and the
following paragraph). The multiplier `(1/m!)` removes the ordering, so the inner
rows `H_{n,j}` have the same property.

This still does **not** give the outer transform:

1. For `j=2`, (T5) is a factorial-weighted (not binomial) convolution of the
   inner rows. The inner rows are not globally interlacing: for `j=2`,
   `H_2 = t` and `H_6 = t(90t^2+830t+517)`. The only root of `H_2` is `0`, while
   the two nonzero roots of `H_6` are both negative (discriminant
   `830^2-4*90*517 = 502780 > 0`, positive coefficients), so `H_2` does not
   interlace `H_6`. The interlacing-sequence hypotheses of the standard
   binomial-transform lemmas are not available.

2. For general `j`, the would-be interlacing family
   `F_n = {H_n} ∪ {binom(n,r) O_{n-r} : 1<=r<=j-1}` of
   `proof-final-transform/NOTES.md` fails at a substep. For `j=3, n=7`,

   ```
   H_7 = t(350t+2801)      roots: -2801/350, 0,
   O_3 = t+12               root: -12.
   ```

   Since `-12 < -2801/350 < 0`, the root of `O_3` is not between the two roots
   of `H_7`; hence `H_7` does **not** interlace `O_3`. But `O_3` enters the
   family member `O_5` through
   `O_5 = H_5 + 5 O_4 + 10 O_3` with a positive coefficient, and `O_5` is a
   member of `F_7`. Thus the pairwise-interlacing criterion cannot establish
   Lemma I there; a weighted separation inequality is necessary. The same
   failure of `H_m` to interlace `O_{m-(j+1)}` (the distance needed when the
   feedback row is expanded once further) occurs for `j=4,5`; see
   `strong2_probe.py`. Direct checks of the whole family `F_n` up to `n=16`
   found no violation (`family_probe2.py`), which is why the gap is at the
   inductive *substep*, not at the final row.

This is exactly the weighted root-separation obstruction already isolated in
`proof-final-wrs/NOTES.md`, now seen to survive even after the inner theorem is
granted and even if one tries to import `j-1` information through the threshold
renewal.

---

## 4. Precise statement of the obstruction to threshold induction

> **Theorem C.** Let `j>=2`. (a) No transformation of the class described in
> Theorems A–B carries `(O_{n,j-1})_{n>=0}` to `(O_{n,j})_{n>=0}`; the sole exact
> relation is the nonlinear renewal (T3)/(T4). (b) Knowledge that
> `(O_{n,j-1})` is real-rooted with interlacing consecutive rows — even augmented
> by the theorem that `(H_{n,j-1})` and `(H_{n,j})` are real-rooted with window
> interlacing — does not imply the interlacing-family criterion needed to run the
> outer convolution: the exact instance `j=3,n=7` above violates it. Consequently
> an induction on `j` cannot be completed by root-preserving linear operations
> alone; it must prove the weighted separation forced by the `(1-t)` factor in
> (T4).

What would suffice (none is established here):

* a renewal-PF theorem for `1/(1-t g_j)` that does not require `g_j` itself to be
  PF/CZDS (the ordinary multiplier and CZDS routes are excluded in
  `proof-final-inner/NOTES.md`, §2.6); or
* an all-index proof of the weighted separation inequality at the level of (T4),
  i.e. control of the product term `Q_m P_{n-j+1-m}` at the roots; or
* a genuinely nonlinear root-preserving transform that moves the `t`-pole from
  `tau_{j-1}(z)` to `tau_j(z)` while preserving interlacing (Theorems A–B rule
  out all linear and Möbius candidates).

---

## 5. Files and verification

* `NOTES.md` — this note.
* `verify_threshold.py` — exact symbolic checks of (T1)–(T5), `j<=7`, `n<=14`; 0 errors.
* `moebius_obstruction.py` — symbolic coefficient comparison for Theorem B, `j<=8`.
* `pole_probe.py` — numerical pole values and failed Möbius fit (illustration only).
* `H_interlace.py`, `family_probe.py`, `family_probe2.py`, `strong_probe.py`,
  `strong2_probe.py`, `H_global.py`, `j_monotone.py`, `roots_probe.py`,
  `explore.py`, `fit.py`, `lp_probe.py` — exact/numerical probes used to locate
  the obstruction; no float result is used as a certificate.

No large caches, installs, model weights, or builds were created. Total size is
well under 1 MB; `__pycache__` is excluded by the repository `.gitignore`.
