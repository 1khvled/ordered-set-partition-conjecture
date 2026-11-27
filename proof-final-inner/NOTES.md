# Inner real-rootedness (H) — structural identities proved, and the exact gap

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-inner` was created or edited.  
Outcome: **(H) is not proved here.** Several new all-index structural identities are
proved below (in particular a *positive* binomial recurrence for the inner array),
and the precise point at which positivity stops being sufficient is isolated.
No PF/Favard route is repeated as a proof attempt.

---

## 0. Setting and conventions

Fix `j >= 2`. Write

```
A_j(z) = sum_{r=1}^{j-1} z^r/r!,      C_j(z) = e^z - sum_{r=0}^{j-1} z^r/r!,
N(z)   = 1 - A_j(z) = 2 - E_{j-1}(z), E_{j-1}(z) = sum_{r=0}^{j-1} z^r/r!,
p(z)   = z^{j-1}/(j-1)! = E_{j-1}(z) - E_{j-2}(z),
g_j(z) = C_j(z)/N(z),
F_j(z,t) = 1/(1 - A_j(z) - t C_j(z)),
H(z,t) = 1/(1 - t g_j(z)) = sum_{n>=0} H_n(t) z^n/n!,
h_{n,m} = [t^m] H_n(t) = n! [z^n] g_j(z)^m.
```

All series are formal power series in `z` with polynomial coefficients in `t`.
We use `T = t d/dt`, `Q(q) = q + (t-1) q'`, and the convention `binom(n,k)=0`
for `k<0` or `k>n`.

Target (H): for every `j>=2, n>=0`, `H_n` has only real nonpositive roots and
`H_n \interl H_{n+1}`. (For `n>=1`, `t=0` is always a root; see Lemma 1.)

---

## 1. Six all-index identities (proved)

### Lemma 1 (basic form and vanishing range). 

`H = t C_j F_j = t \partial_t \log F_j`, and consequently

```
H_0 = 1,    H_n = 0  for  1 <= n <= j-1,    deg_t H_n = floor(n/j)  (n>=j).
```

*Proof.* Since `1-tg = (N-tC_j)/N` and `F_j=1/(N-tC_j)`,

```
H = N/(N-tC_j) = 1 + t C_j/(N-tC_j) = 1 + t C_j F_j = 1 + t \partial_t \log F_j.
```

Since `C_j` has order `j` in `z`, the `n`-th
coefficient of `t C_j F_j` (and hence of `H`) vanishes for `1<=n<=j-1`.
For `t`-degree: `g_j` has order `j`, hence `g_j^m` has order `jm`, so
`h_{n,m}=0` unless `n>=jm`, and `h_{jm,m}=(jm)!/(j!)^m>0`. The identity
`\partial_t\log F_j = C_j F_j` follows from `\partial_t F_j=C_jF_j^2`. ∎

### Lemma 2 (the special differential equation for `g_j`). 

```
(1 - A_j(z)) g_j'(z) = (2 - p(z)) g_j(z) + p(z),    p(z)=z^{j-1}/(j-1)!.   (G)
```

*Proof.* `C_j' = C_j + p` and `A_j' = 1 + A_j - p` (both are termwise
differentiation). From `g_j N = C_j` we get `g_j' N + g_j N' = C_j + p`,
and `N' = -A_j' = N - 2 + p`, hence
`g_j' N = C_j + p - g_j(N-2+p) = g_j N + p - g_j N + 2 g_j - p g_j
= (2-p) g_j + p`. ∎

The useful equivalent form is `2-p = N - N'`, so

```
g_j' - g_j = (p - g_j N')/N.
```

Even better, since `N - C_j = 2 - e^z`, one gets the polynomial-coefficient
form `N^2 g_j' = 2 C_j + p(2 - e^z)`.

### Lemma 3 (linear first-order PDE for `H`). 

```
N(z) \partial_z H  =  t [ 2 + (t-1) p(z) ] \partial_t H  +  t p(z) H,        (PDE)
equivalently  N \partial_z H = (2+(t-1)p) H^2 + (p-2) H.
```

*Proof.* Differentiate `H = 1 + t C_j F_j` and use `F_j'=(2+(t-1)p)F_j^2-F_j`
(the Riccati identity proved in `proof/PROOF_NOTES.md`, §2) together with
`C_j'=C_j+p`. The terms `tC_jF_j` cancel, leaving
`H' = t p F_j + t C_j(2+(t-1)p)F_j^2`. Since `F_j=H/N` and `tC_jF_j=H-1`,
multiplying by `N` gives
`N H' = tp H + (2+(t-1)p)(H^2-H) = (2+(t-1)p)H^2 - (2-p)H`.
Then use `t\partial_t H = H^2-H` (immediate from
`H=1/(1-tg)`) to replace `H^2` by `t\partial_t H + H`, which gives the linear
form (PDE). Both residuals were checked symbolically with SymPy for
`j=2,3,4,5` (`verify_ode.py`, residual 0).

### Lemma 4 (a POSITIVE binomial recurrence for the inner array). 

For every `n>=1` and `m>=1`,

```
h_{n,m} = sum_{k=1}^{j-1} binom(n,k) h_{n-k,m}
        + sum_{q=j}^{n}   binom(n,q) h_{n-q,m-1},                         (R3)
```

with `h_{0,0}=1` and `h_{s,m}=0` for `s<0` or `m<0`. Equivalently, summing
over `m>=1` (and using `H_0=1`, `H_s=0` for `1<=s<j`),

```
H_n = sum_{k=1}^{j-1} binom(n,k) H_{n-k}
    + t sum_{q=j}^{n} binom(n,q) H_{n-q},      for all n >= j.             (R4)
```

*Proof.* Multiply the identity `N g_j^m = C_j g_j^{m-1}` (valid for `m>=1`) and
take `n![z^n]`. Since `N = 1 - sum_{k=1}^{j-1} z^k/k!`, the left side gives
`h_{n,m} - sum_{k=1}^{j-1} binom(n,k) h_{n-k,m}`; since `C_j = sum_{q>=j} z^q/q!`,
the right side gives `sum_{q=j}^n binom(n,q) h_{n-q,m-1}`. This proves (R3).
Summation over `m>=1` yields (R4) for `n>=j` (for `n<j` the `k=n` term would
need the `m>=1` part of `H_0`, which is zero, giving `H_n=0`).
Exact verification for `j=2,...,5`, `n<=18` is in `verify_identities.py`
(0 errors).

### Lemma 5 (decorated associated-Bell interpretation). 

For `n>=1`, `H_n(t)/n!` is the weighted count of **ordered set partitions of
`[n]` whose first block has size at least `j`**, with weight `t` per block of
size `>= j` and weight `1` per block of size `< j`. Equivalently, via the
"chapter" decomposition (each chapter = one big block followed by a string of
small blocks), `h_{n,m}` counts ordered partitions of `[n]` into `m`
chapters, the `m` chapters being ordered; `h_{n,m}/m!` counts the corresponding
*unordered* partitions into `m` chapters.

*Proof.* `H = t C_j F_j`: choose the first block (a single block of size `>=j`,
weight `t`, EGF `tC_j`) and then an arbitrary ordered sequence of blocks
(EGF `F_j`). The chapter form is the same identity read as
`H = sum_{m>=1} t^m g_j^m`, with `g_j = C_j/(1-A_j)` the EGF of a chapter. ∎

### Lemma 6 (the signed recurrence — where positivity is lost). 

Extracting `h_{n,m}` from the linear PDE of Lemma 3 gives, for `n>=1`,

```
h_{n+1,m} = sum_{k=1}^{j-1} binom(n,k) h_{n-k+1,m} + 2m h_{n,m}
          - binom(n,j-1) m ( h_{n-j+1,m} - h_{n-j+1,m-1} ).              (R2)
```

For `j=2` this is `h_{n+1,m}=(n+2m)h_{n,m} - n m (h_{n-1,m}-h_{n-1,m-1})`.
Equivalently, at the polynomial level,

```
H_{n+1} = sum_{k=1}^{j-2} binom(n,k) H_{n-k+1}
        + binom(n,j-1)[ H_{n-j+2} + t Q(H_{n-j+1}) ] + 2T H_n.           (R2')
```

*Proof.* Direct coefficient extraction from (PDE); the `t^2\partial_t` piece
combines with `-T` to give `t Q`. Checked against the exact arrays for
`j=2,...,5`.

---

## 2. What the identities achieve, and what they do not

1. **The inner array has an all-index positive recurrence.** Lemma 4 is new
   relative to `proof-final-structural/NOTES.md`; it is the exact analogue of
   the main family's first-block recurrence, but with the different initial
   data `H_0=1, H_1=...=H_{j-1}=0`. In particular

   ```
   H_n = O_{n,j} - sum_{k=1}^{j-1} binom(n,k) O_{n-k,j}
       = t sum_{q=j}^{n} binom(n,q) O_{n-q,j}.
   ```

2. **Positivity of (R3)/(R4) is not a real-rootedness certificate.** The
   effective kernel of the renewal is `g_j/z^j`, which is not PF
   (`proof-final-structural/NOTES.md`, §2.1), and neither is `1/(1-A_j)`
   for `j>=3`. The general theorem "positive binomial convolution
   recurrence ⇒ PF rows" is false: a kernel with nonnegative coefficients
   need not be a multiplier/PF sequence. The explicit counterexample in
   `proof/PROOF_NOTES.md`, §7 (kernel `t^2 z`, yielding `O_2=1+t^2` with
   complex roots) applies verbatim here. So (R3) cannot by itself
   establish (H).

3. **The linear PDE (PDE) makes the obstruction explicit.** In Lemma 3 the
   coefficient `2+(t-1)p(z)` is responsible for the negative term
   `-binom(n,j-1)m h_{n-j+1,m}` in (R2). One cannot remove it: it is forced
   by `C_j' = C_j + p`, i.e. by the *truncation* of the exponential at order
   `j-1`. For `j=1` there is no truncation (`p=0, N=1`) and the inner family
   degenerates to the Fubini family `1/(1-t(e^z-1))`, real-rooted by
   Comtet's recurrence `O_{n+1}=t(1+D)O_n`; for every `j>=2` the `p`-term is
   present and blocks that argument.

4. **No standard theorem applies.** Brändén–Leite Theorem 4.4 requires a PF
   seed (`g_j` is not PF); Brändén–Saud Maia Leite Theorem 4.5 requires a
   monomial `t x^r` (we have the infinite tail `tC_j`); Chen–Fu–Ruan's
   left-production criterion needs a TN production matrix, which fails for
   every `j` (`proof-final-network/NOTES.md`, §2); the Ahlberg–Alexandersson
   weighted-renewal theorem needs a real-rooted polynomial seed; Simion's
   composition theorem (`Sim84`) covers unweighted compositions of a multiset,
   whereas the chapters here carry the non-PF, non-CZDS weight sequence
   `a_n=n![z^n]g_j` (`1,4,17,86,517,...` for `j=2`); Alexandersson's
   ordered-lift theorem concerns `e^z/(1-tC_k)` and a different multiplier
   window.

5. **Favard fails.** The subsequence with fixed excess `s=n-jm` does not
   satisfy a three-term recurrence with positive coefficients. For `j=2,s=0`
   one does have the closed product
   `h_{2(m+1),m+1}/h_{2m,m}=(m+1)(2m+3)`, but for `s=1` the ratios are
   `20, 63/2, 48, ...`, which are not of Favard form; the earlier
   `recurrence_test.py` probes reached the same obstruction for the main rows.

6. **The complex-zero-decreasing (CZDS) route on the natural chapter weights
   is blocked.** Let `a_n = n![z^n]g_j` (so `a_0=...=a_{j-1}=0, a_j=1`, and
   `a_n>0` for `n>=j`), and let `\Phi(z)=sum_n a_n z^n/n! = g_j(z)`. By the
   Pólya–Schur characterization, a nonnegative sequence is a multiplier
   sequence iff `\Phi` (after removing any initial string of zeros, i.e.
   dividing by a power of `z`) lies in the Laguerre–Pólya class, i.e. is
   entire with only real zeros. But `g_j = C_j/N` has poles at every zero of
   `N(z)=1-A_j(z)` (e.g. `z=1` for `j=2`), hence is not entire. Therefore
   `(a_n)` is not a multiplier sequence, and since every CZDS with nonnegative
   terms is a multiplier sequence (`Craven–Csordás`, Theorem 1.4 and the
   subsequent remark), `(a_n)` is **not** a CZDS. So the "decorated chapter
   weights are CZDS" idea cannot work in the standard basis.

---

## 3. The precise gap

Let `h_{n,m}` be the integer array of Lemma 4. The following statement is
exactly what is missing; (H) is equivalent to it by the Aissen–Schoenberg–
Whitney criterion applied to each row.

> **Missing Lemma (H-array).** For every `j>=2` and `n>=0`, the finite positive
> sequence `(h_{n,0},h_{n,1},...,h_{n,floor(n/j)})` is a Pólya-frequency
> sequence, i.e. every minor of its lower Toeplitz matrix is nonnegative.
> Equivalently, `H_n(t)=sum_m h_{n,m}t^m` has only real nonpositive roots.
> Equivalently, `g_j` is a *renewal-PF* kernel: for every `n`,
> `n![z^n](1-tg_j)^{-1}` is real-rooted, even though `g_j/z^j` is not PF.

What is proved above reduces this to a statement about one kernel:

* By (R3), the row sequence `h_{n,\bullet}` is obtained from `h_{n-1,\bullet}`
  by a positive binomial convolution with the two kernels
  `K_1(k)=[1<=k<j]` and `K_2(k)=[k>=j]`. `K_1+K_2` is the all-ones kernel,
  whose own renewal is the Fubini family (PF); the *split* at `j` is exactly
  what destroys PF of any single kernel description and is not covered by
  convolution closure of PF sequences.

* By (PDE), any proof by rootwise induction must control the signed term
  `-binom(n,j-1)m(h_{n-j+1,m}-h_{n-j+1,m-1})` in (R2). This is the inner
  analogue of the weighted root-separation inequality (WRS) in
  `proof-final-wrs/NOTES.md`, but now at the level of the *coefficient array*
  rather than the polynomial roots. No such inequality is established here.

The exact step that fails in every attempted closure is:

```
positive coefficients in (R3)  =/=>  PF of the row  =/=>  real-rootedness of H_n
```

Because the PF/Favard hypotheses are false, a successful proof needs one of:
(a) a new planar network/LGV model whose boundary minors are the *Toeplitz*
minors of the rows of `(h_{n,m})` (note the network of `proof-final-network`,
§1 gives the reflected matrix and cannot do this); (b) a sector/zero-free-region
theorem for the *reciprocal renewal* `1/(1-tg_j)` that does not require `g_j`
itself to be PF or CZDS (the ordinary multiplier/CZDS routes are excluded by
item 6 above); or (c) the array-level weighted root-separation inequality just
described. No theorem of type (b) applicable to `1/(1-tg_j)` was located or
proved here.

---

## 4. Finite verification (evidence only, not used above)

* `verify_ode.py` — exact symbolic residuals for Lemma 2, Lemma 3 and the
  Riccati form, `j=2,3,4,5`; all zero.
* `verify_identities.py` — confirms (R3) and (R4) with integer/rational
  arithmetic for `j=2,...,5`, `n<=18`; 0 errors.
* `explore_H.py` — prints the exact arrays `h_{n,m}/n!` for `j=2,3,4`, `n<=15`.

These support the identities of §1; they do not constitute a proof of (H).
The real-rootedness probes in `proof-final-structural/H_check.py` and
`H_disc.py` (`j<=5`, `r<=14`) remain finite evidence.

## 5. Storage / cache report

No large caches were created. No model weights, no Mathlib builds, no
database files. New files: this `NOTES.md` plus three small Python scripts
(`explore_H.py`, `verify_identities.py`, `verify_ode.py`), total well under
1 MB. Python `__pycache__` is excluded by the repository `.gitignore`.
