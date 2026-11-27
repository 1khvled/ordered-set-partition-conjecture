# Transform route audit (Alexandersson ordered big-block conjecture)

Date: 2026-11-27
Scope: only `ordered-big-blocks/proof-final-transform-audit` was created/edited.

## 0. Audit verdict

- **No proof of (T) or (H) is obtained here.**
- The quarantined `proof-final-transform/NOTES.md` is **not reliable as written**.
  Its quarantine banner is also inaccurate: it says `H_{n,j}(0)=1` for every `n`,
  but the correct value is `H_{0,j}(0)=1` and `H_{n,j}(0)=0` for `n>=1`
  (equivalently `t | H_{n,j}` for `n>=1`). The note's body claim that
  higher-degree `H_{n,j}` have root `0` is correct.
- Recurrence (T1) is correct. The decomposition (T0) is correct.
- The "sliding window is a common interlacing family" induction is **not proved**.
  Finite checks support it in small ranges, but the all-index step is exactly an
  open aligned common-interlacing / root-separation statement. This is not a
  pairwise-compatibility issue only: pairwise interlacing is insufficient for
  a size >=3 family, and the H-row/O-row boundary is a different-family input.
- The precise remaining obstruction is stated in Section 6.

All computations below were done with exact SymPy rationals / polynomial
arithmetic in `verify_basics.py`, `check_window.py`, `inspect_roots.py`
(evidence only; not used as proof).

## 1. Definitions and checked identities

Fix `j>=2`. Let

```text
A_j(z)=sum_{r=1}^{j-1} z^r/r!,              N_j(z)=1-A_j(z),
C_j(z)=sum_{r>=j} z^r/r!=e^z-sum_{r=0}^{j-1} z^r/r!,
D_j(z)=1/N_j(z),                             g_j(z)=C_j(z)/N_j(z).
```

The ordered big-block EGF is

```text
F_j(z,t)=1/(1-A_j(z)-t C_j(z)) = D_j(z)/(1-t g_j(z)).                (T0)
```

Define the inner and outer rows by

```text
H_j(z,t)=1/(1-t g_j(z))=sum_{n>=0} h_{n,j}(t) z^n,  H_{n,j}(t)=n! h_{n,j}(t),
O_j(z,t)=F_j(z,t)=sum_{n>=0} o_{n,j}(t) z^n,         O_{n,j}(t)=n! o_{n,j}(t).
```

Since `D_j*H_j=F_j`, the convolution is

```text
O_{n,j}(t)=n! sum_{r=0}^{n} d_{n-r}^{(j)} H_{r,j}(t)/r!,  d_k=[z^k]D_j. (C)
```

### 1.1 Basic coefficient facts (proved)

1. `g_j(0)=0` because `C_j` has order `j`. Therefore `H_j(z,0)=1` and
   `H_{0,j}(0)=1`, `H_{n,j}(0)=0` for `n>=1`.
2. `g_j(z)^m` has order `jm`, hence `H_{n,j}=0` for `1<=n<=j-1`, and
   `deg_t H_{n,j}=floor(n/j)` for `n>=j`; moreover `t | H_{n,j}` for `n>=1`.
3. `D_j` has nonnegative coefficients and `d_0=1`, `d_n=sum_{r=1}^{min(n,j-1)} d_{n-r}/r!>0`.
4. `O_{n,j}(t)=n![z^n] sum_{m>=0} t^m C_j^m D_j^{m+1}` has nonnegative
   coefficients, degree `floor(n/j)`, and `O_{n,j}(0)=n! d_n>0`.

Proofs: (1) is immediate from the series definition. (2) is the order of
`g_j^m`. (3) is the recurrence from `D_j=1+A_j D_j`. (4) follows from (C),
nonnegativity of `C_j,D_j`, and the fact that `C_j^m` has order `jm`.
The exact checks in `verify_basics.py` agree for `j=2,3,4`, `n<=12`.

## 2. Recurrence T1 (proved)

From (T0), `(1-A_j)F_j=H_j`. Taking `n![z^n]` gives

```text
O_{n,j}(t)=H_{n,j}(t)+sum_{r=1}^{j-1} binom(n,r) O_{n-r,j}(t),      (T1)
```

with the convention `O_k=0` for `k<0`.

Proof: `[z^n](1-A_j)F_j = o_{n,j} - sum_{r=1}^{min(n,j-1)} o_{n-r,j}/r!`.
Multiplying by `n!` gives `O_{n,j}-sum_{r=1}^{j-1} binom(n,r)O_{n-r,j}`.
This is exactly `[z^n]H_j` multiplied by `n!`, i.e. `H_{n,j}`.  QED.

The same identity is verified exactly in `verify_basics.py` for
`j=2,3,4`, `n<=12`.
## 3. What the quarantined transform note got wrong

The file `proof-final-transform/NOTES.md` is quarantined; that quarantine is
correct. The specific issues are:

1. **Root-at-zero error (banner is wrong; body claim is right).** Its
   quarantine banner says `H_{n,j}(0)=1` for every `n`; this is false for
   `n>=1`. The correct value is `H_{n,j}(0)=0` for `n>=1`, with `H_{0,j}=1`.
   The body's remark that higher-degree `H_{n,j}` have root `0` is correct.
2. **Unproved interlacing-family closure.** Lemma I asserts `F_n` is an
   interlacing family. The note gives only pairwise/finite random-combination
   evidence. Pairwise interlacing does not imply common interlacing for
   families of size >=3.
3. **Unproved consecutive-row step.** Lemma II assumes every member of
   `F_{n+1}` interlaces `P_n` in an aligned way and then invokes a linear-pencil
   fact. Even the linear-pencil fact needs correct degree/order hypotheses;
   the main gap is the aligned interlacing input, which is not proved.
4. **Degree bookkeeping omitted.** Members of the window have degrees
   `q=floor(n/j)` or `q-1`. Any interlacing argument must track this split.

The note should not be cited as a proof or as a conditional proof with a
correct conditional hypothesis.

## 4. Sliding-window interlacing: exact formulation

Let `q=floor(n/j)` and `s=n mod j`, `n>=j`. Define the window

```text
W_n = { H_{n,j} } U { B_{n,r}(t) := binom(n,r) O_{n-r,j}(t) : 1<=r<=j-1 }.
```

Degrees: `deg H_{n,j}=q`, `deg B_{n,r}=q` for `1<=r<=s`, and
`deg B_{n,r}=q-1` for `s<r<=j-1`. All nonzero members have nonnegative
coefficients and therefore all real roots (if any) are <=0.

A rigorous common-interlacing hypothesis for `W_n` is:

> **(CI_n)** There exists a real-rooted polynomial `G_n`, `deg G_n=q+1`,
> with simple real roots, such that every nonzero member of `W_n` interlaces
> `G_n` in the standard mixed-degree sense.

For equal-degree real-rooted families, (CI_n) is equivalent to: every positive
linear combination is real-rooted (the Chudnovsky-Seymour/Fell/Dedieu common
interlacer criterion, used in the MSS interlacing-families setting). For mixed
degrees `q,q-1`, one must use the common interlacer definition directly.

### 4.1 What (CI_n) would give

If (CI_n) holds and every member of `W_n` is real-rooted, then every
nonnegative combination

```text
P = alpha H_{n,j} + sum_{r=1}^{j-1} beta_r B_{n,r},   alpha,beta_r>=0,
```

is real-rooted. In particular (T1) is the all-ones combination, so `O_{n,j}`
is real-rooted.

This implication is the easy direction. It does **not** by itself prove
`O_{n-1,j} interlaces O_{n,j}`. Consecutive-row interlacing needs an aligned
boundary argument, not just membership in a common-interlacing family.

### 4.2 Why pairwise checks are insufficient

For >=3 polynomials, pairwise interlacing does not imply a common interlacer.
A size-3 family can be pairwise interlacing but have no common interlacer when
the root orderings are incompatible. Therefore finite pairwise checks
(including those in `check_window.py`) are evidence only and cannot close (T).

### 4.3 Sliding-window induction attempt and its failure point

A natural induction is:

- Induction hypothesis at `n`: `O_0,...,O_{n-1}` real-rooted, and `W_n`
  satisfies (CI_n).
- Step: (T1) gives `O_n` real-rooted. Then prove `W_{n+1}` satisfies
  (CI_{n+1}) and, if desired, `O_{n-1} interlaces O_n`.

The new window is

```text
W_{n+1}={H_{n+1,j}} U {binom(n+1,1)O_n, binom(n+1,2)O_{n-1}, ..., binom(n+1,j-1)O_{n-j+2}}.
```

The carried-over tail `O_{n-1},...,O_{n-j+2}` is old, but the boundary has
two new members: `O_n` (a positive combination of `W_n`) and `H_{n+1,j}`
(a different family).

The step requires all of the following, none proved here:

- **(A) O_n aligns with the carried-over tail.** `O_n` must interlace
  `O_{n-1},...,O_{n-j+2}` with the same orientation required by `W_{n+1}`.
  This is stronger than "O_n is real-rooted".
- **(B) H_{n+1,j} aligns with the whole O-window.** This couples the open
  inner family (H) to the outer family and requires a cross-family
  root-separation inequality.
- **(C) The degree shift is controlled.** When `n+1` crosses a multiple of
  `j`, degrees in the window change from `q,q-1` to `q+1,q`.
- **(D) Common interlacer, not just pairwise.** Even if every pair in
  `W_{n+1}` interlaces, a common interlacer for the whole window is still
  needed.

Thus the induction does not close. The failure is not a missing pairwise check;
it is a missing aligned common-interlacer propagation theorem.

## 5. Rigorous sign/root observations

These are exact consequences of the definitions; they constrain any proof.

1. **Nonnegative coefficients / nonpositive roots.** `O_{n,j}, H_{n,j}` have
   nonnegative t-coefficients. Any real root is <=0.
2. **Boundary values.** `O_{n,j}(0)=n!d_n>0`; `H_{0,j}=1`; `H_{n,j}(0)=0`
   for `n>=1`.
3. **Top-root separation at 0.** For `n>=j`, `H_n` has root 0 and `O_n` has
   no root at 0. Same-degree interlacing is possible only if O_n's largest
   root lies strictly below 0 and between H_n's two largest roots. This is
   observed in `inspect_roots.py` (e.g. j=2,n=4: H roots [-2.833,0], O_3
   degree 1 root [-0.857]; j=2,n=5: H [-1.075,0], O_4 [-6.922,-0.578]).
4. **No H/O cancellation at 0.** Since O_n(0)>0 and H_n(0)=0, T1 implies
   O_n(0)=sum_{r>=1} binom(n,r)O_{n-r}(0), consistent with (3).
5. **j=2 is not trivial.** For j=2, `O_n=H_n+n O_{n-1}` and `D_2=1/(1-z)`.
   Real-rootedness of O_n still requires a compatibility argument between H_n
   and O_{n-1}; sum of two real-rooted nonnegative-coefficient polynomials is
   not real-rooted in general.

## 6. Exact obstruction

The transform route cannot be closed by the current arguments. The exact open
statement needed by the sliding-window induction is:

> **Gap T-AUDIT (aligned common interlacer for the T1 window).** Fix `j>=2`.
> Assume (H): every `H_{n,j}` is real-rooted and `t|H_{n,j}` for `n>=1`.
> Prove or disprove: for every `n>=j`, the window
> `W_n={H_{n,j}} U {binom(n,r)O_{n-r,j}:1<=r<=j-1}` has a common interlacer
> `G_n`, and these interlacers can be chosen so that `O_n` from (T1) is
> aligned as a carried-over member of `W_{n+1}`.

Equivalent missing content: an all-index root-separation/common-interlacer
propagation theorem for the non-PF T1 transform. It is stronger than pairwise
interlacing, stronger than finite random-combination tests, and not implied by
the known PF/TN obstructions: `g_j` is not PF, `D_j` is not PF for `j>=3`
(`proof-final-structural/NOTES.md`), and the natural network/production routes
are obstructed (`proof-final-network/NOTES.md`).

## 7. What would count as progress

A real proof would need one of the following, proved for all `j,n`:

1. A direct common-interlacer construction for `W_n` with an induction that
   carries alignment across `n -> n+1`.
2. A Wronskian/sign-variation argument proving every positive combination in
   `W_n` is real-rooted, including the `H_n` boundary at 0, without
   pairwise-only handwaves.
3. A different transform theorem that turns (H) into (T) and applies to this
   non-PF kernel.

None of these is established in the repository notes audited here.

## 8. Files in this audit

- `NOTES.md` - this audit.
- `verify_basics.py` - exact SymPy checks for definitions, H_n(0), degrees, T1.
- `check_window.py` - finite pairwise/common-combination evidence only.
- `inspect_roots.py` - root inspection for sign/top-root observations.

No large installs, builds, or caches were created.
