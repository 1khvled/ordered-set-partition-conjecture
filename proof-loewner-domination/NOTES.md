# Loewner domination for the falling-basis Bezout decomposition (j=2)

Date: 2026-11-28
Scope: only `ordered-big-blocks/proof-loewner-domination` created/edited.
No installs/builds. No other agents/models.

## Outcome

**No all-`n` proof of general-degree Loewner domination was obtained. No exact
obstruction in the actual consecutive rows `O_{n,2},O_{n+1,2}` was found.**

Finite exact/numeric evidence (audited below) confirms the domination
`M(O_{n,2},O_{n+1,2}) \succ 0` for `n<=21`, but the normalized margin
`lambda_min(M)/||M||_max` decays rapidly. The generic-TP2 obstruction from
`proof-falling-to-realroots/NOTES.md` remains the exact obstruction to any
proof using only falling-basis TP2; the proved inputs (H), monomial TP2, and
falling-basis TP2 do not by themselves supply a domination certificate for
the actual rows.

---

## 1. Setup and exact decomposition

Let `C(z)=e^z-1-z`, `F(z,t)=1/(1-z-tC(z))`, `O_n(t)=n![z^n]F(z,t)`.
Write the falling-factorial basis `f_i(t)=(t)_i`, `f_0=1`,
`f_{i+1}=(t-i)f_i`, and

```
P=O_{n,2}=sum_i c_{n,i} f_i,  Q=O_{n+1,2}=sum_i c_{n+1,i} f_i,
Delta^{(n)}_{ij}=c_{n,i}c_{n+1,j}-c_{n,j}c_{n+1,i},   i<j.
```

The root-positive Bezout kernel (proof-falling-to-realroots (1),(7)) gives

```
M(P,Q) = sum_{0<=i<j<=D} Delta^{(n)}_{ij} K_{ij},        (1)
K_{ij}(x,y) = (f_j(x)f_i(y)-f_j(y)f_i(x))/(x-y),
```

with `D=max(deg P,deg Q)=floor((n+1)/2)`. The adjacent kernels are PSD rank-1:

```
K_{i,i+1}(x,y)=f_i(x)f_i(y),   K_{i,i+1}=v_i v_i^T.     (2)
```

For `r=j-i>=2`, `K_{i,i+r}` is indefinite (proof-falling-to-realroots ��3).

Split (1) into adjacent and higher-lag parts:

```
M = A_adj + H_lag,
A_adj = sum_{i=0}^{D-1} Delta_{i,i+1} K_{i,i+1}  (PSD),
H_lag = sum_{j-i>=2} Delta_{ij} K_{ij}          (indefinite).
```

The Loewner-domination gap is exactly

```
A_adj + H_lag \succ 0   <=>   M(P,Q) \succ 0.            (3)
```

Equivalently `A_adj \succ -H_lag` in Loewner order. This is the missing
input named in the task; it is equivalent to the consecutive Bezout/Interlacing
target (proof-falling-to-realroots (5)) once real-rootedness/simplicity is
supplied.

---

## 2. Proved inputs used

- Falling-basis TP2: `Delta^{(n)}_{ij}>=0` for all `i<j`, all `n`
  (`proof-j2-tp2-transfer/NOTES.md`, `proof-j2-falling-bezout/NOTES.md`).
- (H) for associated Stirling: proved in `proof-u-operator/NOTES.md` ��7.
- Monomial TP2: proved in `proof-j2-tp2-transfer/NOTES.md`.
- Exact Bezout decomposition (1),(2): proved in
  `proof-falling-to-realroots/NOTES.md` ����1-2.
- Generic obstruction: TP2 weights alone do not imply (3); degree-2 witness
  (14)--(18) in `proof-falling-to-realroots/NOTES.md`.

No new all-`n` theorem is proved from these inputs here.

---

## 3. Degree 2 (exact)

For `D=2` (pairs with `deg<=2`), (1) gives

```
M = [[Delta_01-Delta_02, Delta_02],
     [Delta_02,           Delta_12]],
det M = Delta_01*Delta_12 - Delta_02*(Delta_12+Delta_02).   (4)
```

Thus `M\succ0` iff `Delta_12>0` and (4)>0. This is the same as
proof-falling-to-realroots (11)--(13). For actual rows:

- `n=3` (O3,O4, D=2): `Delta=(138,36,42)`, `det=2988>0`.
- `n=4` (O4,O5, D=2): `Delta=(3984,1920,3084)`, `det=2678976>0`.

Exact residuals are checked in `verify.py`.

---

## 4. Finite evidence for actual rows (not a proof)

`verify.py` computes exact falling-basis coefficients from the associated
Stirling formula, forms `K_ij` by the exact Bezoutian kernel formula, checks
the decomposition residual against the direct monomial Bezoutian, and reports
normalized eigenvalues `(A+A^T)/(2||A||_max)`.

Identity audits (exact, Fraction): `K_{i,i+1}=v_i v_i^T` residual 0 for
`i=0..4`; `K_{0,2}=[[-1,1],[1,0]]`; decomposition residual 0 for `n=2..11`.

Domination finite evidence (eigenvalues of normalized symmetric part):

| n | D | min eig M (norm) | min eig A_adj (norm) | min eig H_lag (norm) |
|---:|---:|---:|---:|---:|
| 2 | 1 | 1.000e+00 | 1.000e+00 | 0.000e+00 |
| 3 | 2 | 2.465e-01 | 3.043e-01 | -1.618e+00 |
| 4 | 2 | 1.905e-01 | 7.741e-01 | -1.618e+00 |
| 5 | 3 | 2.495e-02 | 3.042e-02 | -1.884e+00 |
| 6 | 3 | 1.421e-02 | 6.806e-02 | -1.170e+00 |
| 7 | 4 | 2.693e-03 | 1.076e-03 | -1.128e+00 |
| 8 | 4 | 1.569e-03 | 3.833e-03 | -1.166e+00 |
| 9 | 5 | 2.429e-04 | 3.414e-05 | -1.280e+00 |
| 10 | 5 | 7.760e-05 | 1.575e-04 | -1.504e+00 |
| 11 | 6 | 1.761e-05 | 8.050e-07 | -1.854e+00 |
| 12 | 6 | 5.020e-06 | 3.107e-06 | -1.700e+00 |
| 13 | 7 | 1.214e-06 | 7.407e-09 | -1.621e+00 |
| 14 | 7 | 2.464e-07 | 3.370e-08 | -1.590e+00 |
| 15 | 8 | 5.395e-08 | 5.501e-11 | -1.590e+00 |
| 16 | 8 | 1.256e-08 | 2.829e-10 | -1.609e+00 |
| 17 | 9 | 3.076e-09 | 3.267e-13 | -1.643e+00 |
| 18 | 9 | 6.214e-10 | 1.856e-12 | -1.686e+00 |
| 19 | 10 | 1.274e-10 | 1.565e-15 | -1.738e+00 |
| 20 | 10 | 2.760e-11 | 9.689e-15 | -1.795e+00 |
| 21 | 11 | 6.263e-12 | 6.139e-18 | -1.858e+00 |

All checked `M` are positive definite (finite evidence). The normalized
margin decays rapidly; the adjacent component also becomes nearly singular
relative to the indefinite lag component. This is evidence of a delicate
domination, not a uniform certificate.

---

## 5. Structural gap

### 5.1 Why TP2 alone is insufficient

`proof-falling-to-realroots/NOTES.md` ��5 gives a strict feasible witness
satisfying monomial TP2, falling-basis TP2, row real-rootedness, row PF,
positive coefficients, and no common root, but with `det M = -30 < 0`.
Therefore no argument using only those constraints can prove (3).

### 5.2 What would be needed

A sufficient all-`n` certificate is any of:

1. A direct proof that `M(O_{n,2},O_{n+1,2})\succ0` (equivalent to (3)).
2. A weight inequality on `Delta^{(n)}_{ij}` implying
   `sum_{j-i>=2} Delta_{ij}K_{ij} \preceq sum_i Delta_{i,i+1}v_iv_i^T`
   in Loewner order.
3. A basis/kernel transform making the actual weighted sum a nonnegative
   combination of PSD kernels.

I did not obtain such a certificate from (H), monomial TP2, falling TP2, or
the O_{n,2} recurrence. In particular, (H) is an associated-Stirling ratio
inequality; I found no proved Loewner/positive Cauchy--Binet transfer of (H)
to the indefinite falling-basis Bezout kernels `K_{ij}`.

### 5.3 Obstruction status

No exact obstruction in the actual O_{n,2},O_{n+1,2} rows was found up to
`n=21`. The rapidly shrinking normalized margin means finite evidence is not
strong evidence for large `n`; it is not an obstruction.

---

## 6. Proof-policy audits

### 6.1 Direction audit

Domination is stated as `A_adj + H_lag = M \succ 0` (equivalently
`A_adj \succ -H_lag`). Degree-2 sample `n=4`: `det M = 2678976 > 0`.
Generic counterexample direction: `det M = -30 < 0` (feasible for TP2-only).

### 6.2 Monotonicity audit

No monotonicity claim is used as a proof step. The normalized margin table is
reported as finite evidence only; no monotonicity of that margin is asserted.

### 6.3 Feasible-constraint falsification

The generic TP2-only constraints are falsified by the proof-falling-to-realroots
witness (14)--(18). For the actual O_{n,2} constraints, I attempted no
falsification beyond checking `n<=21`; no witness was found. This is not a
proof.

### 6.4 Boundary/degeneracy audit

`n=0,1`: no Bezoutian/interlacing target needed. `n=2`: D=1, M=[8]>0.
`D=2`: exact formula (4). Zero padding for unequal degree rows is included.
`K_{ii}=0`, `K_{ji}=-K_{ij}`: used in decomposition. All sums finite.

### 6.5 Symbolic residuals

`verify.py` prints exact zero residuals for K_{i,i+1}-v_iv_i^T (i=0..4) and
for the decomposition (1) vs direct monomial Bezoutian (n=2..11). Finite
eigenvalue computations are finite evidence only.

### 6.6 Retraction rule

An initial scratch inspection used an incorrect symmetrization and doubled some K entries; it was discarded before inclusion here. An arithmetic error in the n=3 degree-2 det sample (4536 written instead of 2988) was caught during audit and corrected. No invalid theorem is recorded as proved.

---

## 7. Files

- `NOTES.md` �?" this record.
- `verify.py` �?" exact identity/decomposition residuals and finite domination
  evidence.
- `explore.py`, `explore2.py`, `explore3.py`, `inspect_n5.py` �?" scratch
  probes (not proof).

## 8. Status

- General-degree Loewner domination for actual O_{n,2},O_{n+1,2}: **NOT PROVED**.
- Exact obstruction in actual O_{n,2},O_{n+1,2}: **NOT FOUND** (checked n<=21).
- Equivalent target (consecutive Bezout positivity / interlacing for j=2):
  **still open**.


