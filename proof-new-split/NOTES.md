# A + tB split: audit of the latest higher-model output

Date: 2026-11-28. Scope: `proof-new-split` only. Finite evidence only, not a proof.

## 1. Parameter-free differential induction obstruction — CONFIRMED

With the correct recurrence, for `j=2`:

- `L_5 = 120 + 315t + 54t^2`,
- `t D O_3 = 14t^2 - t`,
- discriminant of `L_5 + λ t D O_3` is `λ^2 - 7350λ + 73305`,
- at `λ = 10` it is `-95 < 0`.

So a parameter-free Obreschkoff/compatible-window induction on the differential
recurrence really does fail at `n=5`. This agrees with the earlier weighted
root-separation obstruction in `proof-weighted/`.

## 2. Derivative-free identity for B — CONFIRMED

With

```
A_n = n O_{n-1} + sum_{r=2}^{j-1} C(n,r) O_{n-r} + C(n-1,j-1) t O_{n-j},
B_n = 2 O'_{n-1} + C(n-1,j-1)(t-1) O'_{n-j},
```

`O_n = A_n + t B_n`. The identity

```
B_n = sum_{r=j}^{n-1} C(n-1,r) (O_{n-r} + O_{n-1-r})
```

is exact. Verified for `j=2..6`, `n <= 30` with zero residual. `B_n` has
nonnegative coefficients for `n >= j+1`, so the spurious positive root of the
`D O_{n-j}` term is gone.

## 3. What actually helps

For `j=2`, exact Sturm checks show `O_{n-1} interlaces t B_n` for `n=6..25`
(no real zero of the Wronskian, negative at zero). By Obreschkoff, proving

```
O_{n-1,2} <= t B_{n,2}    for all n >= 6
```

would immediately prove `O_{n-1,2} <= O_{n,2}`. This is a genuinely cleaner,
middle-sum-free target than the previous consecutive-sign inequality.

For `j>=3`, the higher model's claim that `O_{n-1} interlaces A_n` is **false** in
this normalization: the Wronskian has real zeros starting at `j=3, n=8`, with
many further failures through `j=6, n=25`. So the split is useful for `j=2`,
not directly as an all-`j` package.

## 4. Other items

- The `v = t/(1-t)` recurrence and `h*` reformulation were not checked here;
  they are recorded for a later dedicated attempt, with the warning that the
  positive-root trap survives in the `v`-image unless the A+tB split is used.
- The J-fraction kill-test for the symmetric-matrix track was forwarded to the
  active matrix agent.

## 5. Next concrete target

Prove, for all `n >= 6`, that `O_{n-1,2}` and `t B_{n,2}` interlace:

```
2 (n-1) O'_{n-2}(mu) / O'_{n-1}(mu) < 2 / (1-mu)
```

at every root `mu < 0` of `O_{n-1,2}` (equivalently the exact residue form in
the higher-model message). This is finite-dimensional, derivative-only, and
has no uncontrolled older-tail sum.

Run `verify.py` to reproduce all exact finite checks.
