# Tail attack for Alexandersson's ordered big-block conjecture

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-tail` is created or edited.

## 1. Object under study

Fix `j>=2` and write `n=jm+s`, `0<=s<j`.  The reversed normalized polynomials are

\[
Q_{m,s}(x)=\sum_{\ell=0}^{m}
\frac{c_{jm+s,m-\ell}}{c_{jm+s,m}}
\left(\frac{-x}{D_j m^j}\right)^{\ell},
\qquad
c_{n,k}=\frac{[t^k]O_{n,j}(t)}{n!}.
\]

The left-edge theorem proves that, for every fixed `K`, the first `K` zeros of
`Q_{m,s}` converge to the first `K` positive zeros of

\[
L_s(x)=s!\sum_{\ell\ge0}\frac{(-x)^\ell}{(s+j\ell)!},
\]

and that those zeros are simple.  The right-edge theorem gives a positive
fraction `eta_j>0` of the zeros nearest `0` in the `t`-variable (equivalently
the largest `x`-zeros of `Q_{m,s}`) for all large `m`.

The question for this track is whether these tail facts can be promoted to a
uniform statement for **all** `m` zeros, with an explicit effective bound,
thereby reducing the whole conjecture to the already audited finite range.

## 2. What the tail actually proves

The left-edge argument is a Hurwitz/Rouche argument on a **fixed compact**
set `|x|<=R`.  It yields, for each fixed `K` and each fixed `R` between two
consecutive zeros of `L_s`, an `M(K,R)` such that `Q_{m,s}` has exactly `K`
simple real zeros in `|x|<R` for `m>=M(K,R)`.  The right-edge argument gives
`eta_j n` real interlacing zeros near the other extreme for `m>=N_j`.

Neither statement controls the **middle** zeros:

- number of zeros in the middle band: `m - O(1) - eta_j n = Theta(m)`;
- their location is not fixed; it moves with `m`;
- no fixed compact set contains them.

Thus the tail theorems are, by construction, **edge theorems**.  They do not
touch the bulk of the zero set.

## 3. Precise obstruction: local uniform convergence does not suffice

### 3.1 Quantifier mismatch

For a fixed `K`, local uniform convergence `Q_{m,s}->L_s` on compact sets
gives convergence of the first `K` zeros.  But `K` is fixed, while the degree
is `m`.  To prove real-rootedness of `Q_{m,s}` one needs a statement for
`K=m`, i.e. for zeros that are **not** in any fixed compact set.  The tail
theorem provides no rate or bound uniform in `K=m`; its `M(K,R)` is allowed
to grow with `K`, and indeed must, because the `K`-th zero of `L_s` grows
with `K`.

### 3.2 Coefficient non-uniformity

Write `B_{m,ell}=c_{jm+s,m-ell}/c_{jm+s,m}`.  The tail limit says

\[
\frac{B_{m,\ell}}{(D_j m^j)^\ell}\frac{(s+j\ell)!}{s!}
\longrightarrow 1
\qquad (m\to\infty,\ \ell\ \text{fixed}).
\]

For `ell=alpha m` with `0<alpha<1`, this is false.  Exact computation for
`j=3` gives, for `m=20`:

| `ell` | `B_{m,ell}/((D m^j)^ell s!/(s+j ell)!)` |
|---:|---:|
| 0 | 1 |
| 1 | 1.17 |
| 5 | 12.7 |
| 10 | 1.67e3 |
| 15 | 2.54e5 |
| 19 | 1.95e6 |

(see `coeff_asympt.py`).  The relative error is exponential in `ell`, so the
middle coefficients are not described by `L_s` at all.  Consequently a
Hurwitz argument based on `L_s` cannot see, and cannot control, the middle
zeros.  The normalization `D_j m^j` is tuned to the left edge; it is not a
uniform scaling for the whole polynomial.

### 3.3 The Szegő phenomenon

Local coefficient convergence to an entire function with only real zeros is
not sufficient for eventual real-rootedness of the approximating polynomials.
The classical example is the Taylor polynomials of `e^z`: `e^z` has no zeros,
but its partial sums have complex zeros accumulating on the Szegő curve.
Our `Q_{m,s}` is a triangular array whose degree grows with `m` and whose
coefficients in the middle regime differ exponentially from the limit; the
same logical gap is present.  Finite exact checks up to `n=120` (and the
audited ranges) are evidence, but they do not close a gap that is, in
principle, invisible to fixed-compact Hurwitz arguments.

### 3.4 Consequence for WRS

The weighted root-separation condition (WRS) is a rootwise inequality at
**every** root `mu` of the lagged row `O_{n-j}`.  The tail theorems give
asymptotic locations and derivatives only for `O(1)` largest roots and for
`eta_j n` smallest roots.  The `Theta(m)` middle roots lack uniform
asymptotic locations, so WRS cannot be verified for them from the tail.  Hence
the tail does not yield an all-index WRS proof either.

## 4. What would be needed for a robust tail theorem

A proof of uniform real-rootedness/interlacing from asymptotics would require
one of the following, none of which follows from `L_s`:

1. A **uniform coefficient estimate** `B_{m,ell} = (D_j m^j)^\ell e_\ell
   (1+o(1))` valid for `0<=ell<=m` with an explicit `o(1)` independent of
   `ell`, or at least a uniform log-concavity/ultra-log-concavity bound.
2. A **global root-distribution theorem** (e.g. a saddle-point analysis of
   the bivariate generating function) showing that all roots lie on the real
   axis, not just the edges.
3. An **effective WRS inequality** for the middle roots, using a uniform bound
   on `A_n(mu)/q'(mu)`.

The left/right-edge theorems supply (1) only for `ell=O(1)` or
`ell=m-O(1)`.  The middle regime `ell=Theta(m)` needs a separate
saddle-point/Polya-frequency analysis of the kernel `H(z)` that is not
available in the current literature (the kernel is not PF, see
`literature.md` section 8).

## 5. Numerical evidence (exact, not a proof)

- Exact FLINT/Sturm scan: no non-real roots for `j=3`, `0<=n<=120`
  (`scan_j3.jsonl`).
- Audited ranges in `RESULTS.md`: `j=3..8, n<=45`; `j=5, 141<=n<=800`;
  `j=6..60, 7j<=n<=12j-1`; all zero failures.
- Coefficient asymptotic mismatch as above (`coeff_asympt.py`).

These are finite evidence.  They do not provide an all-`m` proof, and the
tail theorem does not reduce them to a finite check.

## 6. Conclusion

**The asymptotic tail to `L_s` does not suffice to prove WRS, real-rootedness,
or interlacing uniformly for all sufficiently large `m`.**  It proves only
fixed-`K` edge zeros.  The obstruction is not a missing constant but a
structural one: local uniform convergence controls neither the
`Theta(m)` middle zeros nor the coefficients in the `ell=Theta(m)` regime,
and the Szegő phenomenon shows that such control is necessary.  A robust tail
theorem would need a uniform estimate over `0<=ell<=m`, which is not provided
by `L_s` and is not currently known.

Therefore the all-`n` statement remains open; the finite audited range stays
finite evidence rather than a verified tail.
