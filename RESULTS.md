# Result summary

## Full conjecture

**Still open.** No counterexample was found, and no all-parameter proof was
obtained. Finite computations are not presented as a proof.

## Proved identities and reductions

- Correct first-block convolution with `binom(n,r)`.
- Riccati equation and universal differential recurrence.
- Weighted root-separation theorem: WRS implies real-rootedness and
  interlacing with the lagged row.
- Exact obstruction to the naive compatible-window implication.
- Riordan/shear factorization and an exact counterexample to general shear
  preservation of total positivity.
- Bézoutian decomposition and an exact obstruction to positivity from a
  compatible window alone.
- Exact consecutive sign condition: `O_n interlaces O_{n+1}` iff
  `O_{n+1}(mu)/O_n'(mu) < 0` at every simple root `mu` of `O_n`.
- Exact lag-`j` Wronskian target and full-history decomposition of the older tail.
- Inner/outer decomposition, positive coefficient-minor recurrence, threshold
  renewal identities, and bulk saddle-point characterization.
- Exact obstructions for linear/Möbius threshold maps, generic nonlinear renewal,
  crossing involutions, edge-only asymptotics, and naive common-interlacer
  propagation.
- Unconditional proof of the 2-associated-Stirling inequality (H) via a coupled
  interlacing induction for the factorial-weighted rows `proof-u-operator/`.
  Finite exact checks confirm the coupled invariants and (H) through `r=30`.

See `proof/`, `proof-weighted/`, and `proof-matrix/`.

## Audited computations

| Range | Rows | Pairs | Exact Sturm rows | Failures |
|---|---:|---:|---:|---:|
| `j=3..8`, `0 <= n <= 45` | 276 | 270 | 276 | 0 |
| `j=5`, `141 <= n <= 800` | 660 | 660 | 43 | 0 |
| `j=6..60`, `7j <= n <= 12j-1` | 9,075 | 9,075 | 9,075 | 0 |
| lag-Pick Wronskians, `j=2..9`, `n <= min(5j,34)` | 1,026 | — | 1,026 | 0 |
| cross-family lag-Pick, `j=2..8`, `n <= min(5j,30)` | 732 | — | 732 | 0 |
| boundary `W(H_n,O_{n-d})`, `j=2..7`, `n <= min(5j,30)` | 548 | — | 548 | 0 |
| `j=61..92`, `2j <= n <= 5j` (checkpoint, sweep in progress) | 7,376 | 7,344 | 7,376 | 0 |

Independent compact reruns reproduced all 660 `j=5` rows/pairs and all 9,075 wide-range rows/pairs.
The new `j=61..92` figures are a partial checkpoint, not a completed range.
Partial `j=3,4` scans were not audited and are not claimed.

Exact finite checks in `proof-associated-stirling/verify.py` also confirm the
reduction identities and the associated-Stirling ratio inequality through
`m <= 300`, but these are finite evidence only, not an all-index proof.

## Precise remaining target

For consecutive rows, prove the root-sign inequality

\[
\frac{O_{n+1,j}(\mu)}{O_{n,j}'(\mu)}<0
\]

at every root \(\mu\) of \(O_{n,j}\). Equivalently, prove the consecutive
residue/domination inequality in `proof-final-consecutive/NOTES.md`.

For the lag-\(j\) route, prove

\[
-\frac{nO_{n-1,j}(\mu)+2\mu O_{n-1,j}'(\mu)+S_n(\mu)}
        {O_{n-j,j}'(\mu)}
>\binom{n-1}{j-1}\mu(\mu-1),
\]

equivalently \(W(O_{n,j},O_{n-j,j})(t)>0\) on the negative axis. The full-history
identity isolates the missing term as \(W(H_{n,j},O_{n-j,j})\).

## Lean

`lake build` in `lean/` succeeds with no external dependencies. The default
target formalizes the recurrence and verifies sample rows. Experimental
Mathlib EGF/Sturm files are present but unverified and excluded from the
default target.
