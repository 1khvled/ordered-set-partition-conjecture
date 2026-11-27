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

See `proof/`, `proof-weighted/`, and `proof-matrix/`.

## Audited computations

| Range | Rows | Pairs | Exact Sturm rows | Failures |
|---|---:|---:|---:|---:|
| `j=3..8`, `0 <= n <= 45` | 276 | 270 | 276 | 0 |
| `j=5`, `141 <= n <= 800` | 660 | 660 | 43 | 0 |
| `j=6..60`, `7j <= n <= 12j-1` | 9,075 | 9,075 | 9,075 | 0 |

The compact wide rerun independently reproduced all 9,075 rows and pairs.
Additional compact `j=3,4` runs are recorded in `numerical/` when complete.

## Lean

`lake build` in `lean/` succeeds with no external dependencies. The default
target formalizes the recurrence and verifies sample rows. Experimental
Mathlib EGF/Sturm files are present but unverified and excluded from the
default target.
