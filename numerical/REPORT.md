# Audited numerical certificate report

Status: **completed audited ranges; the full conjecture remains open.**

## Definition and independent validation

The corrected first-block recurrence is

```text
O_0 = 1
O_n = sum_{r=1}^{j-1} binom(n,r) O_{n-r}
      + t sum_{r=j}^{n} binom(n,r) O_{n-r}.
```

The block EGF vanishes at zero and `F_j(0,t)=1`. Independent checks confirmed:

- the published `j=3` rows through `n=4`;
- exact rational EGF inversion for `j=3..8, n<=60`;
- explicit ordered-partition enumeration through `n=7`;
- 5,012 Sturm controls against SymPy and an independent rational implementation;
- exact Sturm real-root counts for 276 rows (`j=3..8, n<=45`);
- 270 consecutive pairs by both Arb intervals and independent SymPy rational intervals;
- nine additional exact rational frontier-pair isolations.

## Completed new frontier

| Range | Rows | Pairs | Exact Sturm rows | Failures |
|---|---:|---:|---:|---:|
| `j=5`, `141 <= n <= 800` | 660 | 660 | 43 | 0 |
| `j=6..60`, `7j <= n <= 12j-1` | 9,075 | 9,075 | 9,075 | 0 |

An independent compact fast rerun reproduced all 660 `j=5` rows and pairs with 24 exact Sturm checks and zero failures.

The `j=5` audit includes 61,182 Newton inequalities and 122,892 certified interval comparisons. The wide audit includes 72,600 Newton inequalities and 154,000 interval comparisons. Every wide row received an exact Sturm count.

Partial `j=3,4` scans reached `n=525` and `n=600`, but were not audited and are not claimed.

## Certificate standards

- `ctx.workprec` controls Arb precision.
- An imaginary ball containing zero is never treated as proof of non-reality.
- All-real Arb certificates require pairwise-disjoint root balls.
- Root order uses certified endpoint comparisons, including the final equal-degree endpoint.
- Sturm normalization divides by `abs(content)`.
- Wronskian reality is diagnostic only and never decides interlacing.
- Negative discriminants, Newton failures, non-real balls, or failed interval inequalities trigger exact failure records.

## Disk-safe reproducibility

The original full root-enclosure/discriminant logs were removed after the audit when OpenCode snapshot copies filled the system disk. The scripts and audited summary JSON files were retained. `compact_frontier.py` reruns the completed ranges while storing only SHA-256 hashes, discriminant signs/hashes, exact Sturm counts, and pass/fail interval checks. Its outputs are `compact-j5.jsonl`, `compact-wide.jsonl`, and the corresponding summary JSON files. The `--fast` mode skips full discriminant recomputation and thins the exact-Sturm schedule; it is used for the independent `j=3,4` extension.

Finite checks are not a proof of the conjecture.
