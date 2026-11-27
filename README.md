# Ordered set partitions counted by big blocks

A reproducible exact-arithmetic attack on Per Alexandersson's 2026 conjecture

\[
O_{n,j}(t)=\sum_{\pi}t^{\#\{B\in\pi:|B|\ge j\}},\qquad
\sum_{n\ge0}O_{n,j}(t)\frac{z^n}{n!}=\frac1{1-B_j(z,t)},
\]

where the sum is over **ordered** set partitions and

\[
B_j(z,t)=\sum_{r=1}^{j-1}\frac{z^r}{r!}
+t\left(e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!}\right).
\]

The conjecture asserts that every \(O_{n,j}\) is real-rooted and that
\(O_{n,j}\interl O_{n+1,j}\).

> **Status:** research in progress. Finite exact checks are evidence, not a
> proof of the all-parameter conjecture. This repository will state precisely
> what is proved, what is certified computationally, and what remains open.

## Correct coefficient recurrence

Choosing the first block gives the integer recurrence

\[
O_{0,j}=1,
\qquad
O_{n,j}=
\sum_{r=1}^{j-1}\binom nr O_{n-r,j}
+t\sum_{r=j}^{n}\binom nr O_{n-r,j}.
\]

The binomial coefficient is \(\binom nr\), **not** \(\binom{n-1}{r-1}\):
the latter forces the block containing \(1\) to be first and undercounts
ordered partitions.

For example, for \(j=3\),

```text
O_0 = 1
O_1 = 1
O_2 = 3
O_3 = 12 + t
O_4 = 66 + 9t
O_5 = 450 + 91t
O_6 = 3690 + 973t + 20t^2
```

## Repository layout

- `scan.py`, `interlace.py`, `frontier.py` — exact coefficient generation,
  certified root isolation, Sturm cross-checks, discriminants and Newton tests.
- `guess_recurrence.py` — exact search for short differential recurrences.
- `numerical/` — independent corrected numerical audit and machine-readable logs.
- `proof/` — proof-track identities, candidate recurrences and proof gaps.
- `literature.md` — sourced context and distinction from related ordered lifts.
- `lean/` — Lean 4 definitions and finite Sturm-certificate infrastructure.
- `paper/` — manuscript source and build instructions.

## Python environment

Tested with CPython 3.14.4, `python-flint` 0.9.0 and `gmpy2` 2.3.2.

```powershell
py -3.14 -m pip install python-flint gmpy2 sympy
py -3.14 scan.py --jmin 3 --jmax 20 --nmax 150 --out scan-results.jsonl
py -3.14 interlace.py --jmin 3 --jmax 20 --nmax 150 --wronskian
py -3.14 frontier.py --small-nmin 141 --small-nmax 600 --jmin 6 --jmax 60
```

All final root decisions use rigorous Arb enclosures or exact Sturm arithmetic;
float64 root finders are not used as certificates.

## References

- Per Alexandersson, *Open problems on real-rootedness*, SymCat (2026):
  <https://www.symmetricfunctions.com/realRootedOpenProblems.htm#conjOrderedBigBlocks>
- Public computational audit, Issue #108:
  <https://github.com/1337isnot1337/solved-conjectures/issues/108>
