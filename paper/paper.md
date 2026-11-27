# Ordered Set Partitions Counted by Big Blocks
## Exact Recurrences and Certified Computational Evidence

**Khaled Abdelli**  
27 November 2026

### Abstract

We study Per Alexandersson's ordered big-block conjecture. We derive an integer convolution recurrence, a universal differential recurrence, a Riccati differential equation, a weighted root-separation reduction, and an exact obstruction to a naive compatible-window induction. We describe rigorous exact methods for deciding real-rootedness and consecutive interlacing, and record the independently checked parameter ranges. The distinction between ordering all blocks and ordering only the large blocks is made explicit. Finite certificates are not presented as a proof of the all-parameter conjecture.

## 1. Definitions and conjecture

For `j ≥ 2`, let `π = (B₁ | ⋯ | Bₖ)` be an ordered set partition of `[n]`, and define

> Oₙ,ⱼ(t) = Σπ t^(#{i : |Bᵢ| ≥ j}).

The empty ordered partition contributes O₀,ⱼ(t) = 1. Put

> Aⱼ(z) = Σᵣ₌₁^{j−1} zʳ/r!,  
> Cⱼ(z) = eᶻ − 1 − Aⱼ(z),  
> Bⱼ(z,t) = Aⱼ(z) + t Cⱼ(z).

The labelled-sequence construction gives

> Fⱼ(z,t) = Σₙ≥₀ Oₙ,ⱼ(t) zⁿ/n! = 1/(1 − Bⱼ(z,t)).

**Conjecture (Alexandersson, 2026).** For every `j ≥ 2` and `n ≥ 0`, Oₙ,ⱼ(t) is real-rooted, and Oₙ,ⱼ interlaces Oₙ₊₁,ⱼ.

## 2. Exact coefficient identities

**Proposition.** The coefficients of Oₙ,ⱼ are nonnegative integers and `deg Oₙ,ⱼ = ⌊n/j⌋`. Moreover,

> O₀,ⱼ = 1,  
> Oₙ,ⱼ = Σᵣ₌₁^{j−1} C(n,r) Oₙ₋ᵣ,ⱼ + t Σᵣ₌ⱼⁿ C(n,r) Oₙ₋ᵣ,ⱼ.

**Proof.** Choose the first block B₁. If its size is r, its labels can be selected in C(n,r) ways, and the remaining ordered partition contributes Oₙ₋ᵣ,ⱼ. The first sum handles unweighted blocks of size less than j; the second handles a block contributing one factor of t. Nonnegativity is immediate. A term tᵐ requires at least jm labels, and m = ⌊n/j⌋ is attained by using only large blocks. ∎

For j = 3, the recurrence gives

```text
O_0 = 1
O_1 = 1
O_2 = 3
O_3 = 12 + t
O_4 = 66 + 9t
O_5 = 450 + 91t
O_6 = 3690 + 973t + 20t^2
```

At t = 1, all block weights are equal, so

> Fⱼ(z,1) = 1/(2 − eᶻ),  
> Oₙ,ⱼ(1) = Σₖ₌₀ⁿ k! {n \brace k},

independently of j.

## 3. Riccati and convolution identities

Write

> Hₙ,ⱼ(t) = n![zⁿ]Fⱼ(z,t)² = Σᵣ₌₀ⁿ C(n,r) Oᵣ,ⱼ(t) Oₙ₋ᵣ,ⱼ(t).

Since Fⱼ' = Bⱼ'Fⱼ², coefficient comparison in

> Bⱼ' = 2 + (t − 1) z^{j−1}/(j−1)! − (1 − Bⱼ)

gives the following exact recurrence.

**Theorem.** For n ≥ 0,

> Oₙ₊₁,ⱼ = 2Hₙ,ⱼ − Oₙ,ⱼ + (t − 1) C(n,j−1) Hₙ₋ⱼ₊₁,ⱼ,

where Hₘ,ⱼ = 0 for m < 0.

The final term is signed. An equivalent form exposes nonnegative structural coefficients:

> Fⱼ' = (Bⱼ + Bⱼ' − 1)Fⱼ² + Fⱼ.

Indeed, if bᵣ = 1 for r < j and bᵣ = t for r ≥ j, then the coefficient of zᑫ/q! in Bⱼ + Bⱼ' − 1 is 0 for q = 0, 2 for 1 ≤ q < j−1, 1 + t for q = j−1, and 2t for q ≥ j.

## 4. A universal short recurrence

Put `T = t ∂/∂t`. For every `j ≥ 2`,

> Oₙ,ⱼ = n Oₙ₋₁,ⱼ + 2T Oₙ₋₁,ⱼ + Σᵣ₌₂^{j−1} C(n,r) Oₙ₋ᵣ,ⱼ  
> + C(n−1,j−1) t [1 + (t−1) ∂/∂t] Oₙ₋ⱼ,ⱼ.

Negative-index terms are omitted. For `j = 3`,

> Oₙ,₃ = n Oₙ₋₁,₃ + 2t O'ₙ₋₁,₃ + C(n,2) Oₙ₋₂,₃  
> + t(n−1)(n−2)/2 · [1 + (t−1) ∂/∂t] Oₙ₋₃,₃.

**Proof.** Write `wⱼ = z^{j−1}/(j−1)!` and `Cⱼ = eᶻ − Eⱼ₋₁`. Since `Cⱼ' = Cⱼ + wⱼ` and `∂ₜFⱼ = CⱼFⱼ²`, differentiation of

> t Cⱼ Fⱼ − 2 ∫₀ᶻ T Fⱼ(u,t) du  
> − t/(j−1)! ∫₀ᶻ u^{j−1}[Fⱼ + (t−1)∂ₜFⱼ](u,t) du

shows that it is constant in `z`; its value at zero is zero. Coefficient extraction gives

> t Σᵣ₌ⱼⁿ C(n,r) Oₙ₋ᵣ,ⱼ = 2T Oₙ₋₁,ⱼ  
> + C(n−1,j−1)t[1 + (t−1)∂ₜ]Oₙ₋ⱼ,ⱼ.

Substitution into the first-block recurrence proves the identity. ∎

The same coefficients appear in Chen–Zhu, *A Snail Race Problem* (2026), Theorems 3 and 7.

## 6. A weighted root-separation reduction

Fix `n ≥ j`, and suppose the window `Oₙ₋ⱼ,ⱼ, ..., Oₙ₋₁,ⱼ` is commonly compatible. Put `q = Oₙ₋ⱼ,ⱼ`, `p = Oₙ₋₁,ⱼ`, `Sₙ = Σᵣ₌₂^{j−1} C(n,r)Oₙ₋ᵣ,ⱼ`, and `c = C(n−1,j−1)`. For every root `μ` of `q`, impose

> −[n p(μ) + 2μ p'(μ) + Sₙ(μ)] / q'(μ) > c μ(μ−1). (WRS)

Then `Oₙ,ⱼ` has `deg q + 1` distinct negative real roots and interlaces `q`.

**Proof.** At a root `μ` of `q`,

> t[1 + (t−1)∂ₜ]q evaluated at t=μ equals μ(μ−1)q'(μ).

Thus (WRS) is exactly the assertion `sgn Oₙ,ⱼ(μ) = −sgn q'(μ)`. These signs alternate at the roots of `q`. Since `deg Oₙ,ⱼ = deg q + 1`, the signs at `−∞`, between consecutive roots of `q`, and at `0` force `deg q + 1` distinct roots by the intermediate value theorem. ∎

Mere compatibility is insufficient. For `j = 2, n = 5`, take `p = (t+2)(t+4)` and `q = t+3`. Then `q` interlaces `p`, and also interlaces `R₅(p) = 9t²+42t+40`. Nevertheless, at `μ = −3`, the left side of (WRS) is `5` and the right side is `48`. The resulting test polynomial `17t²+50t+40` has discriminant `−220`. A successful induction must therefore prove the quantitative inequality (WRS), not only a compatible window.

## 7. Related ordered lifts are different

Let Pₙ,ⱼ count ordinary set partitions by their large blocks. Its exponential generating function is e^{Bⱼ}. Ordering only the m large blocks replaces tᵐ by m!tᵐ, and therefore has exponential generating function

> e^{Aⱼ(z)}/(1 − t Cⱼ(z)).

This is **not** 1/(1 − Aⱼ − tCⱼ), which orders all blocks. The distinction is essential when applying theorems for restricted big-block families to the present conjecture.

## 8. Rigorous finite certification

All coefficient generation uses integer arithmetic. Root isolation uses FLINT's integer-polynomial routine with Arb complex balls. A root is certified non-real when its imaginary ball excludes zero; a real-rooted verdict requires every returned root enclosure to lie on the real axis. Consecutive interlacing is checked by comparing rigorous lower and upper endpoints of ordered real-root balls, including every degree-increasing transition.

Independent checks use exact Sturm sequences, exact discriminants, and Newton's ultra-log-concavity inequalities. A negative discriminant is an immediate obstruction; a zero discriminant is not, since repeated real roots are allowed. Floating-point root finders are never used as certificates.

### Audited ranges

All coefficient arithmetic is integral. Real-rootedness and interlacing are certified with rigorous Arb balls; exact Sturm counts, discriminants, and Newton inequalities provide independent checks. No float64 root finder is used. The new scan omits `j = 2` as requested.

| Range | Rows | Consecutive pairs | Exact Sturm rows |
|---|---:|---:|---:|
| `j = 3..8`, `0 ≤ n ≤ 45` | 276 | 270 | 276 |
| `j = 5`, `141 ≤ n ≤ 800` | 660 | 660 | 43 |
| `j = 6..60`, `7j ≤ n ≤ 12j−1` | 9,075 | 9,075 | 9,075 |

No failure was found. The `j = 5` audit checked 61,182 Newton inequalities and 122,892 certified interval comparisons. The wide audit checked 72,600 Newton inequalities and 154,000 interval comparisons. Every wide-range row received an exact Sturm count. A total of 5,012 Sturm controls were cross-checked against SymPy and an independent rational implementation.

Partial `j = 3,4` scans reached `n = 525` and `n = 600`, respectively, but were not audited and are not claimed. Full per-row logs are regenerable from the supplied scripts; compact hash/count certificates are retained.

Wronskian reality is not used as an interlacing criterion. Sampled Wronskians can be nonreal even when the corresponding row pair passes the certified interlacing inequalities.

## 9. Conclusion and remaining work

The all-parameter conjecture is not inferred from finite computation. The proved identities above provide several possible starting points for an interlacing induction, while the certified computations identify the precise unchecked frontier.

## References

1. P. Alexandersson. *Open problems on real-rootedness*. SymCat, 2026. https://www.symmetricfunctions.com/realRootedOpenProblems.htm#conjOrderedBigBlocks
2. Y. Chen and W. Zhu. *A Snail Race Problem*. arXiv:2607.02085, 2026. https://arxiv.org/abs/2607.02085
3. Public computational audit, Issue #108, 2026. https://github.com/1337isnot1337/solved-conjectures/issues/108
