# Conditional real-rootedness lemma

Fix `j ≥ 2` and `n ≥ j`. Write

\[
A_n=nO_{n-1}+2T O_{n-1}+\sum_{r=2}^{j-1}\binom nr O_{n-r},
\qquad T=t\frac{d}{dt},
\]

and

\[
O_n=A_n+\binom{n-1}{j-1}t\bigl[1+(t-1)\partial_t\bigr]O_{n-j}.
\]

All polynomials are understood at the fixed threshold `j`. Assume `O_{n-j}` and `A_n` have simple negative real roots and are compatible in the following precise sense. Put `e = deg O_{n-j}`.

* If `deg A_n = e+1`, with roots `a_1 < ... < a_{e+1}` and roots `q_1 < ... < q_e` of `O_{n-j}`, assume
  `a_1 < q_1 < a_2 < ... < q_e < a_{e+1}`.
* If `deg A_n = e`, with roots `a_1 < ... < a_e`, assume
  `q_1 < a_1 < q_2 < ... < q_e < a_e`.

Then `O_n` has only real roots.

## Proof

The differential term vanishes at every `q_i`, so `O_n(q_i)=A_n(q_i)`. These signs alternate. In the degree-`e+1` case, the sign at `-infinity` is opposite to the sign at `q_1`, and the sign at `q_e` is negative while `O_n(0)=A_n(0)>0`. Thus the intermediate value theorem gives one root left of `q_1`, one between every consecutive pair `q_i,q_{i+1}`, and one right of `q_e`: `e+1` roots.

In the degree-`e` case, `O_n` has degree `e+1` and positive leading coefficient. Its sign at `-infinity` is therefore opposite to `A_n(q_1)`. Again there is one root left of `q_1`, `e-1` between the `q_i`, and one between `q_e` and `0`, because `A_n(q_e)<0<O_n(0)`. Again there are `e+1` roots.

The degree of `O_n` is exactly `e+1`, so these exhaust all roots. ∎

## Computational status of the hypothesis

`lag_probe.py` checks the hypothesis by rigorous Arb root enclosures for `j=2,...,15` and `n<=100`. It holds in every case. This is finite evidence for the compatibility hypothesis, not a proof of it.

If a sliding compatible-window induction proves this hypothesis for all `n,j`, the lemma completes the real-rootedness conjecture. Consecutive interlacing would still require a separate root-order argument.
