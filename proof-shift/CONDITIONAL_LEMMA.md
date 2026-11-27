# Retracted conditional lemma

**Status: REFUTED. Do not use this as a proof.**

An earlier draft claimed that compatibility of `Q = O_{n-j}` with

\[
A_n=nO_{n-1}+2T O_{n-1}+\sum_{r=2}^{j-1}\binom nr O_{n-r}
\]

would imply real-rootedness of `O_n`, because the differential term was
incorrectly assumed to vanish at roots of `Q`.

That assumption is false. At a root `μ` of `Q`,

\[
t[1+(t-1)\partial_t]Q\big|_\mu=\mu(\mu-1)Q'(\mu),
\]

which is generally nonzero.

The weighted-proof agent gives an exact counterexample to the proposed
compatible-window implication in `proof-weighted/NOTES.md`. The correct missing
condition is a weighted root inequality (WRS), not mere compatibility.

This file is retained only as an audit trail.
