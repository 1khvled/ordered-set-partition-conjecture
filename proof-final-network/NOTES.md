# Target-specific network / production-matrix attack: rigorous obstructions

Date: 2026-11-27
Scope: only `proof-final-network` was created or edited (plus two tiny exact probe scripts,
explicitly labelled as exploration, not proof).
Outcome: **no proof obtained**. Three precise, all-index obstructions are proved below. Each
kills one of the three requested routes as stated; §4 records exactly what remains.

---

## 0. Setting

For fixed `j >= 2`,

```
A_j(z) = sum_{r=1}^{j-1} z^r/r!,   C_j(z) = sum_{r>=j} z^r/r!,   D_j(z) = 1/(1-A_j(z)),
a_{n,m} = [t^m] O_{n,j}(t) = n! [z^n] C_j(z)^m D_j(z)^{m+1}.
```

`a_{n,m}` counts ordered set partitions of `[n]` with exactly `m` blocks of size `>= j`;
`a_{n,m} > 0` iff `m = 0` or `n >= j m`. Write `b_{n,m}=a_{n,m}/n!`; row scaling is positive,
so TN is equivalent for `A=(a_{n,m})` and `B=(b_{n,m})`.

---

## 1. The exact path model is correct, but its natural LGV embedding gives the *reflected* matrix

### 1.1 The network

Let `G_j` have vertices `(n,m)` with `n >= j m >= 0`. For `1 <= r <= n` put
`eps(r)=0` if `r<j` and `eps(r)=1` if `r>=j`. There is an edge

```
(n-r, m-eps(r))  -->  (n,m),      weight  binom(n, r),
```

whenever the source vertex is legal (for `r<j` this requires `n >= jm + r`; for `r>=j` it is
automatic). Draw `n` increasing downward and `m` increasing rightward: every edge goes
downward and weakly rightward.

**Lemma 1.** The total weight of paths `(0,0) -> (n,m)` in `G_j` is exactly `a_{n,m}`.

*Proof.* A path is a sequence of block sizes `r_1,...,r_k` summing to `n`, with exactly `m`
entries `>= j`. With `n_i = r_1+...+r_i`, the weight telescopes:

```
prod_i binom(n_i, r_i) = prod_i n_i!/(r_i! n_{i-1}!) = n! / prod_i r_i! .
```

Summing over all such size sequences gives `sum n!/prod r_i!`, which is precisely the number
of ordered set partitions of `[n]` with `m` large blocks, i.e. `a_{n,m}`. ∎

### 1.2 The graph is not planar; planarizing it yields the *reflected* matrix

`G_j` is **not** a planar network in the naive drawing. For `j=2` the two edges
`(2,0) -> (12,1)` and `(4,0) -> (8,1)` have interleaved endpoints (`2<4` but `12>8`)
and share no vertex: in a straight-line drawing they cross without meeting. Hence the
hypothesis of LGV (“vertex-disjoint ⇒ identity matching”) is *false* for the raw digraph.
Indeed for sources `S_1=(2,0), S_2=(4,0)` and sinks `T_1=(12,3), T_2=(12,1)` the two paths
`P_1: S_1 -> (12,1)` (one large block of size 10) and
`P_2: S_2 -> (8,1) -> (12,2) -> (12,3)` (large blocks of sizes 4,4 then a singleton) are
vertex-disjoint but match `S_1 -> T_2, S_2 -> T_1` — a swapped, sign-negative family.

One may planarize by inserting a vertex at every crossing and splitting the two edges there
(assign weights whose products reproduce the original edge weights); path weights are
unchanged and now LGV does apply. Fix `N`, take sources `S_i=(n_i,0)` with
`0 <= n_1 < ... < n_r < N` on the left boundary and sinks `T_j=(N,m^*_j)` with
`N >= j m^*_1 > ... > j m^*_r >= 0` on the bottom boundary, indexed **right-to-left**.
Around the outer face they appear as `S_1,...,S_r,T_r,...,T_1`, so nonintersecting families
force the identity matching. The telescoping of §1.1 gives
`wt(S_i -> T_j) = binom(N,n_i) a_{N-n_i,m^*_j}`, hence

```
det[ binom(N,n_i) a_{N-n_i, m^*_j} ] >= 0        (m^*_1 > ... > m^*_r).      (L)
```

Reindexing columns by *increasing* `m` reverses their order, so

```
det[ binom(N,n_i) a_{N-n_i, m_j} ]_{m_1<...<m_r}
      = (-1)^{binom(r,2)} * (nonnegative).
```

For `r=2` this determinant is `\le 0`, not `\ge 0`.

**Concrete check.** `j=2, N=7, n={2,5}, m={0,1}`: since `a_{5,0}=120, a_{5,1}=311,
a_{2,0}=2, a_{2,1}=1`, the matrix is `21 * [[120,311],[2,1]]`, determinant
`21^2 (120-622) < 0` — exactly the sign (L) predicts.

Row reflection does not repair this: for `j=2`, rows `{2,4}`, cols `{0,1}` of `A` give

```
det[[2,1],[24,45]] = 66 > 0,   but   det[[24,45],[2,1]] = -66 < 0.
```

So the LGV output (L) is a *different* matrix inequality from the target (all minors of
`A >= 0`), and total nonnegativity is not invariant under the reflection `n -> N-n`.

**Consequence.** The weighted path model of §1.1 is a correct exact model, but its natural
planar embedding cannot prove the target. Any successful LGV proof needs a genuinely new
network whose *sources* are indexed by the row indices `n_i` and whose *sinks* are indexed by
the column indices `m_j`; `G_j` has a single source `(0,0)`, so it cannot be that network.

---

## 2. The natural left production matrix fails for **every** `j >= 3` (general exact witness)

With `K_j(z)=C_j(z)/z`, the natural degenerate-Riordan factorization
`b_{n,m} = [z^n] D_j^{m+1} (z K_j)^m` gives the Chen–Fu–Ruan left production matrix

```
Q_{n,0} = d_n := [z^n] D_j,        Q_{n,m} = [z^{n-m}] D_j K_j   (m >= 1).
```

**Lemma 2.** For every `j >= 3`,

```
det Q[ {j, j+1}, {0,1} ] = (1/(j+1)!) * ( (j+2) d_j - (j+1) d_{j+1} ) < 0 .
```

*Proof.* `d_0=d_1=1`. `K_j` begins at degree `j-1` with coefficient `1/j!`, so
`Q_{j,1}=[z^{j-1}]D_jK_j = 1/j!` and
`Q_{j+1,1}=[z^{j}]D_jK_j = 1/j! + 1/(j+1)! = (j+2)/(j+1)!`. Hence the minor equals
`d_j (j+2)/(j+1)! - d_{j+1}/j! = ((j+2)d_j - (j+1)d_{j+1})/(j+1)!`.

It remains to show `(j+1)d_{j+1} > (j+2)d_j`. The `D_j = 1 + A_j D_j` recurrence gives

```
d_{j+1} = d_j + sum_{r=2}^{j-1} d_{j+1-r}/r!,
d_j     = d_{j-1} + sum_{r=2}^{j-1} d_{j-r}/r!   (j >= 3).
```

So it suffices to prove `(j+1) sum_{r=2}^{j-1} d_{j+1-r}/r! > d_j`.
First, `d_n` is nondecreasing: `d_0=d_1=1`, and for every `n >= 2`,
`d_n = d_{n-1} + sum_{r=2}^{j-1} d_{n-r}/r! >= d_{n-1}` (with `d_s=0` for `s<0`). Therefore

```
d_j <= d_{j-1} ( 1 + sum_{r=2}^{j-1} 1/r! ) < d_{j-1} (e-1) < 2 d_{j-1}.
```

On the other hand the `r=2` term alone gives

```
(j+1) sum_{r=2}^{j-1} d_{j+1-r}/r! >= (j+1) d_{j-1}/2 >= 2 d_{j-1}   (j >= 3).
```

Combining, `(j+1) sum >= 2 d_{j-1} > d_j`, as required. ∎

Exact values: `j=3: -1/24`; `j=4: -43/1440`; `j=5: < 0`; etc. For `j=2` the
`{j,j+1}` minor is `+1/6`, but the exact `4x4` witness rows `{2,3,4,5}`, cols `{0,1,2,3}`
with determinant `-1/4320` is recorded in `proof-matrix/NOTES.md` §4.3 and
`literature.md` §8.3. Thus Chen–Fu–Ruan Theorem 1.3 cannot be invoked through the natural
factorization for **any** `j >= 2`; this is an all-index statement, not a finite check.

---

## 3. Ordinary TN of the unsheared Riordan array cannot imply the target shear

Section 2 of `proof-matrix/NOTES.md` gives the exact factorization

```
b_{n,m} = R^{(j)}_{n - j m, m},     R^{(j)} = (D_j, C_j D_j / z^j)  (proper Riordan).
```

**Lemma 3 (exact obstruction).** There exists a proper Riordan array `R` whose production
matrix is TN (hence `R` itself is TN by Chen–Liang–Wang) but whose once-sheared copy
`S_{n,m}=R_{n-2m,m}` has

```
det S[{4,6},{0,1}] = det[[14,3],[132,28]] = 392 - 396 = -4 < 0 .
```

This is proved in `proof-matrix/NOTES.md` §3 by an explicit bidiagonal factorization of the
production matrix; it is an exact finite-dimensional theorem, not a numerical sample.
Consequently no theorem of the form “proper TN Riordan array ⇒ TN after the column shear
`n -> n - j m`” exists, and Lemma 3 applies to *all* proper Riordan arrays, in particular to
`R^{(j)}`. Ordinary TN of `R^{(j)}` is therefore logically insufficient for the
sheared-minor inequality (20) of `proof-matrix/NOTES.md`; a target-specific argument using the
special pair `(D_j, C_j/z^j)` is unavoidable.

---

## 4. Precisely what remains open

The three requested routes are blocked as follows:

1. **Planar network / LGV.** §1 proves the natural weighted path network yields (L), a
   reflected sign-twisted inequality, not the target. Required: a new network with sources
   indexed by `n_i`, sinks by `m_j`, and a planar embedding in which every nonintersecting
   family is the identity matching (or a crossing involution with the same property). No
   such network is known.
2. **TN production matrix.** §2 proves the *natural* left production matrix has negative
   minors for every `j`. Required: a different factorization of `A` into nonnegative
   bidiagonal/TN factors that respects the zero cone `n < jm`. The uniqueness of the left
   production decomposition means this cannot be a reparametrization of `Q`; it must be a
   genuinely different construction.
3. **Sheared-minor lemma.** §3 proves no general shear preservation holds. Required: use the
   special identity `C_j(z) = e^z - (1 + A_j(z))` (equivalently, the Riccati relation
   `F_j' = (2 + (t-1) z^{j-1}/(j-1)!) F_j^2 - F_j`) to exploit the target's structure. The
   signed `(t-1)` term is exactly what prevents a direct positive-convolution induction.

A weaker but still unproved target is row Pólya frequency (real-rootedness of each
`O_{n,j}`); TN of `A` is not implied by it, nor does TN of `A` alone imply it.

---

## 5. Exploration only — not part of the proof

The two scripts `probe.js` and `probe2.js` generate `a_{n,m}` by the exact first-block
recurrence and test determinants with fraction-free (Bareiss) `BigInt` arithmetic. For
`j=2,3,4` and `n<=18` they found **no** negative minor of sizes 2 and 3 exhaustively and
none in 30 000 random size-4 and 10 000 random size-5 samples. This is finite evidence only
and is explicitly *not* used above. No claim of all-index TN is made here.
