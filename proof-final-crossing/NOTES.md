# Sign-reversing crossing involution for the weighted block-size path network: exact obstruction

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-crossing` was created or edited.

## 0. Setting

Fix `j >= 2`.  The network `G_j` has vertices `(n,m)` with `n >= j m >= 0`.  For `1 <= r <= n` put
`eps(r)=0` if `r<j` and `eps(r)=1` if `r>=j`.  There is an edge

```
(n-r, m-eps(r))  -->  (n,m),      weight  binom(n,r),
```

whenever the source vertex is legal.  The total weight of paths `(0,0)->(n,m)` is `a_{n,m}`.
See `proof-final-network/NOTES.md`, Section 1.1.

The target is the nonnegativity of all minors of `A=(a_{n,m})` (equivalently of
`B=(b_{n,m})`, `b_{n,m}=a_{n,m}/n!`).  The natural digraph `G_j` is nonplanar; for
sources `S_i=(n_i,0)` and sinks `T_j=(N,m_j^*)` the LGV determinant is the signed sum
over all families, including vertex-disjoint swapped families.  The question is whether
an explicit sign-reversing involution on the *crossing/swapped* families can cancel them,
leaving only the identity (target) families.

**Outcome.** No such involution exists.  The obstruction is structural and all-index:
crossings between edges of `G_j` never occur at original vertices, and an uncrossing would
have to split an atomic weighted block.  Splitting changes the telescoping weight and the
large-block count `m`, so the resulting objects are no longer paths in `G_j` and do not
correspond to terms of `a_{n,m}`.  The only weight-preserving way to resolve crossings is
planarization, which gives the reflected matrix `(L)` of `proof-final-network/NOTES.md`,
not the target.

## 1. What a valid involution would have to preserve

A path in `G_j` is a sequence of block sizes `r_1,...,r_k` summing to `n`, with exactly
`m` of the `r_i` at least `j`.  Its weight telescopes:

```
prod_i binom(n_i,r_i) = n! / prod_i r_i!,
```

where `n_i=r_1+...+r_i`.  Thus each path corresponds exactly to an ordered partition of
`[n]` into blocks of those sizes, and `m` is the number of large blocks.

For an LGV-style involution to cancel a family and leave the target determinant, any
replacement path must:

1. be a path in the same graph `G_j`;
2. have the same ordered block-size sequence up to the local change, so that the
   telescoping product is unchanged;
3. preserve the number `m` of large blocks, since the column index (sink) is fixed.

In particular, the involution may only reroute at vertices of `G_j`: at a vertex
`(N,M)` the two paths share a state, and swapping their tails preserves both the weight
product and the large-block counts.

## 2. Edge structure lemma (all `j >= 2`)

Let an edge of `G_j` have tail `(n-r, m-eps(r))` and head `(n,m)`.

- **Small block (`1 <= r < j`).**  Then `eps(r)=0`.  The edge is horizontal in the
  `(n,m)`-plane: `m` is constant along it, and `n` increases from `n-r` to `n`.
- **Large block (`r >= j`).**  Then `eps(r)=1`.  The edge raises the large-block count
  by exactly `1` over its length `r`: along the straight segment from
  `(n-r, m-1)` to `(n,m)`, the `m`-coordinate is

  ```
  m(s) = (m-1) + s/r,   0 < s < r.
  ```

**Lemma.** No point in the relative interior of a large edge is a vertex of `G_j`.

*Proof.*  For `0 < s < r`, the value `m(s)` lies strictly between the consecutive
integers `m-1` and `m`.  Every vertex of `G_j` has integer `m`-coordinate.  Hence the
point cannot be a vertex.  ∎

For a small edge, interior points have integer `m` (the same constant), and their
`n`-coordinate may indeed be that of some vertex.  But in `G_j` the edge is *atomic*: it
represents one small block of size `r`.  It does not pass through those intermediate
vertices as states of the path; the path does not "stop" there.  Treating such a point as
a rerouting vertex would split one small block into two or more small blocks, changing
the block-size sequence and hence the telescoping weight.

## 3. Crossings never occur at original vertices

A crossing (in the planar drawing) between two edges can happen in two ways.

**Large-large crossing.**  Both edges are large.  By the Lemma, every interior point of each
has non-integer `m`.  A common point of the two segments would have to have the same
`m`-coordinate on both, hence a non-integer.  It therefore cannot be a vertex of `G_j`.

**Large-small crossing.**  A large edge has `m`-coordinate in `(m-1,m)` for its base level
`m`, while a small edge has constant integer `m_0`.  For them to meet we would need
`m_0 in (m-1,m)`, impossible for integers.  Hence a large edge and a small edge cannot
cross in their interiors at all.

**Small-small crossing.**  Both have constant integer `m`.  If the constants differ, the
edges lie on distinct horizontal lines and do not meet.  If the constants are equal, the
edges overlap in horizontal segments rather than crossing transversally; their common
points are not vertices unless an endpoint coincides, and splitting such an overlap would
again replace one small block by several.

Thus every genuine crossing of `G_j` occurs strictly inside (at least) one large edge, at
a non-vertex point with fractional `m`-coordinate.

## 4. Why this blocks the involution

The standard LGV tail-swap requires a shared *vertex* `v`: two paths
`S_i -> ... -> v -> ... -> T_i` and `S_j -> ... -> v -> ... -> T_j` are replaced by
`S_i -> ... -> v -> ... -> T_j` and `S_j -> ... -> v -> ... -> T_i`.  The weight product
is unchanged because multiplication commutes, and the block counts are unchanged because
the prefixes and suffixes are unchanged as paths.

In `G_j`, crossing families have no shared vertex.  The crossing point is interior to a
large edge, with fractional `m`.  To perform an analogous swap one would have to split that
large edge at the crossing point, producing two "half-edges" whose `m`-coordinates are
fractional.  The resulting objects would not be paths in `G_j`: their intermediate states
are not vertices, and their block-size sequences would no longer be integers.  Consequently
their weights would not telescope to `n!/prod r_i!`, and the large-block count `m` would
not be well-defined.  Hence no weight-preserving, `m`-preserving involution on the actual
path families of `G_j` can cancel the crossings.

The only way to resolve the crossings while preserving weights is to introduce new vertices
at the crossing points and split the edges there, assigning half-edge weights whose
products reproduce the original weights.  That is exactly planarization.  But then LGV on
the planarized graph counts paths that may turn at these artificial vertices; the surviving
nonintersecting families are indexed by the outer-face order `S_1,...,S_r,T_r,...,T_1`,
which yields the reflected inequality `(L)` of `proof-final-network/NOTES.md`:

```
det[ binom(N,n_i) a_{N-n_i, m_j^*} ]_{m_1^* > ... > m_r^*} >= 0,
```

not the target minor of `A` with rows and columns in increasing order.  (Equivalently,
planarization changes the problem to a different, reflected matrix.)

## 5. Exact witness

Take `j=2`.  The two edges

```
E1: (2,0) -> (12,1)   (large block of size 10),
E2: (4,0) -> (8,1)    (large block of size 4)
```

cross.  Solving

```
(n-2)/10 = (n-4)/4
```

gives `n = 16/3`, `m = 1/3`, a non-vertex.  These edges are the ones appearing in the
swapped family recorded in `proof-final-network/NOTES.md` Section 1.2:

```
P1: (2,0) -> (12,1),       P2: (4,0) -> (8,1) -> (12,2) -> (12,3).
```

The crossing point `(16/3, 1/3)` is interior to both `E1` and `E2`; it is not a vertex of
`G_2`.  Any uncrossing must split a large block of size `10` or `4` at a non-integer
position, which is impossible inside the block-size model.

## 6. Conclusion

The obstruction is precise and all-index:

- Crossings in `G_j` occur only in the interiors of large edges (or as overlapping
  small-edge segments), never at original vertices.
- Therefore an uncrossing cannot be performed by tail-swapping at a vertex; it would have
  to split an atomic weighted block, destroying the telescoping weight and the
  large-block count.
- The only weight-preserving resolution is planarization, which does not leave the target
  determinant but the reflected matrix `(L)`.

Hence **no sign-reversing crossing involution on the weighted block-size path network
`G_j` can cancel the swapped families and leave the target nonintersecting determinant.**
A proof of the target minors must use a different network, a different factorization, or
a target-specific argument exploiting `C_j(z)=e^z - (1+A_j(z))`, as recorded in
`proof-final-network/NOTES.md` Section 4.