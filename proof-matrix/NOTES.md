# Ordered big blocks: matrix, stretched-Riordan, and shear notes

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-matrix` was created or edited.

## Status

The original family is

\[
 A_j(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad
 C_j(z)=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!},\qquad
 D_j(z)=\frac1{1-A_j(z)}.
\]

Thus \(B_j(0,t)=0\), \(F_j=1/(1-B_j)\), and

\[
 a_{n,m}=[t^m]O_{n,j}(t)
 =n![z^n]C_j(z)^mD_j(z)^{m+1}.
\]

The support is

\[
 a_{n,m}>0\quad\Longleftrightarrow\quad m=0\ \text{or}\ n\ge jm.
\tag{1}
\]

Consequently literal strict TP is impossible; the natural targets are TN or support-strict TN.

No all-index TN, support-strict TN, or row-PF proof was obtained. The strongest new result below is an **exact obstruction to a general stretched/sheared-Riordan theorem**: a proper Riordan array with a TN production matrix can have a sheared copy with a negative minor. The target therefore needs a genuinely special shear lemma or a different network/factorization.

---

## 1. Literature normalization and Chen–Zhu

Yu Chen and Wenxiang Zhu, *A Snail Race Problem*, arXiv:2607.02085, identify

\[
 r_m(n,k)=a_{n,k}.
\]

Their Theorem 3 is the correct first-block recurrence

\[
 r_m(n,k)=\sum_{i=1}^{m-1}\binom ni r_m(n-i,k)
 +\sum_{i=m}^{n}\binom ni r_m(n-i,k-1),
\tag{2}
\]

and Theorem 7 gives

\[
 \sum_{n\ge0}r_m(n,k)\frac{x^n}{n!}
 =\frac{(e^x-T_{m-1}(x))^k}{(2-T_{m-1}(x))^{k+1}}.
\tag{3}
\]

Here \(2-T_{m-1}=1-A_j\) and \(e^x-T_{m-1}=C_j\), so (3) is exactly the target column EGF. Theorem 10 gives a positive sum over the labels in the large blocks and ordered small-block strings. Section 6, Theorem 15, gives a second formula for \(j=2\) using partial Bell polynomials. None of these results is a minor-positivity theorem.

The normalization-correction discussion in `proof/PROOF_NOTES.md` is not used here. The universal recurrence (10) in that file is valid for the original family.

---

## 2. Exact proper-Riordan/shear factorization

Put

\[
 h_j(z)=\frac{C_j(z)}{z^j}D_j(z),\qquad
 R^{(j)}_{q,m}=[z^q]D_j(z)h_j(z)^m.
\tag{4}
\]

Since \(h_j(0)=1/j!\ne0\), \(R^{(j)}=(D_j,h_j)\) is a proper ordinary Riordan array. The ordinary coefficient matrix of the target is exactly the column shear

\[
 b_{n,m}:= \frac{a_{n,m}}{n!}
 =[z^n]D_j(z)(z^jh_j(z))^m
 =R^{(j)}_{n-jm,m},
\tag{5}
\]

with negative row indices interpreted as zero. Positive row scaling gives

\[
 \det(a_{I,J})=\Bigl(\prod_{i\in I}n_i!\Bigr)\det(b_{I,J}),
\tag{6}
\]

so TN is equivalent for \(A=(a_{n,m})\) and \(B=(b_{n,m})\).

Thus two distinct statements must not be conflated:

1. ordinary TN of the proper Riordan array \(R^{(j)}\);
2. TN of its sheared copy \(S^{(j)}_{n,m}=R^{(j)}_{n-jm,m}\).

Statement 2 is exactly the unresolved all-minor statement for the target. Statement 1 does not imply it in general.

---

## 3. Exact obstruction: TN proper Riordan arrays need not have TN shears

This subsection gives a complete finite-dimensional proof that no theorem of the form “a proper TN Riordan array remains TN after an arbitrary column shear” can hold.

### 3.1 A proper TN Riordan array

Let \(R=(r_{n,k})_{n,k\ge0}\) be defined by the Riordan \(A/Z\) recurrence

\[
 r_{0,0}=1,
\]

\[
 r_{n+1,0}=\sum_{\ell\ge0}z_\ell r_{n,\ell},\qquad
 r_{n+1,k+1}=\sum_{\ell\ge0}a_\ell r_{n,k+\ell},
\tag{7}
\]

where

\[
 (z_0,z_1,z_2,\ldots)=(1,1,0,\ldots),\qquad
 (a_0,a_1,a_2,\ldots)=(1,2,1,0,\ldots).
\tag{8}
\]

Equivalently, \(R=(g,f)\) is the proper Riordan array with \(A\)-sequence \(a\) and \(Z\)-sequence \(z\). Its production matrix is

\[
 J(R)=
 \begin{pmatrix}
 1&0&0&0&0&\cdots\\
 1&1&0&0&0&\cdots\\
 0&2&1&0&0&\cdots\\
 0&1&2&1&0&\cdots\\
 0&0&1&2&1&\cdots\\
 0&0&0&1&2&\ddots\\
 \vdots&\vdots&\vdots&\vdots&\ddots&\ddots
 \end{pmatrix}.
\tag{9}
\]

**Claim.** \(J(R)\) is TN.

**Proof.** Factor \(J(R)=UV\), where \(U,V\) are unit lower bidiagonal matrices. Let their subdiagonals be

\[
 U_{n,n-1}=u_n,\qquad V_{n,n-1}=v_n,
\]

with

\[
 u_1=u_2=0,\qquad u_n=\frac{n-1}{n}\quad(n\ge3),
\tag{10}
\]

and

\[
 v_1=1,\qquad v_2=2,\qquad v_n=\frac{n+1}{n}\quad(n\ge3).
\tag{11}
\]

All entries are nonnegative. The product has first subdiagonal \(u_n+v_n\), namely

\[
 1,2,2,2,\ldots,
\]

and second subdiagonal \(u_nv_{n-1}\), namely

\[
 0,1,1,1,\ldots.
\]

These are exactly the subdiagonals of (9), while both factors have unit diagonal. A nonnegative lower bidiagonal matrix is TN, and the product of TN matrices is TN. Hence \(J(R)\) is TN. \(\square\)

By Chen–Liang–Wang, *Total positivity of Riordan arrays*, arXiv:1601.05637, Theorem 2.1, a proper Riordan array whose production matrix is TN is itself TN. Therefore \(R\) is TN.

The first rows, obtained directly from (7), are

\[
 \begin{array}{c|rrrrrrr}
 n\backslash k&0&1&2&3&4&5&6\\ \hline
 0&1&0&0&0&0&0&0\\
 1&1&1&0&0&0&0&0\\
 2&2&3&1&0&0&0&0\\
 3&5&9&5&1&0&0&0\\
 4&14&28&20&7&1&0&0\\
 5&42&90&75&35&9&1&0\\
 6&132&297&275&154&54&11&1
 \end{array}.
\tag{12}
\]

### 3.2 The shear fails

Define the once-sheared array

\[
 S_{n,m}=R_{n-2m,m}.
\tag{13}
\]

Using (12), the minor with rows \(\{4,6\}\) and columns \(\{0,1\}\) is

\[
 \det S[\{4,6\},\{0,1\}]
 =\det\begin{pmatrix}
 14&3\\
 132&28
 \end{pmatrix}
 =14\cdot28-3\cdot132
 =392-396
 =-4.
\tag{14}
\]

Thus \(R\) is a proper TN Riordan array, but its shear (13) is not TN. This is an exact all-index obstruction to general shear preservation, not a finite numerical test of the ordered-big-block conjecture.

**Consequence for the target.** The factorization (5) reduces the problem to proving a shear theorem for the particular pair

\[
 D_j(z)=\frac1{1-A_j(z)},\qquad
 h_j(z)=\frac{C_j(z)}{z^j}D_j(z),
\]

not to invoking a general stretched-Riordan TP theorem. No such general theorem exists.

---

## 4. What the stretched/compression literature actually proves

### 4.1 Hennessy and He–Sprugnoli

Aoife Hennessy’s thesis, *A Study of Riordan Arrays with Applications to Continued Fractions, Orthogonal Polynomials and Lattice Paths* (2011), defines vertically stretched arrays \((g,xf(x))\), notes that they are noninvertible, and discusses left inversion. It does not state a total-positivity preservation theorem.

Tian-Xiao He and Renzo Sprugnoli, *The Vertical Recursive Relation of Riordan Arrays and Their Matrix Representation*, arXiv:2206.12898, equation (12), gives the vertical recurrence

\[
 d_{n,k}=\sum_{\ell=1}^{n-k+1}f_\ell d_{n-\ell,k-1}.
\tag{15}
\]

This is an algebraic/sequence-characterization result. It is not a TN theorem for vertically stretched arrays.

Paul Barry, *A New Group in the Riordan Family of Matrix Groups: the Sprugnoli Group*, arXiv:2605.16633, gives a production-matrix characterization of the Sprugnoli group and uses stretched arrays in its formulation. It does not assert that stretching an arbitrary TN Riordan array preserves TN; the counterexample (14) shows that such an assertion would be false.

### 4.2 He’s multiple-Riordan compression theorem is not the required shear

Tian-Xiao He, *The Multiple Riordan Group and the Multiple Riordan Type Arrays*, arXiv:2504.04049, Section 6, defines the compression of an \(\ell\)-multiple Riordan array by

\[
 \widehat d_{n,k}=d_{n\ell-(\ell-1)k,k}.
\tag{16}
\]

Theorem 6.4 states that if the compressed series \(\widehat g,\widehat f_1,\ldots,\widehat f_\ell\) are PF formal power series, then the compressed multiple Riordan array is TP.

This is the closest located theorem, but it does **not** apply to (5):

* the row map in (16) is \(n\ell-(\ell-1)m=n+(\ell-1)(n-m)\), whereas the target uses \(n-jm\);
* no fixed \(\ell\) makes these maps equal for all \(n,m\);
* a compressed multiple Riordan array has nonzero diagonal, while the target has \(b_{m,m}=0\) for every \(j\ge2\), \(m>0\).

Thus Theorem 6.4 is a theorem about a different compression, not about arbitrary vertical shears of a proper Riordan array.

### 4.3 Chen–Fu–Ruan

Chen, Fu, and Ruan, *A Unified Planar Network Approach to Total Positivity of Combinatorial Matrices and Real-Rootedness of Polynomials*, arXiv:2512.08369, Theorem 1.3, allows noninvertible lower-triangular arrays. It says that if the left production matrix is TN, then the array is TN, its reversal is TN, and every row polynomial is real-rooted.

That is directly relevant, but the natural left production matrix for the target is not TN. Exact witnesses already recorded in `ordered-big-blocks/literature.md` §8 are:

* for \(j=2\), rows \(\{2,3,4,5\}\) and columns \(\{0,1,2,3\}\) give determinant \(-1/4320\);
* for \(j=3\), rows \(\{3,4\}\) and columns \(\{0,1\}\) give determinant \(-1/24\).

Therefore Theorem 1.3 cannot be applied through the standard natural factorization. A different factorization/network, or a target-specific shear argument, is still required.

---

## 5. Exact combinatorial reductions that remain usable

### 5.1 Weighted path model

There is an exact weighted path representation of \(b_{n,m}\). A state records the current block size \(s\). Consuming one label changes \(s\to s+1\) and has weight \(1/(s+1)\). Closing a block of size \(s<j\) returns to the boundary state without changing the large-block count; closing a block of size \(s\ge j\) returns to the boundary state and increases that count by one.

For a path encoding an ordered partition with block sizes \(r_1,\ldots,r_k\), the edge-weight product is

\[
 \prod_{i=1}^k\prod_{s=1}^{r_i}\frac1{s}
 =\frac1{\prod_i r_i!}.
\]

Multiplication by \(n!\) gives exactly the multinomial number of label assignments. Hence

\[
 a_{n,m}=n!\sum_{\substack{\text{block-size paths of length }n\\
 \text{with }m\text{ large closures}}}
 \prod_i\frac1{|B_i|!}.
\tag{17}
\]

This is a proved path model, but not yet a planar-network proof. The missing step is a planar embedding, or a crossing involution with the correct sign and weight preservation, that realizes arbitrary sheared minors as nonintersecting-path determinants.

### 5.2 Sum over the total number of blocks

Let \(U_{n,k,m}\) count unordered set partitions of \([n]\) with \(k\) total blocks, exactly \(m\) of size at least \(j\). Ordering all \(k\) blocks gives

\[
 a_{n,m}=\sum_{k\ge m}k!\,U_{n,k,m}.
\tag{18}
\]

For fixed \(k\),

\[
 k!U_{n,k,m}
 =\binom km n![z^n]A_j(z)^{k-m}C_j(z)^m.
\tag{19}
\]

Therefore, if every fixed-\(k\) slice in (19) were TN, their sum (18) would be TN. This is an exact reduction. No proof that all the slices (19) are TN was obtained. Positivity of each coefficient is immediate; joint minor positivity is the missing statement.

---

## 6. Precise remaining gap

The matrix route is closed only after one of the following is proved:

1. **Target-specific shear lemma.** For \(R^{(j)}=(D_j,h_j)\), prove
   \[
   \det R^{(j)}[\,\{n_i-jm_\ell\},\{m_\ell\}\,]_{i,\ell=1}^r\ge0
   \tag{20}
   \]
   for all increasing row and column index sets. Ordinary TN of \(R^{(j)}\) is insufficient, as (14) proves.

2. **Alternative planar network.** Embed the weighted path model (17), or another exact model, in a planar network with sources and sinks ordered so that Lindström–Gessel–Viennot gives (20). For support-strict TN, also prove uniqueness of the nonintersecting matching whenever the support admits one.

3. **Alternative production matrix.** Find a factorization whose production matrix is TN. The natural left production matrix fails by the exact minors listed in §4.3.

4. **Row PF only.** A weaker sufficient goal is to prove that every row coefficient sequence is PF. By the Aissen–Edrei–Schoenberg–Whitney theorem this is equivalent to real-rootedness of every \(O_{n,j}\), but it still does not prove all cross-row minors.

The stretched/sheared-Riordan lead is therefore useful for formulating the exact obstruction, but it does not close the proof. The strongest general statement is negative: TN of a proper Riordan array does not imply TN of its shear.

---

## 7. Sources

1. Yu Chen and Wenxiang Zhu, *A Snail Race Problem*, arXiv:2607.02085.  
   https://arxiv.org/abs/2607.02085

2. Aoife Hennessy, *A Study of Riordan Arrays with Applications to Continued Fractions, Orthogonal Polynomials and Lattice Paths* (2011).  
   http://repository.wit.ie/1693/1/AoifeThesis.pdf

3. Tian-Xiao He and Renzo Sprugnoli, *The Vertical Recursive Relation of Riordan Arrays and Their Matrix Representation*, arXiv:2206.12898.  
   https://arxiv.org/abs/2206.12898

4. Paul Barry, *A New Group in the Riordan Family of Matrix Groups: the Sprugnoli Group*, arXiv:2605.16633.  
   https://arxiv.org/abs/2605.16633

5. Tian-Xiao He, *The Multiple Riordan Group and the Multiple Riordan Type Arrays*, especially Section 6 and Theorem 6.4, arXiv:2504.04049.  
   https://arxiv.org/abs/2504.04049

6. Xi Chen, Huyile Liang, and Yi Wang, *Total positivity of Riordan arrays*, especially Theorem 2.1, arXiv:1601.05637.  
   https://arxiv.org/abs/1601.05637

7. Xi Chen, Lang Fu, and Jiajie Ruan, *A Unified Planar Network Approach to Total Positivity of Combinatorial Matrices and Real-Rootedness of Polynomials*, especially Theorems 1.3–1.4, arXiv:2512.08369.  
   https://arxiv.org/abs/2512.08369
