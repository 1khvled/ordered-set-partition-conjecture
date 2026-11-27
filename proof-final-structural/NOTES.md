# Final structural track: decomposition, inner family, and remaining gap

Date: 2026-11-27  
Scope: only `ordered-big-blocks/proof-final-structural` was created or edited.

## 1. Structural decomposition (proved)

Write
\[
A_j(z)=\sum_{r=1}^{j-1}\frac{z^r}{r!},\qquad
C_j(z)=e^z-\sum_{r=0}^{j-1}\frac{z^r}{r!},\qquad
D_j(z)=\frac1{1-A_j(z)},\qquad
g_j(z)=\frac{C_j(z)}{1-A_j(z)}.
\]
Then
\[
\boxed{F_j(z,t)=\frac1{1-A_j-tC_j}=\frac{D_j(z)}{1-t\,g_j(z)}.} \tag{1}
\]
The target row is the convolution
\[
\boxed{O_{n,j}(t)=n!\sum_{r=0}^n d_{n-r}^{(j)}\,H_{r,j}(t),} \tag{2}
\]
where \(d_k=[z^k]D_j\) and \(H_{r,j}(t)=r![z^r]\frac1{1-tg_j(z)}\) is the row of the inner family
\[
\mathcal H_j(z,t)=\frac1{1-t\,g_j(z)}. \tag{3}
\]
Thus the problem splits into two independent structural statements:

- **(H)** every \(H_{r,j}\) is real-rooted (with negative roots) and consecutive rows interlace;
- **(T)** the convolution transform (2) with kernel \(D_j\) preserves real-rootedness and interlacing.

Neither statement follows from the previously recorded routes.

## 2. Inner family \(\mathcal H_j\)

### 2.1 Brändén–Leite route does not explain it

Brändén–Leite Theorem 4.4 would give (H) if \(g_j\) were a Pólya-frequency (PF) power series. It is not:

- For \(j=2\), the shifted kernel \(g_2/z^2\) fails PF at size \(4\) (minor
  rows \(\{1,2,3,4\}\), cols \(\{0,1,2,3\}\) equals \(-13/103680\)).
- For \(j\ge3\), \(g_j/z^j\) already fails PF at size \(2\)
  (e.g. \(j=3\): rows \(\{1,2\}\), cols \(\{0,1\}\) equals \(-19/2880\)).

The constituents are also not PF: \(C_j\) fails PF at size \(3\) for every
\(j\ge2\), and \(D_j=1/(1-A_j)\) fails PF at size \(2\) for every \(j\ge3\).

### 2.2 Nevertheless \(\mathcal H_j\) is real-rooted

Exact computation (SymPy, integer/rational arithmetic, \(j=2,\dots,5\),
\(r\le 14\)) shows:

- every \(H_{r,j}\) has only real roots;
- all discriminants are nonnegative;
- degrees are \(\lfloor r/j\rfloor\).

This is a new, unexplained structural fact: \(g_j\) is not PF, yet the
renewal \(1/(1-tg_j)\) has real-rooted rows. No Favard three-term recurrence
was found for any fixed residue class modulo \(j\): attempts to fit
\(H_{k+1}=(a_k t+b_k)H_k+c_k H_{k-1}\) failed at \(k=2\) for
\(j=2,3\). Thus \(\mathcal H_j\) is not a classical orthogonal sequence in
\(t\).

### 2.3 Hyper-Bessel / multivariate stability

- Truncated \(B_j=A_j+tC_j\) is not real-rooted for any tested \(t>0\), so a
  direct stability of the block EGF does not apply.
- No hyper-Bessel or reciprocal-exponential-polynomial identification was
  found: the degrees grow as \(\lfloor n/j\rfloor\), not \(n\), and the
  coefficients are not hypergeometric in a single index.

## 3. Convolution transform (T)

For \(j=2\), \(D_2=1/(1-z)\) has all coefficients \(1\), so
\(O_{n,2}=n!\sum_{r=0}^n H_{r,2}\). Partial sums of real-rooted polynomials
are not generally real-rooted, yet this one is (verified for \(n\le 40\)).
For \(j\ge3\), \(D_j\) is itself not PF, and the convolution (2) still
produces real-rooted \(O_{n,j}\) (verified for \(j\le 8\), \(n\le 40\)).

No general theorem was located that guarantees (T) from (H) alone. The
kernel \(D_j\) is not a multiplier sequence in the standard basis, and the
transform is not a simple binomial/Eulerian transform with a known
real-rootedness-preservation certificate.

## 4. Precise remaining gap

A complete proof of Alexandersson's conjecture now requires two new lemmas,
neither of which follows from the recorded (and obstructed) classes:

1. **Inner real-rootedness (H).** Prove that for every \(j\ge2\) and
   \(r\ge0\), \(H_{r,j}(t)=r![z^r](1-tg_j(z))^{-1}\) has only negative real
   roots and that consecutive \(r\)'s interlace. The PF hypothesis fails, so
   this needs a new property of \(g_j\) (for instance a
   complex-zero-decreasing or sector-stability property) or a direct
   recurrence with positive coefficients.

2. **Transform preservation (T).** Prove that if a sequence \((H_r)\)
   satisfies (H), then \(n!\sum_{r}d_{n-r}H_r\) is real-rooted and
   interlacing, for the specific (non-PF) kernels \(D_j\).

Equivalently, the original weighted root-separation inequality WRS remains the
single quantitative obstruction at the recurrence level; the decomposition
(1)–(2) shows that WRS is not an artifact of the recurrence's signed term
but reflects the interaction between the non-PF inner kernel \(g_j\) and the
non-PF outer kernel \(D_j\).

## 5. Files

- `explore.py`, `recurrence_test.py`, `pf_check.py`, `pf_check2.py`,
  `pf_check3.py`, `pf_g2.py`, `pf_g3_shift.py`, `stab_test.py`,
  `H_check.py`, `H_deg.py`, `H_disc.py`, `j2_test.py` — exact
  computational probes supporting the statements above.
- `NOTES.md` — this note.
