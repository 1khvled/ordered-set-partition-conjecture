# j=2 falling-basis TP2 transfer and Bezoutian obstruction

Date: 2026-11-28
Scope: only `proof-j2-falling-bezout` created/edited.
No installs/builds. No other agents/models.

## Outcome

**Theorem (falling-basis TP2 transfer).** For every \(n\ge0\) and every
\(i\ge0\),

\[
\boxed{\;
 c_{n,i}c_{n+1,i+1}-c_{n,i+1}c_{n+1,i}\ge0,
 \qquad c_{n,i}=\sum_{k\ge i}{k\brace i}d_{n,k},\;
 d_{n,k}=[t^k]O_{n,2}(1,t).
\;}
\]

More strongly, the falling-basis coefficient matrix
\(C=(c_{n,i})_{n,i\ge0}\) is TP2 in the full \(2\times2\)-minor sense:
for all \(n_1<n_2\), \(i_1<i_2\),

\[
\det C[\{n_1,n_2\},\{i_1,i_2\}]\ge0. \tag{T}
\]

The proof is Cauchy--Binet on the proved monomial TP2 matrix \(d\)
(`proof-j2-tp2-transfer/NOTES.md`) and the classical total nonnegativity of
the Stirling matrix of the second kind \(S_{k,i}={k\brace i}\).

**Bezoutian status (no closure).** The proved falling-basis TP2/ratio
monotonicity does **not** imply positive definiteness of the relevant
consecutive Bezoutian \(M(O_n,O_{n+1})\). This is an exact implication
obstruction, not a small-\(n\) sign failure: `proof-j2-base/NOTES.md` §5.4
gives the explicit feasible witness

\[
 P=10+11t+t^2,\qquad Q=10+22t+100t^2,
\]

whose falling-factorial coefficient ratios are strictly increasing while
\(\operatorname{Disc}_t(Q)=-3516<0\). Therefore monomial TP2 plus (H) does
not, by this route alone, prove Bezoutian positive definiteness. The
\(j=2\) consecutive Bezoutian/interlacing target remains open.

---

## 1. Inputs

### 1.1 Monomial TP2 (proved)

From `proof-j2-tp2-transfer/NOTES.md`, the matrix
\(d_{n,k}=[t^k]O_{n,2}(1,t)\) is TP2 in all \(2\times2\) minors:

\[
\det d[\{n_1,n_2\},\{k_1,k_2\}]
 = d_{n_1,k_1}d_{n_2,k_2}-d_{n_1,k_2}d_{n_2,k_1}\ge0,
 \qquad n_1<n_2,\ k_1<k_2. \tag{D}
\]

The proof uses the degenerate-Riordan left-production factorization and TP2
Cauchy--Binet. No (H) or (17) is used.

### 1.2 Associated-Stirling inputs (context)

The proved inequality (H) (`proof-u-operator/NOTES.md` §7) and the ratio
monotonicity (17) (`proof-j2-base/NOTES.md` §5.1) constrain the
associated-Stirling matrix \(A_{m,k}={m\brace k}_{\ge2}\). They are not used
in the TP2 transfer below. They remain relevant only to alternative routes;
§4 records why they do not close the Bezoutian target.

### 1.3 Falling-basis transform

By `proof-j2-base/NOTES.md` (18),

\[
 O_n(1,t)=\sum_{i\ge0}c_{n,i}(t)_i,
 \qquad c_{n,i}=\sum_{k\ge i}{k\brace i}d_{n,k}, \tag{1}
\]

where \((t)_i=t(t-1)\cdots(t-i+1)\) and \({k\brace i}=S2(k,i)\) is the
Stirling number of the second kind. For fixed \(n\), \(d_{n,k}=0\) unless
\(0\le k\le \lfloor n/2\rfloor\), so every sum in (1) is finite.

In matrix form, with \(S_{k,i}=S2(k,i)\) lower triangular,

\[
 C = D S, \qquad C_{n,i}=c_{n,i},\ D_{n,k}=d_{n,k}. \tag{2}
\]

---

## 2. Total nonnegativity of the Stirling second-kind matrix

**Lemma (S2 is TN).** The infinite lower-triangular matrix
\(S=(S2(k,i))_{k,i\ge0}\) is totally nonnegative: every minor
\(\det S[K,I]\ge0\).

*Justification.* The exponential Riordan array
\([g,f]=[1,e^z-1]\) has entries

\[
 r_{k,i}=\frac{k!}{i!}[z^k]g(z)f(z)^i
        =\frac{k!}{i!}[z^k](e^z-1)^i
        =S2(k,i),
\]

using \(\sum_{k\ge i}S2(k,i)z^k/k!=(e^z-1)^i/i!\). It is proper because
\(f'(0)=1\ne0\). Moreover \(g=1\) is PF\(_\infty\), \(f(z)=e^z-1\) has
compositional inverse \(\bar f(w)=\log(1+w)\), and
\(1/\bar f'(w)=1+w\) is PF\(_\infty\). Therefore Zhu's exponential Riordan
total-positivity theorem (Bao-Xuan Zhu, *Total positivity from the
exponential Riordan arrays*, arXiv:2006.14485v3, SIAM J. Discrete Math. 35
(2021), Thm. 3.6, case \(\lambda=0\)) applies and gives TP\(_\infty\) of
\([1,e^z-1]\). In particular every \(2\times2\) minor is nonnegative. This
is the same Stirling-triangle TN fact recorded in `literature.md` §3
(citing Chen--Fu--Ruan, Cor. 4.3). \(\square\)

We use only the \(2\times2\) consequence:

\[
\det S[\{k_1,k_2\},\{i_1,i_2\}]\ge0,
\qquad k_1<k_2,\ i_1<i_2. \tag{S}
\]

---

## 3. Cauchy--Binet transfer

Fix \(n_1<n_2\) and \(i_1<i_2\). Since \(d_{n,k}=0\) for
\(k>\lfloor n/2\rfloor\), the rows \(n_1,n_2\) of \(D\) are finitely
supported. Thus the minor in (T) is a finite Cauchy--Binet sum:

\[
\det C[\{n_1,n_2\},\{i_1,i_2\}]
 = \sum_{0\le k_1<k_2\le \lfloor n_2/2\rfloor}
   \det D[\{n_1,n_2\},\{k_1,k_2\}]
   \det S[\{k_1,k_2\},\{i_1,i_2\}]. \tag{3}
\]

Every term is a product of two nonnegative numbers by (D) and (S).
Therefore the sum is nonnegative. This proves (T). Taking adjacent rows and
columns gives the requested TP2 inequality. \(\blacksquare\)

---

## 4. Bezoutian obstruction (exact)

### 4.1 What TP2 of \(c\) means

For positive denominators, (T) on adjacent \(n,i\) is equivalent to ratio
monotonicity

\[
\frac{c_{n+1,i+1}}{c_{n,i+1}} \ge \frac{c_{n+1,i}}{c_{n,i}},
\qquad (\text{when }c_{n,i},c_{n,i+1}>0). \tag{4}
\]

Direction audit: multiplying (4) by the positive denominator
\(c_{n,i}c_{n,i+1}\) gives exactly the adjacent minor in (T), same sign.

### 4.2 Ratio monotonicity does not imply Bezoutian PD

`proof-j2-base/NOTES.md` §5.4 gives the exact counterexample

\[
 P(t)=10+11t+t^2,\qquad Q(t)=10+22t+100t^2.
\]

Falling-factorial coefficients are
\(P=10(1)+12(t)_1+1(t)_2\), \(Q=10(1)+122(t)_1+100(t)_2\), with ratios
\(122/12=61/6\) and \(100/1=100\), strictly increasing. However

\[
 \operatorname{Disc}_t(Q)=22^2-4\cdot10\cdot100=484-4000=-3516<0.
\]

Hence \(Q\) is not real-rooted and the Bezoutian/Hermite positivity needed
for consecutive interlacing cannot follow from ratio monotonicity alone.
The counterexample is a polynomial-pair implication obstruction: it shows
that TP2/ratio monotonicity of falling-basis coefficient rows is an
insufficient hypothesis for Bezoutian positive definiteness.

### 4.3 Consequence for the j=2 route

The proved inputs here are monomial TP2 (D), S2 TN (S), hence falling-basis
TP2 (T), and separately (H),(17) for associated Stirling numbers. None of
these, as used here, supplies the missing real-rootedness/Bezoutian
positivity certificate required by `proof-j2-base/NOTES.md` §5.4. In
particular, no all-\(n\) proof of \(M(O_n,O_{n+1})\succ0\) or of the root
velocity (7) is obtained.

---

## 5. Proof-policy audits

### 5.1 Direction audit

Target minor orientation is
\(c_{n,i}c_{n+1,i+1}-c_{n,i+1}c_{n+1,i}\ge0\), the same orientation as
(D). Cauchy--Binet (3) preserves this orientation because it is a
nonnegative sum of products of same-orientation minors.

Explicit feasible sample (exact rationals): \(n=4,i=1\),
\(c_{4,1}=51,c_{4,2}=6,c_{5,1}=421,c_{5,2}=110\), minor
\(51\cdot110-6\cdot421=3084>0\). Another sample \(n=3,i=0\):
\(c_{3,0}=6,c_{3,1}=7,c_{4,0}=24,c_{4,1}=51\), minor
\(6\cdot51-7\cdot24=138>0\).

### 5.2 Monotonicity audit

No monotonicity claim is used to prove (T). The only monotonicity-like
statement is ratio monotonicity (4), which is equivalent to the proved TP2
minor for positive denominators; its direction was audited in §4.1. The
5.4 obstruction is not a reversed monotonicity: it is an insufficient
implication.

### 5.3 Feasible-constraint falsification

The implication "TP2/ratio monotonicity of falling-basis rows => Bezoutian
PD / real-rootedness of the upper row" is falsified by the feasible
polynomial pair in §4.2, which satisfies the ratio condition and violates
Disc>=0. Therefore no proof may deduce Bezoutian PD from (T) alone.

The TP2 transfer itself uses only (D) and (S), both established; a
counterexample to (T) would require a negative Cauchy--Binet term, which is
impossible because every factor minor is nonnegative.

### 5.4 Boundary/degeneracy audit

- \(S2(k,i)=0\) for \(i>k\) and \(S2(k,k)=1\): lower-triangular, no
  division by zero.
- \(d_{n,k}=0\) for \(k>\lfloor n/2\rfloor\): finite support makes (3)
  finite.
- \(k=0/i=0\) columns included: \(d_{n,0}=n!\), \(c_{n,0}=d_{n,0}\).
- Zero minors are allowed (structural zeros from degree/support).
- The Cauchy--Binet sum is finite for each fixed minor; no infinite-sum
  convergence issue.

### 5.5 Symbolic residual audit

Identity (1): \(c_{n,i}=\sum_{k\ge i}{k\brace i}[t^k]O_n(1,t)\) follows from
\((x)^k=\sum_i S2(k,i)(x)_i\) and linearity; this is the standard inverse
Stirling expansion used in `proof-j2-base/NOTES.md` (18). The
Cauchy--Binet identity (3) is exact for finite matrices. `verify.py`
checks (1), (3), and the §4.2 witness with exact `Fraction` arithmetic;
finite checks are audits only.

---

## 6. Files

- `NOTES.md` — this proof/obstruction record.
- `verify.py` — exact finite audit: Cauchy--Binet residuals for (3),
  TP2 samples, and the §4.2 obstruction witness.

## 7. Status

- Falling-basis TP2: **proved for all \(n\)**.
- Bezoutian PD from monomial TP2 + (H): **not proved; implication obstructed
  by §4.2 / j2-base §5.4**.
- j=2 consecutive interlacing/root-velocity: **still open**.
