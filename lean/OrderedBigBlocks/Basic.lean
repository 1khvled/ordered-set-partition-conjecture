import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Data.Nat.Choose.Basic

/-!
# Ordered big-block polynomials

`orderedBigBlock j n` is the polynomial in `t` obtained by choosing the first
block of an ordered set partition of `n` elements.  A block of size `r < j`
has weight `1`, and a block of size `r ≥ j` has weight `t`.

The parameter `j` is positive.  This excludes the degenerate `j = 0` recurrence,
whose second sum would contain `r = 0` and would not be well-founded.
-/

namespace OrderedBigBlocks

open Polynomial
open scoped BigOperators

/-- The weight of a block of size `r`.  Out-of-range binomial coefficients are
zero, so the small-block sum may safely be written as `Icc 1 (j - 1)`. -/
def blockWeight (j : ℕ+) (r : ℕ) : ℤ[X] :=
  if 1 ≤ r ∧ r < (j : ℕ) then 1 else if (j : ℕ) ≤ r then X else 0

/-- The ordered big-block polynomial `O_{n,j}(t)`.  The two sums are exactly
the first-block recurrence

`O_{n,j} = ∑_{r=1}^{j-1} choose n r O_{n-r,j}
          + t ∑_{r=j}^{n} choose n r O_{n-r,j}`.

Here `X` is the polynomial variable `t`, not the EGF variable. -/
def orderedBigBlock (j : ℕ+) : ℕ → ℤ[X]
  | 0 => 1
  | n + 1 =>
      (∑ r ∈ Finset.Icc 1 ((j : ℕ) - 1),
        (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r)) +
      X * (∑ r ∈ Finset.Icc (j : ℕ) (n + 1),
        (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r))
termination_by n => n
decreasing_by
  all_goals
    have hr : 1 ≤ r ∧ r ≤ n + 1 := by
      simpa only [Finset.mem_Icc] using ‹r ∈ _›
    have hj : 1 ≤ (j : ℕ) := j.pos
    omega

@[simp] theorem orderedBigBlock_zero (j : ℕ+) :
    orderedBigBlock j 0 = 1 := rfl

/-- The defining convolution recurrence, stated at a successor index. -/
theorem orderedBigBlock_succ (j : ℕ+) (n : ℕ) :
    orderedBigBlock j (n + 1) =
      (∑ r ∈ Finset.Icc 1 ((j : ℕ) - 1),
        (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r)) +
      X * (∑ r ∈ Finset.Icc (j : ℕ) (n + 1),
        (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r)) := by
  rw [orderedBigBlock]

@[simp] theorem blockWeight_zero (j : ℕ+) : blockWeight j 0 = 0 := by
  simp [blockWeight, j.pos]

theorem blockWeight_of_lt (j : ℕ+) {r : ℕ} (hr : 1 ≤ r) (hrj : r < (j : ℕ)) :
    blockWeight j r = 1 := by
  simp [blockWeight, hr, hrj]

theorem blockWeight_of_le (j : ℕ+) {r : ℕ} (hrj : (j : ℕ) ≤ r) :
    blockWeight j r = X := by
  have hr : 1 ≤ r := j.pos.trans hrj
  simp [blockWeight, hr, not_lt.mpr hrj, hrj]

/-- A single-sum form of the recurrence.  This is useful for coefficient
arguments and makes the zero convention for `Nat.choose` explicit. -/
theorem orderedBigBlock_succ_weight (j : ℕ+) (n : ℕ) :
    orderedBigBlock j (n + 1) =
      ∑ r ∈ Finset.range (n + 2),
        (Nat.choose (n + 1) r : ℤ[X]) * blockWeight j r *
          orderedBigBlock j (n + 1 - r) := by
  let f : ℕ → ℤ[X] := fun r =>
    (Nat.choose (n + 1) r : ℤ[X]) * blockWeight j r *
      orderedBigBlock j (n + 1 - r)
  have hbig : (∑ r ∈ Finset.range ((j : ℕ) + n + 2), f r) =
      (∑ r ∈ Finset.Icc 1 ((j : ℕ) - 1), f r) +
        X * (∑ r ∈ Finset.Icc (j : ℕ) (n + 1), f r) := by
    rw [Finset.sum_range_add_sum_Ico f (j : ℕ) (by omega)]
    have hleft : (∑ r ∈ Finset.range (j : ℕ), f r) =
        ∑ r ∈ Finset.Icc 1 ((j : ℕ) - 1), f r := by
      rw [Finset.sum_range_eq_add_Ico f (by omega : 1 ≤ (j : ℕ))]
      simp only [f, Nat.choose_zero_right, Nat.cast_zero, zero_mul, zero_add]
      apply Finset.sum_congr rfl
      intro r hr
      have hr' : 1 ≤ r ∧ r < (j : ℕ) := by
        simpa only [Finset.mem_Ico] using hr
      rw [blockWeight_of_lt j hr'.1 hr'.2, mul_one]
    have hright : (∑ r ∈ Finset.Ico (j : ℕ) ((j : ℕ) + n + 2), f r) =
        X * ∑ r ∈ Finset.Icc (j : ℕ) (n + 1), f r := by
      have hsub : Finset.Icc (j : ℕ) (n + 1) ⊆
          Finset.Ico (j : ℕ) ((j : ℕ) + n + 2) := by
        intro r hr
        simp only [Finset.mem_Icc] at hr
        simp only [Finset.mem_Ico]
        omega
      have hzero : ∀ r ∈ Finset.Ico (j : ℕ) ((j : ℕ) + n + 2),
          r ∉ Finset.Icc (j : ℕ) (n + 1) → f r = 0 := by
        intro r hr hnot
        have hr' : n + 1 < r := by
          simp only [Finset.mem_Ico] at hr
          simp only [Finset.mem_Icc, not_and_or, not_le] at hnot
          omega
        simp [f, Nat.choose_eq_zero_of_lt hr'.le]
      rw [← Finset.sum_subset hsub hzero, Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro r hr
      have hrj : (j : ℕ) ≤ r := (Finset.mem_Icc.mp hr).1
      simp only [f, blockWeight_of_le j hrj]
      ring
    rw [hleft, hright]
  have hsub : Finset.range (n + 2) ⊆
      Finset.range ((j : ℕ) + n + 2) := by
    intro r hr
    simp only [Finset.mem_range] at hr ⊢
    omega
  have hzero : ∀ r ∈ Finset.range ((j : ℕ) + n + 2),
      r ∉ Finset.range (n + 2) → f r = 0 := by
    intro r hr hnot
    have hr' : n + 1 < r := by
      simp only [Finset.mem_range, not_lt] at hnot ⊢
      omega
    simp [f, Nat.choose_eq_zero_of_lt hr'.le]
  rw [orderedBigBlock_succ]
  change (∑ r ∈ Finset.Icc 1 ((j : ℕ) - 1),
      (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r)) +
      X * (∑ r ∈ Finset.Icc (j : ℕ) (n + 1),
        (Nat.choose (n + 1) r : ℤ[X]) * orderedBigBlock j (n + 1 - r)) = _
  rw [← hbig, ← Finset.sum_subset hsub hzero]
  apply Finset.sum_congr rfl
  intro r hr
  by_cases hrj : r < (j : ℕ)
  · have hr1 : 1 ≤ r := by
      by_contra h
      have : r = 0 := by omega
      subst r
      simp
    rw [blockWeight_of_lt j hr1 hrj]
    ring
  · have hrj' : (j : ℕ) ≤ r := by omega
    rw [blockWeight_of_le j hrj']
    ring

/-! ## Exponential generating functions -/

/-- Cast an integer polynomial to a rational polynomial. -/
def toRat (p : ℤ[X]) : ℚ[X] := p.map (Int.castRingHom ℚ)

/-- The rational version of `blockWeight`. -/
def ratBlockWeight (j : ℕ+) (r : ℕ) : ℚ[X] := (blockWeight j r).map (Int.castRingHom ℚ)

@[simp] theorem toRat_blockWeight (j : ℕ+) (r : ℕ) :
    toRat (blockWeight j r) = ratBlockWeight j r := rfl

/-- The normalized coefficient `O_{n,j}(t) / n!`. -/
noncomputable def egfCoeff (j : ℕ+) (n : ℕ) : ℚ[X] :=
  (Nat.factorial n : ℚ)⁻¹ • toRat (orderedBigBlock j n)

/-- The block EGF `B_j(z,t)`. -/
noncomputable def blockEGF (j : ℕ+) : ℚ[X]⟦X⟧ :=
  PowerSeries.mk fun r => (Nat.factorial r : ℚ)⁻¹ • ratBlockWeight j r

/-- The ordered-partition EGF `∑ O_{n,j}(t) z^n/n!`. -/
noncomputable def orderedEGF (j : ℕ+) : ℚ[X]⟦X⟧ :=
  PowerSeries.mk (egfCoeff j)

/-- The constant coefficient of the ordered EGF is one. -/
@[simp] theorem orderedEGF_coeff_zero (j : ℕ+) :
    PowerSeries.coeff 0 (orderedEGF j) = 1 := by
  simp [orderedEGF, egfCoeff, toRat]

/-- The block EGF has no constant term. -/
@[simp] theorem blockEGF_coeff_zero (j : ℕ+) :
    PowerSeries.coeff 0 (blockEGF j) = 0 := by
  simp [blockEGF, ratBlockWeight]

private theorem factorial_scalar {n r : ℕ} (hr : r ≤ n) :
    (Nat.factorial n : ℚ) * (Nat.factorial r : ℚ)⁻¹ *
      (Nat.factorial (n - r) : ℚ)⁻¹ = (Nat.choose n r : ℚ) := by
  have hnat := Nat.choose_mul_factorial_mul_factorial hr
  have hq : (Nat.choose n r : ℚ) * (Nat.factorial r : ℚ) *
      (Nat.factorial (n - r) : ℚ) = (Nat.factorial n : ℚ) := by
    exact_mod_cast hnat
  rw [← hq]
  field_simp [Nat.factorial_ne_zero]

private theorem factorial_mul_egf_term (j : ℕ+) (n r : ℕ) (hr : r ≤ n + 1) :
    (Nat.factorial (n + 1) : ℚ[X]) *
      (((Nat.factorial r : ℚ)⁻¹ • ratBlockWeight j r) *
        ((Nat.factorial (n + 1 - r) : ℚ)⁻¹ •
          toRat (orderedBigBlock j (n + 1 - r)))) =
      (Nat.choose (n + 1) r : ℚ[X]) * ratBlockWeight j r *
        toRat (orderedBigBlock j (n + 1 - r)) := by
  simp only [smul_eq_mul]
  calc
    _ = ((Nat.factorial (n + 1) : ℚ) * (Nat.factorial r : ℚ)⁻¹ *
          (Nat.factorial (n + 1 - r) : ℚ)⁻¹) *
        (ratBlockWeight j r * toRat (orderedBigBlock j (n + 1 - r))) := by
      push_cast
      ring
    _ = _ := by
      rw [factorial_scalar hr]
      push_cast
      ring

/-- Multiplication by `(n+1)!` turns the EGF convolution into precisely the
integer first-block recurrence. -/
theorem factorial_mul_blockEGF_mul_orderedEGF_coeff (j : ℕ+) (n : ℕ) :
    (Nat.factorial (n + 1) : ℚ[X]) *
      PowerSeries.coeff (n + 1) (blockEGF j * orderedEGF j) =
      toRat (orderedBigBlock j (n + 1)) := by
  rw [PowerSeries.coeff_mul,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ]
  simp only [blockEGF, orderedEGF, PowerSeries.coeff_mk, egfCoeff]
  rw [Finset.mul_sum]
  calc
    _ = ∑ r ∈ Finset.range (n + 2),
        (Nat.choose (n + 1) r : ℚ[X]) * ratBlockWeight j r *
          toRat (orderedBigBlock j (n + 1 - r)) := by
      apply Finset.sum_congr rfl
      intro r hr
      exact factorial_mul_egf_term j n r (by simpa using hr)
    _ = toRat (orderedBigBlock j (n + 1)) := by
      have h := congrArg (Polynomial.map (Int.castRingHom ℚ))
        (orderedBigBlock_succ_weight j n)
      simpa [toRat, ratBlockWeight, Polynomial.map_sum, Polynomial.map_mul,
        Polynomial.map_natCast] using h

/-- Every positive-degree coefficient of `B_j(z,t) O_j(z,t)` vanishes. -/
theorem blockEGF_mul_orderedEGF_coeff_succ (j : ℕ+) (n : ℕ) :
    PowerSeries.coeff (n + 1) (blockEGF j * orderedEGF j) = 0 := by
  have h := factorial_mul_blockEGF_mul_orderedEGF_coeff j n
  have hf : (Nat.factorial (n + 1) : ℚ[X]) ≠ 0 := by
    exact_mod_cast Nat.factorial_ne_zero (n + 1)
  exact (mul_eq_zero.mp h).resolve_left hf

/-- The exact EGF equation `(1 - B_j) O_j = 1`. -/
theorem one_sub_blockEGF_mul_orderedEGF (j : ℕ+) :
    (1 - blockEGF j) * orderedEGF j = 1 := by
  ext n
  cases n with
  | zero => simp
  | succ n =>
      rw [PowerSeries.coeff_one, if_neg (Nat.succ_ne_zero n)]
      simp [sub_mul, blockEGF_mul_orderedEGF_coeff_succ]

end OrderedBigBlocks
