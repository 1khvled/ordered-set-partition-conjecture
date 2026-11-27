import OrderedBigBlocks.Basic
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Algebra.Polynomial.Sturm.Sequence
import Mathlib.Data.List.SignVariations
import Mathlib.Topology.Algebra.Polynomial
import Mathlib.Topology.Order.IntermediateValue

/-!
# Exact Sturm chains and finite real-rootedness certificates

Mathlib currently supplies signed Sturm remainder sequences and sign variations,
but not the analytic Sturm root-counting theorem.  This file therefore keeps the
logical status explicit:

* `computableSturm` is an executable rational Sturm chain and is proved equal to
  Mathlib's `Polynomial.sturmSeq p p.derivative`;
* `sturmDrop` is an exact rational endpoint calculation;
* a `RealRootCertificate` additionally requires a rational sign change in each
  of `natDegree p` pairwise disjoint intervals.  Its soundness proof uses the
  intermediate value theorem and the degree bound, **not** an unproved Sturm
  theorem.

Thus every accepted certificate is a proof that all roots are real.  The Sturm
drop is retained as a exact independent audit of each interval.
-/

namespace OrderedBigBlocks

open Polynomial Set

/-- Fueled, executable construction of a signed Sturm remainder sequence. -/
def sturmChainAux : ℕ → ℚ[X] → ℚ[X] → List ℚ[X]
  | 0, p, _ => if p = 0 then [] else [p]
  | fuel + 1, p, q =>
      if p = 0 then [] else if q = 0 then [p]
      else p :: sturmChainAux fuel q (-p % q)

/-- With `fuel > natDegree q`, the fueled chain is Mathlib's Sturm sequence. -/
theorem sturmChainAux_eq_sturmSeq :
    ∀ (fuel : ℕ) (p q : ℚ[X]), q.natDegree < fuel →
      sturmChainAux fuel p q = p.sturmSeq q := by
  intro fuel
  induction fuel with
  | zero =>
      intro p q h
      omega
  | succ fuel ih =>
      intro p q hq
      by_cases hp : p = 0
      · subst p
        simp [sturmChainAux]
      · by_cases hq0 : q = 0
        · subst q
          simp [sturmChainAux, hp, Polynomial.sturmSeq_zero_right]
        · rw [sturmChainAux, if_neg hp, if_neg hq0,
            Polynomial.sturmSeq_cons hp, ih q (-p % q)]
          by_cases hmod : -p % q = 0
          · simp [hmod]
          · have hdeg := natDegree_lt_natDegree hmod (degree_mod_lt (-p) hq0)
            omega

/-- The complete executable Sturm chain of `p` and its derivative. -/
def computableSturm (p : ℚ[X]) : List ℚ[X] :=
  sturmChainAux (p.derivative.natDegree + 1) p p.derivative

/-- The executable chain agrees with Mathlib's signed remainder sequence. -/
theorem computableSturm_eq_sturmSeq (p : ℚ[X]) :
    computableSturm p = p.sturmSeq p.derivative := by
  exact sturmChainAux_eq_sturmSeq _ p p.derivative (Nat.lt_succ_self _)

/-- Every entry of a Sturm chain is nonzero. -/
theorem computableSturm_ne_zero {p s : ℚ[X]} (hs : s ∈ computableSturm p) : s ≠ 0 := by
  rw [computableSturm_eq_sturmSeq] at hs
  exact Polynomial.ne_zero_of_mem_sturmSeq hs

/-- Exact sign variations of a Sturm chain at a rational point. -/
def sturmVariations (p : ℚ[X]) (x : ℚ) : ℕ :=
  ((computableSturm p).map (fun q => q.eval x)).signVariations

/-- Exact endpoint drop, as an integer so that a negative drop cannot be
silently truncated by natural-number subtraction. -/
def sturmDrop (p : ℚ[X]) (a b : ℚ) : ℤ :=
  (sturmVariations p a : ℤ) - (sturmVariations p b : ℤ)

/-- Cast a rational polynomial to a real polynomial. -/
def toReal (p : ℚ[X]) : ℝ[X] := p.map (Rat.castHom ℝ)

@[simp] theorem toReal_eval_rat (p : ℚ[X]) (a : ℚ) :
    (toReal p).eval (a : ℝ) = (p.eval a : ℝ) := by
  simp [toReal]

@[simp] theorem toReal_natDegree {p : ℚ[X]} (hp : p ≠ 0) :
    (toReal p).natDegree = p.natDegree := by
  exact natDegree_map_eq_of_injective Rat.cast_injective p

theorem toReal_ne_zero {p : ℚ[X]} (hp : p ≠ 0) : toReal p ≠ 0 := by
  exact (Polynomial.map_ne_zero_iff Rat.cast_injective).mpr hp

/-- Opposite endpoint signs give a root in the open rational interval. -/
theorem exists_root_of_sign_change {p : ℚ[X]} {a b : ℚ} (hab : a < b)
    (hsign : p.eval a * p.eval b < 0) :
    ∃ x : ℝ, x ∈ Set.Ioo (a : ℝ) (b : ℝ) ∧ (toReal p).eval x = 0 := by
  have hsignR : (p.eval a : ℝ) * (p.eval b : ℝ) < 0 := by exact_mod_cast hsign
  have hzero : (0 : ℝ) ∈ Set.Ioo ((toReal p).eval (a : ℝ)) ((toReal p).eval (b : ℝ)) := by
    rw [toReal_eval_rat, toReal_eval_rat]
    rcases mul_neg_iff.mp hsignR with h | h
    · exact ⟨h.2, h.1⟩
    · exact ⟨h.1, h.2⟩
  have hcont : ContinuousOn (fun x : ℝ => (toReal p).eval x)
      (Set.Icc (a : ℝ) (b : ℝ)) := (toReal p).continuous.continuousOn
  rcases mul_neg_iff.mp hsignR with h | h
  · obtain ⟨x, hx, hroot⟩ := intermediate_value_Ioo' (show (a : ℝ) ≤ b by exact_mod_cast hab.le)
      hcont (by simpa [toReal_eval_rat] using hzero)
    exact ⟨x, hx, hroot⟩
  · obtain ⟨x, hx, hroot⟩ := intermediate_value_Ioo (show (a : ℝ) ≤ b by exact_mod_cast hab.le)
      hcont (by simpa [toReal_eval_rat] using hzero)
    exact ⟨x, hx, hroot⟩

/-- A finite exact certificate that a nonzero rational polynomial splits over
`ℝ`.  Each interval carries both the Sturm endpoint-drop audit and an actual
sign change. -/
structure RealRootCertificate (p : ℚ[X]) where
  nonzero : p ≠ 0
  intervals : List (ℚ × ℚ)
  length_eq_degree : intervals.length = p.natDegree
  ordered : ∀ I ∈ intervals, I.1 < I.2
  disjoint : intervals.Pairwise (fun I J => I.2 ≤ J.1 ∨ J.2 ≤ I.1)
  sign_change : ∀ I ∈ intervals, p.eval I.1 * p.eval I.2 < 0
  sturm_drop : ∀ I ∈ intervals, sturmDrop p I.1 I.2 = 1

/-- Choose one IVT root in each certified interval. -/
noncomputable def certificateRoot {p : ℚ[X]} (c : RealRootCertificate p)
    (i : Fin c.intervals.length) : ℝ :=
  Classical.choose (exists_root_of_sign_change
    (c.ordered (c.intervals.get i) (List.get_mem _ _))
    (c.sign_change (c.intervals.get i) (List.get_mem _ _)))

/-- The chosen root lies strictly inside its interval and is a root. -/
theorem certificateRoot_spec {p : ℚ[X]} (c : RealRootCertificate p)
    (i : Fin c.intervals.length) :
    certificateRoot c i ∈ Set.Ioo
      ((c.intervals.get i).1 : ℝ) ((c.intervals.get i).2 : ℝ) ∧
      (toReal p).eval (certificateRoot c i) = 0 := by
  exact Classical.choose_spec (exists_root_of_sign_change
    (c.ordered (c.intervals.get i) (List.get_mem _ _))
    (c.sign_change (c.intervals.get i) (List.get_mem _ _)))

/-- Pairwise disjoint open intervals force the chosen roots to be distinct. -/
theorem certificateRoot_injective {p : ℚ[X]} (c : RealRootCertificate p) :
    Function.Injective (certificateRoot c) := by
  intro i j hij
  by_contra hne
  have hroot : certificateRoot c i = certificateRoot c j := hij
  have hi := (certificateRoot_spec c i).1
  have hj := (certificateRoot_spec c j).1
  rcases lt_or_gt_of_ne hne with hlt | hlt
  · have hdisj := c.disjoint.rel_getElem_of_lt i.isLt j.isLt hlt
    have hdisjR : ((c.intervals.get i).2 : ℝ) ≤ ((c.intervals.get j).1 : ℝ) ∨
        ((c.intervals.get j).2 : ℝ) ≤ ((c.intervals.get i).1 : ℝ) := by
      exact_mod_cast hdisj
    rcases hdisjR with h | h <;> linarith
  · have hdisj := c.disjoint.rel_getElem_of_lt j.isLt i.isLt hlt
    have hdisjR : ((c.intervals.get j).2 : ℝ) ≤ ((c.intervals.get i).1 : ℝ) ∨
        ((c.intervals.get i).2 : ℝ) ≤ ((c.intervals.get j).1 : ℝ) := by
      exact_mod_cast hdisj
    rcases hdisjR with h | h <;> linarith

/-- Soundness of the finite certificate: all roots of `p` are real. -/
theorem RealRootCertificate.splits {p : ℚ[X]} (c : RealRootCertificate p) :
    (toReal p).Splits := by
  let S : Finset ℝ := Finset.univ.image (certificateRoot c)
  have hcard : S.card = c.intervals.length := by
    rw [Finset.card_image_of_injective _ (certificateRoot_injective c)]
    simp
  have heval : ∀ x ∈ S, (toReal p).eval x = 0 := by
    intro x hx
    rcases Finset.mem_image.mp hx with ⟨i, _, rfl⟩
    exact (certificateRoot_spec c i).2
  have hp := toReal_ne_zero c.nonzero
  have hdegree : (toReal p).natDegree ≤ S.card := by
    rw [toReal_natDegree c.nonzero, hcard, c.length_eq_degree]
  have hroots : (toReal p).roots = S.val :=
    roots_eq_of_natDegree_le_card_of_ne_zero heval hdegree hp
  apply (Polynomial.splits_iff_card_roots hp).mpr
  rw [hroots, Multiset.card_coe, Finset.card_val, hcard, c.length_eq_degree,
    toReal_natDegree c.nonzero]

/-- All roots of an integer polynomial are real, expressed using the direct
integer-to-real cast. -/
def IntegerRealRooted (p : ℤ[X]) : Prop := (p.map (Int.castRingHom ℝ)).Splits

theorem toRat_toReal (p : ℤ[X]) :
    toReal (toRat p) = p.map (Int.castRingHom ℝ) := by
  ext n
  simp [toReal, toRat, Polynomial.map_map]

/-- Transport a rational finite certificate to an integer polynomial. -/
theorem RealRootCertificate.integerRealRooted {p : ℤ[X]}
    (c : RealRootCertificate (toRat p)) : IntegerRealRooted p := by
  rw [← toRat_toReal]
  exact c.splits

end OrderedBigBlocks
