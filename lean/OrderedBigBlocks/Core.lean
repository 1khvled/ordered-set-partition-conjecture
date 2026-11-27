/-!
A dependency-free executable formalization of the ordered big-block recurrence.
Polynomials are lists of integer coefficients in ascending degree order.
This file does not prove real-rootedness; it provides a small compiling core
that can be checked without downloading Mathlib.
-/

namespace OrderedBigBlocks

/-- Binomial coefficients, defined recursively so no Mathlib import is needed. -/
def choose : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => choose n k + choose n (k + 1)

/-- Addition of coefficient lists. -/
def padd : List Int → List Int → List Int
  | [], ys => ys
  | xs, [] => xs
  | x :: xs, y :: ys => (x + y) :: padd xs ys

/-- Scalar multiplication of a coefficient list. -/
def pscale : Int → List Int → List Int
  | _, [] => []
  | c, x :: xs => c * x :: pscale c xs

/-- Shift a nonempty coefficient list up by one degree. -/
def pshift : List Int → List Int
  | [] => []
  | p => 0 :: p

/-- Multiplication of coefficient lists. -/
def pmul : List Int → List Int → List Int
  | [], _ => []
  | _, [] => []
  | x :: xs, y :: ys =>
      (x * y) :: padd (pscale x ys) (pshift (pmul xs (y :: ys)))

/-- Multiplication by the polynomial variable `t`. -/
def pmulX : List Int → List Int := pshift

/-- Derivative helper; the index is the current exponent. -/
def derivFrom : Nat → List Int → List Int
  | _, [] => []
  | i, x :: xs => (Int.ofNat i * x) :: derivFrom (i + 1) xs

/-- Formal derivative of a coefficient list. -/
def pderiv : List Int → List Int
  | p => derivFrom 1 p

/-- The weight of a block of size `r`. -/
def blockWeight (j r : Nat) : List Int :=
  if 1 ≤ r ∧ r < j then [1]
  else if j ≤ r then [0, 1]
  else [0]

/-- Remove trailing zero coefficients, retaining `[0]` for the zero polynomial. -/
def pnormalize : List Int → List Int
  | [] => [0]
  | p =>
      let r := p.reverse.dropWhile (· == 0)
      if r.isEmpty then [0] else r.reverse

/-- One row of the first-block convolution. -/
def nextRow (j n : Nat) (rows : List (List Int)) : List Int :=
  pnormalize <|
    (List.range (n + 1)).map (fun r =>
      pscale (Int.ofNat (choose n r))
        (pmul (blockWeight j r) (rows.getD (n - r) []))
    ) |>.foldr padd []

/-- Rows `O_{0,j}, ..., O_{n,j}`. -/
def buildRows (j : Nat) : Nat → List (List Int)
  | 0 => [[1]]
  | n + 1 =>
      let rows := buildRows j n
      rows ++ [nextRow j (n + 1) rows]

@[simp] theorem buildRows_zero (j : Nat) : buildRows j 0 = [[1]] := rfl

/-- The published initial rows for `j = 3`, in dense coefficient-list form. -/
example : buildRows 3 6 =
    [[1], [1], [3], [12, 0, 1], [66, 0, 9], [450, 0, 91],
      [3690, 0, 973, 0, 20]] := by
  native_decide

/-- At `t = 1` every block has equal weight, independently of `j`. -/
def rowSum (p : List Int) : Int := p.foldr (· + ·) 0

example : (buildRows 7 6).map rowSum = [1, 1, 3, 13, 75, 541, 4683] := by
  native_decide

end OrderedBigBlocks
