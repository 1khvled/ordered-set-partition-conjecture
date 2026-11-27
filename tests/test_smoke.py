"""Small deterministic regression tests; the large audits live in numerical/."""
import math
import unittest

import sympy as sp
from flint import fmpz_poly

from interlace import check_pair
from scan import nonreal_roots, ordered_rows, sturm_variations


class OrderedBigBlockTests(unittest.TestCase):
    def test_j3_initial_rows(self):
        got = [c for _, c in ordered_rows(3, 6)]
        self.assertEqual(got, [
            [1], [1], [3], [12, 1], [66, 9], [450, 91], [3690, 973, 20]
        ])

    def test_t1_is_ordered_bell(self):
        # At t=1 every block is weighted equally, independently of j.
        bell = [1, 1, 3, 13, 75, 541, 4683, 47293]
        for j in range(2, 9):
            sums = [sum(c) for _, c in ordered_rows(j, len(bell) - 1)]
            self.assertEqual(sums, bell)

    def test_universal_differential_recurrence(self):
        t = sp.Symbol("t")
        for j in range(2, 9):
            rows = [sum(a * t**m for m, a in enumerate(c)) for _, c in ordered_rows(j, 25)]
            for n in range(1, len(rows)):
                rhs = n * rows[n - 1] + 2 * t * sp.diff(rows[n - 1], t)
                for r in range(2, j):
                    rhs += math.comb(n, r) * rows[n - r]
                if n >= j:
                    p = rows[n - j]
                    rhs += math.comb(n - 1, j - 1) * t * (p + (t - 1) * sp.diff(p, t))
                self.assertEqual(sp.expand(rows[n] - rhs), 0, (j, n))

    def test_exact_real_root_count(self):
        self.assertEqual(sturm_variations([-1, 0, 1], minus_inf=True) - sturm_variations([-1, 0, 1]), 2)
        self.assertEqual(sturm_variations([1, 0, 1], minus_inf=True) - sturm_variations([1, 0, 1]), 0)
        self.assertTrue(nonreal_roots([1, 0, 1]))
        self.assertFalse(nonreal_roots([-1, 0, 1]))

    def test_interlacing_example(self):
        # 4x^2-16x+15 has roots 3/2,5/2; the cubic has roots 1,2,3.
        ok, rec = check_pair([15, -16, 4], [-6, 11, -6, 1], 128, wronskian=True)
        self.assertTrue(ok, rec)
        self.assertFalse(rec["common_root"])

    def test_equal_degree_interlacing_checks_last_root(self):
        # f has roots -2,-1; g has roots -3/2,-6/5.  The first gap is valid,
        # but f_2=-1 is not <= g_2=-6/5, so the final inequality fails.
        ok, _ = check_pair([2, 3, 1], [18, 27, 10], 128)
        self.assertFalse(ok)

    def test_nonreal_polynomial_is_rejected(self):
        with self.assertRaises(ValueError):
            check_pair([1, 0, 1], [1, 0, 1], 128)


if __name__ == "__main__":
    unittest.main()
