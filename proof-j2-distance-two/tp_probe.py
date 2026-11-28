"""Check total positivity of the associated-Stirling triangle A_{m,r} and
of coefficient matrices of Q_{m,a}. Finite evidence only."""
import itertools
import sympy as sp
from explore import astirling_table, Q_coeffs

def minors_nonneg(M, maxk=3):
    nrows, ncols = M.shape
    bad = []
    for k in range(1, maxk + 1):
        for rows in itertools.combinations(range(nrows), k):
            for cols in itertools.combinations(range(ncols), k):
                det = sp.Matrix([[M[i, j] for j in cols] for i in rows]).det()
                if det < 0:
                    bad.append((rows, cols, det))
    return bad

def main():
    A = astirling_table(12)
    # build matrix A[m][r] for m=0..11, r=0..5
    M = sp.Matrix([[A[m][r] for r in range(6)] for m in range(12)])
    print("associated-Stirling matrix (m rows 0..11, r cols 0..5):")
    sp.print_latex(M)
    bad = minors_nonneg(M, 3)
    print("negative 1x1/2x2/3x3 minors:", len(bad))
    for b in bad[:10]:
        print("  ", b)

if __name__ == "__main__":
    main()
