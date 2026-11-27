#!/usr/bin/env python3
"""Exact rational certificate that WRS fails for the actual j=2,n=63 row."""
from fractions import Fraction
from flint import fmpq, fmpz_poly
from scan import ordered_rows, sturm_variations


def exact_eval(coeffs, x):
    value = Fraction(0)
    for c in reversed(coeffs):
        value = value*x + c
    return value


def sign_near_center(coeffs, center, radius):
    value = exact_eval(coeffs,center)
    R = max(abs(center-radius),abs(center+radius),Fraction(1))
    derivative_bound = sum(k*abs(coeffs[k])*R**(k-1) for k in range(1,len(coeffs)))
    error = derivative_bound*radius
    if abs(value) > error:
        return (1 if value>0 else -1), value, error
    raise ArithmeticError("sign not certified")


def main():
    rows = {n:c for n,c in ordered_rows(2,63)}
    center = Fraction("-0.301049864539336319020771115593462774725514449468323282248633")
    eps = Fraction(1,10**55)
    a, b = center-eps, center+eps
    q, qp, o63 = rows[61], [k*rows[61][k] for k in range(1,len(rows[61]))], rows[63]
    qa = fmpq(a.numerator,a.denominator)
    qb = fmpq(b.numerator,b.denominator)
    va, vb = sturm_variations(q,qa), sturm_variations(q,qb)
    assert va-vb == 1, (va,vb)
    assert fmpz_poly(q)(qa) != 0 and fmpz_poly(q)(qb) != 0
    sq, qvalue, qerror = sign_near_center(qp,center,eps)
    sp, pvalue, perror = sign_near_center(o63,center,eps)
    assert sq*sp < 0, (sq,sp,qvalue,pvalue)
    print("j=2 n=63 q=O_61")
    print("root_count_in_bracket",va-vb)
    print("a",a)
    print("b",b)
    print("qprime_sign",sq,"value_bits",abs(qvalue.numerator).bit_length(),"error_bits",qerror.numerator.bit_length())
    print("O63_sign",sp,"value_bits",abs(pvalue.numerator).bit_length(),"error_bits",perror.numerator.bit_length())
    print("residue_sign",sq*sp,"=> WRS holds at this certified root")

if __name__ == "__main__":
    main()
