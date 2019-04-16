"""
Polynomial multiplication and integration, for expanding and integrating
factorized polynomials to define interpolatory quadrature rules, particulary
for inner product quadrature.

Created Oct 10, 2013 by Tom Loredo
2019:  Converted to Python 3
"""

# TODO:  Cython version pending; this is currently not used.

from numpy import array, zeros


__all__ = ['poly_mul_2', 'factor_times_poly', 'roots2coefs', 'integrate_poly',
           'indef2def']


def poly_mul_2(a, b):
    """
    Find the polynomial coefficients for the product of two polynomials
    defined by arrays of coefficients.
    """
    pass

def factor_times_poly(u, b):
    """
    Multiply a polynomial factor (x - u) times a polynomial specified by an
    array of coefficients, returning an array of coefficients for the resulting
    polynomial.

    The 2nd polynomial is:  b[0] + x*b[1] + x**2 * b[2] + ...
    """
    coef = zeros(len(b)+1)
    coef[1:] = b
    coef[:-1] += -u*b
    return coef

def roots2coefs(uvals):
    """
    Given a list of roots of a polynomial, return an array of the polynomial
    coefficients for the product of factors (x - u[i]).
    """
    # TODO:  This creates lots of intermediate arrays; rewrite to
    # use scratch space.
    coef = array([-uvals[0], 1.])  #  -u_0 + x
    if len(uvals) == 1:
        return coef
    for u in uvals[1:]:
        coef = factor_times_poly(u, coef)
    return coef

def integrate_poly(coef, a, b):
    """
    Return the definite integral of the polynomial defined by the array of
    coefficients `coef` over the interval [a, b].
    """
    # Array of coefficients of the indefinite integral, dropping the
    # arbitrary constant term (it would have been icoef[0]).
    icoef = array([c/(i+1.) for i, c in enumerate(coef)])
    print(icoef)
    return indef2def(icoef, a, b)

    # This codes indef2def in place:
    # # Evaluate a, b contributions using Horner's algorithm.
    # aval = icoef[-1]
    # for c in icoef[-2::-1]:
    #     aval  = c + a*aval
    # aval *= a  # since const term dropped
    # bval = icoef[-1]
    # for c in icoef[-2::-1]:
    #     bval  = c + b*bval
    # bval *= b  # since const term dropped
    # return bval - aval

def indef2def(indefs, a, b):
    """
    Use the stored polynomial indefinite integral coefficients to calculate
    the definite integral of the polynomial over the interval [a,b].
    """
    # Evaluate a, b contributions using Horner's algorithm.
    aval = indefs[-1]
    for c in indefs[-2::-1]:
        aval  = c + a*aval
    aval *= a  # since const term dropped
    bval = indefs[-1]
    for c in indefs[-2::-1]:
        bval  = c + b*bval
    bval *= b  # since const term dropped
    return bval - aval
