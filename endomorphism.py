# SageMath imports
from sage.all import (
    EllipticCurve,
    PolynomialRing,
    discrete_log,
    ZZ,
    factor
)

import quaternion
import elliptic_curve as ec

# action of square root of -1
def i_action(P, zeta2):
    F = P.base_ring()
    E = P.curve()
    assert zeta2 in F
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    X, Y, Z = P
    return E([-X, zeta2*Y, Z])

# Frobenius endomorphism
def Frobenius(P):
    p = P.base_ring().characteristic()
    E = P.curve()
    X, Y, Z = P
    return E([X**p, Y**p, Z**p])

# return retP s.t. 2*retP = P. Note that retP is over an extension field.
def half_point(P, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    assert F.is_subring(F2)

    E2 = EllipticCurve(F2, [1, 0])
    R = PolynomialRing(F2, name="X")
    X = R.gens()[0]
    if P.is_zero():
        f = X**3 + X
    else:
        f = P[2]*(X**2 - 1)**2 - P[0]*4*(X**3 + X)
    xs = f.roots(multiplicities=False)
    assert len(xs) > 0
    x = xs[0]
    y = (x**3 + x).sqrt()
    retP = E2([x, y])
    if not E(2*retP) == P:
        retP = -retP
    assert E(2*retP) == P
    return retP

# the action of (i + j)/2
def i_j_2_action(P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(i_action(halfP, zeta2) + Frobenius(halfP))

# the action of (1 + ij)/2
def one_ij_2_action(P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(halfP + i_action(Frobenius(halfP), zeta2))

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2
def action(alpha, P, zeta2, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    a, b, c, d = alpha
    ret = a*P
    ret += b*i_action(P, zeta2)
    ret += c*i_j_2_action(P, zeta2, F2)
    ret += d*one_ij_2_action(P, zeta2, F2)
    return ret

# return a generator of E[alpha, ord], where E: y^2 = x^3 + x, P, Q is a basis of E[ord].
def kernel(alpha, basis, ord, zeta2, Fp4):
    assert ZZ(ord).is_prime_power()
    l, e = factor(ord)[0]
    P, Q = basis
    E = P.curve()
    assert Q.curve() == E
    F = E.base_ring()
    p = F.characteristic()
    N = quaternion.norm(alpha, p)
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    assert N % ord == 0
    assert (ord*P).is_zero()
    assert (ord*Q).is_zero()

    R1 = action(alpha, P, zeta2, Fp4)
    R2 = action(alpha, Q, zeta2, Fp4)
    if ec.order(R1, l, e) > ec.order(R2, l, e):
        a = discrete_log(R2, R1, ord, operation='+')
        assert action(alpha, a*P - Q, zeta2, Fp4).is_zero()
        return a*P - Q
    else:
        a = discrete_log(R1, R2, ord, operation='+')
        assert action(alpha, P - a*Q, zeta2, Fp4).is_zero()
        return P - a*Q

