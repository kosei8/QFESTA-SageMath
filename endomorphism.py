# SageMath imports
from sage.all import (
    EllipticCurve,
    PolynomialRing,
    discrete_log,
    ZZ,
    factor,
    matrix,
    identity_matrix,
    vector
)

import elliptic_curve as ec
import utilities_festa as utilities

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

# matrix of the multiplication by alpha w.r.t. basis of N-torsion subgroup
def action_matrix(alpha, basis, N, zeta2, F2):
    P, Q = basis
    aP, aQ = [action(alpha, R, zeta2, F2) for R in [P, Q]]
    a, c = utilities.BiDLP(aP, P, Q, N)
    b, d = utilities.BiDLP(aQ, P, Q, N)
    assert aP == a*P + c*Q
    return matrix([[a, b], [c, d]])

# matrices of the multiplications by i, (i + j)/2, (1 + ij)/2
def action_matrices(basis, N, zeta2, F2):
    one = [1,0,0,0]
    Ms = []
    for alpha in [one[-i:]+one[:-i] for i in range(1, 4)]:
        Ms.append(action_matrix(alpha, basis, N, zeta2, F2))
    return Ms

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2 w.r.t. vector representation by fixed basis
def action_by_matrices(alpha, v, action_matrices):
    Ms = [identity_matrix(2)] + action_matrices
    M = sum(alpha[i]*Ms[i] for i in range(4))
    return M * vector(v)

# return a generator of E[alpha, ord], where E: y^2 = x^3 + x, P, Q is a basis of E[ord].
def kernel(alpha, basis, ord, action_matrices):
    assert ZZ(ord).is_prime_power()
    l, e = factor(ord)[0]
    P, Q = basis

    R1 = action(alpha, [1, 0], action_matrices)
    R2 = action(alpha, [0, 1], action_matrices)

    if ec.order(R1, l, e) > ec.order(R2, l, e):
        a = discrete_log(R2, R1, ord, operation='+')
        return a*P - Q
    else:
        a = discrete_log(R1, R2, ord, operation='+')
        return P - a*Q

