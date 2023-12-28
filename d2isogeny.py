from sage.all import (
    EllipticCurve,
    PolynomialRing
)

import elliptic_curve as ec
import richelot_isogenies as richelot
import utilities_festa as utilities
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

# the images of pairs of points in Rs under a (2^e, 2^e)-isogeny from E1 time E2 with kernel <(P1, P2), (Q1, Q2)>
def D2IsogenyImage(E1, E2, P1, Q1, P2, Q2, e, Rs, strategy, use_theta=False):
    if E1.j_invariant() == E2.j_invariant() == 1728:
        # product of elliptic curves with j-invariant 1728
        cnt = 0
        h = None
        while h == None:
            ret = FromSpecialProdToProd(P1, Q1, P2, Q2, e)
            if ret:
                h, P1, P2, Q1, Q2, phi = ret
                e -= 1
                Rs = [phi(R) for R in Rs]
                cnt += 1
            else:
                break
        assert cnt <= 2

        # product of two elliptic curves
        ret = True
        while not ret == None:
            ret = FromProdToProd(P1, Q1, P2, Q2, e)
            if ret:
                h, P1, P2, Q1, Q2, phi = ret
                e -= 1
                Rs = [phi(R) for R in Rs]

        # Here, we can apply ProdToJac
        # transform to Montgomery curves. ProdToJac requires ((0,0), (0,0)) in the kernel.
        E1, PQ1Rs = ec.WeierstrassToMontgomery(P1.curve(), (2**(e-2)*P1).xy()[0], [P1, Q1] + [R[0] for R in Rs])
        E2, PQ2Rs = ec.WeierstrassToMontgomery(P2.curve(), (2**(e-2)*P2).xy()[0], [P2, Q2] + [R[1] for R in Rs])

        if use_theta:
            kernel = [CouplePoint(PQ1Rs[i], PQ2Rs[i]) for i in range(2)]
            Rs = [CouplePoint(PQ1Rs[i], PQ2Rs[i]) for i in range(2,len(PQ1Rs))]
            phi = EllipticProductIsogenySqrt(kernel, e)
            Rs = [phi(R) for R in Rs]
        else:
            Rs = [CouplePoint(PQ1Rs[i], PQ2Rs[i]) for i in range(2,len(PQ1Rs))]
            if e - 1 in strategy:
                st = strategy[e-1]
            else:
                st = utilities.optimised_strategy(e-1)
            chain = richelot.compute_richelot_chain(PQ1Rs[:2] + PQ2Rs[:2], e, st)
            for phi in chain:
                Rs = [phi((R[0], R[1])) for R in Rs]
        return Rs
    else:
        if use_theta:
            kernel = [CouplePoint(P1, P2), CouplePoint(Q1, Q2)]
            phi = EllipticProductIsogenySqrt(kernel, e)
            Rs = [CouplePoint(R[0], R[1]) for R in Rs]
            Rs = [phi(R) for R in Rs]
        else:
            if e - 1 in strategy:
                st = strategy[e-1]
            else:
                st = utilities.optimised_strategy(e-1)
            chain = richelot.compute_richelot_chain([P1, Q1, P2, Q2], e, st)
            for phi in chain:
                Rs = [phi((R[0], R[1])) for R in Rs]
        return Rs

# (2,2)-isogeny from E1 times E2 with kernel 2**(e-1)*<(P1, P2), (Q1, Q2)>
def FromProdToProd(P1, Q1, P2, Q2, e):
    E1 = P1.curve()
    E2 = P2.curve()
    assert Q1.curve() == E1 and Q2.curve() == E2

    T1 = 2**(e-1)*P1
    T2 = 2**(e-1)*P2
    S1 = 2**(e-1)*Q1
    S2 = 2**(e-1)*Q2

    # the product of 2-isogenies from E
    if T1.is_zero() or T2.is_zero() or S1.is_zero() or S2.is_zero():
        if T1.is_zero():
            T1, S1 = S1, T1
            T2, S2 = S2, T2
        assert not T1.is_zero()
        if not S1.is_zero():
            assert S1 == T1
            S2 -= T2
        assert T2.is_zero() or T2 == S2
        phi1 = E1.isogeny(T1)
        phi2 = E2.isogeny(S2)
        def isogeny(Rs):
            R1, R2 = Rs
            return phi1(R1), phi2(R2)
        imP = isogeny([P1, P2])
        imQ = isogeny([Q1, Q2])
        return [phi1.codomain(), phi2.codomain()], imP[0], imP[1], imQ[0], imQ[1], isogeny
    else:
        return None

# (2,2)-isogeny from the elliptic curve defined by y^2 = x^3 + x
def FromSpecialProdToProd(P1, Q1, P2, Q2, e):
    E = P1.curve()
    assert E == P2.curve() == Q1.curve() == Q2.curve()
    assert E == EllipticCurve(E.base_ring(), [1, 0])

    F = E.base_ring()
    Rx = PolynomialRing(F, name="x")
    x = Rx.gens()[0]

    T1 = 2**(e-1)*P1
    T2 = 2**(e-1)*P2
    S1 = 2**(e-1)*Q1
    S2 = 2**(e-1)*Q2

    # the product of 2-isogenies from E
    if T1.is_zero() or T2.is_zero() or S1.is_zero() or S2.is_zero():
        h, _, _, _, _, phi = FromProdToProd(T1, S1, T2, S2, 1)
        imP = phi([P1, P2])
        imQ = phi([Q1, Q2])
        return h, imP[0], imP[1], imQ[0], imQ[1], phi

    # change T1 to (0, 0) in E
    if S1[0] == 0:
        T1, S1 = S1, T1
        T2, S2 = S2, T2
    elif not T1[0] == 0:
        T1 += S1
        T2 += S2
    assert T1[0] == 0

    if T2[0] == 0:
        if S1 == S2:
            # kernel is {(R, R) | R in E[2]}
            def isogeny(Rs):
                R1, R2 = Rs
                return R1 - R2, R1 + R2
        else:
            # kernel is <((0,0),(0,0)), ((zeta2,0),(-zeta2,0))>, where zeta2^2 = -1
            def isogeny(Rs):
                R1, R2 = Rs
                zeta2 = S1[0]
                assert zeta2**2 == -1
                def iota(R):
                    return E([-R[0], zeta2*R[1], R[2]])
                return R1 + iota(R2), iota(R1) + R2

        imP = isogeny([P1, P2])
        imQ = isogeny([Q1, Q2])
        return None, imP[0], imP[1], imQ[0], imQ[1], isogeny
    else:
        return None
