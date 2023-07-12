from sage.all import (
    EllipticCurve,
    HyperellipticCurve,
    PolynomialRing,
    Matrix
)

import richelot_aux as richelot

# the images of (R1, O) and (S1, O) under a (2^e, 2^e)-isogeny from E1 time E2 with kernel <(P1, P2), (Q1, Q2)>
def D2IsogenyImage(E1, E2, P1, Q1, P2, Q2, e, R1, S1):
    O = E2([0, 1, 0])
    R = (R1, O)
    S = (S1, O)
    if E1.j_invariant() == E2.j_invariant() == 1728:
        chain = []
        h = None

        # product of elliptic curves with j-invariant 1728
        cnt = 0
        while h == None:
            h, P1, P2, Q1, Q2, phi = FromSpecialProdToDim2(P1, Q1, P2, Q2, e)
            e -= 1
            R, S = phi(R), phi(S)
            cnt += 1
        assert cnt <= 2

        # product of two elliptic curves
        while type(h) == list:
            h, P1, P2, Q1, Q2, phi = FromProdToDim2(P1, Q1, P2, Q2, e)
            e -= 1
            R, S = phi(R), phi(S)

        # Jacobian to product
        D11 = P1
        D12 = P2
        D21 = Q1
        D22 = Q2
        next_powers = None
        while e > 1:
            h, D11, D12, D21, D22, phi, next_powers = richelot.FromJacToJac(h, D11, D12, D21, D22, e, powers=next_powers)
            e -= 1
            R, S = phi(R), phi(S)

        G1 = D11
        G2 = D21
        G3, r3 = h.quo_rem(G1 * G2)
        assert r3 == 0

        delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
        assert delta.determinant() == 0
        phi, _ = richelot.FromJacToProd(G1, G2, G3)
        return phi(R), phi(S)
    else:
        P12 = 2**(e-1)*P1
        P22 = 2**(e-1)*P2
        Q12 = 2**(e-1)*Q1
        Q22 = 2**(e-1)*Q2
        R1 = P12 + Q12
        R2 = P22 + Q22
        if P12[0] == P22[0] == 0 or Q12[0] == Q22[0] == 0 or R1[0] == R2[0] == 0:
            print("hoge")
        chain, _ = richelot.Does22ChainSplit(E1, E2, P1, Q1, P2, Q2, e)
        for phi in chain:
            R, S = phi(R), phi(S)
        return R, S

# (2,2)-isogeny from E1 times E2 with kernel 2**(e-1)*<(P1, P2), (Q1, Q2)>
def FromProdToDim2(P1, Q1, P2, Q2, e):
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
        return richelot.FromProdToJac(E1, E2, P1, Q1, P2, Q2, e)

# (2,2)-isogeny from the elliptic curve defined by y^2 = x^3 + x
def FromSpecialProdToDim2(P1, Q1, P2, Q2, e):
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
        h, _, _, _, _, phi = FromProdToDim2(T1, S1, T2, S2, 1)
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
        # change S2 to (0, 0) in E
        if not S2[0] == 0:
            S1 += T1
            S2 += T2

        zeta2 = T2[0]
        assert zeta2**2 == -1
        sign = F(-1)**(not S1[0] == zeta2)
        s1 = 2 * zeta2 / 3
        s1inv = -3 * zeta2 / 2
        t1 = 2 * sign * zeta2 / 3
        t1inv = -3 * sign * zeta2 / 2
        s2 = -zeta2 * sign / 3
        t2 = s2*sign
        a1_t =  -s2 * s1inv
        a2_t = (sign*zeta2 - s2) * s1inv
        a3_t = (-sign*zeta2 - s2) * s1inv
        h = s1 * (x**2 - a1_t) * (x**2 - a2_t) * (x**2 - a3_t)

        def isogeny(Rs):
            J = HyperellipticCurve(h).jacobian()
            R1, R2 = Rs

            if not R1.is_zero():
                xR1, yR1 = R1.xy()
                JR1 = J([s1 * x**2 + s2 - xR1, Rx(yR1 * s1inv)])
            if not R2.is_zero():
                xR2, yR2 = R2.xy()
                JR2 = J([(xR2 - t2) * x**2 - t1, yR2 * x**3 * t1inv])
            if R1.is_zero():
                return JR2
            elif R2.is_zero():
                return JR1
            else:
                return JR1 + JR2

        imP = isogeny((P1, P2))
        imQ = isogeny((Q1, Q2))
        return h, imP[0], imP[1], imQ[0], imQ[1], isogeny
