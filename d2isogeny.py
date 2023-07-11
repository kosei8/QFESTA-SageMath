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
        chain, _ = richelot.Does22ChainSplit(E1, E2, P1, Q1, P2, Q2, e)
        return chain

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
        def isogeny(Ps):
            P1, P2 = Ps
            return phi1(P1), phi2(P2)
        imP = isogeny([P1, P2])
        imQ = isogeny([Q1, Q2])
        return [phi1.codomain(), phi2.codomain()], imP[0], imP[1], imQ[0], imQ[1], isogeny
    else:
        return richelot.FromProdToJac(E1, E2, P1, Q1, P2, Q2, e)

# (2,2)-isogeny from the elliptic curve defined by y^2 = x^3 + x
def FromSpecialProdToDim2(T0, T1, S0, S1, e):
    E = T0.curve()
    assert E == T1.curve() == S0.curve() == S1.curve()
    assert E == EllipticCurve(E.base_ring(), [1, 0])

    F = E.base_ring()
    Rx = PolynomialRing(F, name="x")
    x = Rx.gens()[0]

    T20 = 2**(e-1)*T0
    T21 = 2**(e-1)*T1
    S20 = 2**(e-1)*S0
    S21 = 2**(e-1)*S1

    # the product of 2-isogenies from E
    if T20.is_zero() or T21.is_zero() or S20.is_zero() or S21.is_zero():
        if T20.is_zero() or T21.is_zero():
            if T20.is_zero():
                T20, T21 = T21, T20
                S20, S21 = S21, S20
            assert not T20.is_zero()
            phi1 = E.isogeny(T20)
            phi2 = E.isogeny(S21)
            assert S20 == S21 or S20.is_zero()
            def isogeny(Ps):
                P1, P2 = Ps
                return phi1(P1), phi2(P2)
        else:
            if S20.is_zero():
                T20, T21 = T21, T20
                S20, S21 = S21, S20
            assert not S20.is_zero()
            phi1 = E.isogeny(T21)
            phi2 = E.isogeny(S20)
            assert T20 == T21 or T20.is_zero()
            def isogeny(Ps):
                P1, P2 = Ps
                return phi1(P1), phi2(P2)
        TS0 = isogeny([T0, S0])
        TS1 = isogeny([T1, S1])
        return [phi1.codomain(), phi2.codomain()], TS0[0], TS0[1], TS1[0], TS1[1], isogeny

    # change T20 to (0, 0) in E
    if T21[0] == 0:
        T20, T21 = T21, T20
        S20, S21 = S21, S20
    elif not T20[0] == 0:
        T20 = T20 + T21
        S20 = S20 + S21

    if S20[0] == 0:
        if T21[0] == S21[0]:
            def isogeny(Ps):
                P1, P2 = Ps
                return P1 - P2, P1 + P2
        else:
            def isogeny(Ps):
                P1, P2 = Ps
                zeta2 = T21[0]
                def iota(P):
                    return E([-P[0], zeta2*P[1], P[2]])
                return P1 + iota(P2), iota(P1) + P2

        TS0 = isogeny([T0, S0])
        TS1 = isogeny([T1, S1])
        return None, TS0[0], TS0[1], TS1[0], TS1[1], isogeny
    else:
        # change S21 to (0, 0) in E
        if not S21[0] == 0:
            T21 = T21 + T20
            S21 = S21 + S20

        zeta2 = S20[0]
        sign = F(-1)**(not T21[0] == zeta2)
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

        def isogeny(Ps):
            J = HyperellipticCurve(h).jacobian()
            P1, P2 = Ps

            if not P1.is_zero():
                xP1, yP1 = P1.xy()
                JP1 = J([s1 * x**2 + s2 - xP1, Rx(yP1 * s1inv)])
            if not P2.is_zero():
                xP2, yP2 = P2.xy()
                JP2 = J([(xP2 - t2) * x**2 - t1, yP2 * x**3 * t1inv])
            if P1.is_zero():
                return JP2
            elif P2.is_zero():
                return JP1
            else:
                return JP1 + JP2

        im0 = isogeny((T0, S0))
        im1 = isogeny((T1, S1))
        return h, im0[0], im0[1], im1[0], im1[1], isogeny
