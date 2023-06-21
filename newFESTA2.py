from sage.all import (
    randint,
    ZZ,
    discrete_log
)

import quaternion as quat
import endomorphism as End
import elliptic_curve as ec
import richelot_aux as richelot

# the image of P, Q under a random isogeny of degree N
def NonSmoothRandomIsog(zeta2, Fp4, e, N, P, Q):
    Fp2 = zeta2.base_ring()
    p = Fp2.characteristic()
    E0 = P.curve()
    assert N % 2 == 1
    assert 2**e > N
    assert (p + 1) % 2**e == 0
    assert Fp2.is_subring(Fp4)
    assert ((2**e)*P).is_zero() and ((2**e)*Q).is_zero()
    assert ((2**(e-1))*P).weil_pairing((2**(e-1))*Q, 2) == -1

    # require (0, 0) neq alpha((0, 0)) so that the first (2, 2)-isogeny does not split
    alpha = [1, 0, 0, 0]
    while sum(alpha) % 2 == 1:
        alpha = quat.FullRepresentInteger(N*(2**e - N), p)
        alphaP = End.action(alpha, P, zeta2, Fp4)
        alphaQ = End.action(alpha, Q, zeta2, Fp4)
 
    assert P.weil_pairing(Q, 2**e)**((2**e - N)**2)*alphaP.weil_pairing(alphaQ, 2**e) == 1
    chain, _ = richelot.Does22ChainSplit(E0, E0, (2**e - N)*P, (2**e - N)*Q, alphaP, alphaQ, e)
    O = E0([0,1,0])
    X, Y = (P, O), (Q, O)
    for phi in chain:
        X, Y = phi(X), phi(Y)
    Pd, Qd = X[0], Y[0]
    if not Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N:
        Pd, Qd = X[1], Y[1]
    assert Pd.weil_pairing(Qd, 2**e) == P.weil_pairing(Q, 2**e)**N

    return Pd, Qd

def key_gen(sys_param):
    Fp4, basis2, a, _, _, D1, _, zeta2 = sys_param
    # secret key
    sec_key = 2*randint(0, 2**(a-1)) + 1

    #public key
    P, Q = basis2
    Pd, Qd = NonSmoothRandomIsog(zeta2, Fp4, a, D1, P, Q)
    Pd = sec_key*Pd
    Qd = ZZ(sec_key).inverse_mod(2**a)*Qd
    pub_key = Pd.curve(), [Pd, Qd]

    return sec_key, pub_key

def encrypt(message, sys_param, pub_key):
    Fp4, basis2, a, b, _, _, D2, zeta2 = sys_param
    E1, images = pub_key

    beta = 2*message + 1
    beta_inv = ZZ(beta).inverse_mod(2**a)

    # isogeny from E0 of degree D2
    P, Q = basis2
    P1, Q1 = NonSmoothRandomIsog(zeta2, Fp4, a, D2, P, Q)

    # isogeny from E1 of degree 3^b
    E = E1
    P2, Q2 = images
    for _ in range(b):
        K = ec.point_ord(E, 3, 1)
        phi = E.isogeny(K)
        E = phi.codomain()
        P2, Q2 = phi(P2), phi(Q2)

    return beta*P1, beta_inv*Q1, beta*P2, beta_inv*Q2

def decrypt(ciphertext, sys_param, sec_key, pub_key):
    P1, Q1, P2, Q2 = ciphertext
    E1 = P1.curve()
    E2 = P2.curve()
    _, _, a, b, k, _, D2, _ = sys_param

    EA, tmp = pub_key
    PA, QA = tmp
    PA = D2*ZZ(sec_key).inverse_mod(2**a)*PA
    QA = D2*sec_key*QA

    P1d = 3**b*k*P1
    Q1d = 3**b*k*Q1
    P2d = D2*ZZ(sec_key).inverse_mod(2**a)*P2
    Q2d = D2*sec_key*Q2

    assert P1d.weil_pairing(Q1d, 2**a)*P2d.weil_pairing(Q2d, 2**a) == 1
    chain, _ = richelot.Does22ChainSplit(E1, E2, P1d, Q1d, P2d, Q2d, a)
    O = E2([0,1,0])
    X, Y = (P1, O), (Q1, O)
    for phi in chain:
        X, Y = phi(X), phi(Y)
    R, S = X[0], Y[0]
    if not R.curve().is_isomorphic(EA):
        R, S = X[1], Y[1]
    assert R.curve().is_isomorphic(EA)
    iota = R.curve().isomorphism_to(EA)
    R, S = iota(R), iota(S)
    m = ZZ(discrete_log(R, PA, 2**a, operation='+'))
    md = ZZ(discrete_log(S, QA, 2**a, operation='+'))
    assert (m*md) % 2**a == 1
    if m >= 2**(a-1):
        m = 2**a - m
    return (m - 1)//2



