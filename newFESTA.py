# Sagemath import
from sage.all import (
    randint,
    gcd,
    discrete_log,
    ZZ
)

# local import
import quaternion
import endomorphism as end
import elliptic_curve as ec
import parameter_generate as pgen
import d2isogeny

def key_gen(sys_param):
    Fp4, basis2, basis3, e2, e3, _, D1, D2, zeta2 = sys_param
    # secret key
    sec_key = 2*randint(0, 2**(e2-1)) + 1

    # public key
    p = basis2[0].base_ring().characteristic()
    E = basis2[0].curve()
    alpha = quaternion.FullRepresentInteger(3**e3*D1*D2, p)
    while gcd(gcd(alpha), 6) > 1:
        alpha = quaternion.FullRepresentInteger(3**e3*D1*D2, p)
    K = end.kernel(quaternion.involution(alpha), basis3, 3**e3, zeta2, Fp4)
    phi = E.isogeny(K, algorithm="factored")
    images = [ZZ(3**e3).inverse_mod(2**e2)*phi(end.action(alpha, P, zeta2, Fp4)) for P in basis2]
    images[0] = sec_key*images[0]
    images[1] = ZZ(sec_key).inverse_mod(2**e2)*images[1]
    assert images[0].weil_pairing(images[1], 2**e2) == basis2[0].weil_pairing(basis2[1], 2**e2)**(D1*D2)
    pub_key = (phi.codomain(), images)

    return sec_key, pub_key

def encrypt(message, sys_param, pub_key):
    _, basis2, basis3, e2, e3, e5, _, _, _ = sys_param
    E0 = basis2[0].curve()
    E1, images = pub_key

    beta = 2*message + 1
    beta_inv = ZZ(beta).inverse_mod(2**e2)

    # isogeny from E0 of degree 3^e2
    r, s = 0, 0
    while gcd(r, s) % 3 == 0:
        r, s = [randint(0, 3**e3-1) for _ in range(2)]
    K1 = r*basis3[0] + s*basis3[1]
    phi1 = E0.isogeny(K1, algorithm="factored")
    P1, Q1 = [phi1(P) for P in basis2]

    # isogeny from E1 of degree 5^e5
    u = (-1 + (E0.base_ring()(5).sqrt()))/2
    zeta5 = (u + (-u - 3).sqrt())/2
    assert zeta5**5 == 1
    P2, Q2 = images
    E, P2, Q2 = ec.chain_5radials(E1, P2, Q2, zeta5, e5)

    # for richelot computation
    iota = E.isomorphism_to(E.short_weierstrass_model())
    P2, Q2 = iota(P2), iota(Q2)

    return beta*P1, beta_inv*Q1, beta*P2, beta_inv*Q2

def decrypt(ciphertext, sys_param, sec_key):
    P1, Q1, P2, Q2 = ciphertext
    E1 = P1.curve()
    E2 = P2.curve()
    _, basis2, _, e2, e3, e5, _, D2, _ = sys_param

    P1d = 5**e5*D2*P1
    Q1d = 5**e5*D2*Q1
    P2d = 3**e3*ZZ(sec_key).inverse_mod(2**e2)*P2
    Q2d = 3**e3*sec_key*Q2

    assert P1d.weil_pairing(Q1d, 2**e2)*P2d.weil_pairing(Q2d, 2**e2) == 1
    R, S = ec.basis(E1, 3, e3)
    X, Y = d2isogeny.D2IsogenyImage(E1, E2, P1d, Q1d, P2d, Q2d, e2, R, S)
    Rd, Sd = X[0], Y[0]
    if not Rd.weil_pairing(Sd, 3**e3) == 1:
        Rd, Sd = X[1], Y[1]
    assert Rd.weil_pairing(Sd, 3**e3) == 1
    if ec.order(Rd, 3, e3) > ec.order(Sd, 3, e3):
        k = discrete_log(Sd, Rd, 3**e3, operation='+')
        K = k*R - S
    else:
        k = discrete_log(Rd, Sd, 3**e3, operation='+')
        K = R - k*S
    phi = E1.isogeny(K, algorithm="factored")
    assert phi.codomain().j_invariant() == basis2[0].curve().j_invariant()
    iota = phi.codomain().isomorphism_to(basis2[0].curve())
    phi = iota*phi
    R, S = [phi(P) for P in [P1, Q1]]
    k = ZZ(3**e3).inverse_mod(2**e2)
    m = ZZ(k*discrete_log(R, basis2[0], 2**e2, operation='+')) % 2**e2
    md = ZZ(k*discrete_log(S, basis2[1], 2**e2, operation='+')) % 2**e2
    assert (m*md) % 2**e2 == 1
    if m >= 2**(e2-1):
        m = 2**e2 - m
    return (m - 1)//2


