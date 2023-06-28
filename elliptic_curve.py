# Sagemath import
from sage.all import (
    is_prime,
    PolynomialRing,
    EllipticCurve,
    randint,
    ZZ
)

# return a point on E of order l^e, where l is prime.
def point_ord(E, l, e):
    assert is_prime(l)
    p = E.base_ring().characteristic()
    assert (p + 1) % l**e == 0
    P = (p + 1)//(l**e)*E.random_point()
    while (l**(e-1)*P).is_zero():
        P = (p + 1)//(l**e)*E.random_point()
    assert (l**e*P).is_zero()
    return P

# return a basis of E[l^e], where l is prime.
def basis(E, l, e):
    P = point_ord(E, l, e)
    Q = point_ord(E, l, e)
    while (l**(e-1)*P).weil_pairing(l**(e-1)*Q, l) == 1:
        Q = point_ord(E, l, e)
    return P, Q

def order(P, l, e):
    assert ((l**e)*P).is_zero()
    for k in range(e):
        if P.is_zero():
            return l**k
        P = l*P
    return l**e

# a curve of form y^2 + a1xy + a3y = x^3 isomorphic to E. The image of an order-3 point P is (0, 0).
def curve_for_3radical(E, P, Qs):
    F = E.base_ring()
    R = PolynomialRing(F, 2, names=["x", "y"])
    x, y = R.gens()
    f = E.defining_polynomial()(x=x+P[0], y=y+P[1], z=1)
    a = f.coefficient(x).constant_coefficient()/f.coefficient(y).constant_coefficient()
    g = f(y=y-a*x)
    a1 = g.coefficient(x).coefficient(y)
    a3 = g.coefficient(y).constant_coefficient()
    Ed = EllipticCurve(F, [a1, 0, a3, 0, 0])
    assert E.is_isomorphic(Ed)

    Qds = []
    for Q in Qs:
        Qds.append(Ed([Q[0] - P[0], Q[1] + a*(Q[0] - P[0]) - P[1]]))
    return Ed, Qds

# random 3-isogney by a radical isogeny
def radical_deg3(a1, a3, xs, zeta3):
    p = a1.base_ring().characteristic()
    N = (p**2 - 1)//3
    r = ZZ(3).inverse_mod(N)    # assume 3 || p - 1 
    a = (-a3)**r
    assert a**3 == -a3
    #a = (-a3).nth_root(3) * zeta3**randint(0, 2)
    a *= zeta3**randint(0, 2)
    a1d = -6*a + a1
    a3d = 3*a1*a**2 - a1**2*a + 9*a3

    retXs = []
    for x in xs:
        X = (x**3 + a1*a3*x + a3**2)/x**2
        X -= -a1*a + 3*a**2
        retXs.append(X)

    return a1d, a3d, retXs

# chain of 3-isogenies
def chain_3radials(E, P, Q, zeta3, n):
    if n <= 0:
        return E, P, Q
    F = E.base_ring()
    K = point_ord(E, 3, 1)
    E, points = curve_for_3radical(E, K, [P, Q])
    P, Q = points
    a1, a3 = E.a1(), E.a3()
    xs = [P[0], Q[0], (P + Q)[0]]
    for _ in range(n):
        a1, a3, xs = radical_deg3(a1, a3, xs, zeta3)
    E = EllipticCurve(F, [a1, 0, a3, 0, 0])
    P = E.lift_x(xs[0])
    Q = E.lift_x(xs[1])
    if not (P + Q)[0] == xs[2]:
        Q = -Q
    assert (P + Q)[0] == xs[2]
    return E, P, Q

# a curve of form y^2 + (1-b)xy - by = x^3 - bx^2 isomorphic to E. The image of an order-5 point P is (0, 0).
def curve_for_5radical(E, P, Qs):
    F = E.base_ring()
    R = PolynomialRing(F, 2, names=["x", "y"])
    x, y = R.gens()
    f = R(E.defining_polynomial()(x=x+P[0], y=y+P[1], z=1))
    a = f.coefficient(x).constant_coefficient()/f.coefficient(y).constant_coefficient()
    g = f(y=y-a*x)
    u = -g.coefficient(y).constant_coefficient()/g.coefficient(x**2)
    g = g(x=u**2*x, y=u**3*y)/u**6
    b = g.coefficient(x**2)
    Ed = EllipticCurve(F, [1 - b, -b, -b, 0, 0])
    assert E.is_isomorphic(Ed)

    Qds = []
    for Q in Qs:
        Qds.append(Ed([(Q[0] - P[0])/u**2, (Q[1] + a*(Q[0] - P[0]) - P[1])/u**3]))

    return Ed, Qds

# random 5-isogney by a radical isogeny
def radical_deg5(b, xs, zeta5):
    p = b.base_ring().characteristic()
    N = (p**2 - 1)//5
    r = ZZ(5).inverse_mod(N)    # assume 5 || p - 1 
    a = b**r
    assert a**5 == b
    a *= zeta5**randint(0, 4)
    bd = a*(a**4 + 3*a**3 + 4*a**2 + 2*a + 1)/(a**4 - 2*a**3 + 4*a**2 - 3*a + 1)

    retXs = []
    for x in xs:
        X = (b**3*x**3 + b**4*x - 3*b**3*x**2 + 3*b**2*x**3 - 2*b*x**4 + x**5 + b**4 - 3*b**3*x + 3*b**2*x**2 - b*x**3)/(b**2*x**2 - 2*b*x**3 + x**4)
        X -= 5*a**4 + (b - 3)*a**3 + (b + 2)*a**2 + (2*b - 1)*a - 2*b
        u = -((-5*b + 15)*a**4 + (-b**2 + 6*b - 9)*a**3 + (b**2 - 21*b + 4)*a**2 + (-4*b**2 + 29*b - 1)*a - 25*b)/((b - 3)*a**4 + (b + 2)*a**3 + (2*b - 1)*a**2 + (3*b + 1)*a + 5*b)
        X /= u**2
        retXs.append(X)

    return bd, retXs

# chain of 5-isogenies
def chain_5radials(E, P, Q, zeta5, n):
    if n <= 0:
        return E, P, Q
    F = E.base_ring()
    K = point_ord(E, 5, 1)
    E, points = curve_for_5radical(E, K, [P, Q])
    P, Q = points
    b = -E.a3()
    xs = [P[0], Q[0], (P + Q)[0]]
    for _ in range(n):
        b, xs = radical_deg5(b, xs, zeta5)
    E = EllipticCurve(F, [1-b, -b, -b, 0, 0])
    P = E.lift_x(xs[0])
    Q = E.lift_x(xs[1])
    if not (P + Q)[0] == xs[2]:
        Q = -Q
    assert (P + Q)[0] == xs[2]
    return E, P, Q
