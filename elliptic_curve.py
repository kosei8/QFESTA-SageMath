# Sagemath import
from sage.all import (
    is_prime
)

# return a point on E of order l^e, where l is prime.
def point_ord(E, l, e):
    assert is_prime(l)
    p = E.base_ring().characteristic()
    assert E.order() == (p + 1)**2
    assert (p + 1) % l**e == 0
    P = (p + 1)//(l**e)*E.random_point()
    while (l**(e-1)*P).is_zero():
        P = (p + 1)//(l**e)*E.random_point()
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