# SageMath imports
from sage.all import (
    log,
    xgcd,
    gcd,
    floor,
    ceil,
    GF,
    is_prime,
    ZZ,
    PolynomialRing,
    Primes
)

# prime of form 2^e2*3^e3*5*f - 1
def calc_p(lam, modified_coeff):
    if modified_coeff:
        e2 = ceil(modified_coeff*lam)
    else:
        e2 = 5*lam
    e3 = ceil(2*lam*log(2, 3))
    e5 = ceil(2*lam*log(2, 5)) 
    f = 1
    p = 2**e2*3**e3*5*f - 1
    while not is_prime(p):
        f += 1
        if gcd(f, 2*3*5) == 1:
            p = 2**e2*3**e3*5*f - 1
    return p, e2, e3, e5, f

# return D1, D2 > 0 s.t. D1d1 + D2d2 = n
def calcD1D2(d1, d2, n):
    g, a, b = xgcd(d1, d2)
    assert g == 1
    c = floor(n*(b - a)/(d1 + d2))
    D1 = n*a + c*d2
    D2 = n*b - c*d1
    while gcd(D1, 2*a) > 1:
        D1 += d2
        D2 -= d1
    assert D1*d1 + D2*d2 == n and D1 > 0 and D2 > 0
    return D1, D2

# return GF(p^4), GF(p^2), square root of -1 in GF(p^2)
def calcFields(p):
    assert p % 4 == 3
    R = PolynomialRing(GF(p), name="x")
    x = R.gens()[0]
    Fp2 = GF(p**2, modulus=x**2+1, name="i")
    z = Fp2.random_element()
    while z.is_square():
        z = Fp2.random_element()
    t = ZZ(z + z**p)
    n = ZZ(z**(p+1))
    R = PolynomialRing(ZZ, name="x")
    x = R.gens()[0]
    Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
    Fp2 = Fp4.subfield(2)
    i = Fp2(-1).sqrt()
    return Fp4, Fp2, i

# generate system parameter
def params(lam, modified_coeff=None):
    p, e2, e3, e5, f = calc_p(lam, modified_coeff)
    D1, D2 = calcD1D2(3**e3, 5**e5, 2**e2)
    Fp4, Fp2, i = calcFields(p)
    return p, e2, e3, e5, f, D1, D2, Fp4, Fp2, i

# whether n is (small factors) * (prime >= 2^lam)
def IsNearPrime(n, lam):
    for p in Primes()[:100]:
        while n % p == 0:
            n //= p
    return n > 2**lam and is_prime(n)

# generate system parameter
def SysParam(lam):
    a = lam
    b = ceil(log(2, 3)*lam)
    k = 1
    D1 = 2**a - k*3**b
    D2 = 2**a + k*3**b
    while not(IsNearPrime(D1, lam) and IsNearPrime(D2, lam)):
        if 2**a > 3**b*(k+2):
            k += 2
        else:
            a += 1
            k = 1
        D1 = 2**a - k*3**b
        D2 = 2**a + k*3**b
    
    f = 1
    while (not is_prime(2**(2*a)*3*f - 1)) or f % 3 == 0:
        # f should be odd for using the Tate pairing
        f += 2
    return 2*a, 2*b, f, k, D1, D2

def BoundedFactor(n, B):
    d = 1
    for l in Primes()[:10000]:
        while n % l == 0:
            n //= l
            d *= l
            if d > B:
                return d
    return None

def SysParam2(lam):
    a = lam
    b1 = ceil(log(2, 3)*lam)
    b2 = 2*b1
    d1 = 2**a - 3**b1
    d2 = 2**(2*a) + 2**a * 3**b1 + 3**b2
    while d1 < 2**lam or d2 < 2**(2*lam):
        a += 1
        d1 = 2**a - 3**b1
        d2 = 2**(2*a) + 2**a * 3**b1 + 3**b2
    f = 1
    while (not is_prime(2**(3*a)*3*f - 1)) or f % 3 == 0:
        # f should be odd for using the Tate pairing
        f += 2
    assert d1*d2 + 3**(b1 + b2) == 2**(3*a)
    return 3*a, b1, b2, f, d1, d2
