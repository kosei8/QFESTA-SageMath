# SageMath imports
from sage.all import (
    ceil,
    log,
    Primes,
    is_prime
)

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
    while not is_prime(2**(2*a)*3*f - 1):
        f += 1
    return 2*a, 2*b, f, k, D1, D2
