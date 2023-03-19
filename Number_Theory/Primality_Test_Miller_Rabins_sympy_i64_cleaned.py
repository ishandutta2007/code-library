# 4x Faster than Hardcoded
import time
import random
from collections import defaultdict
from functools import reduce
import random
import math
from bisect import bisect


def _test(n, base, s, t):
    """Miller-Rabin strong pseudoprime test for one base.
    Return False if n is definitely composite, True if n is
    probably prime, with a probability greater than 3/4.

    """
    b = pow(base, t, n)
    if b == 1 or b == n - 1:
        return True
    else:
        for j in range(1, s):
            b = pow(b, 2, n)
            if b == n - 1:
                return True
            if b == 1:
                return False
    return False


powers = [1 << _ for _ in range(300)]


def python_bitcount(n):
    """Calculate bit size of the nonnegative integer n."""
    bc = bisect(powers, n)
    if bc != 300:
        return bc
    bc = int(math.log(n, 2)) - 4
    return bc + bctable[n >> bc]


mpmath_bitcount = python_bitcount


def bitcount(n):
    """Return smallest integer, b, such that |n|/2**b < 1."""
    return mpmath_bitcount(abs(int(n)))


SYMPY_INTS = (int,)


def as_int(n):
    return int(n)


small_trailing = [0] * 256
for j in range(1, 8):
    small_trailing[1 << j :: 1 << (j + 1)] = [j] * (1 << (7 - j))


def trailing(n):
    """Count the number of trailing zero digits in the binary
    representation of n, i.e. determine the largest power of 2
    that divides n.

    Examples
    ========

    >>> from sympy import trailing
    >>> trailing(128)
    7
    >>> trailing(63)
    0
    """
    n = abs(int(n))
    if not n:
        return 0
    low_byte = n & 0xFF
    if low_byte:
        return small_trailing[low_byte]

    z = bitcount(n) - 1
    if isinstance(z, SYMPY_INTS):
        if n == 1 << z:
            return z

    if z < 300:
        t = 8
        n >>= 8
        while not n & 0xFF:
            n >>= 8
            t += 8
        return t + small_trailing[n & 0xFF]
    t = 0
    p = 8
    while not n & 1:
        while not n & ((1 << p) - 1):
            n >>= p
            t += p
            p *= 2
        p //= 2
    return t


def mr(n, bases):
    """Perform a Miller-Rabin strong pseudoprime test on n using a
    given list of bases/witnesses.

    References
    ==========
    .. [1] Richard Crandall & Carl Pomerance (2005), "Prime Numbers:
           A Computational Perspective", Springer, 2nd edition, 135-138

    A list of thresholds and the bases they require are here:
    https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants

    Examples
    ========

    >>> from sympy.ntheory.primetest import mr
    >>> mr(1373651, [2, 3])
    False
    >>> mr(479001599, [31, 73])
    True
    """
    n = as_int(n)
    if n < 2:
        return False
    s = trailing(n - 1)
    t = n >> s
    for base in bases:
        if base >= n:
            base %= n
        if base >= 2:
            if not _test(n, base, s, t):
                return False
    return True


def isprime(n):
    """
    Test if n is a prime number (True) or not (False). For n < 2^64 the
    answer is definitive; larger n values have a small probability of actually
    being pseudoprimes.

    Negative numbers (e.g. -2) are not considered prime.

    The first step is looking for trivial factors, which if found enables
    a quick return.  Next, if the sieve is large enough, use bisection search
    on the sieve.  For small numbers, a set of deterministic Miller-Rabin
    tests are performed with bases that are known to have no counterexamples
    in their range.  Finally if the number is larger than 2^64, a strong
    BPSW test is performed.  While this is a probable prime test and we
    believe counterexamples exist, there are no known counterexamples.

    Examples
    ========

    >>> from sympy.ntheory import isprime
    >>> isprime(13)
    True
    >>> isprime(13.0)
    False
    >>> isprime(15)
    False

    Notes
    =====

    This routine is intended only for integer input, not numerical
    expressions which may represent numbers. Floats are also
    rejected as input because they represent numbers of limited
    precision. While it is tempting to permit 7.0 to represent an
    integer there are errors that may "pass silently" if this is
    allowed:

    >>> from sympy import Float, S
    >>> int(1e3) == 1e3 == 10**3
    True
    >>> int(1e23) == 1e23
    True
    >>> int(1e23) == 10**23
    False

    >>> near_int = 1 + S(1)/10**19
    >>> near_int == int(near_int)
    False
    >>> n = Float(near_int, 10)
    >>> n == int(n)
    True
    >>> n = Float(near_int, 20)
    >>> n == int(n)
    False

    See Also
    ========

    sympy.ntheory.generate.primerange : Generates all primes in a given range
    sympy.ntheory.generate.primepi : Return the number of primes less than or equal to n
    sympy.ntheory.generate.prime : Return the nth prime

    References
    ==========
    - https://en.wikipedia.org/wiki/Strong_pseudoprime
    - "Lucas Pseudoprimes", Baillie and Wagstaff, 1980.
      http://mpqs.free.fr/LucasPseudoprimes.pdf
    - https://en.wikipedia.org/wiki/Baillie-PSW_primality_test
    """
    try:
        n = as_int(n)
    except ValueError:
        return False
    if n in [2, 3, 5]:
        return True
    if n < 2 or (n % 2) == 0 or (n % 3) == 0 or (n % 5) == 0:
        return False
    if n < 49:
        return True
    if (
        (n % 7) == 0
        or (n % 11) == 0
        or (n % 13) == 0
        or (n % 17) == 0
        or (n % 19) == 0
        or (n % 23) == 0
        or (n % 29) == 0
        or (n % 31) == 0
        or (n % 37) == 0
        or (n % 41) == 0
        or (n % 43) == 0
        or (n % 47) == 0
    ):
        return False
    if n < 2809:
        return True
    if n < 31417:
        return pow(2, n, n) == 2 and n not in [7957, 8321, 13747, 18721, 19951, 23377]
    if n < 341531:
        return mr(n, [9345883071009581737])
    if n < 885594169:
        return mr(n, [725270293939359937, 3569819667048198375])
    if n < 350269456337:
        return mr(n, [4230279247111683200, 14694767155120705706, 16641139526367750375])
    if n < 55245642489451:
        return mr(n, [2, 141889084524735, 1199124725622454117, 11096072698276303650])
    if n < 7999252175582851:
        return mr(
            n,
            [
                2,
                4130806001517,
                149795463772692060,
                186635894390467037,
                3967304179347715805,
            ],
        )
    if n < 585226005592931977:
        return mr(
            n,
            [
                2,
                123635709730000,
                9233062284813009,
                43835965440333360,
                761179012939631437,
                1263739024124850375,
            ],
        )
    if n < 18446744073709551616:
        return mr(n, [2, 325, 9375, 28178, 450775, 9780504, 1795265022])
    if n < 318665857834031151167461:
        return mr(n, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37])
    if n < 3317044064679887385961981:
        return mr(n, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41])

    return mr(n, [2]) and is_strong_lucas_prp(n)


st = time.time()
r = 10 ** 9
isprime(r)
r = 10 ** 18
isprime(r)
for i in range(10000):
    r = random.randint(10 ** 9, 2 * 10 ** 9)
    isprime(r)
for i in range(100000):
    r = random.randint(10 ** 18, 2 * 10 ** 18)
    isprime(r)
print(f"{time.time() - st} sec")
