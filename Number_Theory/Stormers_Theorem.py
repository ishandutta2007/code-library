from math import sqrt
from time import time
import pprint as pp

inf = float("inf")


def isqrt(n):
    return int(sqrt(n))


try:
    from gmpy2 import mpz, isqrt

    mpzv, inttypes = 2, (int, type(mpz(1)))
except ImportError:
    mpz, mpzv, inttypes = int, 0, (int,)


def issmooth(n, ps):
    for p in ps:
        while n % p == 0:
            n //= p
        if n == 1:
            return True
    return False


def simplepell(D, bail=inf):  # copied from labmath (https://pypi.org/project/labmath/)
    """
    Generates the positive solutions of x**2 - D * y**2 == 1.  We use some
    optimizations specific to this case of the Pell equation that makes this
    more efficient than calling pell(D,1)[0].  Note that this function is not
    equivalent to calling pell(D,1)[0]: pell() is concerned with the general
    equation, which may or may not have trivial solutions, and as such yields
    all non-negative solutions, whereas this function is concerned only with the
    simple Pell equation, which always has an infinite family of positive
    solutions generated from a single primitive solution and always has the
    trivial solution (1,0); since that trivial solution always exists for this
    function's scope, we omit it from the output.
    Input:
        D -- an integer, assumed to be positive
        bail -- yield no solutions whose x-coordinate is > this number.
                Default == inf.
    Output: sequence of 2-tuples of positive integers
    Examples:
    >>> list(islice(simplepell(2), 6))
    [(3, 2), (17, 12), (99, 70), (577, 408), (3363, 2378), (19601, 13860)]
    >>> list(islice(simplepell(3), 7))
    [(2, 1), (7, 4), (26, 15), (97, 56), (362, 209), (1351, 780), (5042, 2911)]
    >>> next(simplepell(61))
    (1766319049, 226153980)
    >>> next(simplepell(661))
    (16421658242965910275055840472270471049, 638728478116949861246791167518480580)
    """
    d = isqrt(D)
    i, B0, G0, P, Q = False, 1, d, d, D - d * d
    if Q == 0:
        return
    B1 = a = (2 * d) // Q
    G1 = a * d + 1
    while Q != 1:
        P = a * Q - P
        Q = (D - P**2) // Q
        a = (P + d) // Q
        i, B1, B0, G1, G0 = not i, a * B1 + B0, B1, a * G1 + G0, G1
        if G0 > bail:
            return
    x, y = a, b = (G0, B0) if i else (G0**2 + D * B0**2, 2 * G0 * B0)
    while x <= bail:
        yield (x, y)
        x, y = x * a + y * b * D, y * a + x * b


def sqfrgen(ps):  # copied from labmath (https://pypi.org/project/labmath/)
    """
    Generates squarefree products of elements of ps.
    Input: ps -- indexable iterable of primes
    Output: sequence of integers
    Examples:
    >>> sorted(filter(lambda x: x < 100, sqfrgen(list(primegen(12)))))
    [1, 2, 3, 5, 6, 7, 10, 11, 14, 15, 21, 22, 30, 33, 35, 42, 55, 66, 70, 77]
    """
    if len(ps) == 0:
        yield 1
        return
    for n in sqfrgen(ps[1:]):
        yield n
        yield n * ps[0]


def stormer2(
    ps, *ps2, abc=None, procs=1
):  # derived from labmath (https://pypi.org/project/labmath/)
    """
    For any given set ps of prime numbers, there are only finitely many pairs of
    consecutive integers (n,n+2) that are both ps-smooth.  Stormer's theorem
    provides a method to find them all.  We implement Lehmer's simplification
    of that method.  It is worth noting that the time to complete this iteration
    for the first n primes appears to scale superexponentially in n, while
    iterating hamming() over the nth-prime-smooth numbers up to max(stormer2(1st
    n primes)) appears to scale singly exponentially; however, max(stormer2(ps))
    cannot yet be computed without actually executing the Stormer-Lehmer
    algorithm.
    Let S be a set of primes, let x and x+2 be S-smooth, and let T be the
    product of the elements of S.  Then on the abc conjecture we have
    x+2 < k * rad(2 * x * (x+2)) ** d < k * T**d.  This enables a major speedup.
    Input:
        ps -- indexable iterable whose elements are assumed to be prime
        abc -- Assume an effective abc conjecture of the form
               c < abc[0] * rad(a*b*c)**abc[1].
               Default == None; i.e., make no assumptions.
        procs -- Use this many processes to work in parallel.  Default == 1.
                 When procs == 1, we don't invoke multiprocessing an instead
                 work directly in the main thread, because multiprocessing has
                 significant overhead.
    Output: finite sequence of pairs of integers
    Example:
    """
    if isinstance(ps, inttypes):
        ps = [ps] + list(ps2)
    pl = [mpz(x) for x in set(ps)]
    k = max(3, (max(pl) + 1) // 2)
    bail = 2 + 2 * abc[0] * prod(pl) ** abc[1] if abc else inf
    if procs == 1:
        for sqfr in sqfrgen(pl):
            for n, (x, y) in enumerate(simplepell(sqfr, bail)):
                if n >= k:
                    break
                # We now check that we have found a smooth pair.  We don't outsource to factorint since we only need to divide
                # out the small primes and check whether anything remains --- we don't need to get stuck factoring RSA numbers.
                if issmooth(x - 1, pl):
                    if issmooth(x + 1, pl):
                        yield (x - 1, x + 1)
                elif n == 0:
                    break
                # Pell solutions have some nice divisibility properties that allow us to build a sieve to knock out subsequent
                # solutions if one of them turns out to be rough.  In the simplest case, if the fundamental solution is rough,
                # then we can skip all subsequent solutions from that Pell equation.  This is a major speedup.
                # See https://projecteuclid.org/download/pdf_1/euclid.ijm/1256067456 page 11/67.
        return


start = time()

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]  # , 41, 43 , 41, 47]
solutions = list(stormer2(primes))
solutions.sort()
only_even_by2solutions = [
    (num1 // 2, num2 // 2)
    for num1, num2 in solutions
    if num1 % 2 == 0 and num2 % 2 == 0
]
# print(solutions)
pp.pprint(only_even_by2solutions)
ns_only = [pair[0] for pair in only_even_by2solutions]
print("max({}-smooth) = {} sum(ns) = {}".format(primes[-1], max(ns_only), sum(ns_only)))
print(time() - start, "secs")
