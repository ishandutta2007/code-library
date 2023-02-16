import logging

import math
import cmath
import operator
import itertools
import functools
import fractions
import numbers
import random
import bitarray

# credits: https://github.com/goulu/Goulib


def identity(x):
    """Do nothing and return the variable untouched"""
    return x


def compress(iterable, key=identity, buffer=None):
    """
    generates (item,count) pairs by counting the number of consecutive items in iterable)
    :param iterable: iterable, possibly infinite
    :param key: optional function defining which elements are considered equal
    :param buffer: optional integer. if defined, iterable is sorted with this buffer
    """
    key = key or identity

    if buffer:
        iterable = sorted_iterable(iterable, key, buffer)

    prev, count = None, 0
    for item in iterable:
        if count and key(item) == key(prev):
            count += 1
        else:
            if prev is not None:  # to skip initial junk
                yield prev, count
            count = 1
            prev = item
    if count:
        yield prev, count


def nextprime(n):
    """Determines, with some semblance of efficiency, the least prime number strictly greater than n."""
    # from https://pypi.python.org/pypi/primefac
    if n < 2:
        return 2
    if n == 2:
        return 3
    n = (n + 1) | 1  # first odd larger than n
    m = n % 6
    if m == 3:
        if is_prime(n + 2):
            return n + 2
        n += 4
    elif m == 5:
        if is_prime(n):
            return n
        n += 2
    for m in itertools.count(n, 6):
        if is_prime(m):
            return m
        if is_prime(m + 4):
            return m + 4


def primes_gen(start=2, stop=None):
    """generate prime numbers from start"""
    if start == 1:
        yield 1  # if we asked for it explicitly
        start = 2
    if stop is None or stop > start:
        n = start - 1  # to include start if it is prime
        while True:
            n = nextprime(n)
            if (stop is None) or n <= stop:
                yield n
            else:
                break
    else:  # backwards
        n = start + 1
        while True:
            n = prevprime(n)
            if n and n >= stop:
                yield n
            else:
                break


def prime_factors(num, start=2):
    """generates all prime factors (ordered) of num"""
    for p in primes_gen(start):
        if num == 1:
            break
        if is_prime(num):  # because it's fast
            yield num
            break
        if p > num:
            break
        while num % p == 0:
            yield p
            num = num // p


def factorize(n):
    """find the prime factors of n along with their frequencies. Example:
    >>> factor(786456)
    [(2,3), (3,3), (11,1), (331,1)]
    """

    if n == 1:  # allows to make many things quite simpler...
        return [(1, 1)]
    return compress(prime_factors(n))


class Sieve:
    # should be derived from bitarray but ...
    # https://github.com/ilanschnell/bitarray/issues/69
    # TODO: simplify when solved
    def __init__(self, f, init):
        self._ = bitarray.bitarray(init)
        self.f = f

    def __len__(self):
        return len(self._)

    def __getitem__(self, index):
        return self._[index]

    def __call__(self, n):
        self.resize(n)
        return (i for i, v in enumerate(self._) if v)

    def resize(self, n):
        l = len(self) - 1
        if n <= l:
            return
        n = int(n)  # to tolerate n=1E9, which is float
        self._.extend([True] * (n - l))
        for i in self(n):
            if i == 2:
                i2, s = 4, 2
            else:
                i2, s = self.f(i)
            if i2 > n:
                break
            self._[i2::s] = False  # bitarray([False]*int((n-i2)/s+1))


def erathostene(n):
    return n * n, 2 * n


# array of bool indicating primality
_sieve = Sieve(erathostene, [False, False, True])


def sieve(n, oneisprime=False):
    """prime numbers from 2 to a prime < n"""
    res = _sieve(n)
    if oneisprime:
        res = itertools.chain([1], res)
    return list(res)


_primes = sieve(1000)  # primes up to 1000
_primes_set = set(_primes)  # to speed us primality tests below


def is_prime(n, oneisprime=False, tb=(3, 5, 7, 11), eb=(2,), mrb=None):
    """main primality test.
    :param n: int number to test
    :param oneisprime: bool True if 1 should be considered prime (it was, a long time ago)
    :param tb: trial division basis
    :param eb: Euler's test basis
    :param mrb: Miller-Rabin basis, automatic if None
    :see: https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
    It’s an implementation of the BPSW test (Baillie-Pomerance-Selfridge-Wagstaff)
    with some prefiltes for speed and is deterministic for all numbers less than 2^64
    Iin fact, while infinitely many false positives are conjectured to exist,
    no false positives are currently known.
    The prefilters consist of trial division against 2 and the elements of the tuple tb,
    checking whether n is square, and Euler’s primality test to the bases in the tuple eb.
    If the number is less than 3825123056546413051, we use the Miller-Rabin test
    on a set of bases for which the test is known to be deterministic over this range.
    """
    # https://pypi.python.org/pypi/primefac

    if n <= 0:
        return False

    if n == 1:
        return oneisprime

    if n < len(_sieve):
        return _sieve[n]

    if n in _primes_set:
        return True

    if any(n % p == 0 for p in tb):
        return False

    if is_square(n):
        return False  # it's quick ...

    if not is_prime_euler(n):
        return False  # Euler's test

    s, d = pfactor(n)
    if not sprp(n, 2, s, d):
        return False
    if n < 2047:
        return True

    # BPSW has two phases: SPRP with base 2 and SLPRP.
    # We just did the SPRP; now we do the SLPRP
    if n >= 3825123056546413051:
        d = 5
        while True:
            if gcd(d, n) > 1:
                p, q = 0, 0
                break
            if jacobi(d, n) == -1:
                p, q = 1, (1 - d) // 4
                break
            d = -d - 2 * d // abs(d)
        if p == 0:
            return n == d
        s, t = pfactor(n + 2)
        u, v, u2, v2, m = 1, p, 1, p, t // 2
        k = q
        while m > 0:
            u2, v2, q = (u2 * v2) % n, (v2 * v2 - 2 * q) % n, (q * q) % n
            if m % 2 == 1:
                u, v = u2 * v + u * v2, v2 * v + u2 * u * d
                if u % 2 == 1:
                    u += n
                if v % 2 == 1:
                    v += n
                u, v, k = (u // 2) % n, (v // 2) % n, (q * k) % n
            m //= 2
        if (u == 0) or (v == 0):
            return True
        for i in range(1, s):
            v, k = (v * v - 2 * k) % n, (k * k) % n
            if v == 0:
                return True
        return False

    # Miller-Rabin
    if not mrb:
        if n < 1373653:
            mrb = [3]
        elif n < 25326001:
            mrb = [3, 5]
        elif n < 3215031751:
            mrb = [3, 5, 7]
        elif n < 2152302898747:
            mrb = [3, 5, 7, 11]
        elif n < 3474749660383:
            mrb = [3, 5, 6, 11, 13]
        elif n < 341550071728321:
            # This number is also a false positive for primes(19+1).
            mrb = [3, 5, 7, 11, 13, 17]
        elif n < 3825123056546413051:
            # Also a false positive for primes(31+1).
            mrb = [3, 5, 7, 11, 13, 17, 19, 23]
    return all(sprp(n, b, s, d) for b in mrb)


def kempner(n):
    """ "Kempner function, also called Smarandache function
    :return: int smallest positive integer m such that n divides m!.
    :param n: int
    :see: https://en.wikipedia.org/wiki/Kempner_function
    :see: http://mathworld.wolfram.com/SmarandacheFunction.html
    """
    if n == 1:
        return 1
    if is_prime(n):
        return n

    # @decorators.memoize
    def _np(n, p):
        # n^p . use https://codereview.stackexchange.com/a/129868/37671
        k = 0
        while p > n:
            k += n
            p -= n + 1
            t = k
            while t % n != 0:
                t = t // n
                p -= 1
        p = max(0, p)

        return (k + p) * n

    return max(_np(f, p) for f, p in factorize(n))


print(kempner(10))
