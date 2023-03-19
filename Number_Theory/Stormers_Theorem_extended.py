#! /usr/bin/env python3

# This program generates various sequences associated with Stormer's theorem, including:
# A002071
# A002072
# A117581 - 1
# A117582
# A117583
# A145604
# A145606
# A175607 - 1
# A181471

from itertools import count, islice, takewhile, compress
from multiprocessing import Pool
from math import prod, isqrt
from time import time
from sys import argv

inf = float("inf")

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


def primegen(limit=inf):  # copied from labmath (https://pypi.org/project/labmath/)
    """
    Generates primes < limit lazily via a segmented sieve of Eratosthenes.
    Memory use is, possibly up to some logarithmic factors, O(sqrt(p)),
    where p is the most recently yielded prime.
    Input: limit -- a number (default = inf)
    Output: sequence of integers
    Examples:
    >>> list(islice(primegen(), 20))
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    >>> list(primegen(73))
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    """
    # We don't sieve 2, so we ought to be able to get sigificant savings by halving the length of the sieve.
    # But the tiny extra computation involved in that seems to exceed the savings.
    yield from takewhile(
        lambda x: x < limit, (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    )
    pl, pg = [3, 5, 7], primegen()
    for p in pl:
        next(pg)
    n = next(pg)
    nn = n * n
    while True:
        n = next(pg)
        ll, nn = nn, n * n
        sl = nn - ll
        sieve = bytearray([True]) * sl
        for p in pl:
            k = (-ll) % p
            sieve[k::p] = bytearray([False]) * ((sl - k) // p + 1)
        if nn > limit:
            break  # TODO bringing this condition up to the while statement may increase speed.
        yield from (ll + k for k in compress(range(0, sl, 2), sieve[::2]))
        pl.append(n)
    yield from takewhile(
        lambda x: x < limit, (ll + k for k in compress(range(0, sl, 2), sieve[::2]))
    )


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
        Q = (D - P ** 2) // Q
        a = (P + d) // Q
        i, B1, B0, G1, G0 = not i, a * B1 + B0, B1, a * G1 + G0, G1
        if G0 > bail:
            return
    x, y = a, b = (G0, B0) if i else (G0 ** 2 + D * B0 ** 2, 2 * G0 * B0)
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


def stormerfactory(data):  # derived from labmath (https://pypi.org/project/labmath/)
    sqfr, pl, bail, k = data
    ans = []
    for n, (x, y) in enumerate(simplepell(sqfr, bail)):
        if n >= k:
            break
        # We now check that we have found a smooth pair.  We don't outsource to factorint since we only need to divide
        # out the small primes and check whether anything remains --- we don't need to get stuck factoring RSA numbers.
        if issmooth(x - 1, pl):
            if issmooth(x + 1, pl):
                ans.append((x - 1, x + 1))
        elif n == 0:
            break
        # Pell solutions have some nice divisibility properties that allow us to build a sieve to knock out subsequent
        # solutions if one of them turns out to be rough.  In the simplest case, if the fundamental solution is rough,
        # then we can skip all subsequent solutions from that Pell equation.  This is a major speedup.
        # See https://projecteuclid.org/download/pdf_1/euclid.ijm/1256067456 page 11/67.
    return ans


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
    with Pool(procs) as procpool:
        chunksize = max(1, 2 ** (len(pl) // 2))
        for ans in procpool.imap_unordered(
            stormerfactory,
            ((sqfr, pl, bail, k) for sqfr in sqfrgen(pl)),
            chunksize=chunksize,
        ):
            yield from ans


# Now we parse the arguments.

print("Invoke this program as")
print()
print("    ./stormer.py [abc=x,y] [p=z]")
print()
print("The order of the arguments may be swapped.")
print("If abc=x,y is present, then we assume an abc conjecture c < x * rad(abc)^y.")
print("If p=z is present, then we use z processes.  By default, z = 1.")
print()
print("------------------------")
print()

try:
    abc, pc = None, 1
    for arg in argv[1:]:
        a, b = arg.split("=")
        if a == "abc":
            x, y = b.split(",")
            abc = (int(x), int(y))
        if a == "p":
            pc = int(b)
except:
    print("Error while parsing arguments.")
    exit()

# The arguments are parsed.  Now tell the user how they were interpreted.

print(
    "Assuming no abc conjectures."
    if abc == None
    else "Assuming an abc conjecture: c < %d * rad(abc)^%d." % abc
)
print("Using 1 process." if pc == 1 else "Using %d processes." % pc)
print(
    "Using native arithmetic (gmpy2 not found)."
    if len(inttypes) == 1
    else "Using gmpy2."
)
print()

# Now we begin computation.

# Note that every smooth pair (x, x+1) yields a smooth pair (x, x+2) upon doubling,
#  and every even smooth pair (x, x+2) yields a smooth pair (x, x+1) upon halving.

primes = []
for n, p in enumerate(primegen(), start=1):
    primes.append(p)
    oldtime = time()
    (
        total,
        maxsol,
        totalnew,
        maxnew,
        total1,
        max1,
        total1new,
        max1new,
        total2,
        max2,
        total2new,
        max2new,
        total1sq,
        total1tri,
    ) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    for x, y in stormer2(primes, abc=abc, procs=pc):
        assert x + 2 == y, (x, y, primes, abc)
        assert issmooth(x, primes), (x, y, primes, abc)
        assert issmooth(y, primes), (x, y, primes, abc)
        # total     = A#?         = number of x such that x and x+2 are both primes-smooth.
        # maxsol    = A#?         = maximum   x such that x and x+2 are both primes-smooth.
        # totalnew  = A181471     = number of x such that x and x+2 are both primes-smooth and          one is divisible by p.
        # maxnew    = A175607 - 1 = largest   x such that x and x+2 are both primes-smooth and          one is divisible by p.
        # total1    = A002071     = number of x such that x and x+1 are both primes-smooth.
        #                           number of x such that x and x+2 are both primes-smooth and even.
        # max1      = A002072     = largest   x such that x and x+1 are both primes-smooth.  Equivalently,
        #                           largest   x such that x and x+2 are both primes-smooth and even.
        # total1new = A145604     = number of x such that x and x+1 are both primes-smooth and          one is divisible by p.  Equivalently,
        #                           number of x such that x and x+2 are both primes-smooth and even and one is divisible by p.
        # max1new   = A145606     = largest   x such that x and x+1 are both primes-smooth and          one is divisible by p.  Equivalently,
        #                           largest   x such that x and x+2 are both primes-smooth and even and one is divisible by p.
        # total2    = A#?         = number of x such that x and x+2 are both primes-smooth and odd.  Equivalently,
        #                           number of x such that x and x+2 are both primes-smooth and (x, x+2) does not arise from an (x, x+1) pair.
        # max2      = A#?         = largest   x such that x and x+2 are both primes-smooth and odd.  Equivalently,
        #                           largest   x such that x and x+2 are both primes-smooth and (x, x+2) does not arise from an (x, x+1) pair.
        # total2new = A#?         = number of x such that x and x+2 are both primes-smooth and odd and  one is divisible by p.  Equivalently,
        #                           number of x such that x and x+2 are both primes-smooth and (x, x+2) does not arise from an (x, x+1) pair.
        # max2new   = A#?         = largest   x such that x and x+2 are both primes-smooth and odd and  one is divisible by p.  Equivalently,
        #                           number of x such that x and x+2 are both primes-smooth and (x, x+2) does not arise from an (x, x+1) pair.
        # total1sq  = A117582     = number of x such that x and x+1 are both primes-smooth and x+1 is square.
        # total1tri = A117583     = number of x such that x and x+1 are both primes-smooth and x+1 is triangular.
        total += 1
        if x > maxsol:
            maxsol = x
        if x % 2 == y % 2 == 0:
            total1 += 1
            x2, y2 = x // 2, y // 2
            if x2 > max1:
                max1 = x2
            if y2 == isqrt(y2) ** 2:
                total1sq += 1
            yroot = isqrt(y)
            if y == yroot * (yroot + 1):
                total1tri += 1
        else:
            assert x % 2 == 1 or y % 2 == 1, (x, y, primes)
            total2 += 1
            if x > max2:
                max2 = x
        if x % p == 0 or y % p == 0:
            totalnew += 1
            if x > maxnew:
                maxnew = x
            if x % 2 == y % 2 == 0:
                total1new += 1
                if x // 2 > max1new:
                    max1new = x // 2
            else:
                assert x % 2 == 1 or y % 2 == 1, (x, y, primes)
                total2new += 1
                if x > max2new:
                    max2new = x
        print(
            "\b" * 1000,
            n,
            primes[-1],
            total,
            maxsol,
            totalnew,
            maxnew,
            total1,
            max1,
            total1new,
            max1new,
            total2,
            max2,
            total2new,
            max2new,
            total1sq,
            total1tri,
            end="   ",
            flush=True,
        )
    print(time() - oldtime)


"""
Running this code on my computer generated the following data.
Note that the results for n >= 20 are conditional on the effective abc conjecture c < rad(abc)**2.

                                                                                    A117581 - 1
   A000040                                A181471  A175607 - 1              A002071 A002072                  A145604   A145606                                                                           A117582  A117583
n  pmax    total maxsol                   totalnew maxnew                   total1  max1                     total1new max1new                  total2 max2                total2new max2new             total1sq total1tri Assumptions Time
1  2       1     2                        1        2                        1       1                        1         1                        0      0                   0         0                   0        0         none        9.870529174804688e-05
2  3       5     16                       4        16                       4       8                        3         8                        1      1                   1         1                   2        1         none        0.000141143798828125
3  5       13    160                      8        160                      10      80                       6         80                       3      25                  2         25                  5        3         none        0.0003108978271484375
4  7       29    8748                     16       8748                     23      4374                     13        4374                     6      243                 3         243                 10       7         none        0.0008060932159423828
5  11      49    19600                    20       19600                    40      9800                     17        9800                     9      243                 3         75                  15       9         none        0.0018792152404785156
6  13      83    246400                   34       246400                   68      123200                   28        123200                   15     1573                6         1573                24       16        none        0.004002571105957031
7  17      130   672280                   47       672280                   108     336140                   40        336140                   22     2023                7         2023                34       22        none        0.01859140396118164
8  19      202   23718420                 72       23718420                 167     11859210                 59        11859210                 35     3969                13        3969                46       29        none        0.018877267837524414
9  23      297   23718420                 95       10285000                 241     11859210                 74        5142500                  56     1447873             21        1447873             57       35        none        0.0306851863861084
10 29      423   354365440                126      354365440                345     177182720                104       177182720                78     26578123            22        26578123            74       39        none        0.06382274627685547
11 31      591   3222617398               168      3222617398               482     1611308699               137       1611308699               109    287080365           31        287080365           90       50        none        0.17664265632629395
12 37      799   9447152317               208      9447152317               653     3463199999               171       3463199999               146    9447152317          37        9447152317          114      57        none        0.6038854122161865
13 41      1061  127855050750             262      127855050750             869     63927525375              216       63927525375              192    9447152317          46        2470954913          141      68        none        2.3266193866729736
14 43      1404  842277599278             343      842277599278             1153    421138799639             284       421138799639             251    9447152317          59        1181631185          174      84        none        10.393783569335938
15 47      1837  2218993446250            433      2218993446250            1502    1109496723125            349       1109496723125            335    9447152317          84        4768304959          208      100       none        56.3782160282135
16 53      2344  2907159732048            507      2907159732048            1930    1453579866024            428       1453579866024            414    9447152317          79        907922169           244      112       none        382.6211712360382
17 59      2978  41257182408960           634      41257182408960           2454    20628591204480           524       20628591204480           524    122187528125        110       122187528125        287      127       none        3632.322960615158
18 61      3777  63774701665792           799      63774701665792           3106    31887350832896           652       31887350832896           671    122187528125        147       37517343435         334      151       none        55125.35413146019
19 67      4753  63774701665792           976      25640240468750           3896    31887350832896           790       12820120234375           857    932784765625        186       932784765625        387      167       none        968867.9461529255

19 67      4753  63774701665792           976      25640240468750           3896    31887350832896           790       12820120234375           857    932784765625        186       932784765625        387      167       abc=(1,2)   19.819456577301025
20 71      5899  238178082107392          1146     238178082107392          4839    119089041053696          943       119089041053696          1060   932784765625        203       430188203733        433      186       abc=(1,2)   58.99822473526001
21 73      7338  4573663454608288         1439     4573663454608288         6040    2286831727304144         1201      2286831727304144         1298   2704807796873       238       2704807796873       494      204       abc=(1,2)   125.97661519050598
22 79      9036  19182937474703818750     1698     19182937474703818750     7441    9591468737351909375      1401      9591468737351909375      1595   6774628298373       297       6774628298373       561      230       abc=(1,2)   280.43266677856445
23 83      11118 19182937474703818750     2082     34903240221563712        9179    9591468737351909375      1738      17451620110781856        1939   473599589105797     344       473599589105797     635      253       abc=(1,2)   612.6705558300018
24 89      13489 19182937474703818750     2371     332110803172167360       11134   9591468737351909375      1955      166055401586083680       2355   473599589105797     416       211089142289023     709      275       abc=(1,2)   1440.6109733581543
25 97      16223 19182937474703818750     2734     99913980938200000        13374   9591468737351909375      2240      49956990469100000        2849   473599589105797     494       151377426693205     787      304       abc=(1,2)   2827.535365343094
26 101     19583 19182937474703818750     3360     8216517931479010998      16167   9591468737351909375      2793      4108258965739505499      3416   30945502301275315   567       30945502301275315   871      327       abc=(1,2)   1839.5827414989471 (p=5)
27 103     23605 38632316754147847668000  4022     38632316754147847668000  19507   19316158377073923834000  3340      19316158377073923834000  4098   30945502301275315   682       2728575616834375    973      354       abc=(1,2)   3429.167060136795 (p=5)
28 107     28242 38632316754147847668000  4637     773079686222382448       23367   19316158377073923834000  3860      386539843111191224       4875   30945502301275315   777       15562980136200913   1096     394       abc=(1,2)   6616.486833810806 (p=5)
29 109     33817 38632316754147847668000  5575     181101212761682433220    27949   19316158377073923834000  4582      90550606380841216610     5868   31584322190711673   993       31584322190711673   1227     432       abc=(1,2)   13841.13929772377 (p=5)
30 113     40241 38632316754147847668000  6424     410284126426376207278    33233   19316158377073923834000  5284      205142063213188103639    7008   95281435950754443   1140      95281435950754443   1355     462       abc=(1,2)   29441.672268152237 (p=5)
31 127     47509 38632316754147847668000  7268     106469590255765459648    39283   19316158377073923834000  6050      53234795127882729824     8226   767308615386551775  1218      767308615386551775  1498     496       abc=(1,2)   60842.86489534378 (p=5)
32 131     55860 38632316754147847668000  8351     8228608891233272032062   46166   19316158377073923834000  6883      4114304445616636016031   9694   6845433178730663473 1468      6845433178730663473 1619     529       abc=(1,2)   124984.86140799522 (p=5)
33 137     65521 248451871690466638878346 9661     248451871690466638878346 54150   124225935845233319439173 7984      124225935845233319439173 11371  6845433178730663473 1677      3425684783380442539 1773     551       abc=(1,2)   263680.98630714417 (p=5)
"""
