#Credits: Jatin Yadav

from math import log
import math, bisect
from fractions import gcd
import sys, time, os
from random import randint
from multiprocessing import Pool, Process
from array import *

def sqrt_mod_p(n, p):
    n %= p
    if n <= 1 : return n
    if pow(n, (p - 1) / 2, p) != 1:
        return -1
    z, Q, M = 2, p - 1, 0
    while pow(z, (p - 1) / 2, p) == 1: z += 1
    while Q % 2 == 0: Q, M = Q / 2, M + 1
    c, t, R = pow(z, Q, p), pow(n, Q, p), pow(n, (Q + 1) / 2, p)
    while t > 1:
        t2 = (t * t) % p
        for i in xrange(1, M):
            if t2 == 1:
                b = pow(c, 1 << (M - i - 1), p)
                M = i
                c = (b * b) % p
                t = (t * c) % p
                R = (R * b) % p
                break
            t2 = (t2 * t2) % p
    return R

def get_factors(n):
    global factors
    factors = array('l', [0]*n)
    primes = array('l')
    for i in xrange(2, n):
        if factors[i] == 0:
            for j in xrange(i, n, i):
                factors[j] = i
            primes.append(i)
    return primes
MAXP =  3 * 10 ** 5
small_primes = get_factors(MAXP)

def ceil_kth_root(n, k):
    lo, hi = 1, int(10 * (n ** (1. / k)))
    while lo < hi:
        mid = (lo + hi) / 2
        if (mid ** k) >= n: hi = mid
        else: lo = mid + 1
    return lo

def getFirst(p, L, r):
    md = L % p
    if md <= r:
        return L + r - md
    else :
        return L + r + p - md

def inv(x, mod):
    return pow(x, mod - 2, mod)

def hensel_lift(pe, p, x, n):
    k = (n - x ** 2) / pe
    k = (k * inv( (2 * x) % p, p)) % p
    if k < 0: k += p
    return k * pe + x

# https://rosettacode.org/wiki/Miller%E2%80%93Rabin_primality_test
def _try_composite(a, d, n, s):
    if pow(a, d, n) == 1:
        return False
    for i in xrange(s):
        if pow(a, 2**i * d, n) == n-1:
            return False
    return True # n  is definitely composite

def is_prime(n, _precision_for_huge_n=16):
    if n < 2: return False
    if n <= _known_primes[-1]:
        if n in _known_primes:
            return True
    for p in _known_primes:
        if n % p == 0:
            return False
    d, s = n - 1, 0
    while not d % 2:
        d, s = d >> 1, s + 1
    # Returns exact according to http://primes.utm.edu/prove/prove2_3.html
    if n < 1373653: 
        return not any(_try_composite(a, d, n, s) for a in (2, 3))
    if n < 25326001: 
        return not any(_try_composite(a, d, n, s) for a in (2, 3, 5))
    if n < 118670087467: 
        if n == 3215031751: 
            return False
        return not any(_try_composite(a, d, n, s) for a in (2, 3, 5, 7))
    if n < 2152302898747:
        return not any(_try_composite(a, d, n, s) for a in (2, 3, 5, 7, 11))
    if n < 3474749660383: 
        return not any(_try_composite(a, d, n, s) for a in (2, 3, 5, 7, 11, 13))
    if n < 341550071728321:
        return not any(_try_composite(a, d, n, s) for a in (2, 3, 5, 7, 11, 13, 17))
    # otherwise
    return not any(_try_composite(a, d, n, s) 
                for a in _known_primes[:_precision_for_huge_n])


_known_primes = [2, 3]
_known_primes += [x for x in xrange(5, 1000, 2) if is_prime(x)]
 
def modinv(a, b):
    lr, r = abs(a), abs(b)
    x, lx, y, ly = 0, 1, 1, 0
    while r:
        lr, (q, r) = r, divmod(lr, r)
        x, lx = lx - q * x, x
        y, ly = ly - q * y, y
    
    if lr != 1:
        raise ValueError
    return (lx % b + b) % b

def nextprime(n):
    while not is_prime(n):
        n += 1
    return n

def get_primes(n):
    if n <= MAXP:
        return small_primes[:bisect.bisect_right(small_primes, n)]
    notPrime = array('B', [0] * n)
    primes = array('l')
    for i in xrange(2, n):
        if notPrime[i] == 0:
            for j in xrange(i, n, i):
                notPrime[j] = 1
            primes.append(i)
    return primes
print_process_wise_progress = False
if len(sys.argv) > 3 and sys.argv[3] in ["True", "1"]:
    print_process_wise_progress = True

class Parallel_Factorizer:
    def __init__(self, p):
        self.num_processes = p
        self.extra = 0
    # # For large integers only
    def parallel_mpqs(self, n):
        primes, primepowers, residues, logprime = self.primes, self.primepowers, self.residues, self.logprime
        MAX, basis, x, _p, _e = primes[-2], [], [], [], []
        M = 2 * MAX
        a, b, Q, iter, req = 0, 0, ceil_kth_root( ceil_kth_root(2 * n, 2) / M, 2), 0, len(primes) / self.num_processes + self.extra
        tried, success, intermediate = 0, 0, 0
        t1, t2, t3 = 0, 0, 0
        enter_time, REQ = time.time(), req
        q = Q
        while req > 0:
            find_enter = time.time()
            iter += 1
            if self.num_processes == 1:
                q = q + 1
            else:
                q = randint(Q, 2 * Q)
            while req > 0:
                q = nextprime(q + 1)
                b = sqrt_mod_p(n, q)
                if b != -1:
                    a = q * q
                    b = hensel_lift(q, q, b, n)
                    if n % q == 0: return [[-1], [q], [], []]
                    invq = modinv(q, n)
                    break
            sieve_start = time.time()
            t3 += sieve_start - find_enter
            log_estimate = array('B', [0]*(2 * M + 1))
            for i in xrange(len(primepowers)):
                pe, r, lg = primepowers[i], residues[i], logprime[i]
                inva = modinv(a, pe)
                R = inva * (r - b) % pe
                if R < 0: R += pe
                for k in xrange( M * (1 - pe) + getFirst(pe, M * (pe - 1), R), 2 * M + 1, pe):
                    log_estimate[k] += lg
                if not (pe & (pe - 1)):
                    continue
                R = inva * (-r - b) % pe
                if R < 0: R += pe
                for k in xrange( M * (1 - pe) + getFirst(pe, M * (pe - 1), R), 2 * M + 1, pe):
                    log_estimate[k] += lg
            sieve_end = time.time()
            t1 += sieve_end - sieve_start
            threshold, logMAX = int(log(n / a, 2) - log(MAX, 2) - 0.5), log(MAX) / log(2) / 1.4
            for k in xrange(0, 2 * M + 1):
                if log_estimate[k] < threshold or req <= 0: continue
                intermediate += 1
                _v =(a * (k - M) + b) ** 2 - n
                v = abs(_v / a)
                if log_estimate[k] < int(log(v + 1, 2) - logMAX): continue
                tried += 1
                num, nonzeroP, nonzeroE = 0, array('l'), array('B')
                for j in xrange(len(primes) - 1):
                    p, Ej = primes[j], 0
                    num <<= 1
                    while v % p == 0: v, Ej = v / p, Ej + 1
                    if Ej != 0:
                        nonzeroP.append(j)
                        nonzeroE.append(Ej)
                    num += Ej % 2
                num <<= 1
                if _v < 0: num += 1
                if v == 1:
                    success += 1
                    basis.append(num)
                    x.append( invq * (a * (k - M) + b) % n)
                    _p.append(nonzeroP)
                    _e.append(nonzeroE)
                    req -= 1
            t2 += time.time() - sieve_end
            if iter % 100 == 0 and print_process_wise_progress:
                print "Process ", os.getpid(), ":\n"
                print iter, "polynomials used"
                print intermediate, "values of f(x) considered for factoring"
                print tried, "values of f(x) actually factored" 
                print "Found ", success, "congruences, ", req, "yet to be found"
                print "Total time till now ->", time.time() - enter_time
                print "Time taken to find new prime q ->", t3
                print "Time taken to find inverses and sieve ->", t1
                print "Time taken in factoring ->", t2
                print "Expected time remaining = ", ((time.time() - enter_time) * req) / (REQ - req), "seconds\n\n................................................\n\n"
                
        return [basis, x, _p, _e]
    def factor(self, n):
        while True:
            self.extra  += 10
            ret = self.factor_utility(n)
            if ret != None:
                return ret

    def factor_utility(self, n):
        if n < MAXP: return factors[n]
        for p in small_primes:
            if n % p == 0:
                return p
        enter_time = time.time()
        init_B = 2 * int(math.exp( 0.5 * math.sqrt(0.707 * log( n ) * log( log( n ) ) ) ) ) 
        init_primes, primes, primepowers, residues, logprime = get_primes(init_B), array('l'), array('l'), array('l'), array('B')
        MAX = init_primes[-1]
        for p in init_primes:
            lg = int(log(p) / log(2) + 0.5)
            y = sqrt_mod_p(n, p)
            if y != -1:
                primes.append(p)
                if p == 2:
                    for e in xrange(1, 100):
                        pe = p**e
                        if pe > MAX:
                            break
                        for j in xrange(0, pe):
                            if (j * j - n) % (pe) == 0:
                                primepowers.append(pe)
                                residues.append(j)
                                logprime.append(lg)
                    continue
                pe = p
                while pe <= MAX:
                    primepowers.append(pe)
                    residues.append(y)
                    logprime.append(lg)
                    y = hensel_lift(pe, p, y, n)
                    pe *= p
        primes.append(1)
        self.primes = primes
        self.primepowers = primepowers
        self.residues = residues
        self.logprime = logprime
        print "Factorizing ", len(str(n)), "digit number"
        print "Size of factor base -> ", len(primes)
        print "Collecting data in parallel using", self.num_processes, "processes"
        if self.num_processes == 1:
            data = [self.parallel_mpqs(n)]
        else:
            pool = Pool(self.num_processes)
            data = pool.map(self.parallel_mpqs, [n] * self.num_processes)

        basis, masks, x, _p, _e = [], [], [], [], []
        print "Time taken to gather data -> ", time.time() - enter_time, "seconds"
        print "Solving equations ..."
        eq_begin = time.time()
        for pid in xrange(len(data)):
            D = data[pid]
            if D[0][0] == -1:
                return D[1][0]
            for j in xrange(len(D[0])):
                num = D[0][j]
                _x = D[1][j]
                fact_p = D[2][j]
                fact_e = D[3][j]
                mask = 1 << len(x)
                E = [0] * len(primes)
                for k in xrange(len(fact_p)):
                    E[fact_p[k]] = fact_e[k]
                for pos in xrange(len(x)):
                    if basis[pos] & -basis[pos] & num:
                        num ^= basis[pos]
                        mask ^= masks[pos]
                if num == 0:
                    X, Y = _x, 1
                    for pos in xrange(len(x)):
                        if (mask >> pos & 1) == 1:
                            X = (X * x[pos]) % n
                            for j in xrange(len(_p[pos])): E[_p[pos][j]] += _e[pos][j]
                    for j in xrange(len(primes)):
                        Y = (Y * pow(primes[j], E[j] / 2, n)) % n
                    assert((X**2 - Y**2) % n == 0)
                    ret = gcd(X + Y, n)
                    if ret != 1 and ret != n:
                        print "Time taken to solve equations -> ", time.time() - eq_begin, "seconds"
                        return ret
                    else:
                        continue
                for pos in xrange(len(x)):
                    if num & -num & basis[pos]:
                        basis[pos] ^= num
                        masks[pos] ^= mask
                masks.append(mask)
                basis.append(num)
                x.append(_x)
                _p.append(fact_p)
                _e.append(fact_e)
def factor(n):
    num_processes = 1
    if n > 10**20 and len(sys.argv) > 2: num_processes = int(sys.argv[2])
    return Parallel_Factorizer(num_processes).factor(n)

def factorize(n):
    D = {}
    if n == 1:
        return D
    if is_prime(n):
        if n in D:
            D[n] += 1
        else: D[n] = 1
        return D
    for k in xrange(2, int(log(n) / log(2)) + 1):
        t = ceil_kth_root(n, k)
        if t ** k == n:
            d = factorize(t)
            for i in d:
                if i in D:
                    D[i] += d[i] * k
                else:
                    D[i] = d[i] * k
            return D
    
    x = factor(n)
    D = factorize(n / x)
    d = factorize(x)

    for i in d:
        if D.has_key(i):
            D[i] += d[i]
        else:
            D[i] = d[i]
    return D
n = int(sys.argv[1])
ans = str(n)+" = "
if n == 0:
    sys.exit()
F = factorize(n)
keys = []
for i in F: keys.append(i)
keys.sort()
for i in keys:
    ans += str(i) +"^" + str(F[i]) + " * "
print ans[:-3]
