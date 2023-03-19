# smallest positive integer m such that n divides m!.
# But this doesnt just use sieve . It uses the trick similar to legendre's prime counting

from math import sqrt
from itertools import count


def prime_sieve(N):
    is_prime = [1] * (N + 1)
    is_prime[0] = 0
    v = isqrt(N)
    for p in range(2, v + 1):
        if not is_prime[p]:
            continue
        for k in range(p * p, N + 1, p):
            is_prime[k] = 0
    return [p for p in range(2, N + 1) if is_prime[p]]


def isqrt(n):
    x = int(sqrt(n * (1 + 1e-14)))
    while True:
        y = (x + n // x) >> 1
        if y >= x:
            return x
        x = y


def icbrt(n):
    if n <= 0:
        return 0
    x = int(n ** (1.0 / 3.0) * (1 + 1e-12))
    while True:
        y = (2 * x + n // (x * x)) // 3
        if y >= x:
            return x
        x = y


def tabulate_all_prime_sum(N):
    def T(n):
        return n * (n + 1) // 2 - 1

    if N <= 1:
        return [0, 0], [0, 0]

    v = isqrt(N)

    smalls = [T(i) for i in range(v + 1)]
    larges = [0 if i == 0 else T(N // i) for i in range(v + 1)]

    for p in range(2, v + 1):
        if smalls[p - 1] == smalls[p]:
            continue
        p_sum = smalls[p - 1]
        q = p * p
        end = min(v, N // q)
        for i in range(1, end + 1):
            d = i * p
            if d <= v:
                larges[i] -= (larges[d] - p_sum) * p
            else:
                larges[i] -= (smalls[N // d] - p_sum) * p
        for i in range(v, q - 1, -1):
            smalls[i] -= (smalls[i // p] - p_sum) * p
    return smalls, larges


N = 10 ** 8


def rec(n, beg, s, primes):
    ret = s
    for pi in range(beg, len(primes)):
        p = primes[pi]
        if p > n:
            break
        if p > s and p * p > n:
            ret += larges[N // n] if n > sqrtN else smalls[n]
            ret -= smalls[p - 1 if p <= sqrtN else sqrtN]
            break
        q = 1
        for e in count(1):
            q *= p
            if q > n:
                break
            ret += rec(n // q, pi + 1, max(s, ss[pi][e]), primes)
    return ret


sqrtN = isqrt(N)
smalls, larges = tabulate_all_prime_sum(N)
primes = prime_sieve(sqrtN)
primes += [sqrtN + 1]  # dummy

ans = 0
ss = []
for p in primes:
    q = p
    c, t, e = 0, 0, 1
    seq = [0]
    while q <= N:
        while c < e:
            t += p
            s = t
            while s % p == 0:
                s //= p
                c += 1
        seq += [t]
        q *= p
        e += 1
    ss += [seq]
ans += rec(N, 0, 0, primes)
print(ans)
