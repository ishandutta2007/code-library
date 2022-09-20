def meissel(n=2e6):
    primes = list(primerange(1, int(n ** (1 / 2)) + 1))
    cr_primes = [x for x in primes if x < n ** (1 / 3)]

    def prime_sum(n):
        ps = list(primerange(1, n + 1))
        return sum(ps)

    @lru_cache(maxsize=32)
    def sigma(x, a):
        if a == 1:
            return (x // 2) ** 2 if x % 2 == 0 else (x // 2 + 1) ** 2
        else:
            return sigma(x, a - 1) - primes[a - 1] * sigma(x // primes[a - 1], a - 1)

    def total(x, a):
        res, index = 0, a
        for p in primes[a:]:
            res += p * (prime_sum(x // p) - sum(primes[:index]))
            index += 1
        return res

    a = len(cr_primes)
    ans = sigma(n, a) + sum(cr_primes) - total(n, a) - 1
    return ans
