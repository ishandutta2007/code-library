import sympy

# prime_pow_delta[i] = 0 if i is not a power of prime
# prime_pow_delta[i] = p^k - p^(k-1) if i is kth power of prime p


def prime_power_delta(sqrt_N):
    primes = list(sympy.primerange(sqrt_N))
    pis = [0] * (sqrt_N + 1)
    prime_pow_delta = [0] * (sqrt_N + 1)
    for p in primes:
        pis[p] = 1
        q = p
        p2 = p * p
        while q <= sqrt_N:
            prime_pow_delta[q] = q if q % p2 != 0 else q - q // p
            q *= p
    return prime_pow_delta


# lpp_sum[i] = sum of largest powers of prime under i
# lpp_sum[10] = 2^3 + 3^2 + 5 + 7
# lpp_sum[20] = 2^4 + 3^2 + 5 + 7 + 11 + 13 + 17 + 19
# lpp_sum[i] = 2^[log(i,2)] + 3^[log(i,3)] + 5^[log(i,5)] + ... + p_l^[log(i,p_l)] (where p_l is the largest prime under i)


def prime_power_delta_cumulitive(prime_pow_delta, sqrt_N):
    lpp_sum = [prime_pow_delta[0], prime_pow_delta[1]]
    for i in range(2, sqrt_N):
        lpp_sum.append(lpp_sum[i - 1] + prime_pow_delta[i])
    return lpp_sum


N = 20
prime_pow_delta = prime_power_delta(N)
print(
    "prime power delta ie p^k-p^(k-1)                               = ", prime_pow_delta
)
ppdc = prime_power_delta_cumulitive(prime_pow_delta, N)
print("largest prime power sum p1^k1+p2^k2+... (such that p2^k2 <= i) = ", ppdc)

# This concept was used in https://projecteuler.net/problem=347
