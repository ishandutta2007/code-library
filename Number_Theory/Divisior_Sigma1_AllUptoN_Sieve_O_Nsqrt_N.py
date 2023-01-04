from math import sqrt

# Sieve method is used to compute for all

MAX = 1000000


def sum_proper_divisors(n):
    if n <= 1:
        return 0
    divisors = [1]
    for i in range(2, int(sqrt(n)) + 1):
        if n % i == 0:
            if i == n // i:
                divisors += [i]
            else:
                divisors += [i, n // i]
    # print(n, divisors)
    return sum(divisors)


sigma1 = [0] * (MAX + 1)
sigma1[0] = 0
sigma1[1] = 1
for i in range(2, MAX + 1):
    sigma1[i] = sum_proper_divisors(i)
