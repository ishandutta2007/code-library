import math

# sum of divisors = d(1)+d(2)+d(3)+...+d(n)
# It can be computed in O(sqrt(n)) without having ro compute same for any of 1 to N-1
# Cumulative sum of divisors c[i]=d(i)+C[i-1]
# For a diven n i can compute sum of all divisors


def isqrt(n):
    if n > 0:
        x = 1 << (n.bit_length() + 1 >> 1)
        while True:
            y = (x + n // x) >> 1
            if y >= x:
                return x
            x = y
    elif n == 0:
        return 0


def A006218(n):
    return 2 * sum(n // k for k in range(1, isqrt(n) + 1)) - isqrt(n) ** 2


MAX = 100
for n in range(1, MAX + 1):
    print(n, "=>", A006218(n) - n)
