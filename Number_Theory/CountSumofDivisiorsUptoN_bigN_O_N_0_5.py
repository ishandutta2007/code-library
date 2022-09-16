import math


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


n = 10
for n in range(1, 11):
    print(A006218(n) - n)
