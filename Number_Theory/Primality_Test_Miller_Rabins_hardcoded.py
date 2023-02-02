import time
import random


def mr_pass(a, s, d, n):
    x = pow(a, d, n)
    if x == 1:
        return True
    for i in range(s - 1):
        if x == n - 1:
            return True
        x = (x * x) % n
    return x == n - 1


def is_prime_miller(n):
    if n in [
        2,
        3,
        5,
        7,
        11,
        13,
        17,
        19,
        23,
        29,
        31,
        37,
        41,
        47,
        53,
        59,
        61,
        67,
        71,
        79,
        83,
        89,
        97,
        101,
        103,
        107,
        109,
        113,
    ]:
        return True
    if n < 2047:
        prime = [2, 3]
    if n < 1373653:
        prime = [2, 3]
    if n < 9080191:
        prime = [31, 73]
    if n < 25326001:
        prime = [2, 3, 5]
    if n < 3215031751:
        prime = [2, 3, 5, 7]
    if n < 4759123141:
        prime = [2, 7, 61]
    if n < 1122004669633:
        prime = [2, 13, 23, 1662803]
    if n < 2152302898747:
        prime = [2, 3, 5, 7, 11]
    if n < 3474749660383:
        prime = [2, 3, 5, 7, 11, 13]
    if n < 341550071728321:
        prime = [2, 3, 5, 7, 11, 13, 17]
    if n < 3825123056546413051:
        prime = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    if n < 18446744073709551616:
        prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    if n < 318665857834031151167461:
        prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    if n < 3317044064679887385961981:
        prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    else:
        prime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67, 71]
    d = n - 1
    s = 0
    while d % 2 == 0:
        d >>= 1
        s += 1
    for repeat in range(len(prime)):
        a = prime[repeat]
        if not mr_pass(a, s, d, n):
            return False
    return True


print(17, is_prime_miller(17))
print(100, is_prime_miller(100))

st = time.time()
for i in range(100000):
    r = random.randint(10**18, 2 * 10**18)
    is_prime_miller(r)
print(f"{time.time() - st} sec")
