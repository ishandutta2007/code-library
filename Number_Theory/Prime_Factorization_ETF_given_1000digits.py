import random
import sys


def gcd(a, b):
    while b > 0:
        a, b = b, a % b
    return a


def rabinMiller(num):
    s = num - 1
    t = 0
    while s % 2 == 0:
        s = s >> 1
        t += 1
    for trials in range(5):
        a = random.randrange(2, num - 1)
        v = pow(a, s, num)
        if v != 1:
            i = 0
            while v != (num - 1):
                if i == t - 1:
                    return False
                else:
                    i = i + 1
                    v = (v ** 2) % num
    return True


def prime(n):
    if n < 2:
        return False
    elif n < 4:
        return True
    elif n & 1 == 0:
        return False
    else:
        return rabinMiller(n)


def fact(n, phi):
    while phi & 1 == 0:
        phi >>= 1
    for base in range(1, 1000, 2):
        a = pow(base, phi, n)
        if a == 1:
            continue
        while a != n - 1:
            b = a
            a = a * a % n
            if a == 1:
                return gcd(b + 1, n)
            if a == b:
                break
    assert 0


def split(n, p):
    if prime(n):
        return set([n])
    y = fact(n, p)
    a, b = split(y, p), split(n // y, p)
    return set.union(a, b)


t = int(input())
for _ in range(t):
    n, p = map(int, input().split())
    ans = []
    if n & 1 == 0:
        p <<= 1
        c = 0
        while n & 1 == 0:
            n, p, c = n >> 1, p >> 1, c + 1
        ans.append((2, c))

    while n != 1:
        g = gcd(n, p)
        num, den = n // g, p // g
        out = split(num, p)
        for x in out:
            p = (p * x) // (x - 1)
            c = 0
            while n % x == 0:
                n, p, c = n // x, p // x, c + 1
            ans.append((x, c))
    print("total prime factors", len(ans))
    ans.sort()
    for x in ans:
        print("prime=", x[0], ":", x[1], "times")
