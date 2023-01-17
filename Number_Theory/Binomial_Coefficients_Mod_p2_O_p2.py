# ncr modulo p^2 where p is 10^7
# n and r are also 64 bit
# Time complexity O(p^2)
# For O(p) check the BinomialCoefficients_p2.cpp


def countps(n):
    pp = p * p
    return n // p + n // pp


def factorial_stripped(n, p):
    pp = p * p
    f_stripped = 1
    for i in range(1, n + 1):
        v = i
        while v % p == 0:
            v = v // p
        f_stripped = f_stripped * v % pp
        if f_stripped % p == 0:
            print(n, f_stripped, v, i)
            raise Exception("WTF")
    return f_stripped


def binom_mod_p2(n, r, p):
    pp = p * p
    nr = n - r
    excessp = countps(n) - countps(r) - countps(nr)

    if excessp >= 2:
        ans = 0
    else:
        try:
            f_n = factorial_stripped(n, p)
            f_r = factorial_stripped(r, p)
            f_nr = factorial_stripped(nr, p)
            f_r_inv = pow(f_r, -1, pp)
            f_nr_inv = pow(f_nr, -1, pp)
            ans = f_n
            ans = ans * f_r_inv % pp
            ans = ans * f_nr_inv % pp
            ans = ans * (p**excessp) % pp
        except Exception as e:
            print(e)
            print(n, r)
            print(f_n)
            print(f_r)
            print(f_nr)
            print(f_r_inv)
            print(f_nr_inv)
            raise e

    return ans


p = 10**4 + 7
n = 99930048
r = n // 2

print("binom_mod_p2(", n, ",", r, ",", p, ")=", binom_mod_p2(n, r, p))
