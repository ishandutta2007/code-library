import sympy
import bisect
import random
import math
import time

# Without any precomputation it can answer 20 queries in 1 sec
# but if n is small (n<10^7) and we are allowed to prcompute all mobius then it can answer 10^4 queries in ~1 sec

def T(n):
    return n * (n + 1) // 2


def sum_of_coprimes_of_n_till_k(n, k):
    divs = list(sympy.divisors(n))
    ans = 0
    for d in divs:
        ans += sympy.mobius(d) * d * T(k//d)
    return ans


start_time = time.time()
print(
    sum_of_coprimes_of_n_till_k(
        2 * 3 * 5 * 7 * 11 * 13 * 17 * 19, 2 * 3 * 5 * 7 * 11 * 13
    )
)
print((time.time() - start_time))
