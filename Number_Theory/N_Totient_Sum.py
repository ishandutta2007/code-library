# https://math.stackexchange.com/questions/4661467/calculate-sum-k-1n-k-cdot-varphik?noredirect=1&lq=1
# O(n^{3/4}) (for sum k * u(k)) + O(n^{1/2}) (for sum k * phi(k) using O(1) access to sum k * u(k))

# i*phi(i)
# i_phi_i_good(10^4)= 202643891472849
# time: 2.7748780250549316 432
# =====
# i_phi_i_good(10^5)= 202642368741515819
# time: 8.449854850769043 1368
# =====
# i_phi_i_good(10^7)= 202642380629476099463
# time: 25.690521955490112 4324
# =====
# i_phi_i_good(10^8)= 202642367994273571457613
# time: 80.64575815200806 13676
# =====
# i_phi_i_good(10^10.999999995739572)= 202642361326330276710977268231599
# time: 5015.025189876556 632435
# =====

import math
import sympy
from functools import lru_cache


def i_phi_i_brute(n):
    s = 0
    iphi_sigma = []
    # A002618
    for i in range(1, n + 1):
        v = i * sympy.totient(i)
        # print(v, end=",")
        s += v
        iphi_sigma.append(s)

    # A011755
    # print()
    # print(iphi_sigma)
    return iphi_sigma[-1]


def S(n):
    return n * (n + 1) * (2 * n + 1) // 6


import time
import math

MAX_N = 10 ** 11
MAX_NSQ = 10 ** 6
NN = None
N_SQ = None


def T(n):
    return n * (n + 1) // 2


snu = [None] * MAX_NSQ
snu_large_inv = [None] * MAX_NSQ

g_sum_k_mu_ctr = 0


@lru_cache(maxsize=None)
def sum_k_mu(x):
    global snu
    global snu_large_inv
    if x == 1:
        snu[x] = 1
        return snu[x]
    if x <= N_SQ:
        if snu[x] is not None:
            return snu[x]
    else:
        if snu_large_inv[NN // x] is not None:
            return snu_large_inv[NN // x]
    global g_sum_k_mu_ctr
    g_sum_k_mu_ctr += 1
    # print("called sum_k_mu with x=", x)

    x_sr = int(x ** 0.5)
    ans = 1

    for d in range(2, x_sr + 1):
        ans -= d * sum_k_mu(x // d)

    st = x // x_sr - 1
    for m in range(st, 0, -1):
        ans -= (T(x // m) - T(x // (m + 1))) * sum_k_mu(m)

    if x <= N_SQ:
        snu[x] = ans
        return snu[x]
    else:
        snu_large_inv[NN // x] = ans
        return snu_large_inv[NN // x]


def i_phi_i_good(n):
    global NN
    global N_SQ
    global snu_large_inv
    n_sr = int(n ** 0.5)
    s = 0
    for d in range(1, n_sr + 1):
        s += sympy.mobius(d) * d * S(n // d)
    # print("s after 1st half = ", s)

    for v in range(1, n // (n_sr + 1) + 1):
        l = n // (v + 1)
        NN = l
        N_SQ = int(l ** 0.5)
        snu_large_inv = [None] * MAX_NSQ
        s_l = sum_k_mu(l)

        r = n // v
        NN = r
        N_SQ = int(r ** 0.5)
        snu_large_inv = [None] * MAX_NSQ
        s_r = sum_k_mu(r)

        sum_k_mu_range = s_r - s_l
        s += S(v) * sum_k_mu_range

    return s


ns = [
    # 10,
    100,
    # 1000,
    # 10 ** 4,
    # 10 ** 5,
    # 10 ** 6,
    # 10 ** 7,
    # 10 ** 8,
    # 10 ** 9,
    # 10 ** 10,
    # 10 ** 11,
    99999999019
]
for n in ns:
    if n <= 10 ** 5:
        start = time.time()
        print("i_phi_i_brute({})=".format(n), i_phi_i_brute(n))
        print("brute time", time.time() - start)
    start = time.time()
    g_sum_k_mu_ctr = 0
    print("i_phi_i_good(10^{})=".format(math.log10(n)), i_phi_i_good(n))
    print("time:", time.time() - start, g_sum_k_mu_ctr)
    print("=====")
