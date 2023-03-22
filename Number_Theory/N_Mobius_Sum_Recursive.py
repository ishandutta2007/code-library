#  https://math.stackexchange.com/questions/4661944/calculate-sum-k-1n-k-cdot-muk/4662128#4662128
#  \sum k*phi(k)
#  O(n^{3/4}) Time and O(n^{1/2}) Space
#  It is to be noted recursive version is of same speed as iterative if not a tad faster as the recusion depth is 1 so doesn't make any difference

# 10^4.0  :  -209614
# t: 0.009712934494018555
# 10^5.0  :  -4283169
# t: 0.025476932525634766
# 10^6.0  :  226244412
# t: 0.13820314407348633
# 10^7.0  :  10564681565
# t: 0.7327280044555664
# 10^8.0  :  180523857298
# t: 4.295886039733887
# 10^9.0  :  -732495361483
# t: 23.655582189559937
# 10^10.0  :  -332900112839490
# t: 138.67993998527527
# 10^11.0  :  -8439137382134677
# t: 776.7232780456543
# Assuming python as source of truth I found c++ long long has conflict on 10^10

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


def sum_k_mu(x):
    if x == 1:
        snu[x] = 1
        return snu[x]
    if x <= N_SQ:
        if snu[x] is not None:
            return snu[x]
    else:
        if snu_large_inv[NN // x] is not None:
            return snu_large_inv[NN // x]

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


for i in range(1, 30 + 1):
    snu_large_inv = [None] * MAX_NSQ
    NN = i
    N_SQ = int(i ** 0.5)
    print(i, " : ", sum_k_mu(i))

ns = [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9, 10 ** 10, 10 ** 11]
for n in ns:
    start_time = time.time()
    snu_large_inv = [None] * MAX_NSQ
    NN = n
    N_SQ = int(n ** 0.5)
    ret = sum_k_mu(n)
    print("10^{}".format(math.log10(n)), " : ", ret)
    print("t:", time.time() - start_time)
