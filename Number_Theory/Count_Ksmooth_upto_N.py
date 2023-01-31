# n = 1e+10
# 0.4 sec in pypy2.7 and 5 sec in python 3 (ideone)

def count_k_smooth(k, N): 
    assert k ** 2 <= N < (k + 1) ** 2

    lo = [i - 1 for i in range(k + 1)]
    hi = [0] + [N // i - 1 for i in range(1, k + 1)]

    tot = N
    for p in range(2, k + 1): 
        if lo[p] == lo[p - 1]: continue

        tot -= p # all multiples of p less than p * p is not

        p_cnt = lo[p - 1]
        q = p * p
        end = min(k, N // q)
        for i in range(1, end + 1): 
            d = i * p
            if d <= k: hi[i] -= hi[d] - p_cnt
            else: hi[i] -= lo[N // d] - p_cnt
        for i in range(k, q - 1, -1): 
            lo[i] -= lo[i // p] - p_cnt

    for k in range(1, k): 
        tot -= k * (hi[k] - hi[k + 1])

    tot -= k * (hi[k] - lo[k]) # correction: for primes in (k, N // k]
    return tot

print(count_k_smooth(10 ** 5,10 ** 10))

