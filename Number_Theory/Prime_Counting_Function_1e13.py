# number of primes
# Thanks for pyranine
# URL: https://judge.yosupo.jp/submission/31819
def Prime_Pi(N):
    """number of primes less than or equal to N
    N: int
    """

    if N < 2:
        return 0
    v = int(N**0.5) + 1
    smalls = [i // 2 for i in range(1, v + 1)]
    smalls[1] = 0
    s = v // 2
    roughs = [2 * i + 1 for i in range(s)]
    larges = [(N // roughs[i] + 1) // 2 for i in range(s)]
    skip = [False] * v

    pc = 0
    for p in range(3, v):
        if smalls[p] <= smalls[p - 1]:
            continue

        q = p * p
        pc += 1
        if q * q > N:
            break
        skip[p] = True
        for i in range(q, v, 2 * p):
            skip[i] = True

        ns = 0
        for k in range(s):
            i = roughs[k]
            if skip[i]:
                continue
            d = i * p
            larges[ns] = (
                larges[k] - (larges[smalls[d] - pc] if d < v else smalls[N // d]) + pc
            )
            roughs[ns] = i
            ns += 1
        s = ns
        for j in range((v - 1) // p, p - 1, -1):
            c = smalls[j] - pc
            e = min((j + 1) * p, v)
            for i in range(j * p, e):
                smalls[i] -= c

    for k in range(1, s):
        m = N // roughs[k]
        s = widths[k] - (pc + k - 1)
        for l in range(1, k):
            p = roughs[l]
            if p * p > m:
                break
            s -= smalls[m // p] - (pc + l - 1)
        larges[0] -= s

    return larges[0]
