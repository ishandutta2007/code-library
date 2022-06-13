import sys


def setup(a, b):
    global m, r, n, K, w, theta, s, N
    m = 10
    r = 5

    n = 1 << m
    K = r * n
    w = 1 << (r << 1)
    theta = 1 << r
    s = (K - m - 1) // 2
    N = s * n
    if a.bit_length() + b.bit_length() > N:
        print("Input is too large", file = sys.stderr)
        sys.exit()
    print("m = ", m)
    print("r = ", r)
    print("n = ", n)
    print("K = ", K)
    print("w = ", w)
    print("theta = ", theta)
    print("s = ", s)
    print("N = ", N)


def normalize(x):
    if x < 0:
        x += (1 << K) + 1
    return x % ((1 << K) + 1)


def mul(x, y):
    return normalize(x * y)


def add(x, y):
    return normalize(x + y)


def sub(x, y):
    return normalize(x - y)


def left_shift(x, nbits):
    nbits = nbits % (2 * K)
    if nbits >= K:
        return normalize(-(x << (nbits - K)))
    return normalize(x << nbits)


def bit_rev(x, l):
    return sum(1 << (l - 1 - i) for i in range(l) if (x >> i) & 1)


def bit_reverse(a):
    l = (len(a) - 1).bit_length()
    pos = [bit_rev(n, l) for n in range(len(a))]
    return [a[i] for i in pos]


def fft_loop(a, l, pos, depth):
    l1 = l >> 1
    if l > 1:
        fft_loop(a, l1,  pos, depth + 1)
        fft_loop(a, l1,  pos + l1, depth + 1)
    for i in range(l1):
        ide = i + pos
        ido = i + pos + l1
        # w ^ (i * 2 ^ depth) = 2 ^ (2 * r * i * 2 ^ depth) = 2 ^ (r * i * 2 ^
        # (depth + 1))
        term = mul(left_shift(1, r * i << (depth + 1)), a[ido])
        a[ido] = sub(a[ide], term)
        a[ide] = add(a[ide], term)


def ifft_loop(a, l, pos, depth):
    l1 = l >> 1
    if l > 1:
        ifft_loop(a, l1,  pos, depth + 1)
        ifft_loop(a, l1,  pos + l1, depth + 1)
    for i in range(l1):
        ide = i + pos
        ido = i + pos + l1
        # w ^ (-i * 2 ^ depth) = w ^ ((n - 1) * i * 2 ^ depth)
        # = 2 ^ ((n - 1) * 2 * r * i * 2 ^ depth)
        # = 2 ^ ((n - 1) * r * i * 2 ^ (depth + 1))
        # (depth + 1))
        term = mul(left_shift(1, (n - 1) * r * i << (depth + 1)), a[ido])
        a[ido] = sub(a[ide], term)
        a[ide] = add(a[ide], term)


def ifft(a):
    a_copy = bit_reverse(a[:] + [0] * (n - len(a)))
    ifft_loop(a_copy, len(a_copy),  0, 0)
    # n = 2 ^ m => 1/n = 2 ^ -m = 2 ^ (2K - m)
    return [left_shift(i, 2 * K - m) for i in a_copy]



def fft(a):
    a_copy = bit_reverse(a[:] + [0] * (n - len(a)))
    fft_loop(a_copy, len(a_copy),  0, 0)
    return a_copy


def gen_data():
    return 1403333232314433, 1231213323345552


def split_data(a):
    ret = []
    while a > 0:
        ret.append(a & int('1' * s, 2))
        a >>= s
    return ret


def combine_data(arr):
    ret = 0
    arr.reverse()
    for i in arr:
        ret <<= s
        ret += i
    arr.reverse()
    return ret


def acyclic_leverage(a):
    # a[i] * theta ^ i = a[i] * 2 ^ (r * i) = a[i] << (r * i) mod= a[i] << ((r *
    # i) % (2k))
    return [left_shift(a[i], r * i) for i in range(len(a))]


def acyclic_recover(a):
    # theta ^ (2n) mod= 1 => theta ^ (-k) mod= theta ^ (2n - k) = 2 ^ (r * (2n -
    # k))
    return [left_shift(a[i], r * (2 * n - i)) for i in range(len(a))]


def main():
    a, b = gen_data()
    setup(a,b)
    a_arr = split_data(a)
    b_arr = split_data(b)
    a_acyc = acyclic_leverage(a_arr)
    b_acyc = acyclic_leverage(b_arr)
    a_fft = fft(a_acyc)
    b_fft = fft(b_acyc)
    c_fft = []
    for i in range(n):
        c_fft.append(mul(a_fft[i], b_fft[i]))
    c_acyc = ifft(c_fft)
    c_arr = acyclic_recover(c_acyc)
    c = combine_data(c_arr)
    print("a = %d" % a)
    print("b = %d" % b)
    print("a * b = %d" % (a * b))
    print("algorithm result = %d" % c)
    # print("a array = ", a_arr)
    # print("b array = ", b_arr)
    # print("a acyclic = ", a_acyc)
    # print("b acyclic = ", b_acyc)
    # print("a fft = ", a_fft)
    # print("b fft = ", b_fft)
    # print("c fft = ", c_fft)
    # print("c acyclic = ", c_acyc)
    # print("c array = ", c_arr)


main()
