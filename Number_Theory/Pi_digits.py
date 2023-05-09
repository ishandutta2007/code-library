import time
import os
import math
from random import randint
from bisect import bisect

COS_SIN_CACHE_PREC = 400
COS_SIN_CACHE_STEP = 8
EXP_SERIES_U_CUTOFF = 1500
EXP_COSH_CUTOFF = 600

_1_800 = 1 << 800
_1_600 = 1 << 600
_1_400 = 1 << 400
_1_200 = 1 << 200
_1_100 = 1 << 100
_1_50 = 1 << 50

CHUD_A = 13591409
CHUD_B = 545140134
CHUD_C = 640320
CHUD_D = 12

fzero = (0, 0, 0, 0)
fnzero = (1, 0, 0, 0)
fone = (0, 1, 0, 1)
fhalf = (0, 1, -1, 1)
intern = lambda x: x

round_nearest = intern("n")
round_floor = intern("f")
round_ceiling = intern("c")
round_up = intern("u")
round_down = intern("d")
round_fast = round_down

small_trailing = [0] * 256
for j in range(1, 8):
    small_trailing[1 << j :: 1 << (j + 1)] = [j] * (1 << (7 - j))


def trailing(n):
    n = abs(int(n))
    if not n:
        return 0
    low_byte = n & 0xFF
    if low_byte:
        return small_trailing[low_byte]
    t = 8
    n >>= 8
    return t + small_trailing[n & 0xFF]


shifts_down = {
    round_floor: (1, 0),
    round_ceiling: (0, 1),
    round_down: (1, 1),
    round_up: (0, 0),
}
trailtable = [trailing(n) for n in range(256)]


def normalize(sign, man, exp, bc, prec, rnd):
    n = bc - prec
    if n > 0:
        if shifts_down[rnd][sign]:
            man >>= n
        exp += n
        bc = prec
    if not man & 1:
        t = trailtable[int(man & 255)]
        if not t:
            while not man & 255:
                man >>= 8
                exp += 8
                bc -= 8
            t = trailtable[int(man & 255)]
        man >>= t
        exp += t
        bc -= t
    return sign, man, exp, bc


def normalize1(sign, man, exp, bc, prec, rnd):
    if bc <= prec:
        return sign, man, exp, bc
    n = bc - prec
    if shifts_down[rnd][sign]:
        man >>= n
    exp += n
    bc = prec
    if not man & 1:
        t = trailtable[int(man & 255)]
        if not t:
            while not man & 255:
                man >>= 8
                exp += 8
                bc -= 8
            t = trailtable[int(man & 255)]
        man >>= t
        exp += t
        bc -= t
    return sign, man, exp, bc


# def to_int(s):
#     sign, man, exp, bc = s
#     return man >> (-exp)


powers = [1 << _ for _ in range(300)]


def bitcount(n):
    bc = bisect(powers, n)
    if bc != 300:
        return bc
    bc = int(math.log(n, 2)) - 4
    return bc + bctable[n >> bc]


bctable = [bitcount(n) for n in range(1024)]


# def from_man_exp(man, exp, prec=None, rnd=round_fast):
#     man = int(man)
#     sign = 0
#     if man < 0:
#         sign = 1
#         man = -man
#     if man < 1024:
#         bc = bctable[int(man)]
#     else:
#         bc = bitcount(man)
#     if not prec:
#         if not man:
#             return fzero
#         if not man & 1:
#             if man & 2:
#                 return (sign, man >> 1, exp + 1, bc - 1)
#             t = trailtable[int(man & 255)]
#             if not t:
#                 while not man & 255:
#                     man >>= 8
#                     exp += 8
#                     bc -= 8
#                 t = trailtable[int(man & 255)]
#             man >>= t
#             exp += t
#             bc -= t
#         return (sign, man, exp, bc)
#     return normalize(sign, man, exp, bc, prec, rnd)


# def def_mpf_constant(fixed):
#     def f(prec, rnd=round_fast):
#         wp = prec + 20
#         v = fixed(wp)
#         return normalize(0, v, -wp, bitcount(v), prec, rnd)

#     f.__doc__ = fixed.__doc__
#     return f


def bs_chudnovsky(a, b, level):
    if b - a == 1:
        g = int((6 * b - 5) * (2 * b - 1) * (6 * b - 1))
        p = b**3 * CHUD_C**3 // 24
        q = (-1) ** b * g * (CHUD_A + CHUD_B * b)
    else:
        mid = (a + b) // 2
        g1, p1, q1 = bs_chudnovsky(a, mid, level + 1)
        g2, p2, q2 = bs_chudnovsky(mid, b, level + 1)
        p = p1 * p2
        g = g1 * g2
        q = q1 * p2 + q2 * g1
    return g, p, q


# def constant_memo(f):
#     f.memo_prec = -1
#     f.memo_val = None

#     def g(prec, **kwargs):
#         memo_prec = f.memo_prec
#         if prec <= memo_prec:
#             return f.memo_val >> (memo_prec - prec)
#         newprec = int(prec * 1.05 + 10)
#         f.memo_val = f(newprec, **kwargs)
#         f.memo_prec = newprec
#         return f.memo_val >> (newprec - prec)

#     g.__name__ = f.__name__
#     g.__doc__ = f.__doc__
#     return g


@constant_memo
def pi_fixed(prec):
    N = int(prec / 3.3219280948 / 14.181647462 + 2)
    g, p, q = bs_chudnovsky(0, N, 0)
    sqrtC = isqrt_fast(CHUD_C << (2 * prec))
    v = p * CHUD_C * sqrtC // ((q + CHUD_A * p) * CHUD_D)
    return v


# mpf_pi = def_mpf_constant(pi_fixed)

# int_cache = dict((n, from_man_exp(n, 0)) for n in range(-10, 257))


# def from_int(n, prec=0, rnd=round_fast):
#     if not prec:
#         if n in int_cache:
#             return int_cache[n]
#     return from_man_exp(n, 0, prec, rnd)


def mpf_mul(s, t, prec=0, rnd=round_fast):
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    sign = ssign ^ tsign
    man = sman * tman
    if man:
        bc = sbc + tbc - 1
        bc += int(man >> bc)
        if prec:
            return normalize1(sign, man, sexp + texp, bc, prec, rnd)
        else:
            return (sign, man, sexp + texp, bc)
    s_special = (not sman) and sexp
    t_special = (not tman) and texp
    return fzero


def mpf_div(s, t, prec, rnd=round_fast):
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    sign = ssign ^ tsign
    if tman == 1:
        return normalize1(sign, sman, sexp - texp, sbc, prec, rnd)
    extra = prec - sbc + tbc + 5
    if extra < 5:
        extra = 5
    quot, rem = divmod(sman << extra, tman)
    if rem:
        quot = (quot << 1) + 1
        extra += 1
        return normalize1(sign, quot, sexp - texp - extra, bitcount(quot), prec, rnd)
    return normalize(sign, quot, sexp - texp - extra, bitcount(quot), prec, rnd)


def mpf_add(s, t, prec=0, rnd=round_fast, _sub=0):
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    tsign ^= _sub
    if sman and tman:
        offset = sexp - texp
        if offset:
            if offset > 0:
                if offset > 100 and prec:
                    delta = sbc + sexp - tbc - texp
                if ssign == tsign:
                    man = tman + (sman << offset)
                else:
                    man = (sman << offset) - tman
                    ssign = 0
                bc = bitcount(man)
                return normalize1(ssign, man, texp, bc, prec or bc, rnd)
            elif offset < 0:
                if offset < -100 and prec:
                    delta = tbc + texp - sbc - sexp
                if ssign == tsign:
                    man = sman + (tman << -offset)
                else:
                    if tsign:
                        man = sman - (tman << -offset)
                    ssign = 0
                bc = bitcount(man)
                return normalize1(ssign, man, sexp, bc, prec or bc, rnd)
        if ssign == tsign:
            man = tman + sman
        else:
            man = sman - tman
            ssign = 0
        bc = bitcount(man)
        return normalize(ssign, man, texp, bc, prec or bc, rnd)
    if not sman:
        if tman:
            return normalize1(tsign, tman, texp, tbc, prec or tbc, rnd)
    return normalize1(ssign, sman, sexp, sbc, prec or sbc, rnd)


def mpf_sub(s, t, prec=0, rnd=round_fast):
    return mpf_add(s, t, prec, rnd, 1)


def from_rational(p, q, prec, rnd=round_fast):
    return mpf_div(from_int(p), from_int(q), prec, rnd)


def giant_steps(start, target, n=2):
    L = [target]
    while L[-1] > start * n:
        L = L + [L[-1] // n + 2]
    return L[::-1]


def isqrt_small(x):
    if x < _1_800:
        r = int(x**0.5 * 1.00000000000001) + 1
    while 1:
        y = (r + x // r) >> 1
        if y >= r:
            return r
        r = y


def isqrt_fast(x):
    if x < _1_800:
        y = int(x**0.5)
        if x >= _1_100:
            y = (y + x // y) >> 1
            if x >= _1_200:
                y = (y + x // y) >> 1
                if x >= _1_400:
                    y = (y + x // y) >> 1
        return y
    bc = bitcount(x)
    guard_bits = 10
    x <<= 2 * guard_bits
    bc += 2 * guard_bits
    bc += bc & 1
    hbc = bc // 2
    startprec = min(50, hbc)
    r = int(2.0 ** (2 * startprec) * (x >> (bc - 2 * startprec)) ** -0.5)
    pp = startprec
    for p in giant_steps(startprec, hbc):
        r2 = (r * r) >> (2 * pp - p)
        xr2 = ((x >> (bc - p)) * r2) >> p
        r = (r * ((3 << p) - xr2)) >> (pp + 1)
        pp = p
    return (r * (x >> hbc)) >> (p + guard_bits)


# def isqrt(x):
#     return sqrtrem(x)[0]


# def sqrtrem(x):
#     if x < _1_600:
#         y = isqrt_small(x)
#         return y, x - y * y
#     y = isqrt_fast(x) + 1
#     rem = x - y * y
#     while rem < 0:
#         y -= 1
#         rem += 1 + 2 * y
#     return y, rem


# def mpf_sqrt(s, prec, rnd=round_fast):
#     sign, man, exp, bc = s
#     if exp & 1:
#         exp -= 1
#         man <<= 1
#         bc += 1
#     elif man == 1:
#         return normalize1(sign, man, exp // 2, bc, prec, rnd)
#     shift = max(4, 2 * prec - bc + 4)
#     shift += shift & 1
#     if rnd in "fd":
#         man = isqrt(man << shift)
#     return from_man_exp(man, (exp - shift) // 2, prec, rnd)


# cos_sin_cache = {}


# def cos_sin_basecase(x, prec):
#     if prec > COS_SIN_CACHE_PREC:
#         return exponential_series(x, prec, 2)
#     precs = prec - COS_SIN_CACHE_STEP
#     t = x >> precs
#     n = int(t)
#     if n not in cos_sin_cache:
#         w = t << (10 + COS_SIN_CACHE_PREC - COS_SIN_CACHE_STEP)
#         cos_t, sin_t = exponential_series(w, 10 + COS_SIN_CACHE_PREC, 2)
#         cos_sin_cache[n] = (cos_t >> 10), (sin_t >> 10)
#     cos_t, sin_t = cos_sin_cache[n]
#     offset = COS_SIN_CACHE_PREC - prec
#     cos_t >>= offset
#     sin_t >>= offset
#     x -= t << precs
#     cos = 1 << prec
#     sin = x
#     k = 2
#     a = -((x * x) >> prec)
#     while a:
#         a //= k
#         cos += a
#         k += 1
#         a = (a * x) >> prec
#         a //= k
#         sin += a
#         k += 1
#         a = -((a * x) >> prec)
#     return ((cos * cos_t - sin * sin_t) >> prec), ((sin * cos_t + cos * sin_t) >> prec)


# def mod_pi2(man, exp, mag, wp):
#     if mag > 0:
#         i = 0
#         while 1:
#             cancellation_prec = 20 << i
#             wpmod = wp + mag + cancellation_prec
#             pi2 = pi_fixed(wpmod - 1)
#             pi4 = pi2 >> 1
#             offset = wpmod + exp
#             if offset >= 0:
#                 t = man << offset
#             n, y = divmod(t, pi2)
#             if y > pi4:
#                 small = pi2 - y
#             else:
#                 small = y
#             if small >> (wp + mag - 10):
#                 n = int(n)
#                 t = y >> mag
#                 wp = wpmod - mag
#                 break
#             i += 1
#     else:
#         wp += -mag
#         offset = exp + wp
#         if offset >= 0:
#             t = man << offset
#         n = 0
#     return t, n, wp


# def mpf_cos_sin(x, prec, rnd=round_fast, which=0):
#     sign, man, exp, bc = x
#     mag = bc + exp
#     wp = prec + 10
#     t, n, wp = mod_pi2(man, exp, mag, wp)
#     c, s = cos_sin_basecase(t, wp)
#     m = n & 3
#     if m == 1:
#         c, s = -s, c
#     elif m == 2:
#         c, s = -c, -s
#     elif m == 3:
#         c, s = s, -c
#     if which == 1:
#         return from_man_exp(c, -wp, prec, rnd)
#     if which == 2:
#         return from_man_exp(s, -wp, prec, rnd)


# def mpf_cos(x, prec, rnd=round_fast):
#     return mpf_cos_sin(x, prec, rnd, 1)


# def mpf_sin(x, prec, rnd=round_fast):
#     return mpf_cos_sin(x, prec, rnd, 2)


# def is_quad_residue(a, p):
#     a, p = int(a), int(p)
#     if a >= p or a < 0:
#         a = a % p
#     if a < 2 or p < 3:
#         return True
#     return pow(a, (p - 1) // 2, p) == 1


# def legendre_symbol(a, p):
#     a, p = int(a), int(p)
#     a = a % p
#     if pow(a, (p - 1) // 2, p) == 1:
#         return 1
#     return -1


# def igcdex(a, b):
#     x_sign = 1
#     y_sign = 1
#     x, y, r, s = 1, 0, 0, 1
#     while b:
#         (c, q) = (a % b, a // b)
#         (a, b, r, s, x, y) = (b, c, x - q * r, y - q * s, r, s)
#     return (x * x_sign, y * y_sign, a)


# def jacobi_symbol(a, n):
#     a %= n
#     result = 1
#     while a != 0:
#         while a % 2 == 0:
#             a /= 2
#             n_mod_8 = n % 8
#             if n_mod_8 in (3, 5):
#                 result = -result
#         a, n = n, a
#         if a % 4 == 3 and n % 4 == 3:
#             result = -result
#         a %= n
#     return result


# def _sqrt_mod_tonelli_shanks(a, p):
#     s = trailing(p - 1)
#     t = p >> s
#     while 1:
#         d = randint(2, p - 1)
#         r = legendre_symbol(d, p)
#         if r == -1:
#             break
#     A = pow(a, t, p)
#     D = pow(d, t, p)
#     m = 0
#     for i in range(s):
#         adm = A * pow(D, m, p) % p
#         adm = pow(adm, 2 ** (s - 1 - i), p)
#         if adm % p == p - 1:
#             m += 2**i
#     x = pow(a, (t + 1) // 2, p) * pow(D, m // 2, p) % p
#     return x


# def sqrt_mod_prime_power(a, p, k):
#     pk = p**k
#     a = a % pk
#     if k == 1:
#         if p % 4 == 3:
#             res = pow(a, (p + 1) // 4, p)
#         elif p % 8 == 5:
#             sign = pow(a, (p - 1) // 4, p)
#             if sign == 1:
#                 res = pow(a, (p + 3) // 8, p)
#             else:
#                 b = pow(4 * a, (p - 5) // 8, p)
#                 x = (2 * a * b) % p
#                 if pow(x, 2, p) == a:
#                     res = x
#         else:
#             res = _sqrt_mod_tonelli_shanks(a, p)
#         return sorted([res, p - res])
#     if k > 1:
#         if p == 2:
#             rv = [1, 3, 5, 7]
#             n = 3
#             res = []
#             for r in rv:
#                 nx = n
#                 while nx < k:
#                     r1 = (r**2 - a) >> nx
#                     if r1 % 2:
#                         r = r + (1 << (nx - 1))
#                     nx += 1
#                 if r not in res:
#                     res.append(r)
#                 x = r + (1 << (k - 1))
#                 if x < (1 << nx) and x not in res:
#                     if (x**2 - a) % pk == 0:
#                         res.append(x)
#             return res
#         rv = sqrt_mod_prime_power(a, p, 1)
#         r = rv[0]
#         fr = r**2 - a
#         n = 1
#         px = p
#         while 1:
#             n1 = n
#             n1 *= 2
#             if n1 > k:
#                 break
#             n = n1
#             px = px**2
#             frinv = igcdex(2 * r, px)[0]
#             r = (r - fr * frinv) % px
#             fr = r**2 - a
#         if n < k:
#             px = p**k
#             frinv = igcdex(2 * r, px)[0]
#             r = (r - fr * frinv) % px
#         return [r, px - r]


# def _pre():
#     maxn = 10**5
#     global _factor
#     global _totient
#     _factor = [0] * maxn
#     _totient = [1] * maxn
#     lim = int(maxn**0.5) + 5
#     for i in range(2, lim):
#         if _factor[i] == 0:
#             for j in range(i * i, maxn, i):
#                 if _factor[j] == 0:
#                     _factor[j] = i
#     for i in range(2, maxn):
#         if _factor[i] == 0:
#             _factor[i] = i
#             _totient[i] = i - 1
#             continue
#         x = _factor[i]
#         y = i // x
#         if y % x == 0:
#             _totient[i] = _totient[y] * x
#         else:
#             _totient[i] = _totient[y] * (x - 1)


# def _a(n, k, prec):
#     if k == 1:
#         return fone
#     k1 = k
#     e = 0
#     p = _factor[k]
#     while k1 % p == 0:
#         k1 //= p
#         e += 1
#     k2 = k // k1
#     v = 1 - 24 * n
#     pi = mpf_pi(prec)
#     if k1 == 1:
#         if p == 2:
#             mod = 8 * k
#             v = mod + v % mod
#             v = (v * pow(9, k - 1, mod)) % mod
#             m = sqrt_mod_prime_power(v, 2, e + 3)[0]
#             arg = mpf_div(mpf_mul(from_int(4 * m), pi, prec), from_int(mod), prec)
#             return mpf_mul(
#                 mpf_mul(
#                     from_int((-1) ** e * jacobi_symbol(m - 1, m)),
#                     mpf_sqrt(from_int(k), prec),
#                     prec,
#                 ),
#                 mpf_sin(arg, prec),
#                 prec,
#             )
#         if p == 3:
#             mod = 3 * k
#             v = mod + v % mod
#             if e > 1:
#                 v = (v * pow(64, k // 3 - 1, mod)) % mod
#             m = sqrt_mod_prime_power(v, 3, e + 1)[0]
#             arg = mpf_div(mpf_mul(from_int(4 * m), pi, prec), from_int(mod), prec)
#             return mpf_mul(
#                 mpf_mul(
#                     from_int(2 * (-1) ** (e + 1) * legendre_symbol(m, 3)),
#                     mpf_sqrt(from_int(k // 3), prec),
#                     prec,
#                 ),
#                 mpf_sin(arg, prec),
#                 prec,
#             )
#         v = k + v % k
#         if v % p == 0:
#             if e == 1:
#                 return mpf_mul(
#                     from_int(jacobi_symbol(3, k)), mpf_sqrt(from_int(k), prec), prec
#                 )
#             return fzero
#         if not is_quad_residue(v, p):
#             return fzero
#         _phi = p ** (e - 1) * (p - 1)
#         v = v * pow(576, _phi - 1, k)
#         m = sqrt_mod_prime_power(v, p, e)[0]
#         arg = mpf_div(mpf_mul(from_int(4 * m), pi, prec), from_int(k), prec)
#         return mpf_mul(
#             mpf_mul(
#                 from_int(2 * jacobi_symbol(3, k)), mpf_sqrt(from_int(k), prec), prec
#             ),
#             mpf_cos(arg, prec),
#             prec,
#         )
#     if p != 2 or e >= 3:
#         d1, d2 = math.gcd(k1, 24), math.gcd(k2, 24)
#         e = 24 // (d1 * d2)
#         n1 = (
#             (d2 * e * n + (k2**2 - 1) // d1)
#             * pow(e * k2 * k2 * d2, _totient[k1] - 1, k1)
#         ) % k1
#         n2 = (
#             (d1 * e * n + (k1**2 - 1) // d2)
#             * pow(e * k1 * k1 * d1, _totient[k2] - 1, k2)
#         ) % k2
#         return mpf_mul(_a(n1, k1, prec), _a(n2, k2, prec), prec)
#     if e == 2:
#         n1 = ((8 * n + 5) * pow(128, _totient[k1] - 1, k1)) % k1
#         n2 = (4 + ((n - 2 - (k1**2 - 1) // 8) * (k1**2)) % 4) % 4
#         return mpf_mul(mpf_mul(from_int(-1), _a(n1, k1, prec), prec), _a(n2, k2, prec))
#     n1 = ((8 * n + 1) * pow(32, _totient[k1] - 1, k1)) % k1
#     n2 = (2 + (n - (k1**2 - 1) // 8) % 2) % 2
#     return mpf_mul(_a(n1, k1, prec), _a(n2, k2, prec), prec)


# def mpf_shift(s, n):
#     sign, man, exp, bc = s
#     return sign, man, exp + n, bc


# def bsp_acot(q, a, b, hyperbolic):
#     if b - a == 1:
#         a1 = int(2 * a + 3)
#         if hyperbolic or a & 1:
#             return 1, a1 * q**2, a1
#     m = (a + b) // 2
#     p1, q1, r1 = bsp_acot(q, a, m, hyperbolic)
#     p2, q2, r2 = bsp_acot(q, m, b, hyperbolic)
#     return q2 * p1 + r1 * p2, q1 * q2, r1 * r2


# def acot_fixed(a, prec, hyperbolic):
#     N = int(0.35 * prec / math.log(a) + 20)
#     p, q, r = bsp_acot(a, 0, N, hyperbolic)
#     return ((p + q) << prec) // (q * a)


# def machin(coefs, prec, hyperbolic=False):
#     extraprec = 10
#     s = 0
#     for a, b in coefs:
#         s += int(a) * acot_fixed(int(b), prec + extraprec, hyperbolic)
#     return s >> extraprec


# @constant_memo
# def ln2_fixed(prec):
#     return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)


# def exponential_series(x, prec, type=0):
#     sign = 0
#     r = int(0.5 * prec**0.5)
#     xmag = bitcount(x) - prec
#     r = max(0, xmag + r)
#     extra = 10 + 2 * max(r, -xmag)
#     wp = prec + extra
#     x <<= extra - r
#     one = 1 << wp
#     alt = type == 2
#     if prec < EXP_SERIES_U_CUTOFF:
#         x2 = a = (x * x) >> wp
#         x4 = (x2 * x2) >> wp
#         s0 = s1 = 0
#         k = 2
#         while a:
#             a //= (k - 1) * k
#             s0 += a
#             k += 2
#             a //= (k - 1) * k
#             s1 += a
#             k += 2
#             a = (a * x4) >> wp
#         s1 = (x2 * s1) >> wp
#         if alt:
#             c = s1 - s0 + one
#         else:
#             c = s1 + s0 + one
#     else:
#         u = int(0.3 * prec**0.35)
#         x2 = a = (x * x) >> wp
#         xpowers = [one, x2]
#         for i in range(1, u):
#             xpowers.append((xpowers[-1] * x2) >> wp)
#         sums = [0] * u
#         k = 2
#         while a:
#             for i in range(u):
#                 a //= (k - 1) * k
#                 if alt and k & 2:
#                     sums[i] -= a
#                 else:
#                     sums[i] += a
#                 k += 2
#             a = (a * xpowers[-1]) >> wp
#         for i in range(1, u):
#             sums[i] = (sums[i] * xpowers[i]) >> wp
#         c = sum(sums) + one
#     if type == 0:
#         s = isqrt_fast(c * c - (one << wp))
#         v = c + s
#         for i in range(r):
#             v = (v * v) >> wp
#         return v >> extra
#     else:
#         pshift = wp - 1
#         for i in range(r):
#             c = ((c * c) >> pshift) - one
#         s = isqrt_fast(abs((one << wp) - c * c))
#         return (c >> extra), (s >> extra)


# def exp_basecase(x, prec):
#     if prec > EXP_COSH_CUTOFF:
#         return exponential_series(x, prec, 0)
#     r = int(prec**0.5)
#     prec += r
#     s0 = s1 = 1 << prec
#     k = 2
#     a = x2 = (x * x) >> prec
#     while a:
#         a //= k
#         s0 += a
#         k += 1
#         a //= k
#         s1 += a
#         k += 1
#         a = (a * x2) >> prec
#     s1 = (s1 * x) >> prec
#     s = s0 + s1
#     u = r
#     while r:
#         s = (s * s) >> prec
#         r -= 1
#     return s >> u


# def mpf_exp(x, prec, rnd=round_fast):
#     sign, man, exp, bc = x
#     if man:
#         mag = bc + exp
#         wp = prec + 14
#         if mag > 1:
#             wpmod = wp + mag
#             offset = exp + wpmod
#             t = man >> (-offset)
#             lg2 = ln2_fixed(wpmod)
#             n, t = divmod(t, lg2)
#             n = int(n)
#             t >>= mag
#         man = exp_basecase(t, wp)
#         return from_man_exp(man, n - wp, prec, rnd)


# def exp_expneg_basecase(x, prec):
#     if prec > EXP_COSH_CUTOFF:
#         cosh, sinh = exponential_series(x, prec, 1)
#         return cosh + sinh, cosh - sinh
#     a = exp_basecase(x, prec)
#     b = (1 << (prec + prec)) // a
#     return a, b


# def mpf_cosh_sinh(x, prec, rnd=round_fast, tanh=0):
#     sign, man, exp, bc = x
#     mag = exp + bc
#     wp = prec + 14
#     if mag > 10:
#         if 3 * (1 << (mag - 1)) > wp:
#             c = s = mpf_shift(mpf_exp(x, prec, rnd), -1)
#             return c, s
#     if mag > 1:
#         wpmod = wp + mag
#         offset = exp + wpmod
#         t = man >> (-offset)
#         lg2 = ln2_fixed(wpmod)
#         n, t = divmod(t, lg2)
#         n = int(n)
#         t >>= mag
#     else:
#         offset = exp + wp
#         t = man >> (-offset)
#         n = 0
#     a, b = exp_expneg_basecase(t, wp)
#     cosh = a + (b >> (2 * n))
#     sinh = a - (b >> (2 * n))
#     cosh = from_man_exp(cosh, n - wp - 1, prec, rnd)
#     sinh = from_man_exp(sinh, n - wp - 1, prec, rnd)
#     return cosh, sinh


# def _d(n, j, prec, sq23pi, sqrt8):
#     j = from_int(j)
#     pi = mpf_pi(prec)
#     a = mpf_div(sq23pi, j, prec)
#     b = mpf_sub(from_int(n), from_rational(1, 24, prec), prec)
#     c = mpf_sqrt(b, prec)
#     ch, sh = mpf_cosh_sinh(mpf_mul(a, c), prec)
#     D = mpf_div(mpf_sqrt(j, prec), mpf_mul(mpf_mul(sqrt8, b), pi), prec)
#     E = mpf_sub(mpf_mul(a, ch), mpf_div(sh, c, prec), prec)
#     return mpf_mul(D, E)



def A068425(n):
    return pi_fixed(i)

for i in range(1,101):
    print(A068425(i))

for i in range(999,1010):
    print(A068425(i))
