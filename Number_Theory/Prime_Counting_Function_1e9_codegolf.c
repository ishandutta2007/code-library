#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//Source https://codegolf.stackexchange.com/a/74372
// Meissel-Lehmer - Fast till 32 bit breaks beyond

#define cache_size 16384
#define Phi_prec_max (47 * a)

#define bit(k) (1ULL << ((k) & 63))
#define word(k) sieve[(k) >> 6]
#define sbit(k) ((word(k >> 1) >> (k >> 1)) & 1)
#define ones(k) (~0ULL >> (64 - (k)))
#define m2(k) ((k + 1) / 2)
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define ns(t) (1000000000 * t.tv_sec + t.tv_nsec)
#define popcnt __builtin_popcountll

#define mask_build(i, p, o, m) mask |= m << i, i += o, i -= p * (i >= p)
#define Phi_prec_bytes ((m2(Phi_prec_max) + 1) * sizeof(int16_t))
#define Phi_prec(i, j) Phi_prec_pointer[(j) * (m2(Phi_prec_max) + 1) + (i)]
#define Phi_6_next ((i / 1155) * 480 + Phi_5[i % 1155] - Phi_5[(i + 6) / 13])
#define Phi_6_upd_1() t = Phi_6_next, i += 1, *(l++) = t
#define Phi_6_upd_2() t = Phi_6_next, i += 2, *(l++) = t, *(l++) = t
#define Phi_6_upd_3() t = Phi_6_next, i += 3, *(l++) = t, *(l++) = t, *(l++) = t

typedef unsigned __int128 uint128_t;
struct timespec then, now;
uint64_t a, primes[4648] = { 2, 3, 5, 7, 11, 13, 17, 19 }, *primes_fastdiv;
uint16_t *Phi_6, *Phi_prec_pointer;

static inline uint64_t Phi_6_mod(uint64_t y)
{
    if (y < 30030)
        return Phi_6[m2(y)];
    else
        return (y / 30030) * 5760 + Phi_6[m2(y % 30030)];
}

static inline uint64_t fastdiv(uint64_t dividend, uint64_t fast_divisor)
{
    return ((uint128_t) dividend * fast_divisor) >> 64;
}

uint64_t Phi(uint64_t y, uint64_t c)
{
    uint64_t *d = primes_fastdiv, i = 0, r = Phi_6_mod(y), t = y / 17;

    r -= Phi_6_mod(t), t = y / 19;

    while (i < c && t > Phi_prec_max) r -= Phi(t, i++), t = fastdiv(y, *(d++));

    while (i < c && t) r -= Phi_prec(m2(t), i++), t = fastdiv(y, *(d++));

    return r;
}

uint64_t Phi_small(uint64_t y, uint64_t c)
{
    if (!c--) return y;

    return Phi_small(y, c) - Phi_small(y / primes[c], c);
}

uint64_t pi_small(uint64_t y)
{
    uint64_t i, r = 0;

    for (i = 0; i < 8; i++) r += (primes[i] <= y);

    for (i = 21; i <= y; i += 2)
        r += i % 3 && i % 5 && i % 7 && i % 11 && i % 13 && i % 17 && i % 19;

    return r;
}

int output(int result)
{
    clock_gettime(CLOCK_REALTIME, &now);
    printf("pi(x) = %9d    real time:%9ld ns\n", result , ns(now) - ns(then));

    return 0;
}

uint64_t Lehmer(uint64_t x)
{
    uint64_t b, i, j, k, limit, mask, P2, *p, start, t = 8;
    uint64_t root2 = sqrt(x), root3 = pow(x, 1./3), top = x / root3 + 1;
    uint64_t halftop = m2(top), *sieve, sieve_length = (halftop + 63) / 64;
    uint64_t i3 = 1, i5 = 2, i7 = 3, i11 = 5, i13 = 6, i17 = 8, i19 = 9;
    uint16_t Phi_3[] = { 0, 1, 1, 1, 2, 2, 3, 4, 4, 5, 6, 6, 7, 7, 7, 8 };
    uint16_t *l, *m, Phi_4[106], Phi_5[1156];

    clock_gettime(CLOCK_REALTIME, &then);

    sieve = malloc(sieve_length * sizeof(int64_t));

    if (x < 529) return output(pi_small(x));

    for (i = 0; i < sieve_length; i++)
    {
        mask  = 0;

        mask_build( i3,  3,  2, 0x9249249249249249ULL);
        mask_build( i5,  5,  1, 0x1084210842108421ULL);
        mask_build( i7,  7,  6, 0x8102040810204081ULL);
        mask_build(i11, 11,  2, 0x0080100200400801ULL);
        mask_build(i13, 13,  1, 0x0010008004002001ULL);
        mask_build(i17, 17,  4, 0x0008000400020001ULL);
        mask_build(i19, 19, 12, 0x0200004000080001ULL);

        sieve[i] = ~mask;
    }

    limit = min(halftop, 8 * cache_size);

    for (i = 21; i < root3 + 1; i += 2)
        if (sbit(i))
            for (primes[t++] = i, j = i * i / 2; j < limit; j += i)
                word(j) &= ~bit(j);

    a = t;

    for (i = (root3 + 1) | 1; i < root2 + 1; i += 2)
        if (sbit(i)) primes[t++] = i;

    b = t;

    while (limit < halftop)
    {
        start = 2 * limit + 1, limit = min(halftop, limit + 8 * cache_size);

        for (p = &primes[8]; p < &primes[a]; p++)
            for (j = max(start / *p | 1, *p) * *p / 2; j < limit; j += *p)
                word(j) &= ~bit(j);
    }

    P2 = (a - b) * (a + b - 1) / 2;

    for (i = m2(root2); b --> a; P2 += t, i = limit)
    {
        limit = m2(x / primes[b]), j = limit & ~63;

        if (i < j)
        {
            t += popcnt((word(i)) >> (i & 63)), i = (i | 63) + 1;

            while (i < j) t += popcnt(word(i)), i += 64;

            if (i < limit) t += popcnt(word(i) & ones(limit - i));
        }
        else if (i < limit) t += popcnt((word(i) >> (i & 63)) & ones(limit - i));
    }

    if (a < 7) return output(Phi_small(x, a) + a - 1 - P2);

    a -= 7, Phi_6 = malloc(a * Phi_prec_bytes + 15016 * sizeof(int16_t));
    Phi_prec_pointer = &Phi_6[15016];

    for (i = 0; i <= 105; i++)
        Phi_4[i] = (i / 15) * 8 + Phi_3[i % 15] - Phi_3[(i + 3) / 7];

    for (i = 0; i <= 1155; i++)
        Phi_5[i] = (i / 105) * 48 + Phi_4[i % 105] - Phi_4[(i + 5) / 11];

    for (i = 1, l = Phi_6, *l++ = 0; i <= 15015; )
    {
        Phi_6_upd_3(); Phi_6_upd_2(); Phi_6_upd_1(); Phi_6_upd_2();
        Phi_6_upd_1(); Phi_6_upd_2(); Phi_6_upd_3(); Phi_6_upd_1();
    }

    for (i = 0; i <= m2(Phi_prec_max); i++)
        Phi_prec(i, 0) = Phi_6[i] - Phi_6[(i + 8) / 17];

    for (j = 1, p = &primes[7]; j < a; j++, p++)
    {
        i = 1, memcpy(&Phi_prec(0, j), &Phi_prec(0, j - 1), Phi_prec_bytes);
        l = &Phi_prec(*p / 2 + 1, j), m = &Phi_prec(m2(Phi_prec_max), j) - *p;

        while (l <= m)
            for (k = 0, t = Phi_prec(i++, j - 1); k < *p; k++) *(l++) -= t;

        t = Phi_prec(i++, j - 1);

        while (l <= m + *p) *(l++) -= t;
    }

    primes_fastdiv = malloc(a * sizeof(int64_t));

    for (i = 0, p = &primes[8]; i < a; i++, p++)
    {
        t = 96 - __builtin_clzll(*p);
        primes_fastdiv[i] = (bit(t) / *p + 1) << (64 - t);
    }

    return output(Phi(x, a) + a + 6 - P2);
}

int main() {
  uint64_t n =20*1e8;
    uint64_t res = Lehmer(n);
    printf("%d = result: %lld\n", n, res);
  return 0;
}