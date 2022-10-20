#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>

#define UINT128_C(u) ((__uint128_t)u)
const __uint128_t UINT128_MAX = UINT128_C(UINT64_MAX) << 64 | UINT64_MAX;

void pu64(__uint64_t u) { printf("%" PRIu64, u); }
void pu640(__uint64_t u) { printf("%019" PRIu64, u); }

#define D19_ UINT64_C(10000000000000000000)
const __uint128_t d19_ = D19_;
const __uint128_t d38_ = UINT128_C(D19_) * D19_;

void pu128(__uint128_t u) {
  if (u < d19_)
    pu64(u);
  else if (u < d38_) {
    pu64(u / d19_);
    pu640(u % d19_);
  } else {
    pu64(u / d38_);
    u %= d38_;
    pu640(u / d19_);
    pu640(u % d19_);
  }
}

// from https://stackoverflow.com/a/45637368/5674289
#include <errno.h>
typedef __uint128_t u128;

static int strdigit__(char c) {
  /* This is ASCII / UTF-8 specific, would not work for EBCDIC */
  return (c >= '0' && c <= '9') ? c - '0' : (c >= 'a' && c <= 'z')
                                                ? c - 'a' + 10
                                                : (c >= 'A' && c <= 'Z')
                                                      ? c - 'A' + 10
                                                      : 255;
}

static u128 strtou128__(const char *p, char **endp, int base) {
  u128 v = 0;
  int digit;

  if (base == 0) { /* handle octal and hexadecimal syntax */
    base = 10;
    if (*p == '0') {
      base = 8;
      if ((p[1] == 'x' || p[1] == 'X') && strdigit__(p[2]) < 16) {
        p += 2;
        base = 16;
      }
    }
  }
  if (base < 2 || base > 36) {
    errno = EINVAL;
  } else if ((digit = strdigit__(*p)) < base) {
    v = digit;
    /* convert to unsigned 128 bit with overflow control */
    while ((digit = strdigit__(*++p)) < base) {
      u128 v0 = v;
      v = v * base + digit;
      if (v < v0) {
        v = ~(u128)0;
        errno = ERANGE;
      }
    }
    if (endp) {
      *endp = (char *)p;
    }
  }
  return v;
}
// from https://stackoverflow.com/a/45637368/5674289

__uint128_t set_128(unsigned long h, unsigned long l) {
  __uint128_t r = h;
  return (r << 64) + l;
}

typedef __uint128_t uint256_t[2];

void set_256(uint256_t d, __uint128_t l, __uint128_t h) {
  d[0] = l;
  d[1] = h;
}
void add_128(uint256_t d, uint256_t x, __uint128_t a) {
  d[0] = x[0] + a;
  d[1] = (d[0] < a) ? x[1] + 1 : x[1];
}
void add_256(uint256_t d, uint256_t x, uint256_t a) {
  add_128(d, x, a[0]);
  d[1] += a[1];
}

void shl_256(uint256_t d, long s) {
  d[1] <<= s;
  d[1] += (d[0] >> (128 - s));
  d[0] <<= s;
}

void sqr_128(uint256_t d, __uint128_t x) {
  if (x <= UINT64_MAX) {
    d[0] = x * x;
    d[1] = 0;
  } else {
    __uint128_t l = x & UINT64_MAX;
    __uint128_t h = x >> 64;

    uint256_t yz, w;
    set_256(d, l * l, 0);
    set_256(yz, l * h, 0);
    shl_256(yz, 64 + 1);
    set_256(w, 0, h * h);
    add_256(d, d, yz);
    add_256(d, d, w);
  }
}

void q(long l) { printf("%016lx ", l); }
void p(long l) { printf("%lx", l); }
void P(__uint128_t i) {
  if (i > UINT64_MAX)
    p(i >> 64);
  q(i & 0xFFFFFFFFFFFFFFFF);
}
void _P_256(uint256_t x) {
  if (x[1] > 0)
    P(x[1]);
  P(x[0]);
}
#define P_256(l)                                                               \
  {                                                                            \
    printf("%s: ", #l);                                                        \
    _P_256(l);                                                                 \
  }

#define min(x, y) ((x < y) ? (x) : (y))

int clz(__uint128_t u) {
  unsigned long long h = u >> 64;
  return (h != 0) ? __builtin_clzll(h) : 64 + __builtin_clzll(u);
}

__uint128_t mod_256(uint256_t x, __uint128_t n) {
  if (x[1] == 0)
    return x[0] % n;
  else {
    __uint128_t r = x[1] % n;
    int F = clz(n);
    int R = clz(r);
    for (int i = 0; i < 128; ++i) {
      if (R > F + 1) {
        int h = min(R - (F + 1), 128 - i);
        r <<= h;
        R -= h;
        i += (h - 1);
        continue;
      }
      r <<= 1;
      if (r >= n) {
        r -= n;
        R = clz(r);
      }
    }
    r += (x[0] % n);
    if (r >= n)
      r -= n;

    return r;
  }
}

__uint128_t a, n;
__uint128_t f(__uint128_t x) {
  uint256_t p, s;
  sqr_128(p, x);
  add_128(s, p, a);
  return mod_256(s, n);
}

inline int ctz(__uint128_t u) {
  unsigned long long h = u;
  return (h != 0) ? __builtin_ctzll(h) : 64 + __builtin_ctzll(u >> 64);
}

unsigned long long _gcd(unsigned long long u, unsigned long long v) {
  for (;;) {
    if (u > v) {
      unsigned long long a = u;
      u = v;
      v = a;
    }

    v -= u;

    if (v == 0)
      return u;

    v >>= __builtin_ctzll(v);
  }
}

__uint128_t gcd(__uint128_t u, __uint128_t v) {
  if (u == 0) {
    return v;
  } else if (v == 0) {
    return u;
  }

  int i = ctz(u);
  u >>= i;
  int j = ctz(v);
  v >>= j;
  int k = (i < j) ? i : j;

  for (;;) {
    if (u > v) {
      __uint128_t a = u;
      u = v;
      v = a;
    }

    if (v <= UINT64_MAX)
      return _gcd(u, v) << k;

    v -= u;

    if (v == 0)
      return u << k;

    v >>= ctz(v);
  }
}

#define dif(x, y) ((x) < (y) ? (y) - (x) : (x) - (y))

long c;
__uint128_t pollard_rho_1st(__uint128_t n) {
  __uint128_t g, x, y;
  c = 0;
  for (;; ++a) {
    x = 2;
    y = 2;
    g = 1;
    for (; g == 1;) {
      //      printf("rho("); P(n); printf(","); P(a); printf(")\n");
      ++c;
      x = f(x);
      y = f(f(y));
      g = gcd(n, dif(x, y));
      //      P(x); printf(" "); P(y); printf(" "); P(g); printf("\n");
    }
    //    printf("found: "); P(g); printf("\n");
    if (g < n)
      break;
  }

  return g;
}

int main(int argc, char *argv[]) {
  struct timespec t0, t1;
  __uint128_t g, p, q;

  assert(argc == 3);
  p = strtou128__(argv[1], NULL, 0);
  q = strtou128__(argv[2], NULL, 0);
  n = p * q;
  a = set_128(0, 1);

  clock_gettime(CLOCK_REALTIME, &t0);
  g = pollard_rho_1st(n);
  clock_gettime(CLOCK_REALTIME, &t1);

  pu128(n);
  printf(": ");
  pu128(g);
  printf(" ");
  pu128(n / g);
  printf("\n [%ld] %ldns\n", c,
         (t1.tv_sec - t0.tv_sec) * 1000000000 + t1.tv_nsec - t0.tv_nsec);

  return 0;
}