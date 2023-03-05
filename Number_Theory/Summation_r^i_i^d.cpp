#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
using namespace std;
typedef int i32;
typedef unsigned u32;
typedef long long i64;
typedef unsigned long long u64;
// wheel-sieve
namespace prime_sieve {
const int max_primes = 160000;
const int sieve_span = 1 << 22;
const int sieve_words = sieve_span >> 7;
const int wheel_size = 3 * 5 * 7 * 11 * 13;

u64 mask[64]; // mask[i] = 1ull << i;
bool built_mask;
void build_mask() {
  if (built_mask)
    return;
  built_mask = true;
  mask[0] = 1ull;
  for (int i = 1; i < 64; ++i)
    mask[i] = mask[i - 1] << 1ull;
}

int primes[max_primes], mcnt;
u64 sieve[sieve_words];
u64 pattern[wheel_size];
void mark(u64 *s, int o) { s[o >> 6] |= mask[o & 63]; }
void unmark(u64 *s, int o) { s[o >> 6] &= ~mask[o & 63]; }
int test(u64 *s, int o) { return (s[o >> 6] & mask[o & 63]) == 0; }

int all[6000000], pcnt; // pi(10^8) = 5761455

bool inside(int k) { return k < pcnt; }
int get_prime(int k) { return all[k < pcnt ? k : pcnt - 1]; }

bool pre_sieved;

void pre_sieve() {
  if (pre_sieved)
    return;
  pre_sieved = true;
  for (int i = 0; i < (1048576 >> 7); ++i)
    sieve[i] = 0;
  for (int i = 3; i < 1024; i += 2)
    if (test(sieve, i >> 1))
      for (int j = (i * i) >> 1; j < 1048576; j += i)
        mark(sieve, j);
  mcnt = 0;
  for (int i = 8; i < 1048576; ++i)
    if (test(sieve, i))
      primes[mcnt++] = (i << 1) + 1;
  // printf("m is %d\n", m);
  for (int i = 0; i < wheel_size; ++i)
    pattern[i] = 0;
  for (int i = 1; i < wheel_size * 64; i += 3)
    mark(pattern, i);
  for (int i = 2; i < wheel_size * 64; i += 5)
    mark(pattern, i);
  for (int i = 3; i < wheel_size * 64; i += 7)
    mark(pattern, i);
  for (int i = 5; i < wheel_size * 64; i += 11)
    mark(pattern, i);
  for (int i = 6; i < wheel_size * 64; i += 13)
    mark(pattern, i);
}

void update_sieve(int base) {
  int o = base % wheel_size;
  o = (o + ((o * 105) & 127) * wheel_size) >> 7;
  for (int i = 0, k; i < sieve_words; i += k, o = 0) {
    k = min(wheel_size - o, sieve_words - i);
    memcpy(sieve + i, pattern + o, sizeof(*pattern) * k);
  }
  if (base == 0) { // mark 1 as not prime, and mark 3, 5, 7, 11, and 13 as prime
    sieve[0] |= mask[0];
    sieve[0] &= ~(mask[1] | mask[2] | mask[3] | mask[5] | mask[6]);
  }
  for (int i = 0; i < mcnt; ++i) {
    i64 j = primes[i] * primes[i];
    if (j > base + sieve_span - 1)
      break;
    if (j > base)
      j = (j - base) >> 1;
    else {
      j = primes[i] - base % primes[i];
      if ((j & 1) == 0)
        j += primes[i];
      j >>= 1;
    }
    while (j < sieve_span >> 1) {
      mark(sieve, j);
      j += primes[i];
    }
  }
}

// sieve [base, min(base+span, lim)]
void segment_sieve_with_lower(int base, int lower, int upper) {
  update_sieve(base);
  int u = min(base + sieve_span, upper);
  for (int i = 0; i < sieve_words; ++i) {
    u64 o = ~sieve[i];
    while (o) {
      int p = __builtin_ctzll(o);
      i64 u = base + (i << 7) + (p << 1) + 1;
      if (u >= upper)
        break;
      if (u >= lower)
        all[pcnt++] = u;
      o -= o & ((~o) + 1);
    }
  }
}

void segment_sieve(int base, int upper) {
  update_sieve(base);
  int u = min(base + sieve_span, upper);
  for (int i = 0; i < sieve_words; ++i) {
    u64 o = ~sieve[i];
    while (o) {
      int p = __builtin_ctzll(o);
      i64 u = base + (i << 7) + (p << 1) + 1;
      if (u >= upper)
        break;
      all[pcnt++] = u;
      o -= o & ((~o) + 1);
    }
  }
}

// sieve from [0, n)
void fast_sieve(int lim) {
  build_mask();
  pre_sieve();
  pcnt = 0;
  all[pcnt++] = 2;
  int i = 0, now = 1;
  for (int base = 0; base < lim; base += sieve_span)
    segment_sieve(base, lim);
}
// sieve from [lo, hi)
void segmented_sieve(int lo, int hi) {
  build_mask();
  pre_sieve();
  pcnt = 0;
  if (lo <= 2)
    all[pcnt++] = 2;
  int base = (lo / sieve_span) * sieve_span;
  segment_sieve_with_lower(base, lo, hi), base += sieve_span;
  while (base < hi)
    segment_sieve(base, hi), base += sieve_span;
}
} // namespace prime_sieve

const int mod = 998244353;
inline int trim(i64 x) {
  return (x >= 0) ? (x % mod) : (mod - (((-x) % mod) ? ((-x) % mod) : mod));
}
inline int trim_pow(i64 x) {
  return x % (mod - 1);
} // assume x > 0, phi(mod) = mod - 1
inline int add(int a, int b) {
  return (a >= mod - b) ? (a - mod + b) : (a + b);
}
inline int sub(int a, int b) { return (a >= b) ? (a - b) : (a + mod - b); }
inline int neg(int a) { return sub(0, a); }
inline int mul(int a, int b) { return (1ll * a * b) % mod; }
inline int quick_pow(int x, i64 p) {
  p = trim_pow(p), x = trim(x);
  int ans = 1;
  while (p) {
    if (p & 1)
      ans = mul(ans, x);
    x = mul(x, x);
    p >>= 1;
  }
  return ans;
}
inline int inv(int x) { return quick_pow(x, mod - 2); }
inline int divi(int a, int b) { return mul(a, inv(b)); }

const int lim = 100000000;
const int maxn = lim + 100;
// save memory by the cost of time.
namespace Binomial {
int invfac_[maxn];
int fac_sz;
void init(int n) {
  if (fac_sz >= n + 10)
    return;
  fac_sz = n + 10;
  n += 9;
  invfac_[0] = 1;
  for (int i = 1; i <= n; ++i)
    invfac_[i] = mul(invfac_[i - 1], i);
  invfac_[n] = inv(invfac_[n]);
  for (int i = n - 1; i > 0; i--)
    invfac_[i] = mul(invfac_[i + 1], i + 1);
}
inline int fac(int i) { return inv(invfac_[i]); }
inline int finv(int i) { return invfac_[i]; }
inline int inv(int i) { return mul(fac(i - 1), invfac_[i]); }
} // namespace Binomial
// get rid of g, pd
int f_sz;
int f[maxn];
int dp[maxn];

// given y(x=0)...y(k) , return y(x)  deg(y) <= k+1
int lagrange_interpolation(const int *y, int y_sz, i64 x) {
  int N = (int)y_sz - 1;
  if (x <= N)
    return y[x];
  int ret = 0;
  dp[0] = 1;
  int a = trim(x);
  for (int i = 0; i < N; i++)
    dp[i + 1] = mul(dp[i], a), a = sub(a, 1);
  int tmp = a;
  for (int i = N; i > 0; i--)
    dp[i - 1] = mul(dp[i - 1], tmp), a = add(a, 1), tmp = mul(tmp, a);
  for (int i = 0; i <= N; i++) {
    int tmp =
        mul(mul(y[i], dp[i]), mul(Binomial::finv(i), Binomial::finv(N - i)));
    ret = ((N - i) & 1) ? sub(ret, tmp) : add(ret, tmp);
  }
  return ret;
}
// given f(0)...f(k) (deg(f) = k)
// return \sum_{i=0...n-1} a^i f(i)
// save memory, so g is replaced by f (f is modified here)
int sum_of_exp(int *f, int f_sz, int fac_f_sz, int a, i64 n) {
  if (n == 0)
    return 0;
  if (a == 0)
    return f[0];
  if (a == 1) {
    int tmp = f[0], tmp2 = 0; // g[0] = 0; g[i] = \sum_{k=0}^{i-1} f[k]
    f[0] = 0;
    for (int i = 1; i < (int)f_sz + 1; i++)
      tmp2 = f[i], f[i] = add(f[i - 1], tmp), tmp = tmp2;
    // g[i] = add(g[i - 1], f[i - 1]);
    return lagrange_interpolation(f, f_sz + 1, n);
  }
  int K = f_sz - 1;
  int buf = 1;
  for (int i = 0; i < (int)f_sz; i++)
    f[i] = mul(f[i], buf), buf = mul(buf, a); // g[i] = mul(f[i], buf)
  for (int i = 1; i < (int)f_sz; i++)
    f[i] = add(f[i], f[i - 1]); // g[i] = add(g[i], g[i - 1])
  int c = 0, buf2 = 1;
  // C(K + 1, i) = mul(fac[K + 1], mul(invfac[K + 1 - i], invfac[i]))
  for (int i = 0; i <= K; i++) // g[K - i]
    c = add(
        c, mul(mul(fac_f_sz, mul(Binomial::finv(K + 1 - i), Binomial::finv(i))),
               mul(buf2, f[K - i]))),
    buf2 = mul(buf2, neg(a));

  c = divi(c, quick_pow(add(neg(a), 1), K + 1));

  int buf3 = 1, ia = inv(a);
  for (int i = 0; i < (int)f_sz; i++)
    f[i] = mul(sub(f[i], c), buf3), buf3 = mul(buf3, ia);
  int tn = lagrange_interpolation(f, f_sz, n - 1);
  return add(mul(tn, quick_pow(a, n - 1)), c);
}
// given f(0)...f(k) (deg(f) = k)
// return \sum_{i=0...infty} a^i f(i)  f(i) \in (-1, 1), modulo prime.
// save memory, so g is replaced by f (f is modified here)
int sum_of_exp_limit(int *f, int f_sz, int fac_f_sz, int a) {
  if (a == 0)
    return f[0];
  int K = f_sz - 1;
  int buf = 1;
  for (int i = 0; i < (int)f_sz; i++)
    f[i] = mul(f[i], buf), buf = mul(buf, a); // g[i] = mul(f[i], buf)
  for (int i = 1; i < (int)f_sz; i++)
    f[i] = add(f[i], f[i - 1]); // g[i] = add(g[i], g[i - 1])
  int c = 0, buf2 = 1;
  // C(K + 1, i) = mul(fac[K + 1], mul(invfac[K + 1 - i], invfac[i]))
  for (int i = 0; i <= K; i++) // g[K - i]
    c = add(
        c, mul(mul(fac_f_sz, mul(Binomial::finv(K + 1 - i), Binomial::finv(i))),
               mul(buf2, f[K - i]))),
    buf2 = mul(buf2, neg(a));
  c = divi(c, quick_pow(sub(1, a), K + 1));
  return c;
}
void exp_enamurate(int *f, int &f_sz, int p, int n) {
  f_sz = n + 1;
  memset(f, 0, f_sz * sizeof(int));
  if (!p) {
    f[0] = 1;
    return;
  }
  f[1] = 1;
  for (int i = 0; prime_sieve::inside(i) && prime_sieve::get_prime(i) <= n; ++i)
    f[prime_sieve::get_prime(i)] = quick_pow(prime_sieve::get_prime(i), p);
  for (int i = 2; i <= n; ++i)
    for (int j = 0;
         prime_sieve::inside(j) && (prime_sieve::get_prime(j) <= i) &&
         (1ll * i * prime_sieve::get_prime(j) <= n);
         ++j) {
      f[i * prime_sieve::get_prime(j)] =
          mul(f[i], f[prime_sieve::get_prime(j)]);
      if (!(i % prime_sieve::get_prime(j)))
        break;
    }
}
// \sum_{i=0}^{n-1}(r^i)*(i^d)
int sum_of_exp2(int d, int r, long long n) {
  exp_enamurate(f, f_sz, d, d);
  return sum_of_exp(f, f_sz, Binomial::fac(f_sz), r, n);
}
// \sum_{i=0}^{+\infty}(r^i)*(i^d)  r \in (-1, 1), modulo prime.
int sum_of_exp_limit2(int d, int r) {
  exp_enamurate(f, f_sz, d, d);
  return sum_of_exp_limit(f, f_sz, Binomial::fac(f_sz), r);
}
long long r, d, n;
// #define FILE_TEST
// #define TIME_TEST
int main() {
#ifdef FILE_TEST
  freopen("MOON4-testfile-input.txt", "r", stdin);
// freopen("MOON4-testfile-output.txt", "w", stdout);
#endif

#ifdef TIME_TEST
  clock_t start, finish;
  start = clock();
#endif
  prime_sieve::segmented_sieve(1, 10000000);
  Binomial::init(10000000);
  scanf("%lld%lld%lld", &r, &d, &n);
  printf("%d", sum_of_exp2(d, trim(r), n));
#ifdef TIME_TEST
  finish = clock();
  fprintf(stderr, "time : %.3f s\n",
          ((double)(finish - start)) / CLOCKS_PER_SEC);
#endif
}
// CYMario
// https://judge.yosupo.jp/problem/sum_of_exponential_times_polynomial
