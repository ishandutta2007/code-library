/***
 * Complexity: Roughly ~O(n^(2/3))
 * Credits: sgtlaugh
 ***/
// n = 1e+10 : 455055057(time: 0.83s)
// n = 1e+11 : 4118094746(time: 0.85s)
// n = 1e+12 : 37607404117(time: 1.05s)
// n = 1e+13 : 346029867905(time: 2.86s)
#include <bits/stdc++.h>
using namespace std;
using i128 = __int128;

std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}
/// Magic constants, optimized to answer prime counting queries for n=10^13 but
/// can be tweaked
/// Reduce MAXN and MAXM to sacrifice time for memory

const int MAXN = 50;
const int MAXM = 2000010;
const int MAXV = 20000010;

const auto fast_div = [](const uint64_t &a, const uint32_t &b) -> uint64_t {
  return double(a) / b + 1e-9;
};

int pi[MAXV];
uint64_t pi_sum[MAXV], dp[MAXN][MAXM];

vector<int> primes;
bitset<MAXV> is_prime;

inline uint64_t func(int p, int k) {
  if (k == 0)
    return 1;
  else if (k == 1)
    return p;
}

inline i128 accfunc(uint64_t n, int k) {
  if (k == 0)
    return n;
  else if (k == 1)
    return (i128)n * (i128)(n + 1) / 2;
}

void sieve() {
  is_prime[2] = true;
  for (int i = 3; i < MAXV; i += 2)
    is_prime[i] = true;

  for (int i = 3; i * i < MAXV; i += 2) {
    for (int j = i * i; is_prime[i] && j < MAXV; j += (i << 1)) {
      is_prime[j] = false;
    }
  }

  primes.push_back(-1);
  for (int i = 1; i < MAXV; i++) {
    pi[i] = pi[i - 1], pi_sum[i] = pi_sum[i - 1];
    if (is_prime[i]) {
      primes.push_back(i);
      pi[i]++, pi_sum[i] += func(i, 1); // i;
    }
  }
}

void gen() {
  sieve();
  for (int i = 0; i < MAXM; i++)
    dp[0][i] = (uint64_t)accfunc(i, 1); // i * (i + 1) / 2;
  for (int i = 1; i < MAXN; i++) {
    for (int j = 1; j < MAXM; j++) {
      dp[i][j] = dp[i - 1][j] - dp[i - 1][fast_div(j, primes[i])] *
                                    func(primes[i], 1); // primes[i];
    }
  }
}

typedef uint64_t int128;

template <typename T> T phi(T m, int n) {
  if (!n)
    return (T)accfunc(m, 1); // m * (m + 1) / 2;
  if (n < MAXN && m < MAXM)
    return dp[n][m];
  if (m < MAXV && (uint64_t)primes[n] * primes[n] >= m)
    return pi_sum[m] - pi_sum[primes[n]] + 1;
  return phi(m, n - 1) - phi((T)fast_div(m, primes[n]), n - 1) * primes[n];
}

template <typename T> T lehmer(T n, int k) {
  if (n < MAXV)
    return pi_sum[n];

  int s = sqrt(0.5 + n), c = cbrt(0.5 + n);
  T res = phi(n, pi[c]) + pi_sum[c] - 1;

  for (int i = pi[c] + 1; i <= pi[s]; i++) {
    T w = lehmer(fast_div(n, primes[i]), k) - pi_sum[primes[i] - 1];
    res -= w * primes[i];
  }

  return res;
}

__int128 prime_sum(long long n) {
  if (n <= UINT_MAX)
    return lehmer((uint64_t)n, 1);
  return lehmer((__int128)n, 1);
}

__int128 prime_pi(long long n) {
  if (n <= UINT_MAX)
    return lehmer((uint64_t)n, 0);
  return lehmer((__int128)n, 0);
}

int main() {
  auto start = clock();
  gen();
  printf("Pre-process time = %0.3f\n\n",
         (clock() - start) / (double)CLOCKS_PER_SEC); /// 0.940
  start = clock();

  assert(prime_sum(1000) == 76127);
  assert(prime_sum(1e6) == 37550402023);
  assert(prime_sum(1e7) == 3203324994356);
  assert(prime_sum(1e8) == 279209790387276);
  assert(prime_sum(1e9) == 24739512092254535LL);
  assert(prime_sum(1e10) == 2220822432581729238LL);
  assert(prime_sum(UINT_MAX) == 425649736193687430LL);
  assert(prime_sum(1e11) ==
         (__int128)603698 * 333721625289043LL); /// 201467077743744681014
  assert(prime_sum(1e12) ==
         (__int128)15929208151LL * 1157344946327LL); /// 18435588552550705911377
  assert(prime_sum(1e13) ==
         (__int128)10166702 *
             167138413556114797LL); /// 1699246443377779418889494

  printf("\nCalculation time = %0.3f\n",
         (clock() - start) / (double)CLOCKS_PER_SEC); /// 2.130
  return 0;
}
