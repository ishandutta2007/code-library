/***
 *
 * Prime counting function in sublinear time with the Meissel-Lehmer algorithm
 *
 * The function lehmer(n) returns the number of primes not exceeding n
 * Complexity: Roughly ~O(n^(2/3))
 *
 ***/

// n = 1e+09 : 50847513(time: 0.28s)
// n = 1e+10 : 455052511(time: 0.29s)
// n = 1e+11 : 4118054813(time: 0.29s)
// n = 1e+12 : 37607912018(time: 0.39s)
// n = 1e+13 : 346065536839(time: 1.27s)
#include <bits/stdc++.h>

using namespace std;

/// Magic constants, optimized to answer prime counting queries for n=10^12 but
/// can be tweaked

const int MAXV = 20000010;

const int MAXP = 7;
const int MAXN = 50;
const int MAXM =
    2 * 3 * 7 * 5 * 11 * 13 * 17; /// Product of the first MAXP primes

const auto fast_div = [](const long long &a, const int &b) -> long long {
  return double(a) / b + 1e-9;
};

vector<int> primes;
bitset<MAXV> is_prime;
int prod[MAXP], pi[MAXV], dp[MAXN][MAXM];

void sieve() {
  is_prime[2] = true;
  for (int i = 3; i < MAXV; i += 2)
    is_prime[i] = true;

  for (int i = 3; i * i < MAXV; i += 2) {
    for (int j = i * i; is_prime[i] && j < MAXV; j += (i << 1)) {
      is_prime[j] = false;
    }
  }

  for (int i = 1; i < MAXV; i++) {
    pi[i] = pi[i - 1] + is_prime[i];
    if (is_prime[i])
      primes.push_back(i);
  }
}

void gen() {
  int i, j;
  assert(MAXN >= MAXP);

  sieve();
  for (prod[0] = primes[0], i = 1; i < MAXP; i++) {
    prod[i] = prod[i - 1] * primes[i];
  }

  for (i = 0; i < MAXM; i++)
    dp[0][i] = i;
  for (i = 1; i < MAXN; i++) {
    for (j = 1; j < MAXM; j++) {
      dp[i][j] = dp[i - 1][j] - dp[i - 1][fast_div(j, primes[i - 1])];
    }
  }
}

uint64_t phi(long long m, int n) {
  if (!n)
    return m;
  if (n < MAXN && m < MAXM)
    return dp[n][m];
  if (n < MAXP)
    return dp[n][m % prod[n - 1]] +
           fast_div(m, prod[n - 1]) * dp[n][prod[n - 1]];

  long long p = primes[n - 1];
  if (m < MAXV && p * p >= m)
    return pi[m] - n + 1;
  if (p * p * p < m || m >= MAXV)
    return phi(m, n - 1) - phi(fast_div(m, p), n - 1);

  int lim = pi[(int)sqrt(0.5 + m)];
  uint64_t res = pi[m] - (lim + n - 2) * (lim - n + 1) / 2;
  for (int i = n; i < lim; i++) {
    res += pi[fast_div(m, primes[i])];
  }

  return res;
}

uint64_t lehmer(long long n) {
  if (n < MAXV)
    return pi[n];

  int s = sqrt(0.5 + n), c = cbrt(0.5 + n);
  uint64_t res = phi(n, pi[c]) + pi[c] - 1;
  for (int i = pi[c]; i < pi[s]; i++) {
    res -= lehmer(fast_div(n, primes[i])) - i;
  }

  return res;
}

int main() {
  auto start = clock();
  gen();
  printf("Pre-process time = %0.3f\n\n",
         (clock() - start) / (double)CLOCKS_PER_SEC); /// 0.237

  start = clock();

  assert(lehmer(7) == 4);
  assert(lehmer(1000) == 168);
  assert(lehmer(1000000) == 78498);
  assert(lehmer(1e7) == 664579);
  assert(lehmer(1e8) == 5761455);
  assert(lehmer(1e9) == 50847534);
  assert(lehmer(1e10) == 455052511);
  assert(lehmer(1e11) == 4118054813);
  assert(lehmer(1e12) == 37607912018LL);
  assert(lehmer(1e13) == 346065536839LL);

  printf("\nCalculation time = %0.3f\n",
         (clock() - start) / (double)CLOCKS_PER_SEC); /// 0.997
  return 0;
}