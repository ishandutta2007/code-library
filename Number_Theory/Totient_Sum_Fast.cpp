// This can be achieved in time complexity O(N^(2/3)) by presieving values of
// phi. We use that sum of PHI(d) floor(N/d)^2 over 1 <= d <= N is the same as
// what we need except it counts each pair where i != j twice (but we can easily
// correct for this).

// O(N^(2/3)) precompute O(1) per query
#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// Prime modulus for this problem
#define MOD 998244353
// Inverse of 2 modulo MOD
#define INV2 499122177
#define N 100000000000 // 10^11
// We must have D * I <= N < (D + 1) * I.
// Optimal performance is achieved with I approximately N^(1/3).
#define I 4641
#define D 21547080

// Lowest prime divisor of each integer. Necessary for linear sieve.
int32_t LP[D + 1];
// Prime numbers below D. Will be much smaller than D + 1 in actuality.
int32_t PRIMES[D + 1];
// Euler totient function below D.
int32_t PHI[D + 1];
// Summatory totient function.
int32_t SUM_PHI[D + 1];
// S[i] will be the sum of the euler totient function up to N / i.
int32_t S[I + 1];

// Multiply 64-bits integers modulo MOD.
int64_t mul_mod(int64_t a, int64_t b) {
  return ((int64_t)a * (int64_t)b) % MOD;
}

// nth triangular number modulo MOD.
int32_t T_mod(int64_t n) {
  return mul_mod(mul_mod(n % MOD, (n + 1) % MOD), INV2);
}

int main() {
  PHI[1] = 1;
  int32_t nprimes = 0;
  for (int32_t k = 2; k <= D; ++k) {
    if (LP[k] == 0) {
      PHI[k] = k - 1;
      LP[k] = k;
      PRIMES[nprimes++] = k;
    }

    for (int64_t j = 0;
         PRIMES[j] < LP[k] && (int64_t)PRIMES[j] * (int64_t)k <= (int64_t)D;
         ++j) {
      LP[PRIMES[j] * k] = PRIMES[j];
      PHI[PRIMES[j] * k] = (PRIMES[j] - 1) * PHI[k];
    }

    if ((int64_t)LP[k] * (int64_t)k <= (int64_t)D) {
      LP[LP[k] * k] = LP[k];
      PHI[LP[k] * k] = LP[k] * PHI[k];
    }
  }
  for (int32_t d = 1; d <= D; ++d) {
    SUM_PHI[d] = (SUM_PHI[d - 1] + PHI[d]) % MOD;
    // cout << "etf_sum(" << d << ")=" << SUM_PHI[d] << endl;  // == etf(N/i)
  }

  S[I] = SUM_PHI[D];

  // Use recurrence SUM_PHI(N) = T(N) - sum_d >=2 SUM_PHI(N/d) to compute S
  // values.
  for (int64_t i = I - 1; i > 0; --i) {
    const int64_t x_i = N / i;
    S[i] = T_mod(x_i);

    int64_t d = 2;
    for (; d * i < I; ++d) {
      S[i] = (S[i] - S[d * i] + MOD) % MOD;
    }

    for (; d * d <= x_i; ++d) {
      S[i] = (S[i] - SUM_PHI[x_i / d] + MOD) % MOD;
    }
    for (int64_t k = x_i / d; k > 0; --k) {
      ll m = mul_mod((x_i / k - x_i / (k + 1)) % MOD, SUM_PHI[k]);
      S[i] = (S[i] - m + MOD) % MOD;
    }
    // cout << "etf_sum(" << N << "/" << i << "=" << N/i << ")=" << S[i] <<
    // endl;  // == etf(N/i)
  }

  ll ii = 3456;
  if (ii < I) cout << "etf_sum(" << N << "/" << ii << "=" << N / ii << ")=" << S[ii]  << endl;
  else cout << "etf_sum(" << N << "/" << ii << "=" << N / ii << ")=" << SUM_PHI[N / ii] << endl;
  cout << clock() / (double)CLOCKS_PER_SEC << endl;
}
