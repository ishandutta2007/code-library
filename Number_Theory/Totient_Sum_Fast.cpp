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

// Least prime factor of each integer. Necessary for linear sieve.
int32_t lpf[D + 1];
// Prime numbers below D. Will be much smaller than D + 1 in actuality.
int32_t primes[D + 1];
// Euler totient function below D.
int32_t phi[D + 1];
// Summatory totient function.
int32_t phi_sum[D + 1];
// large_phi_sum_compressed[i] will be the sum of the euler totient function up to N / i.
int32_t large_phi_sum_compressed[I + 1];

// Multiply 64-bits integers modulo MOD.
ll mul_mod(ll a, ll b) {
  return ((ll)a * (ll)b) % MOD;
}

// nth triangular number modulo MOD.
int32_t T_mod(ll n) {
  return mul_mod(mul_mod(n % MOD, (n + 1) % MOD), INV2);
}

int main() {
  phi[1] = 1;
  int32_t nprimes = 0;
  for (int32_t k = 2; k <= D; ++k) {
    if (lpf[k] == 0) {
      phi[k] = k - 1;
      lpf[k] = k;
      primes[nprimes++] = k;
    }

    for (ll j = 0;
         primes[j] < lpf[k] && (ll)primes[j] * (ll)k <= (ll)D;
         ++j) {
      lpf[primes[j] * k] = primes[j];
      phi[primes[j] * k] = (primes[j] - 1) * phi[k];
    }

    if ((ll)lpf[k] * (ll)k <= (ll)D) {
      lpf[lpf[k] * k] = lpf[k];
      phi[lpf[k] * k] = lpf[k] * phi[k];
    }
  }
  for (int32_t d = 1; d <= D; ++d) {
    phi_sum[d] = (phi_sum[d - 1] + phi[d]) % MOD;
    // cout << "etf_sum(" << d << ")=" << phi_sum[d] << endl;
  }

  large_phi_sum_compressed[I] = phi_sum[D];

  // Use recurrence phi_sum(N) = T(N) - sum_d >=2 phi_sum(N/d) to compute large_phi_sum_compressed
  // values.
  cout << "Now computing all larges, amortised O(N^(2/3)) for all query or O(1) per query" << endl;
  for (ll i = I - 1; i > 0; --i) {
    const ll x_i = N / i;
    large_phi_sum_compressed[i] = T_mod(x_i);

    ll d = 2;
    for (; d * i < I; ++d) {
      large_phi_sum_compressed[i] = (large_phi_sum_compressed[i] - large_phi_sum_compressed[d * i] + MOD) % MOD;
    }

    for (; d * d <= x_i; ++d) {
      large_phi_sum_compressed[i] = (large_phi_sum_compressed[i] - phi_sum[x_i / d] + MOD) % MOD;
    }
    for (ll k = x_i / d; k > 0; --k) {
      ll m = mul_mod((x_i / k - x_i / (k + 1)) % MOD, phi_sum[k]);
      large_phi_sum_compressed[i] = (large_phi_sum_compressed[i] - m + MOD) % MOD;
    }
    // cout << "etf_sum(" << N << "/" << i << "=" << N/i << ")=" << large_phi_sum_compressed[i] <<
    // endl;  // == etf_sum(N/i)
  }

  ll ii = 3456;
  if (ii < I) cout << "etf_sum(" << N << "/" << ii << "=" << N / ii << ")=" << large_phi_sum_compressed[ii]  << endl;
  else cout << "etf_sum(" << N << "/" << ii << "=" << N / ii << ")=" << phi_sum[N / ii] << endl;
  cout << clock() / (double)CLOCKS_PER_SEC << endl;
}
