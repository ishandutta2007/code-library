#include "../../headers/headers.h"
#include "../CRT/CRT.cpp"
#include "../primalityChecks/millerRabin/millerRabin.cpp"
#include "../primalityChecks/sieveEratosthenes/sieve.cpp"

// Modular computation of nCr using lucas theorem, granville theorem and CRT

ll num; // Set num to the corresponding mod for the nCr calculations
umap<ll, int> MOD; // MOD[P]=V_p(mod)
umap<ll, vector<ll>> FMOD; // n! mod p if MOD[p]=1 else the product of all i mod
                           // P^MOD[P], where 1<=i<=n and (i,p)=1
umap<ll, vector<ll>> invFMOD; // the inverse of FMOD[n] in the corresponding MOD

void preCompute() {
  // Factor mod->MOD
  vi primes = sieve(num);
  ll m = num;
  for (auto p : primes) {
    if (p * p > m)
      break;
    while (m % p == 0) {
      MOD[p]++;
      if ((m /= p) == 1)
        goto next;
    }
  }
  if (m > 1)
    MOD[m] = 1;
next:
  // Compute FMOD and invFMOD
  for (auto p : MOD) {
    int m = pow(p.first, p.second); // p^V_p(n)
    FMOD[p.first].assign(m, 1);
    invFMOD[p.first].assign(m, 1);
    repx(i, 2, FMOD[p.first].size()) {
      if (i % p.first == 0 and p.second > 1)
        FMOD[p.first][i] = FMOD[p.first][i - 1];
      else
        FMOD[p.first][i] = mul(FMOD[p.first][i - 1], i, FMOD[p.first].size());

      // Compute using Euler's theorem i.e. a^phi(m)=1 mod m with (a.m)=1
      invFMOD[p.first][i] =
          fastPow(FMOD[p.first][i], m / p.first * (p.first - 1) - 1, m);
    }
  }
}

// Compute nCr using Granville's theorem (prime powers)
// Auxiliary functions

// V_p(n!) using Legendre's theorem
int V(ll n, int p) {
  int e = 0;
  while ((n /= p) > 0)
    e += n;
  return e;
}

//
ll f(ll n, ll p) {
  ll m = pow(p, MOD[p]);
  int e = n / m;
  return mul(fastPow(FMOD[p][m - 1], e, m), FMOD[p][n % m], m);
}
ll F(ll n, ll p) {
  ll m = pow(p, MOD[p]);
  ll ans = 1;
  do {
    ans = mul(ans, f(n, p), m);
  } while ((n /= p) > 0);
  return ans;
}
// Granville theorem
ll granville(ll n, ll r, int p) {
  int e = V(n, p) - V(n - r, p) - V(r, p);
  ll m = pow(p, MOD[p]);
  if (e >= MOD[p])
    return 0;
  ll ans = fastPow(p, e, m);
  ans = mul(ans, F(n, p), m);
  ans = mul(ans, fastPow(F(r, p), pow(p, MOD[p] - 1) * (p - 1) - 1, m), m);
  ans = mul(ans, fastPow(F(n - r, p), pow(p, MOD[p] - 1) * (p - 1) - 1, m), m);
  return ans;
}

// Compute nCr using Lucas theorem (primes)
ll lucas(ll n, ll r, int p) {
  // Trivial cases
  if (r > n or r < 0)
    return 0;
  if (r == 0 or n == r)
    return 1;
  if (r == 1 or r == n - 1)
    return n % p;
  // Base case
  if (n < p and r < p) {
    ll ans = mul(invFMOD[p][r], invFMOD[p][n - r], p); // 1/(r!(n-r)!) mod p
    ans = mul(ans, FMOD[p][n], p);                     // n!/(r!(n-r!)) mod p
    return ans;
  }
  ll ans = lucas(n / p, r / p, p); // Recursion
  ans = mul(ans, lucas(n % p, r % p, p), p); // False recursion
  return ans;
}

// Given the prime decomposition of mod;
ll nCr(ll n, ll r) {
  // Trivial cases
  if (n < r or r < 0)
    return 0;
  if (r == 0 or r == n)
    return 1;
  if (r == 1 or r == n - 1)
    return (n % num);
  // Non-trivial cases
  ll ans = 0;
  ll mod = 1;
  for (auto p : MOD) {
    ll temp = pow(p.first, p.second);
    if (p.second > 1) {
      ans = CRT(ans, mod, granville(n, r, p.first), temp).first;
    } else {
      ans = CRT(ans, mod, lucas(n, r, p.first), temp).first;
    }
    mod *= temp;
  }
  return ans;
}

//https://github.com/ignaciohermosillacornejo/apuntes_icpc/blob/234221acf8ceb6ddb0097fd6f94cea47b2923883/math/nCr/nCr.cpp
//https://github.com/pavel-zeman/CodeChef/blob/a8a6b7638d3af4ed9fe7f64843a86e218882194a/CodeChef/src/y2017/m05/challenge/LongSandwich.cpp
// https://people.math.rochester.edu/faculty/doug/otherpapers/granville.pdf
