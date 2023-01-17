// Computes C(n, k) mod p^e
// Constraints:
// - n, k fit within the long long range
// - 1 ≤ k ≤ n
// - p*e ≤ 10^6
// - p primes

#include <iostream>
#include <vector>

using namespace std;

long long power(long long x, long long y, long long m) {
  long long ret = 1;
  for (; y; y >>= 1) {
    if (y & 1)
      ret = (__int128)ret * x % m;
    x = (__int128)x * x % m;
  }

  return ret;
}

namespace Binomial {
long long p = 1, pe = 1;
int e = 0;

// p^k
vector<long long> pPower;

// First kind
vector<vector<long long>> stirling;

// (p-1)!^k
vector<long long> pstirling;

// factorials and its inverse, but without the factor p
vector<long long> fact, ifact;
// p-adic valuation of factorial
vector<int> pfact;

vector<long long> lagrangeCoeff;

vector<long long> prods;

// Computes (np+1)(np+2) ... (np+m) % p^e
// with m < p
long long rising_factorial(long long n, long long m) {
  long long ret = 0;
  long long pn = 1;
  int ep = (e < p ? e : p);
  for (int j = 0; j < ep; ++j) {
    ret = ((__int128)stirling[m + 1][j + 1] * pn + ret) % pe;
    pn = (__int128)pn * p * n % pe;
  }

  return ret;
}

// Compute inverse of n mod p^e
inline long long inv(long long n) { return power(n, pe / p * (p - 1) - 1, pe); }

long long lagrangeInterpolate(long long ndp) {
  int len = e * 2 - 1;
  long long ptot = 0;

  int pfactorsnum[len];
  long long prenum[len], sufnum[len];

  for (int i = 0; i < len; ++i) {
    long long num = ndp - i;
    pfactorsnum[i] = 0;
    while (num % p == 0) {
      num /= p;
      ++pfactorsnum[i];
    }
    ptot += pfactorsnum[i];
    prenum[i] = sufnum[i] = num;

    if (i > 0)
      prenum[i] = (__int128)prenum[i - 1] * prenum[i] % pe;
  }

  for (int i = len - 2; i >= 0; --i) {
    sufnum[i] = (__int128)sufnum[i] * sufnum[i + 1] % pe;
  }

  long long sum = 0;
  for (int j = 0; j < len; ++j) {
    int j2 = len - 1 - j;
    long long pfactor = ptot - pfactorsnum[j] - pfact[j] - pfact[j2];

    if (pfactor >= e)
      continue;

    long long numerator = (__int128)(j > 0 ? prenum[j - 1] : 1) *
                          (j < len - 1 ? sufnum[j + 1] : 1) % pe;

    sum =
        ((__int128)numerator * lagrangeCoeff[j] % pe * pPower[pfactor] + sum) %
        pe;
  }

  return sum;
}

// n! % p^e, but p, 2p, ... is not multiplied
long long factorial_without_p(long long n) {
  long long ndp = n / p;
  long long ret = rising_factorial(ndp, n % p);

  if (ndp <= e * 2 - 2) {
    return (__int128)ret * prods[ndp] % pe;
  }

  ret = (__int128)ret * lagrangeInterpolate(ndp) % pe;
  return ret;
}

// n! but without the p factors
long long factorial(long long n) {
  long long ret = 1;

  while (n > 0) {
    ret = (__int128)ret * factorial_without_p(n) % pe;
    n /= p;
  }

  return ret;
}

// Legendre formula to compute p-adic valuation of n!
long long legendre(long long n) {
  long long ret = 0;

  do {
    ret += n /= p;
  } while (n >= p);

  return ret;
}

// O(p*min(e, p) + e^2)
void setMod(long long _p, int _e) {
  p = _p, e = _e;

  pPower.resize(e);
  pe = 1;
  for (int i = 0; i < e; ++i) {
    pPower[i] = pe;
    pe *= p;
  }

  int ep = (e < p ? e : p);
  stirling = vector<vector<long long>>(p + 1, vector<long long>(ep + 1, 0));
  stirling[0][0] = 1;
  for (int i = 1; i <= p; ++i) {
    for (int j = 1; j <= ep; ++j) {
      stirling[i][j] =
          ((__int128)(i - 1) * stirling[i - 1][j] + stirling[i - 1][j - 1]) %
          pe;
    }
  }

  prods.resize(e * 2 - 1);
  long long prod = 1;
  long long invStirling = inv(stirling[p][1]);

  for (int i = 0; i <= e * 2 - 3; ++i) {
    prods[i] = prod;
    prod = (__int128)prod * rising_factorial(i, p - 1) % pe * invStirling % pe;
  }
  prods[e * 2 - 2] = prod;

  pstirling.resize(e);
  pstirling[0] = 1;
  for (int i = 1; i < e; ++i)
    pstirling[i] = (__int128)pstirling[i - 1] * stirling[p][1] % pe;

  int len = e * 2 - 1;
  fact.resize(len);
  ifact.resize(len);
  pfact.resize(len);
  fact[0] = 1LL;
  pfact[0] = 0;
  for (int i = 1; i < len; ++i) {
    long long num = i;
    pfact[i] = 0;
    while (num % p == 0) {
      ++pfact[i];
      num /= p;
    }

    ifact[i - 1] = num;
    pfact[i] += pfact[i - 1];
    fact[i] = (__int128)fact[i - 1] * num % pe;
  }

  ifact[len - 1] = inv(fact[len - 1]);
  for (int i = len - 2; i >= 0; --i) {
    ifact[i] = (__int128)ifact[i + 1] * ifact[i] % pe;
  }

  lagrangeCoeff.resize(len);
  for (int i = 0; i < len; ++i) {
    long long denominator = (__int128)ifact[i] * ifact[len - 1 - i] % pe;
    if (len - 1 - i & 1)
      denominator = pe - denominator;

    lagrangeCoeff[i] = (__int128)denominator * prods[i] % pe;
  }
}

// O(log(n) * e)
long long binomial_coefficient(long long n, long long k) {
  long long binom_padic = legendre(n) - legendre(k) - legendre(n - k);
  if (binom_padic >= e)
    return 0LL;

  return (__int128)pPower[binom_padic] * factorial(n) % pe *
         inv((__int128)factorial(k) * factorial(n - k) % pe) % pe *
         pstirling[binom_padic] % pe;
}
}

void test() {
  vector<vector<long long>> C(100, vector<long long>(100, 0));
  C[0][0] = 1;

  int mod = 3 * 3 * 3 * 3 * 3;

  Binomial::setMod(3, 5);

  for (int i = 1; i < 100; ++i) {
    C[i][0] = 1;
    for (int j = 1; j <= i; ++j) {
      C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
      assert(C[i][j] == Binomial::binomial_coefficient(i, j));
    }
  }
}

int main() {
  // Compute binom mod 2^60
  Binomial::setMod(2, 60);
  cout << Binomial::binomial_coefficient(1234567890123456789LL,
                                         123456789012345678LL)
       << endl;
  cout << Binomial::binomial_coefficient(1e18, 5e17) << endl;
  // ~3.5 ms (including setMod)

  // Compute binom mod 999979^3
  Binomial::setMod(999979, 3);
  cout << Binomial::binomial_coefficient(1234567890123456789LL,
                                         123456789012345678LL)
       << endl;
  cout << Binomial::binomial_coefficient(1e18, 5e17) << endl;
  // ~270 ms (including setMod)
  return 0;
}
