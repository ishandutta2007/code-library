#include "bits/stdc++.h"
using namespace std;
using ll = long long;
using int64 = long long;
using uint64 = long long;
#define MOD 1000000000000000009
// n = 1e+10(time: 0.336641s)
// n = 1e+11(time: 1.64541s)
// n = 1e+12(time: 8.17151s)

// mod should be not greater than 2^62 (about 4e18)
// return a * b % mod when mod is less than 2^31
inline uint64 mul_mod(uint64 a, uint64 b, uint64 mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  if (mod < int(1e9))
    return a * b % mod;
  uint64 k = (uint64)((long double)a * b / mod);
  uint64 res = a * b - k * mod;
  if ((int64)res < 0)
    res += mod;
  return res;
}

int64 pow_mod(int64 a, int64 n, int64 m) {
  int64 res = 1;
  for (a %= m; n; n >>= 1) {
    if (n & 1)
      res = mul_mod(res, a, m);
    a = mul_mod(a, a, m);
  }
  return res;
}

inline int64 add_mod(int64 x, int64 y, int64 mod) { return (x + y) % mod; }

inline int64 sub_mod(int64 x, int64 y, int64 mod) {
  return (x - y + mod) % mod;
}

inline int64 sqr(int64 x) { return x * x; }
inline int64 cub(int64 x) { return x * x * x; }

const int SZ = 100000000, MN = 10000;
int pi[SZ], pl[SZ], m;
void sieve() {
  m = 0;
  for (int i = 2; i < SZ; ++i)
    pi[i] = 1;
  for (int i = 2; i < SZ; ++i) {
    if (pi[i])
      pl[m++] = i;
    for (int j = 0; j < m && pl[j] * i < SZ; ++j) {
      pi[pl[j] * i] = 0;
      if (i % pl[j] == 0)
        break;
    }
  }
  for (int i = 2; i < SZ; ++i)
    pi[i] += pi[i - 1];
}
std::map<int64, int64> cache;
int64 phi(int64 x, int64 a) {
  if (a == 1 || !x)
    return (x + 1) / 2;
  int64 &r = cache[(x << 10) + a];
  if (r)
    return r;
  return r = phi(x, a - 1) - phi(x / pl[a - 1], a - 1);
}
int64 get_pi(int64 n) {
  if (n < SZ)
    return pi[n];
  int64 a = get_pi(pow(n, .25));
  int64 b = get_pi(sqrt(n));
  int64 c = get_pi(pow(n, 1. / 3));
  int64 r = phi(n, a) + (b + a - 2) * (b - a + 1) / 2;
  for (int i = a + 1; i <= b; ++i) {
    int64 w = n / pl[i - 1];
    r -= get_pi(w);
    if (i <= c) {
      int upp = get_pi(sqrt(w));
      for (int j = i; j <= upp; ++j)
        r += j - 1 - get_pi(w / pl[j - 1]);
    }
  }
  return r;
}

// return the sum of p^k for all p <= m, where m is in form floor(n / i)
// for m <= sqrt{n}, stored in ssum[m]; for m > sqrt{n} stored in lsum[n / m]
// note: if you need all correct value of ssum and lsum, please remove `mark`
// and make `delta` always be 1
std::pair<std::vector<int64>, std::vector<int64>> prime_count(int64 n, int64 k,
                                                              int64 mod) {
  auto pow_sum = [](int64 n, int64 k, int64 mod) {
    if (k == 0)
      return n;
    if (k == 1)
      return n * (n + 1) / 2 % mod;
  };
  const int64 v = static_cast<int64>(sqrt(n));
  std::vector<int64> ssum(v + 1), lsum(v + 1);
  // std::vector<bool> mark(v + 1);
  for (int i = 1; i <= v; ++i) {
    ssum[i] = pow_sum(i, k, mod) - 1;
    lsum[i] = pow_sum(n / i, k, mod) - 1;
  }
  for (int64 p = 2; p <= v; ++p) {
    if (ssum[p] == ssum[p - 1])
      continue;
    int64 psum = ssum[p - 1], q = p * p, ed = std::min(v, n / q);
    int64 pk = pow_mod(p, k, mod);
    int delta = (p & 1) + 1;
    for (int i = 1; i <= ed; i += delta) {
      // if (!mark[i])
      {
        int64 d = i * p;
        if (d <= v) {
          lsum[i] =
              sub_mod(lsum[i], sub_mod(lsum[d], psum, mod) * pk % mod, mod);
        } else {
          lsum[i] =
              sub_mod(lsum[i], sub_mod(ssum[n / d], psum, mod) * pk % mod, mod);
        }
      }
    }
    // for (int64 i = q; i <= ed; i += p * delta) mark[i] = true;
    for (int64 i = v; i >= q; --i) {
      ssum[i] =
          sub_mod(ssum[i], sub_mod(ssum[i / p], psum, mod) * pk % mod, mod);
    }
  }
  return {std::move(ssum), std::move(lsum)};
}

int main() {
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);
  cout.tie(NULL);

  vector<ll> ns = {1, 10, 100, 1000, (ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    cout << "n = " << (double)n;
    auto res = prime_count(n, 2, MOD);
    vector<ll> res_ssum = res.first, res_lsum = res.second;
    cout << "result: " << res_lsum.back() - 1;
    cout << "(time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s)"
         << endl;
  }

  return 0;
}
