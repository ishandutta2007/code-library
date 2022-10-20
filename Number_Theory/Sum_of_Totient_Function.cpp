#include <cstdio>
#include <cassert>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <vector>
#include <functional>

using namespace std;

using i64 = long long;
using u32 = unsigned;
using u64 = unsigned long long;
using i128 = __int128_t;

using u8 = unsigned char;

int isqrt(i64 N) {
  int x = sqrtl(N);
  return x;
}

int main() {
  i64 N;
  while (~scanf("%lld", &N)) {
    const int v = isqrt(N);
    vector<int> primes;
    vector<i64> s0(v + 1), s1(v + 1), l0(v + 1);
    vector<i128> l1(v + 1);

    // phi
    auto f = [&](int p, int e) -> i64 {
      i64 ret = p - 1;
      while (e > 1)
        --e, ret *= p;
      return ret;
    };

    auto divide = [](i64 n, i64 d) { return double(n) / d; };

    // sum f(p)
    for (int i = 1; i <= v; ++i)
      s0[i] = i - 1, s1[i] = i64(i) * (i + 1) / 2 - 1;
    for (int i = 1; i <= v; ++i)
      l0[i] = N / i - 1, l1[i] = i128(N / i) * (N / i + 1) / 2 - 1;
    for (int p = 2; p <= v; ++p)
      if (s0[p] > s0[p - 1]) {
        primes.push_back(p);
        i64 q = i64(p) * p, M = N / p, t0 = s0[p - 1], t1 = s1[p - 1];
        int t = v / p, u = min<i64>(v, N / q);
        for (int i = 1; i <= t; ++i)
          l0[i] -= (l0[i * p] - t0), l1[i] -= (l1[i * p] - t1) * p;
        for (int i = t + 1; i <= u; ++i)
          l0[i] -= (s0[divide(M, i)] - t0),
              l1[i] -= (s1[divide(M, i)] - t1) * p;
        for (int i = v; i >= q; --i)
          s0[i] -= (s0[divide(i, p)] - t0),
              s1[i] -= (s1[divide(i, p)] - t1) * p;
      }
    for (int i = 1; i <= v; ++i)
      s1[i] -= s0[i];
    for (int i = 1; i <= v; ++i)
      l1[i] -= l0[i];

    // sum g(n) (g(p^e) := f(p)^e)
    for (auto it = primes.rbegin(); it != primes.rend(); ++it) {
      int p = *it;
      i64 q = i64(p) * p, M = N / p, s = s1[p - 1];
      int t = v / p, u = min<i64>(v, N / q);
      for (i64 i = q; i <= v; ++i)
        s1[i] += (s1[divide(i, p)] - s) * f(p, 1);
      for (int i = u; i > t; --i)
        l1[i] += (s1[divide(M, i)] - s) * f(p, 1);
      for (int i = t; i >= 1; --i)
        l1[i] += (l1[i * p] - s) * f(p, 1);
    }
    for (int i = 1; i <= v; ++i)
      s1[i] += 1;
    for (int i = 1; i <= v; ++i)
      l1[i] += 1;

    // sum f(n)
    function<i128(i64, int, i128)> rec = [&](i64 n, size_t beg,
                                             i64 coeff) -> i128 {
      if (!coeff)
        return 0;
      i128 ret = i128(coeff) * (n > v ? l1[divide(N, n)] : s1[n]);
      for (size_t i = beg; i < primes.size(); ++i) {
        i64 p = primes[i], q = i64(p) * p;
        if (q > n)
          break;
        i64 nn = divide(n, q);
        for (int e = 2; nn > 0; nn = divide(nn, p), ++e) {
          ret += rec(nn, i + 1, coeff * (f(p, e) - f(p, 1) * f(p, e - 1)));
        }
      }
      return ret;
    };
    const int mod = 998244353;
    printf("%lld\n", i64(rec(N, 0, 1) % mod));
  }
  return 0;
}
// https://judge.yosupo.jp/problem/sum_of_totient_function
// min_25