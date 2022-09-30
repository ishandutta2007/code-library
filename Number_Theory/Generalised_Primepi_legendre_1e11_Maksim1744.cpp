#include <bits/stdc++.h>
using ll = long long;
using i128 = __int128;
using namespace std;
// credits:  Maksim1744
// Legendre
// O(n^(3/4))
// Prime Counting timings:
// n = 1e+10 : 455052511(time: 0.283723s)
// n = 1e+11 : 4118054813(time: 1.40211s)
// n = 1e+12 : 37607912018(time: 7.02708s)

// Prime summing timings:
// n = 1e+10 : 2220822432581729238(time: 0.322144s)
// n = 1e+11 : 201467077743744681014(time: 1.63643s)
// n = 1e+12 : 18435588552550705911377(time: 8.39142s)

i128 func(ll n, ll k) {
  // should be multiplicative (func(ab) = func(a)func(b)),
  // the result will be sum func(p) over all primes
  if (k == 0)
    return 1;
  if (k == 1)
    return (i128)n;
  if (k == 2)
    return (i128)n * n;
};

i128 accfunc(ll n, ll k) {
  // should return sum_{i=1..n} func(i)
  if (k == 0)
    return n;
  if (k == 1) // can handle upto 63 bit n
    return ((i128)n * n + n) / 2;
  if (k == 2) // can handle upto 42 bit n
    return (i128)n * (n + 1) * (2 * n + 1) / 6;
};

i128 primepi_generalised(ll n, ll k) {
  vector<ll> v;
  v.reserve((int)sqrt(n) * 2 + 20);
  ll sq;
  {
    ll k = 1;
    for (; k * k <= n; ++k) {
      v.push_back(k);
    }
    --k;
    sq = k;
    if (k * k == n)
      --k;
    for (; k >= 1; --k) {
      v.push_back(n / k);
    }
  }
  vector<i128> s(v.size());
  for (int i = 0; i < s.size(); ++i)
    s[i] = accfunc(v[i], k) - 1;
  auto geti = [&](ll x) {
    if (x <= sq)
      return (int)x - 1;
    else
      return (int)(v.size() - (n / x));
  };
  for (ll p = 2; p * p <= n; ++p) {
    if (s[p - 1] != s[p - 2]) {
      ll sp = s[p - 2];
      ll p2 = p * p;
      for (int i = (int)v.size() - 1; i >= 0; --i) {
        if (v[i] < p2) {
          break;
        }
        s[i] -= (s[geti(v[i] / p)] - sp) * func(p, k);
      }
    }
  }

  return s[s.size() - 1];
}

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

int main() {
  // for (ll i = 1; i <= 10; i++)
  // printf("%lld %lld %lld %lld\n", i, primepi_generalised(i, 0),
  //        primepi_generalised(i, 1), primepi_generalised(i, 2));
  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    i128 res = primepi_generalised(n, 0);
    cout << "n = " << n << " : " << res
         << "(time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s)"
         << endl;
  }

  return 0;
}
