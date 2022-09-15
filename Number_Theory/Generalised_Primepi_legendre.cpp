#include <bits/stdc++.h>
using ll = long long;
using namespace std;
// Legendre
// O(n^(3/4))
// n = 1e+10 : 2220822432581729238(time: 0.194799s)
// n = 1e+11 : -1447107067060386762(time: 1.02261s)
// n = 1e+12 : -1932149121990928815(time: 5.73414s)

ll f(ll n, ll k) {
  // should be multiplicative (f(ab) = f(a)f(b)),
  // the result will be sum f(p) over all primes
  if (k == 0)
    return 1;
  if (k == 1)
    return (ll)n;
  if (k == 2)
    return (ll)n * n;
};

ll pref(ll n, ll k) {
  // should return sum_{i=1..n} f(i)
  if (k == 0)
    return n;
  if (k == 1)
    return ((ll)n * n + n) / 2;
  if (k == 2)
    return (ll)n * (n + 1) * (2 * n + 1) / 6;
};

ll sum_prime_pow_k(ll n, ll k) {
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
  vector<ll> s(v.size());
  for (int i = 0; i < s.size(); ++i)
    s[i] = pref(v[i], k) - 1;
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
        s[i] -= (s[geti(v[i] / p)] - sp) * f(p, k);
      }
    }
  }

  return s[s.size() - 1];
}

int main() {
  // for (ll i = 1; i <= 10; i++)
  // printf("%lld %lld %lld %lld\n", i, sum_prime_pow_k(i, 0),
  //        sum_prime_pow_k(i, 1), sum_prime_pow_k(i, 2));
  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    ll res = sum_prime_pow_k(n, 1);
    cout << "n = " << n << " : " << res
         << "(time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s)"
         << endl;
  }

  return 0;
}
