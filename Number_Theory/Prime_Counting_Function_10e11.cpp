#include "bits/stdc++.h"

using namespace std;
// n = 1e+10
// result: 455052511
// time: 0.146883s

// n = 1e+11
// result: 4118054813
// time: 0.762825s

// n = 1e+12
// result: 37607912018
// time: 3.89003s
using ll = long long;
ll count_primes(ll n) {
  vector<ll> v;
  for (ll k = 1; k * k <= n; ++k) {
    v.push_back(n / k);
    v.push_back(k);
  }
  sort(v.begin(), v.end());
  v.erase(unique(v.begin(), v.end()), v.end());

  // return i such that v[i] = x
  // since v[i] = i + 1 for i <= sqrt(n) and v[v.size() - i] = n / i for i <=
  // sqrt(n),
  // we can calculate index in O(1)
  ll sq = sqrt(n);
  auto geti = [&](ll x) {
    if (x <= sq)
      return (int)x - 1;
    else
      return (int)(v.size() - (n / x));
  };

  vector<ll> dp(v.size());

  // S(n, 0) = n
  for (int i = 0; i < v.size(); ++i)
    dp[i] = v[i];

  int a = 0;
  for (ll p = 2; p * p <= n; ++p) {
    // this condition is true for primes
    if (dp[geti(p)] != dp[geti(p - 1)]) {
      ++a;
      for (int i = (int)v.size() - 1; i >= 0; --i) {
        if (v[i] < p * p)
          break;
        dp[i] -= dp[geti(v[i] / p)] - a;
      }
    }
  }

  return dp[geti(n)] - 1;
}

int main() {
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);
  cout.tie(NULL);

  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    cout << "n = " << (double)n << endl;
    ll res = count_primes(n);
    cout << "result: " << res << endl;
    cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s"
         << endl;
    cout << endl;
  }

  return 0;
}
