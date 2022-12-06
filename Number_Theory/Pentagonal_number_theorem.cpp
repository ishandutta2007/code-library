#include <bits/stdc++.h>
#include <atcoder/all>

using namespace std;
using namespace atcoder;

using ll = long long;
using mint = modint;

const int mod = 200003;

mint fact[mod];
mint ifact[mod];

// binom(n, k) mod 200003
mint binom(ll n, ll k) {
  if (k < 0 or k > n)
    return 0;
  mint res = 1;
  while (n) {
    ll n0 = n % mod;
    ll k0 = k % mod;
    if (n0 < k0)
      return 0;
    res *= fact[n0] * ifact[k0] * ifact[n0 - k0];
    n /= mod;
    k /= mod;
  }
  return res;
}

int main() {
  mint::set_mod(mod);
  fact[0] = 1;
  for (int i = 1; i < mod; i++)
    fact[i] = fact[i - 1] * i;
  ifact[mod - 1] = fact[mod - 1].inv();
  for (int i = mod - 1; i >= 1; i--)
    ifact[i - 1] = ifact[i] * i;

  ll n, m;
  cin >> n >> m;

  mint ans = 0;
  for (ll i = 0; i * (3 * i + 1) / 2 <= m - n; i++) {
    ll j = i * (3 * i + 1) / 2;
    ans += binom(m + n - j - 1, 2 * n - 1) * (i % 2 ? -1 : 1);
  }
  for (ll i = 1; i * (3 * i - 1) / 2 <= m - n; i++) {
    ll j = i * (3 * i - 1) / 2;
    ans += binom(m + n - j - 1, 2 * n - 1) * (i % 2 ? -1 : 1);
  }
  cout << ans.val() << endl;
}
// https://atcoder.jp/contests/abc279/editorial/5314
