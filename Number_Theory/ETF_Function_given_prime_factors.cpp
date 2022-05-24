#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

ll mulmod(ll a, ll b, ll m) {
  ll x = 0, y = a % m;
  while (b > 0) {
    if (b % 2 == 1) {
      x = (x + y) % m;
    }
    y = (y * 2) % m;
    b /= 2;
  }
  return x % m;
}

ll euler_totient(ll n, vector<ll> prime_factors) {
  ll etf = n;
  for (int i = 0; i < prime_factors.size(); i++) {
    etf = etf / prime_factors[i];
  }
  etf = etf % MOD;
  for (int i = 0; i < prime_factors.size(); i++) {
    etf = mulmod(etf, prime_factors[i] - 1, MOD);
  }
  // cout << "etf of " << n << "=" << etf << endl;
  return etf;
}
