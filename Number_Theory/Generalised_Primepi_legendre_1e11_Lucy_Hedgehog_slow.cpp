#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128;
// credits: Lucy_Hedgehog
// https://github.com/Bodigrim/arithmoi/issues/60
// Legendre
// O(n^(3/4))

// Prime Counting Times:
// n = 1e+10 : 470299869(time: 1.33058s)
// n = 10+11 : 99998937295(time: 2.45336s)
// n = 10+12 : 999987237212(time: 6.91859s)

// Prime Summing Times:
// n = 1e+10 :  7413165246951120256(time: 0.843584s)
// n = 1e+11 :  8746113789707973631(time: 1.74779s)
// n = 1e+12 :  8247104618519365631(time: 5.27918s)

using ull = unsigned long long;
using namespace std;
vector<ull> V;

unordered_map<ull, ull> S;

ll f(ll n) {
  // return 1;
  return n;
};

ll g(ll n) {
  // return n;
  return n * (n + 1) / 2;
};

ull prime_sum(ull n) {
  ull r = (ull)sqrt(n);
  V.reserve(r * 2 + 20);
  for (int i = 1; i <= r; i++)
    V.push_back(n / i);
  for (int i = V[V.size() - 1] - 1; i > 0; i--)
    V.push_back(i);
  for (int i = 0; i < V.size(); i++)
    S[V[i]] = g(V[i]) - 1;

  for (int p = 2; p <= r; p++)
    if (S[p] > S[p - 1]) {
      ull sp = S[p - 1];
      for (int i = 0; i < V.size() && V[i] >= p * p; i++)
        S[V[i]] -= f(p) * (S[V[i] / p] - sp);
    }
  return S[n];
}

int main() {
  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    ll res = prime_sum(n);
    cout << "n = " << n << " : " << res
         << "(time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s)"
         << endl;
  }
  return 0;
}
