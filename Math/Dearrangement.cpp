#include <bits/stdc++.h>
#pragma optimize("Ofast")
using namespace std;
typedef long long ll;
#define rep(i, n) for (ll i = 0; (i) < (n); (i)++)
#define _rep(i, n) for (ll i = 1; (i) <= (n); (i)++)

#ifdef VinceBlack
#define _debug(...) fprintf(stderr, __VA_ARGS__);
#else
#define _debug(...)
#endif

#define MAXN 30

const ll N = MAXN;
const ll M = MAXN;
ll binom[N][M];
void init_binom() {
  _rep(i, N) binom[i - 1][0] = 1;
  _rep(i, N - 1) _rep(j, M - 1) binom[i][j] =
      binom[i - 1][j] + binom[i - 1][j - 1];
}

ll derangement[MAXN];
void init_derangement() {
  derangement[0] = 1;
  derangement[1] = 0, derangement[2] = 1;
#ifdef mod
  for (ll i = 3; i < MAXN; i++)
    derangement[i] =
        (i - 1) * ((derangement[i - 1] + derangement[i - 2]) % mod) % mod;
#else
  for (ll i = 3; i < MAXN, i <= 20; i++)
    derangement[i] = (i - 1) * (derangement[i - 1] + derangement[i - 2]);
// MAXIMUM 20
#endif
}

signed main() {
#ifdef Yee_172
  freopen("../in.txt", "r", stdin);
#endif
  init_binom();
  init_derangement();
  ll n, res;
  while (~scanf("%lld", &n) && n) {
    res = 0;
    for (ll i = 0; 2 * i <= n; i++)
      res += derangement[i] * binom[n][i];
    printf("%lld\n", res);
  }
  return 0;
}

// HDU_2068.