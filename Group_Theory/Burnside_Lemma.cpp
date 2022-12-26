#include <cstdio>
#include <cstring>
#include <algorithm>
#define ll long long
using namespace std;
const ll N = 1e7 + 10, XJQ = 1e9 + 7, inv6 = 833333345000000041ll;
ll T, n, m, k, ans, P, cnt, num;
ll pri[N / 10], mu[N], p[110], c[110];
bool v[N];
void init() {
  mu[1] = 1;
  for (ll i = 2; i < N; i++) {
    if (!v[i])
      pri[++num] = i, mu[i] = -1;
    for (ll j = 1; j <= num && i * pri[j] < N; j++) {
      v[i * pri[j]] = 1;
      dans
    }
    for (ll i = 1; i < N; i++)
      mu[i] += mu[i - 1];
    return;
  }
  ll ksc(ll a, ll b, ll p) {
    a %= p;
    b %= p;
    ll tmp = (long double)a * b / p;
    long double ans = a * b - tmp * p;
    if (ans >= p)
      ans -= p;
    else if (ans < 0)
      ans += p;
    return ans;
  }
  ll power(ll x, ll b, ll p) {
    ll ans = 1;
    while (b) {
      if (b & 1)
        ans = ksc(ans, x, p);
      x = ksc(x, x, p);
      b >>= 1;
    }
    return ans;
  }
  void Get_M(ll n) {
    m = 2;
    for (ll l = 1, r; l <= n; l = r + 1) {
      r = n / (n / l);
      ll k = n / l;
      (m += ksc(ksc(ksc(k, k, P), k + 3, P), mu[r] - mu[l - 1], P)) %= P;
    }
    m = ksc(m, inv6, P);
    return;
  }
  ll f(ll n) {
    return power(m - 1, n, P) + ((n & 1) ? (P - m + 1) : (m - 1)) % P;
  }
  void dfs(ll dep, ll x, ll val) {
    if (dep > cnt) {
      (ans += ksc(f(n / x), val, P)) %= P;
      return;
    }
    dfs(dep + 1, x, val);
    val *= p[dep] - 1;
    x *= p[dep];
    dfs(dep + 1, x, val);
    for (ll i = 2; i <= c[dep]; i++)
      val = val * p[dep], x *= p[dep], dfs(dep + 1, x, val);
    return;
  }
  signed main() {
    scanf("%lld", &T);
    init();
    while (T--) {
      scanf("%lld%lld", &n, &k);
      if (n % XJQ != 0)
        P = XJQ;
      else
        P = XJQ * XJQ;
      Get_M(k);
      // printf("%lld",m);
      ll x = n;
      cnt = ans = 0;
      for (ll i = 1; i <= num && pri[i] * pri[i] <= x; i++) {
        if (x % pri[i] == 0) {
          p[++cnt] = pri[i];
          c[cnt] = 0;
          while (x % pri[i] == 0)
            c[cnt]++, x /= pri[i];
        }
      }
      if (x > 1)
        p[++cnt] = x, c[cnt] = 1;
      dfs(1, 1, 1);
      if (n % XJQ == 0)
        printf("%lld\n", ans / XJQ * power(n / XJQ, XJQ - 2, XJQ) % XJQ);
      else
        printf("%lld\n", ans * power(n, XJQ - 2, XJQ) % XJQ);
    }
    return 0;
  }

// https://www.luogu.com.cn/problem/P3307
// https://www.fatalerrors.org/a/p3307-sdoi2013-necklace-burnside-lemma-mobius-inversion-characteristic-equation.html#:~:text=Now%20consider%20the%20scheme%20of%20necklace.%20According%20to,the%20number%20of%20permutation%20rings%20with%20%28i%29%20%29.
// https://www.geeksforgeeks.org/orbit-counting-theorem-or-burnsides-lemma/
// https://crypto.stanford.edu/~blynn/polya/burnside.html
