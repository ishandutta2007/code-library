#include "iostream"
#include "algorithm"
#include "cstring"
#include "cstdio"
#include "cmath"
#include "vector"
#include "map"
#include "set"
#include "queue"
using namespace std;
#define MAXN 200006
//#define int long long
#define rep(i, a, b) for (int i = (a), i##end = (b); i <= i##end; ++i)
#define per(i, a, b) for (int i = (a), i##end = (b); i >= i##end; --i)
#define pii pair<int, int>
#define fi first
#define se second
#define mp make_pair
#define pb push_back
#define eb emplace_back
#define vi vector<int>
#define all(x) (x).begin(), (x).end()
#define mem(a) memset(a, 0, sizeof a)
typedef long long ll;
const int P = 1e9 + 7;
const int iv2 = P + 1 >> 1, iv3 = (P + 1) / 3;
ll n;

namespace Min_25 {
ll n;
int B;
int pri[MAXN], en, sp[MAXN], sp2[MAXN];
void sieve() {
  B = sqrt(n + 0.5);
  rep(i, 2, B) {
    if (!pri[i])
      pri[++en] = i, sp[en] = (sp[en - 1] + i) % P,
      sp2[en] = (sp2[en - 1] + i * 1ll * i) % P;
    for (int j = 1; j <= en && i * pri[j] <= B; ++j) {
      pri[i * pri[j]] = 1;
      if (i % pri[j] == 0)
        break;
    }
  }
}
int cc(ll x) {
  x %= P;
  return (x * (x + 1) / 2 + P - 1) % P;
}
int cc2(ll x) {
  x %= P;
  return (x * (x + 1) / 2 % P * (2 * x + 1) % P * iv3 + P - 1) % P;
}
ll g1[MAXN], g2[MAXN], A[MAXN];
int A1[MAXN], A2[MAXN];
int po(ll x) { return x <= B ? A1[x] : A2[n / x]; }
void getG() {
  int tot = 0;
  for (ll l = 1, r; l <= n; l = r + 1) {
    r = n / (n / l);
    ++tot;
    A[tot] = n / l;
    if (A[tot] <= B)
      A1[A[tot]] = tot;
    else
      A2[n / A[tot]] = tot;
    g1[tot] = cc(A[tot]), g2[tot] = cc2(A[tot]);
  }
  rep(i, 1, en) for (int j = 1; pri[i] * 1ll * pri[i] <= A[j]; ++j) {
    g1[j] = (g1[j] + P -
             pri[i] * 1ll * (g1[po(A[j] / pri[i])] + P - sp[i - 1]) % P) %
            P;
    g2[j] = (g2[j] + P -
             pri[i] * 1ll * pri[i] % P *
                 (g2[po(A[j] / pri[i])] + P - sp2[i - 1]) % P) %
            P;
  }
}

int g(int x) { return (g2[x] < g1[x] ? g2[x] - g1[x] + P : g2[x] - g1[x]); }
int f(ll x) {
  x %= P;
  return x * (x - 1) % P;
}

int S(ll n, int k) {
  if (n <= pri[k])
    return 0;
  int re = (1ll * g(po(n)) + P - sp2[k] + sp[k]) % P;
  for (int i = k + 1; pri[i] * 1ll * pri[i] <= n && i <= en; ++i) {
    ll cur = pri[i];
    for (int e = 1; cur <= n; ++e) {
      re = (re + f(cur) * 1ll * (S(n / cur, i) + (e != 1))) % P;
      cur *= pri[i];
    }
  }
  return re;
}

int solve(ll x) {
  n = x;
  sieve();
  getG();
  return (S(n, 0) + 1) % P;
}
}

void solve() {
  cin >> n;
  cout << Min_25::solve(n) << endl;
}

signed main() {
  //    int T;cin >> T;while( T-- ) solve();
  solve();
}