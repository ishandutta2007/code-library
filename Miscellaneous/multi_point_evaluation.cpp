#include <bits/stdc++.h>

#define reg register
#define pr std::pair<int, int>
#define fi first
#define se second
#define FIN(s) freopen(s, "r", stdin)
#define FOUT(s) freopen(s, "w", stdout)
#define debug(...) fprintf(stderr, __VA_ARGS__)
#define rep(i, l, r) for (int i = l; i <= r; ++i)
#define lep(i, l, r) for (int i = l; i < r; ++i)
#define irep(i, r, l) for (int i = r; i >= l; --i)
#define ilep(i, r, l) for (int i = r; i > l; --i)
#define Rep(i, n) rep(i, 1, n)
#define Lep(i, n) lep(i, 1, n)
#define IRep(i, n) irep(i, n, 1)
#define ILep(i, n) ilep(i, n, 1)
typedef long long ll;
typedef long double ld;

namespace modular {
const int MOD = 786433;
inline int add(int x, int y) { return (x += y) >= MOD ? x -= MOD : x; }
inline void inc(int &x, int y) { (x += y) >= MOD ? x -= MOD : 0; }
inline int mul(int x, int y) { return 1LL * x * y % MOD; }
inline int qpow(int x, int y) {
  int ans = 1;
  for (; y; y >>= 1, x = mul(x, x))
    if (y & 1)
      ans = mul(ans, x);
  return ans;
}
}; // namespace modular
using namespace modular;

namespace Base {
template <typename Tp> inline Tp input() {
  Tp x = 0, y = 1;
  char c = getchar();
  while ((c < '0' || '9' < c) && c != EOF) {
    if (c == '-')
      y = -1;
    c = getchar();
  }
  if (c == EOF)
    return 0;
  while ('0' <= c && c <= '9')
    x = x * 10 + c - '0', c = getchar();
  return x *= y;
}
template <typename Tp> inline void read(Tp &x) { x = input<Tp>(); }
template <typename Tp> inline void chmax(Tp &x, Tp y) { x < y ? x = y : 0; }
template <typename Tp> inline void chmin(Tp &x, Tp y) { x > y ? x = y : 0; }
}; // namespace Base
using namespace Base;
/*----------------------------------------------------------------------------*/

#define G 10
#define MAX_N 1000007

int N, M, Q;
int a[MAX_N], b[MAX_N], q[MAX_N], rev[MAX_N], ans[MAX_N], mp[MAX_N];

void NTT(int *A, int lim) {
  lep(i, 0, lim) if (i < rev[i]) std::swap(A[i], A[rev[i]]);
  for (int l = 1; l < lim; l <<= 1) {
    int gn = qpow(G, (MOD - 1) / (l << 1));
    for (int i = 0, L = l << 1; i < lim; i += L) {
      for (int j = i, x, y, g = 1; j < i + l; ++j, g = mul(g, gn)) {
        x = A[j], y = mul(A[j + l], g);
        A[j] = add(x, y), A[j + l] = add(x, MOD - y);
      }
    }
  }
}

void solve() {
  M = 1 << 18;
  int cur = 1;
  lep(i, 0, MOD - 1) {
    mp[cur] = i;
    cur = mul(cur, G);
  }
  Rep(i, Q) q[i] = !q[i] ? -1 : mp[q[i]];
  lep(i, 0, M) rev[i] = (rev[i >> 1] >> 1) | (i & 1) << 17;
  lep(i, 0, N) b[i] = a[i];
  NTT(b, M);
  Rep(i, Q) if (q[i] % 3 == 0) ans[i] = b[q[i] / 3];
  int bs = 1;
  lep(i, 0, N) b[i] = mul(a[i], bs), bs = mul(bs, G);
  lep(i, N, M) b[i] = 0;
  NTT(b, M);
  Rep(i, Q) if (q[i] % 3 == 1) ans[i] = b[(q[i] - 1) / 3];
  int g2 = mul(G, G);
  bs = 1;
  lep(i, 0, N) b[i] = mul(a[i], bs), bs = mul(bs, g2);
  lep(i, N, M) b[i] = 0;
  NTT(b, M);
  Rep(i, Q) if (q[i] % 3 == 2) ans[i] = b[(q[i] - 2) / 3];
  Rep(i, Q) printf("%d\n", q[i] == -1 ? a[0] : ans[i]);
}

int main() {
#ifndef ONLINE_JUDGE
  FIN("a.in");
  FOUT("a.out");
#endif
  read(N);
  N++;
  lep(i, 0, N) read(a[i]);
  read(Q);
  Rep(i, Q) read(q[i]);
  solve();
  return 0;
}

// https://www.codechef.com/problems/POLYEVAL
// multi point evaluation is O(n*log^2(n))
// interpolation is O(n*log^3(n))