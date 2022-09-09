#include <bits/stdc++.h>

using namespace std;

#define GODSPEED                                                               \
  ios::sync_with_stdio(0);                                                     \
  cin.tie(0);                                                                  \
  cout.tie(0);                                                                 \
  cout << fixed;                                                               \
  cout << setprecision(15);
#define f first
#define s second
#define newl cout << "\n";
#define pb push_back
#define mset(a, x) memset(a, x, sizeof(a))
#define debv(a)                                                                \
  for (auto it : a)                                                            \
    cout << it << " ";                                                         \
  newl;
#define deb1(a) cout << a << "\n";
#define deb2(a, b) cout << a << " " << b << "\n";
#define deb3(a, b, c) cout << a << " " << b << " " << c << "\n";
#define deb4(a, b, c, d) cout << a << " " << b << " " << c << " " << d << "\n";
#define uniq(a) a.resize(unique(a.begin(), a.end()) - a.begin());
#define all(a) a.begin(), a.end()

typedef long long ll;
typedef long double ld;

const ll N = 1e5 + 5;
const ll mod = 1e9 + 7;
const ll INF = 0x7f7f7f7f7f7f7f7f;
const int INFi = 0x7f7f7f7f;

struct gauss {
  int b, size;
  vector<int> basis;

  gauss() {}
  gauss(int b) : b(b), size(0) { basis.assign(b + 1, 0); }

  void insert(long long x) {
    for (int i = b; i >= 0; i--) {
      if (((x >> i) & 1) == 0)
        continue;
      if (!basis[i]) {
        basis[i] = x;
        size++;
        break;
      } else
        x ^= basis[i];
    }
  }

  bool is_part(long long x) {
    for (int i = b; i >= 0; i--) {
      if (((x >> i) & 1) == 0)
        continue;
      x ^= basis[i];
    }
    return x == 0;
  }
};

int n, q, a[N], l[N], x[N];

ll fastmod(ll x, ll y = mod - 2) {
  ll res = 1;
  x %= mod;
  while (y > 0) {
    if (y & 1)
      res = res * x % mod;
    x = x * x % mod;
    y /= 2;
  }
  return res;
}

void solve() {
  gauss g(19);

  cin >> n >> q;
  for (int i = 1; i <= n; i++)
    cin >> a[i];
  for (int i = 1; i <= q; i++)
    cin >> l[i] >> x[i];

  vector<int> ind(q), ans(q + 1);
  iota(all(ind), 1);
  sort(all(ind), [&](int i, int j) { return l[i] < l[j]; });

  int i = 1;
  for (auto it : ind) {
    while (i <= l[it]) {
      g.insert(a[i++]);
    }

    ans[it] = g.is_part(x[it]) * fastmod(2, i - 1 - g.size);
  }

  for (int i = 1; i <= q; i++)
    deb1(ans[i])
}

int main() {
  GODSPEED;
  int test = 1;
  // cin >> test;
  for (int i = 1; i <= test; i++) {
    solve();
  }
}

// https://codeforces.com/contest/959/problem/F
// https://www.youtube.com/watch?v=DBW7I6VwJgk
// https://codeforces.com/blog/entry/100066