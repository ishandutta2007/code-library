// Retired?
#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using ull = long unsigned long;
using ld = double long;

int n, k;

template <class T, class D, class E, class G, class U, class F = plus<D>>
struct tree_dp {
  vector<T> a;
  D e;
  vector<D> up, dn, rt;
  F f;
  G g;
  U u;
  vector<pair<int, E> *> p;
  vector<vector<pair<int, E>>> h;

  // a : vertex array - typed T
  // b : array of edges - typed E
  // g : (D, T) -> D - extends the sum of subtrees with the root vertex
  // u : (D, E) -> D - extends a subtree with an edge
  // e : neutral element of D
  // f : (D, D) -> D - adds together two (extended) subtrees
  tree_dp(const vector<T> &a, const vector<tuple<int, int, E>> &b, G g, U u,
          D e = D(), F f = F())
      : a(a), e(e), up(a.size()), dn(a.size()), f(f), g(g), u(u), p(a.size()),
        h(a.size()) {
    for (auto[x, y, z] : b) {
      h[x].emplace_back(y, z);
      h[y].emplace_back(x, z);
    }
    dfs1(0, 0);
    dfs2(0);
    rt.resize(a.size());
    rt[0] = dn[0];
    for (int x = 1; x < (int)a.size(); x++) {
      D d = u(up[x], p[x]->second);
      for (auto &hh : h[x]) {
        auto & [ y, w ] = hh;
        if (y != p[x]->first)
          d = f(d, u(dn[y], w));
      }
      rt[x] = g(d, a[x]);
    }
  }

  D operator()(int x) { return rt[x]; }

  void dfs1(int x, int o) {
    D d = e;
    for (auto &hh : h[x]) {
      auto & [ y, w ] = hh;
      if (y == o) {
        p[x] = &hh;
        continue;
      }
      dfs1(y, x);
      d = f(d, u(dn[y], w));
    }
    dn[x] = g(d, a[x]);
  }

  void dfs2(int x) {
    D d = x ? u(up[x], p[x]->second) : e, t = e;
    int m = h[x].size() - !!x;
    vector<D> pref(m, e);
    int i = 0;
    for (auto &hh : h[x]) {
      auto & [ y, w ] = hh;
      if (p[x] && p[x]->first == y)
        continue;
      pref[i++] = t;
      t = f(t, u(dn[y], w));
    }
    for (int j = (int)h[x].size() - 1; j >= 0; j--) {
      auto &hh = h[x][j];
      auto & [ y, w ] = hh;
      if (p[x] && p[x]->first == y)
        continue;
      i--;
      up[y] = g(f(pref[i], d), a[x]);
      d = f(d, u(dn[y], w));
    }
    for (auto &hh : h[x]) {
      auto & [ y, w ] = hh;
      if (p[x] && p[x]->first == y)
        continue;
      dfs2(y);
    }
  }
};

struct value {
  int a[40];

  value operator+(const value &b) const {
    value t;
    for (int i = 0; i < 2 * k; i++) {
      t.a[i] = a[i] ^ b.a[i];
    }
    return t;
  }
};

value extend_edge(value v, int) {
  rotate(v.a, v.a + 2 * k - 1, v.a + 2 * k);
  return v;
}

value extend_root(value v, int x) {
  v.a[0] ^= x;
  return v;
}

const value neut = []() {
  value t;
  memset(t.a, 0, sizeof(t.a));
  return t;
}();

int main() {
  ios::sync_with_stdio(!cin.tie(0));

  cin >> n >> k;
  vector<int> a(n);
  vector<tuple<int, int, int>> gr(n - 1);
  for (int i = 0; i < n - 1; i++) {
    int x, y;
    cin >> x >> y;
    gr[i] = {x - 1, y - 1, 1};
  }
  for (int i = 0; i < n; i++)
    cin >> a[i];

  tree_dp dp(a, gr, extend_root, extend_edge, neut, plus<value>());
  for (int i = 0; i < n; i++) {
    auto v = dp(i);
    cout << !!accumulate(v.a + k, v.a + k + k, 0, bit_xor<int>())
         << " \n"[i == n - 1];
  }
}
// https://codeforces.com/contest/1498/submission/126109184
