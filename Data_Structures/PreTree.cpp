#include <bits/stdc++.h>

using namespace std;

struct PreTree {  
  struct Node { int l = -1, r = -1; };
  
  int n;
  vector<Node> T;

  PreTree(int n) : n(n), T(n) {}

  int merge(int u, int v, int b, int e) {
    if (u == -1) return v;
    if (v == -1) return u;
    if (u > v) swap(u, v);

    int m = (b + e) / 2;
      
    T[u].l = merge(T[u].l, T[v].l, b, m);
    T[u].r = merge(T[u].r, T[v].r, m, e);
    T[v].l = T[v].r = -1;
    
    if (v < m)
      T[u].l = merge(T[u].l, v, b, m);
    else
      T[u].r = merge(T[u].r, v, m, e);
    return u;
  }

  int get(int u, int pos) {
    if (u == -1 || u > pos) return -1;
    int ret = get(T[u].r, pos);
    if (ret == -1)
      ret = get(T[u].l, pos);
    if (ret == -1)
      ret = u;
    return ret;
  } 

  int Merge(int u, int v) { 
    return merge(u, v, 0, n);
  }

  int Get(int u, int pos) {
    return get(u, pos);
  }
};

int main() {
  ios_base::sync_with_stdio(false);
  cin.tie(0);
  
  int n; cin >> n;
  vector<vector<int>> graph(n);
  for (int i = 1; i < n; ++i) {
    int u, v; cin >> u >> v; --u; --v;
    graph[u].push_back(v);
    graph[v].push_back(u);
  }
  int q; cin >> q;
  vector<array<vector<int>, 3>> qrs(n);
  for (int i = 0; i < q; ++i) {
    int typ, v; cin >> typ >> v; --typ; --v;
    qrs[v][typ].push_back(i);
  }
  PreTree pt(q);
  set<int> filled;
  vector<int> ans(q, -1);
  function<int(int, int)> dfs = [&](int v, int par) {
    for (auto i : qrs[v][0])
      filled.insert(i);
    
    int ret = -1;
    for (auto i : qrs[v][1])
      ret = pt.Merge(ret, i);
    for (auto ch : graph[v]) {
      if (ch == par) continue;
      ret = pt.Merge(ret, dfs(ch, v));
    }

    for (auto i : qrs[v][2]) {
      auto it = filled.upper_bound(i);
      int last_empty = pt.Get(ret, i);
      if (it == filled.begin()) ans[i] = 0;
      else ans[i] = (last_empty < *prev(it));
    }

    for (auto i : qrs[v][0])
      filled.erase(i);
    return ret;
  };
  dfs(0, -1);

  for (int i = 0; i < q; ++i)
    if (ans[i] != -1)
      cout << ans[i] << '\n';

  return 0;
}

// https://codeforces.com/contest/343/problem/D
// https://codeforces.com/blog/entry/83170
