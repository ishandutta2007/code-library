#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef long double ld;

namespace bo {

typedef complex<ld> cd;

ld x_coord, x_pinpoint;

struct pt {
  ld x, y;
};

struct seg {
  pt a, b;
  int id;
  bool closed_a, closed_b;

  ld ycut() const {
    if (a.x == b.x)
      return a.y;
    return (a.x * b.y - a.y * b.x + (a.y - b.y) * x_coord) / (a.x - b.x);
  }

  bool operator<(const seg &t) const {
    ld u = ycut();
    ld v = t.ycut();
    if (u == v)
      return id < t.id;
    else
      return u < v;
  }
};

ld rerr(ld x, ld y) { return abs(x - y) / max(ld(1), y); }

struct event {
  ld t;
  int i, j, tag;

  bool operator<(const event &b) const {
    if (rerr(t, b.t) < 1e-14)
      return tag < b.tag;
    return t < b.t;
  }
};

struct solver {

  int n;
  vector<seg> a;

  solver(const vector<seg> &arg) : n(arg.size()), a(arg) {
    for (int i = 0; i < n; i++) {
      a[i].id = i;
    }
  }

  set<seg> ms;
  multiset<event> pq;
  string alive;

  void enlist(seg u, seg v) {
    auto a = u.a, b = u.b, c = v.a, d = v.b;
    ld xp =
        (a.y * b.x * c.x - a.x * b.y * c.x - a.y * b.x * d.x + a.x * b.y * d.x -
         a.x * c.y * d.x + b.x * c.y * d.x + a.x * c.x * d.y - b.x * c.x * d.y);
    ld xq = a.y * c.x - b.y * c.x - a.x * c.y + b.x * c.y - a.y * d.x +
            b.y * d.x + a.x * d.y - b.x * d.y;
    if (xq != 0) {
      ld x = xp / xq;
      if (x >= x_pinpoint)
        pq.insert({x, u.id, v.id, 0});
    }
  }

  bool seg_cut(int id) {
    swap(x_coord, x_pinpoint);
    auto it = ms.upper_bound(a[id]);
    ld mre = 1;
    if (it != ms.end()) {
      mre = min(mre, rerr(a[id].ycut(), it->ycut()));
    }
    if (it != ms.begin()) {
      --it;
      mre = min(mre, rerr(a[id].ycut(), it->ycut()));
    }
    swap(x_coord, x_pinpoint);
    return mre < 1e-14;
  }

  void seg_insert(int id) {
    auto it = ms.insert(a[id]).first;
    alive[id] = 1;
    if (it != ms.begin())
      enlist(*prev(it), *it);
    if (next(it) != ms.end())
      enlist(*it, *next(it));
  }

  void seg_erase(int id) {
    auto it = ms.find(a[id]);
    if (it != ms.begin() && next(it) != ms.end())
      enlist(*prev(it), *next(it));
    ms.erase(it);
    alive[id] = 0;
  }

  ld solve() {
    for (int i = 0; i < n; i++) {
      if (a[i].a.x == a[i].b.x) {
        pq.insert({a[i].a.x, i, -3, 0});
      } else {
        pq.insert({a[i].a.x, i, -1, a[i].closed_a ? -1 : 1});
        pq.insert({a[i].b.x, i, -2, a[i].closed_b ? 1 : -1});
      }
    }

    alive = string(n, 0);

    while (pq.size()) {
      auto ev = *pq.begin();
      // cerr << ev.t << ' ' << ev.i << ' ' << ev.j << ' ' << ev.tag << '\n';
      x_pinpoint = ev.t;
      x_coord = ev.t + ev.tag * 1e-12;
      pq.erase(pq.begin());
      if (ev.j == -1) {
        if (seg_cut(ev.i))
          return ev.t;
        seg_insert(ev.i);
      } else if (ev.j == -2) {
        seg_erase(ev.i);
      } else if (ev.j == -3) {
        if (seg_cut(ev.i))
          return ev.t;
      } else if (alive[ev.i] && alive[ev.j]) {
        return ev.t;
      }
    }

    return 1e50;
  }
};

} // end namespace bentley_ottmann

namespace tree {

int n;
basic_string<int> e[100005];
vector<int> p, h, d, r, pos;

int dfs(int x) {
  int sz = 1, best = 0;
  for (int y : e[x]) {
    if (y == p[x])
      continue;
    p[y] = x;
    d[y] = d[x] + 1;
    int z = dfs(y);
    if (z > best) {
      best = z;
      h[x] = y;
    }
    sz += z;
  }
  return sz;
}

// indeksirano od 0, hld radi za grane, pazi na nekomutativne op!
void init() {
  p = h = d = r = pos = vector<int>(n, -1);
  d[0] = 0;
  dfs(0);
  int k = 0;
  for (int i = 0; i < n; i++) {
    if (i == 0 || h[p[i]] != i) {
      vector<int> put;
      for (int j = i; j != -1; j = h[j])
        put.push_back(j);
      reverse(put.begin(), put.end());
      for (int x : put)
        pos[x] = k++, r[x] = i;
    }
  }
}

template <class Op> void path(int x, int y, Op op) {
  while (r[x] != r[y]) {
    if (d[r[x]] > d[r[y]]) {
      op(pos[x], pos[r[x]], -1);
      x = p[r[x]];
    } else {
      op(pos[y], pos[r[y]], 1);
      y = p[r[y]];
    }
  }

  if (d[x] > d[y]) {
    op(pos[x], pos[y], -2);
  } else if (d[y] > d[x]) {
    op(pos[y], pos[x], 2);
  } else {
    op(pos[x], pos[x], 0);
  }
}
}

int main() {
  ios_base::sync_with_stdio(false);
  cin.tie(nullptr);
  cout.tie(nullptr);
  cerr.tie(nullptr);
  cout << setprecision(20) << fixed;
  cerr << setprecision(20) << fixed;

  int m;
  cin >> tree::n >> m;
  for (int i = 1; i < tree::n; i++) {
    int u, v;
    cin >> u >> v;
    u--, v--;
    tree::e[u] += v;
    tree::e[v] += u;
  }

  tree::init();

  using bo::seg;
  using bo::pt;
  using bo::solver;

  vector<seg> a;

  while (m--) {
    int t, c, v, u;
    cin >> t >> c >> v >> u;

    v--, u--;

    vector<pair<int, int>> g, h, gg, hh, mid;

    auto op = [&](int x, int y, int sgn) {
      if (sgn == -2)
        gg.emplace_back(x, y);
      else if (sgn == -1)
        g.emplace_back(x, y);
      else if (sgn == 0)
        mid.emplace_back(x, y);
      else if (sgn == 1)
        h.emplace_back(x, y);
      else
        hh.emplace_back(x, y);
    };

    tree::path(v, u, op);

    reverse(h.begin(), h.end());

    ld tm = t;

    // cerr << "szvec " << g.size() << ' ' << gg.size() << ' ' << mid.size();
    // cerr << ' ' << hh.size() << ' ' << h.size() << '\n';

    for (auto[x, y] : g) {
      ld t2 = tm + ld(y - x + 1) / c;
      ld xx = x;
      ld yy = y + 1;
      a.push_back({{tm, xx}, {t2, yy}, 0, true, false});
      tm = t2;
    }

    for (auto[x, y] : gg) {
      ld t2 = tm + ld(y - x) / c;
      ld xx = x;
      ld yy = y;
      a.push_back({{tm, xx}, {t2, yy}, 0, true, true});
      tm = t2;
    }

    for (auto[x, y] : mid) {
      ld xx = x;
      a.push_back({{tm, xx}, {tm, xx}, 0, true, true});
    }

    for (auto[x, y] : hh) {
      ld t2 = tm + ld(y - x) / c;
      ld yy = x;
      ld xx = y;
      a.push_back({{tm, xx}, {t2, yy}, 0, true, true});
      tm = t2;
    }

    for (auto[x, y] : h) {
      ld t2 = tm + ld(y - x + 1) / c;
      ld yy = x;
      ld xx = y + 1;
      a.push_back({{tm, xx}, {t2, yy}, 0, false, true});
      tm = t2;
    }
  }

  // for (auto s : a)
  // cerr << s.a.x << ' ' << s.a.y << " --- " << s.b.x << ' ' << s.b.y << '\n';

  solver sv(a);
  ld t = sv.solve();

  // check for duplicates
  set<array<ld, 2>> sd;
  for (auto s : a) {
    if (s.a.x != s.b.x)
      continue;
    auto p = sd.insert({s.a.x, s.a.y});
    if (!p.second)
      t = min(t, s.a.x);
  }

  if (t > 1e20)
    cout << "-1\n";
  else
    cout << t << "\n";
}

/// https://codeforces.com/problemset/problem/704/E