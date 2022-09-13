#include <bits/stdc++.h>

using namespace std;


using int64 = long long;
const int mod = 1000003;

const int64 infll = (1LL << 62) - 1;
const int inf = (1 << 30) - 1;

struct IoSetup {
  IoSetup() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(10);
    cerr << fixed << setprecision(10);
  }
} iosetup;


template< typename T1, typename T2 >
ostream &operator<<(ostream &os, const pair< T1, T2 > &p) {
  os << p.first << " " << p.second;
  return os;
}

template< typename T1, typename T2 >
istream &operator>>(istream &is, pair< T1, T2 > &p) {
  is >> p.first >> p.second;
  return is;
}

template< typename T >
ostream &operator<<(ostream &os, const vector< T > &v) {
  for(int i = 0; i < (int) v.size(); i++) {
    os << v[i] << (i + 1 != v.size() ? " " : "");
  }
  return os;
}

template< typename T >
istream &operator>>(istream &is, vector< T > &v) {
  for(T &in : v) is >> in;
  return is;
}

template< typename T1, typename T2 >
inline bool chmax(T1 &a, T2 b) { return a < b && (a = b, true); }

template< typename T1, typename T2 >
inline bool chmin(T1 &a, T2 b) { return a > b && (a = b, true); }

template< typename T = int64 >
vector< T > make_v(size_t a) {
  return vector< T >(a);
}

template< typename T, typename... Ts >
auto make_v(size_t a, Ts... ts) {
  return vector< decltype(make_v< T >(ts...)) >(a, make_v< T >(ts...));
}

template< typename T, typename V >
typename enable_if< is_class< T >::value == 0 >::type fill_v(T &t, const V &v) {
  t = v;
}

template< typename T, typename V >
typename enable_if< is_class< T >::value != 0 >::type fill_v(T &t, const V &v) {
  for(auto &e : t) fill_v(e, v);
}

template< typename F >
struct FixPoint : F {
  FixPoint(F &&f) : F(forward< F >(f)) {}

  template< typename... Args >
  decltype(auto) operator()(Args &&... args) const {
    return F::operator()(*this, forward< Args >(args)...);
  }
};

template< typename F >
inline decltype(auto) MFP(F &&f) {
  return FixPoint< F >{forward< F >(f)};
}

template< typename T >
struct edge {
  int src, to;
  T cost;

  edge(int to, T cost) : src(-1), to(to), cost(cost) {}

  edge(int src, int to, T cost) : src(src), to(to), cost(cost) {}

  edge &operator=(const int &x) {
    to = x;
    return *this;
  }

  operator int() const { return to; }
};

template< typename T >
using Edges = vector< edge< T > >;
template< typename T >
using WeightedGraph = vector< Edges< T > >;
using UnWeightedGraph = vector< vector< int > >;
template< typename T >
using Matrix = vector< vector< T > >;

template< typename G >
struct CentroidDecomposition {
  const G &g;
  vector< int > sub;
  vector< vector< int > > belong;
  vector< bool > v;

  CentroidDecomposition(const G &g) : g(g), sub(g.size()), v(g.size()), belong(g.size()) {}

  inline int build_dfs(int idx, int par) {
    sub[idx] = 1;
    for(auto &to : g[idx]) {
      if(to == par || v[to]) continue;
      sub[idx] += build_dfs(to, idx);
    }
    return sub[idx];
  }

  inline int search_centroid(int idx, int par, const int mid) {
    for(auto &to : g[idx]) {
      if(to == par || v[to]) continue;
      if(sub[to] > mid) return search_centroid(to, idx, mid);
    }
    return idx;
  }

  inline void belong_dfs(int idx, int par, int centroid) {
    belong[idx].emplace_back(centroid);
    for(auto &to : g[idx]) {
      if(to == par || v[to]) continue;
      belong_dfs(to, idx, centroid);
    }
  }

  inline int build(UnWeightedGraph &t, int idx) {
    int centroid = search_centroid(idx, -1, build_dfs(idx, -1) / 2);
    v[centroid] = true;
    // belong_dfs(centroid, -1, centroid);
    for(auto &to : g[centroid]) {
      if(!v[to]) t[centroid].emplace_back(build(t, to));
    }
    v[centroid] = false;
    return centroid;
  }

  inline int build(UnWeightedGraph &t) {
    t.resize(g.size());
    return build(t, 0);
  }
};


template< class T >
struct ArrayPool {
  vector< T > pool;
  vector< T * > stock;
  int ptr;

  ArrayPool(int sz) : pool(sz), stock(sz) {}

  inline T *alloc() { return stock[--ptr]; }

  inline void free(T *t) { stock[ptr++] = t; }

  void clear() {
    ptr = (int) pool.size();
    for(int i = 0; i < pool.size(); i++) stock[i] = &pool[i];
  }
};

template< typename T, T x_low, T x_high, T id >
struct DynamicLiChaoTree {

  struct Line {
    T a, b;

    Line(T a, T b) : a(a), b(b) {}

    inline T get(T x) const { return a * x + b; }
  };

  struct Node {
    Line x;
    Node *l, *r;

    Node(const Line &x) : x{x}, l{nullptr}, r{nullptr} {}
  };

  Node *root;

  DynamicLiChaoTree() : root{nullptr} {}

  Node *add_line(Node *t, Line &x, const T &l, const T &r, const T &x_l, const T &x_r) {
    if(!t) return new Node(x);

    T t_l = t->x.get(l), t_r = t->x.get(r);

    if(t_l >= x_l && t_r >= x_r) {
      return t;
    } else if(t_l <= x_l && t_r <= x_r) {
      t->x = x;
      return t;
    } else {
      T m = (l + r) / 2;
      if(m == r) --m;
      T t_m = t->x.get(m), x_m = x.get(m);
      if(t_m < x_m) {
        swap(t->x, x);
        if(x_l <= t_l) t->l = add_line(t->l, x, l, m, t_l, t_m);
        else t->r = add_line(t->r, x, m + 1, r, t_m + x.a, t_r);
      } else {
        if(t_l <= x_l) t->l = add_line(t->l, x, l, m, x_l, x_m);
        else t->r = add_line(t->r, x, m + 1, r, x_m + x.a, x_r);
      }
      return t;
    }
  }

  void add_line(const T &a, const T &b) {
    Line x(a, b);
    root = add_line(root, x, x_low, x_high, x.get(x_low), x.get(x_high));
  }

  T query(const Node *t, const T &l, const T &r, const T &x) const {
    if(!t) return id;
    if(l == r) return t->x.get(x);
    T m = (l + r) / 2;
    if(m == r) --m;
    if(x <= m) return max(t->x.get(x), query(t->l, l, m, x));
    else return max(t->x.get(x), query(t->r, m + 1, r, x));
  }

  T query(const T &x) const {
    return query(root, x_low, x_high, x);
  }
};

int main() {
  int N;
  cin >> N;
  UnWeightedGraph g(N);
  for(int i = 1; i < N; i++) {
    int a, b;
    cin >> a >> b;
    --a, --b;
    g[a].emplace_back(b);
    g[b].emplace_back(a);
  }
  vector< int > A(N);
  cin >> A;
  CentroidDecomposition< UnWeightedGraph > cpd(g);
  UnWeightedGraph t;
  int root = cpd.build(t);


  vector< int > used(N);


  const int64 LIM = 150001;
  DynamicLiChaoTree< int64, 0, LIM, -infll > lct;
  auto down_path = MFP([&](auto down_path, int idx, int par, int64 ans, int64 sz, int64 all) -> int64 {
    ++sz;
    ans += A[idx] * sz;
    all += A[idx];
    bool have = false;
    int64 ret = ans;
    for(auto &to : g[idx]) {
      if(to == par || used[to]) continue;
      chmax(ret, down_path(to, idx, ans, sz, all));
      have = true;
    }
    if(!have) lct.add_line(all, ans);
    return ret;
  });

  auto up_path = MFP([&](auto up_path, int idx, int par, int64 ans, int64 sum, int64 sz) -> int64 {
    sum += A[idx];
    ans += sum;
    sz++;
    int64 ret = -infll;
    bool have = false;
    for(auto &to : g[idx]) {
      if(to == par || used[to]) continue;
      chmax(ret, up_path(to, idx, ans, sum, sz));
      have = true;
    }
    if(!have) chmax(ret, max(lct.query(sz), 0LL) + ans);
    return ret;
  });

  auto rec = MFP([&](auto rec, int idx) -> int64 {
    used[idx] = true;
    int64 ret = A[idx];
    lct.root = nullptr;
    for(auto &to : g[idx]) {
      if(used[to]) continue;
      chmax(ret, up_path(to, idx, A[idx], A[idx], 1));
      chmax(ret, down_path(to, idx, 0, 0, 0));
    }
    reverse(begin(g[idx]), end(g[idx]));
    lct.root = nullptr;
    for(auto &to : g[idx]) {
      if(used[to]) continue;
      chmax(ret, up_path(to, idx, A[idx], A[idx], 1));
      chmax(ret, down_path(to, idx, 0, 0, 0));
    }
    for(auto &to : t[idx]) chmax(ret, rec(to));
    used[idx] = false;
    return ret;
  });
  cout << rec(root) << "\n";
}
//https://codeforces.com/problemset/problem/1303/G
