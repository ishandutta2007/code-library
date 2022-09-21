// Credits: naayan
#include <bits/stdc++.h>
#define whlie while
#define pb push_back
#define eb emplace_back
#define fi first
#define se second
#define rep(i, N) for (int i = 0; i < (N); i++)
#define repr(i, N) for (int i = (N)-1; i >= 0; i--)
#define rep1(i, N) for (int i = 1; i <= (N); i++)
#define repr1(i, N) for (int i = (N); i > 0; i--)
#define each(x, v) for (auto &x : v)
#define all(v) (v).begin(), (v).end()
#define sz(v) ((int)(v).size())
#define vrep(v, it) for (auto it = v.begin(); it != v.end(); it++)
#define vrepr(v, it) for (auto it = v.rbegin(); it != v.rend(); it++)
#define ini(...)                                                               \
  int __VA_ARGS__;                                                             \
  in(__VA_ARGS__)
#define inl(...)                                                               \
  ll __VA_ARGS__;                                                              \
  in(__VA_ARGS__)
#define ins(...)                                                               \
  string __VA_ARGS__;                                                          \
  in(__VA_ARGS__)
using namespace std;
void solve();
using ll = long long;
using vl = vector<ll>;
using vi = vector<int>;
using vvi = vector<vector<int>>;
constexpr int inf = 1001001001;
constexpr ll infLL = (1LL << 61) - 1;
struct IoSetupNya {
  IoSetupNya() {
    cin.tie(nullptr);
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(15);
    cerr << fixed << setprecision(7);
  }
} iosetupnya;
template <typename T, typename U> inline bool amin(T &x, U y) {
  return (y < x) ? (x = y, true) : false;
}
template <typename T, typename U> inline bool amax(T &x, U y) {
  return (x < y) ? (x = y, true) : false;
}
template <typename T, typename U>
ostream &operator<<(ostream &os, const pair<T, U> &p) {
  os << p.first << " " << p.second;
  return os;
}
template <typename T, typename U>
istream &operator>>(istream &is, pair<T, U> &p) {
  is >> p.first >> p.second;
  return is;
}
template <typename T> ostream &operator<<(ostream &os, const vector<T> &v) {
  int s = (int)v.size();
  rep(i, s) os << (i ? " " : "") << v[i];
  return os;
}
template <typename T> istream &operator>>(istream &is, vector<T> &v) {
  for (auto &x : v)
    is >> x;
  return is;
}
void in() {}
template <typename T, class... U> void in(T &t, U &... u) {
  cin >> t;
  in(u...);
}
void out() { cout << "\n"; }
template <typename T, class... U> void out(const T &t, const U &... u) {
  cout << t;
  if (sizeof...(u))
    cout << " ";
  out(u...);
}
template <typename T> void die(T x) {
  out(x);
  exit(0);
}
#ifdef NyaanDebug
#include "NyaanDebug.h"
#define trc(...)                                                               \
  do {                                                                         \
    cerr << #__VA_ARGS__ << " = ";                                             \
    dbg_out(__VA_ARGS__);                                                      \
  } while (0)
#define trca(v, N)                                                             \
  do {                                                                         \
    cerr << #v << " = ";                                                       \
    array_out(v, N);                                                           \
    cout << endl;                                                              \
  } while (0)
#else
#define trc(...)
#define trca(...)
int main() { solve(); }
#endif

using P = pair<int, int>;
using vp = vector<P>;
constexpr int MOD = /** 1000000007; //*/ 998244353;
////////////////

template <typename T, T INF> struct DynamicLiChaoTree {

  struct Line {
    T slope, intercept;
    Line() : slope(0), intercept(INF) {}
    Line(T slope, T intercept) : slope(slope), intercept(intercept) {}
    inline T get(T x) const { return slope * x + intercept; }
    inline bool over(const Line &other, const T &x) const {
      return get(x) < other.get(x);
    }
  };

  // remind 1-indexed
  T xmin, xmax, _size;
  unordered_map<T, Line> seg;

  // [l , r]Li Chao Tree
  DynamicLiChaoTree(T xmin, T xmax) : xmin(xmin), xmax(xmax) {
    _size = 1;
    while (_size < xmax - xmin + 1)
      _size <<= 1;
  }

  //
  void update(T a, T b, T left, T right, T seg_idx) {
    Line line(a, b);
    while (1) {
      int mid = (left + right) >> 1;
      bool l_over = line.over(seg[seg_idx], min(xmax, left + xmin));
      bool r_over = line.over(seg[seg_idx], min(xmax, right - 1 + xmin));
      if (l_over == r_over) {
        if (l_over)
          swap(seg[seg_idx], line);
        return;
      }
      bool m_over = line.over(seg[seg_idx], min(xmax, mid + xmin));
      if (m_over)
        swap(seg[seg_idx], line);
      if (l_over != m_over)
        seg_idx = (seg_idx << 1), right = mid;
      else
        seg_idx = (seg_idx << 1) | 1, left = mid;
    }
  }
  //
  void update(T a, T b, T seg_idx) {
    T left, right;
    T upper_bit = 63 - __builtin_clzll(seg_idx);
    left = (_size >> upper_bit) * (seg_idx - (1 << upper_bit));
    right = left + (_size >> upper_bit);
    update(a, b, left, right, seg_idx);
  }
  /*
  void update(T a,T b,T left,T right){
    T lower_bit = (left | right) & (- (left | right) );
    T seg_idx = ( _size + left ) / lower_bit;
    update(a , b , left , right , seg_idx);
  }
  */
  // y = ax + b
  void update(T a, T b) { update(a, b, 0, _size, 1); }

  // x in [left , right]y = ax + b
  void update_line_segment(T a, T b, T left, T right) {
    left -= xmin - _size, right -= xmin - _size - 1;
    for (; left < right; left >>= 1, right >>= 1) {
      if (left & 1)
        update(a, b, left++);
      if (right & 1)
        update(a, b, --right);
    }
  }

  // x
  T query(T x) {
    T left = 0, right = _size, seg_idx = 1, idx = x - xmin, ret = INF;
    while (1) {
      T cur = seg[seg_idx].get(x);
      // if(cur == INF) break;
      ret = min(ret, cur);
      if (left + 1 >= right)
        break;
      int mid = (left + right) >> 1;
      if (idx < mid)
        seg_idx = (seg_idx << 1), right = mid;
      else
        seg_idx = (seg_idx << 1) | 1, left = mid;
    }
    return ret;
  }

  void merge(DynamicLiChaoTree<T, INF> &other) {
    if (seg.size() < other.seg.size())
      seg.swap(other.seg);
    for (auto &x : other.seg) {
      if (x.second.intercept == INF)
        continue;
      update(x.second.slope, x.second.intercept, x.first);
    }
  }
};

/*
void solve(){
  ini(N);
  vl a(N) , x(N) , y(N); in(a , x , y);
  vl dp(N , 0);
  DynamicLiChaoTree<ll , infLL> lichao(0 , 100000);
  lichao.update_line_segment(1 , y[0] - x[0] , x[0] , 100000);
  lichao.update_line_segment(-1, y[0] + x[0] , 0 , x[0]);
  each(x , lichao.seg) trc(x.first, x.second.slope , x.second.intercept);
  rep(i,N){
    if(i){
      lichao.update_line_segment(1 , dp[i-1] + y[i] - x[i] , x[i] , 100000);
      lichao.update_line_segment(-1, dp[i-1] + y[i] + x[i] , 0 , x[i]);
    }
    dp[i] = lichao.query(a[i]);
  }
  trc(dp);
  out(dp[N - 1]);
}
*/

int N;
vvi g;

template <typename G> G convert_rooted_tree(const G &g, int r = 0) {
  int N = (int)g.size();
  G rg(N);
  vector<int> v(N);
  v[r] = 1;
  queue<int> que;
  que.emplace(r);
  while (!que.empty()) {
    auto p = que.front();
    que.pop();
    for (auto &to : g[p]) {
      if (v[to] == 0) {
        v[to] = 1;
        que.emplace(to);
        rg[p].emplace_back(to);
      }
    }
  }
  return rg;
}

vl dp, a, b;

void merge(DynamicLiChaoTree<ll, infLL> &lichao, int cur) {
  if (g[cur].size() == 0) {
    // dp[cur]=0
    lichao.update(b[cur], dp[cur]);
    return;
  }
  for (int i = 0; i < sz(g[cur]); i++) {
    // g[cur][i]
    DynamicLiChaoTree<ll, infLL> lichao2(-100000, 100000);
    merge(lichao2, g[cur][i]);
    //
    lichao.merge(lichao2);
  }
  // dp[cur]
  dp[cur] = lichao.query(a[cur]);
  lichao.update(b[cur], dp[cur]);
}

void solve() {
  in(N);
  a.resize(N);
  b.resize(N);
  in(a, b);
  g.resize(N);
  rep(i, N - 1) {
    ini(x, y);
    x--;
    y--;
    g[x].pb(y);
    g[y].pb(x);
  }
  vvi rg = convert_rooted_tree(g);
  g.swap(rg);
  dp.resize(N, 0);
  DynamicLiChaoTree<ll, infLL> lichao(-100000, 100000);
  merge(lichao, 0);
  out(dp);
}
// https://codeforces.com/problemset/problem/932/F
