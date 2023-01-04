#include <bits/stdc++.h>
using namespace std;

#define ll long long
#define MOD 998244353

// written by pauloamed

template <int MAXN, class T, T (*BinaryOp)(T, T)> struct segTreeLazy {
  T segt[3 * MAXN]; // segTree
  T lp[3 * MAXN];   // lazy prop

  T NEUTRAL_LP;
  T NEUTRAL_VAL;
  int n;

  void init(int _n, T neutralVal, T neutralLp, vector<T> *vals = NULL) {
    n = _n;
    NEUTRAL_VAL = neutralVal;
    NEUTRAL_LP = neutralLp;
    build(0, n - 1, 0, vals);
  }

  /*************************************************************/
  /*********************    CHANGE ME  *************************/
  /*************************************************************/
  void lazyProp(int node, T val, int l, int r) {
    assert(l == r);
    // if(l == r){
    segt[node] = val; // ATTENTION
    // }
    // if(l != r){
    //     lp[(node<<1)+1] = val;  // ATTENTION
    //     lp[(node<<1)+2] = val;  // ATTENTION
    // }
  }
  /*************************************************************/

  void maybeUpdtLazy(int l, int r, int node) {
    if (lp[node] != NEUTRAL_LP) {
      lazyProp(node, lp[node], l, r);
      lp[node] = NEUTRAL_LP;
    }
  }

  void build(int l, int r, int node, vector<T> *initVals) {
    lp[node] = NEUTRAL_LP;
    if (l == r) {
      if (initVals)
        segt[node] = (*initVals)[l];
      else
        segt[node] = NEUTRAL_VAL;
    } else {
      int mid = (l + r) / 2;
      build(l, mid, (node << 1) + 1, initVals);
      build(mid + 1, r, (node << 1) + 2, initVals);
      segt[node] = BinaryOp(segt[(node << 1) + 1], segt[(node << 1) + 2]);
    }
  }

  void _updateRange(int l, int r, int node, int lu, int ru, T val) {
    // maybeUpdtLazy(l,r,node);
    if (lu > r || l > ru)
      return;
    if (lu <= l && r <= ru) {
      lazyProp(node, val, l, r);
    } else {
      int mid = (l + r) / 2;
      _updateRange(l, mid, (node << 1) + 1, lu, ru, val);
      _updateRange(mid + 1, r, (node << 1) + 2, lu, ru, val);
      segt[node] = BinaryOp(segt[(node << 1) + 1], segt[(node << 1) + 2]);
    }
  }

  void update(int lu, int ru, T val) { _updateRange(0, n - 1, 0, lu, ru, val); }

  T _queryRange(int l, int r, int node, int lq, int rq) {
    // maybeUpdtLazy(l, r, node);
    if (rq < l || lq > r)
      return NEUTRAL_VAL;
    else if (l >= lq && r <= rq)
      return segt[node];
    else {
      int mid = (l + r) / 2; // mid
      return BinaryOp(_queryRange(l, mid, (node << 1) + 1, lq, rq),
                      _queryRange(mid + 1, r, (node << 1) + 2, lq, rq));
    }
  }

  T query(int lq, int rq) { return _queryRange(0, n - 1, 0, lq, rq); }
};

pair<int, int> mergeFunc(pair<int, int> a, pair<int, int> b) {
  if (a.first == b.first) {
    return make_pair(a.first, b.second + a.second);
  } else {
    if (b.first < a.first)
      return b;
    else
      return a;
  }
}

template <int MAXN, bool VALS_IN_EDGES, class T, T (*BinaryOp)(T, T)>
struct HLD {
  int n;
  vector<int> v[MAXN];
  int par[MAXN];
  int root[MAXN];
  int depth[MAXN];
  int sz[MAXN];
  int time;
  int pos[MAXN];
  vector<int> rpos; // rpos not used, but could be useful
  segTreeLazy<MAXN, T, BinaryOp> segTree;

  void addEdge(int x, int y) {
    v[x].push_back(y);
    v[y].push_back(x);
  }

  void dfsSz(int x) {
    sz[x] = 1;
    for (int i = 0; i < v[x].size(); ++i) {
      par[v[x][i]] = x;
      depth[v[x][i]] = depth[x] + 1;
      v[v[x][i]].erase(find(v[v[x][i]].begin(), v[v[x][i]].end(), x));
      dfsSz(v[x][i]);
      sz[x] += sz[v[x][i]];
      if (sz[v[x][i]] > sz[v[x][0]])
        swap(v[x][i], v[x][0]);
    }
  }

  void dfsHld(int x) {
    pos[x] = time++;
    rpos.push_back(x);
    for (auto y : v[x]) {
      if (y == v[x][0])
        root[y] = root[x];
      else
        root[y] = y;
      dfsHld(y);
    }
  }

  void init(int _n, int rootNode) {
    n = _n;
    par[rootNode] = depth[rootNode] = time = 0;
    dfsSz(rootNode);
    root[rootNode] = rootNode;
    dfsHld(rootNode);
  }

  int lca(int x, int y) {
    while (root[x] != root[y]) {
      if (depth[root[x]] > depth[root[y]])
        swap(x, y);
      y = par[root[y]];
    }
    return depth[x] < depth[y] ? x : y;
  }

  /// int dist(int x, int y) { // # edges on path
  /// 	return depth[x]+depth[y]-2*depth[lca(x,y)]; }

  template <class Op> void processPath(int x, int y, Op op, bool skipTop) {
    while (root[x] != root[y]) {
      if (depth[root[x]] > depth[root[y]])
        swap(x, y);
      op(pos[root[y]], pos[y]);
      y = par[root[y]];
    }
    if (depth[x] > depth[y])
      swap(x, y);
    op(pos[x] + VALS_IN_EDGES + skipTop, pos[y]); // assert valid range
  }

  void modifyPath(int x, int y, T v, bool skipTop = false) {
    processPath(x, y, [this, &v](int l, int r) { segTree.update(l, r, v); },
                skipTop);
  }

  T queryPath(int x, int y, bool skipTop = false) {
    T res = segTree.NEUTRAL_VAL;
    processPath(x, y, [this, &res](int l, int r) {
      // watch out if operation is not commutative
      //
      res = BinaryOp(segTree.query(l, r), res);
    }, skipTop);
    return res;
  }

  void modifySubtree(int x, int v) {
    segTree.upd(pos[x] + VALS_IN_EDGES, pos[x] + sz[x] - 1, v);
  }
};

typedef pair<ll, ll> line;
typedef pair<line, line> BothLines;

line mergeLine(line x, line y) {
  int a = x.first;
  int b = x.second;

  int c = y.first;
  int d = y.second;

  return make_pair((1LL * a * c) % MOD, (1LL * a * d + b) % MOD);
}

BothLines merge(BothLines a, BothLines b) {

  auto x = mergeLine(a.first, b.first);
  auto y = mergeLine(b.second, a.second);
  return make_pair(x, y);
}

HLD<200010, false, BothLines, merge> hld;
int main() {
  cin.tie(NULL);
  iostream::sync_with_stdio(false);
  int n, q;
  cin >> n >> q;

  vector<line> v(n);
  for (auto &x : v) {
    cin >> x.first >> x.second;
  }

  for (int i = 0; i < n - 1; ++i) {
    int a, b;
    cin >> a >> b;
    hld.addEdge(a, b);
  }

  hld.init(n, 0);
  const line neut = make_pair(1LL, 0LL);
  hld.segTree.init(n, make_pair(neut, neut), BothLines());

  for (int i = 0; i < n; ++i) {
    hld.modifyPath(i, i, make_pair(v[i], v[i]));
  }

  while (q--) {
    int id;
    cin >> id;
    if (id == 0) {
      int i, a, b;
      cin >> i >> a >> b;
      auto x = make_pair(a, b);
      hld.modifyPath(i, i, make_pair(x, x));
    } else {
      int u, v, x;
      cin >> u >> v >> x;
      int lca = hld.lca(u, v);
      auto fs = hld.queryPath(u, lca);
      auto gs = hld.queryPath(v, lca, true);
      auto f = mergeLine(gs.second, fs.first);
      int a = (1LL * x * f.first) % MOD;
      int b = f.second;
      cout << (0LL + a + b) % MOD << endl;
    }
  }
}
