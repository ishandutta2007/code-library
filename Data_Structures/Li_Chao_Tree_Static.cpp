
struct Node {
  Line line;
  Node *lc = nullptr;
  Node *rc = nullptr;
  void clear() {
    lc = nullptr;
    rc = nullptr;
    line = Line();
  }
};
auto *root = new Node[80000];
int sz;
class LiChao {
public:
  void InsertLineKnowingly(Node *&n, int tl, int tr, Line x) {
    if (n == nullptr)
      n = new Node();
    if (n->line.get(tl) < x.get(tl))
      swap(n->line, x);
    if (n->line.get(tr) >= x.get(tr))
      return;
    if (tl == tr)
      return;
    int mid = (tl + tr) / 2;
    if (n->line.get(mid) > x.get(mid)) {
      InsertLineKnowingly(n->rc, mid + 1, tr, x);
    } else {
      swap(n->line, x);
      InsertLineKnowingly(n->lc, tl, mid, x);
    }
  }

  void InsertLine(Node *&n, int tl, int tr, int l, int r, Line x) {
    if (tr < l || r < tl || tl > tr || l > r)
      return;
    if (n == nullptr)
      n = new Node();
    if (l <= tl && tr <= r)
      return InsertLineKnowingly(n, tl, tr, x);
    int mid = (tl + tr) / 2;
    InsertLine(n->lc, tl, mid, l, r, x);
    InsertLine(n->rc, mid + 1, tr, l, r, x);
  }

  int Query(Node *&n, int tl, int tr, int x) {
    if (n == nullptr)
      return -INF;
    if (tl == tr)
      return n->line.get(x);
    int res = n->line.get(x);
    int mid = (tl + tr) / 2;
    if (x <= mid) {
      res = max(res, Query(n->lc, tl, mid, x));
    } else {
      res = max(res, Query(n->rc, mid + 1, tr, x));
    }
    return res;
  }

  void InsertLine(int l, int r, Line x) {
    return InsertLine(root, 0, sz - 1, l, r, x);
  }

  int Query(int x) { return Query(root, 0, sz - 1, x); }

  void build(int n) { forn(i, 0, 4 * n) root[i].clear(); }
};
// https://github.com/ishandutta2007/codeforces/blob/a6901ca26964c035dee28818fe525b015d93e0b4/vercingetorix/normal/1175/G.cpp
// https://codeforces.com/problemset/problem/1175/G