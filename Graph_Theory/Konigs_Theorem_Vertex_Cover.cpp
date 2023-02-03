struct match_1 {
  int nx, ny, m, adj[maxm], next[maxm], last[maxn], matx[maxn], maty[maxn],
      que[maxn], level[maxn], run[maxn];
  void init(int _nx, int _ny) {
    nx = _nx;
    ny = _ny;
    m = 0;
    for (int i = 1; i <= nx; i++)
      matx[i] = -1, last[i] = -1;
    for (int i = 1; i <= ny; i++)
      maty[i] = -1;
  }
  void add(int u, int v) {
    adj[m] = v;
    next[m] = last[u];
    last[u] = m++;
  }
  bool bfs() {
    bool found = 0;
    int size = 0;
    for (int x = 1; x <= nx; x++) {
      if (matx[x] != -1)
        level[x] = -1;
      else {
        level[x] = 0;
        que[size++] = x;
      }
    }
    for (int i = 0; i < size; i++) {
      for (int x = que[i], e = last[x], y; e != -1; e = next[e]) {
        y = adj[e];
        if (maty[y] == -1)
          found = 1;
        else if (level[maty[y]] == -1) {
          level[maty[y]] = level[x] + 1;
          que[size++] = maty[y];
        }
      }
    }
    return found;
  }
  int dfs(int x) {
    for (int &e = run[x], y; e != -1; e = next[e]) {
      y = adj[e];
      if (maty[y] == -1 || (level[maty[y]] == level[x] + 1 && dfs(maty[y]))) {
        maty[y] = x;
        matx[x] = y;
        return 1;
      }
    }
    return 0;
  }
  int maxMatching() {
    int total = 0;
    while (bfs()) {
      for (int x = 1; x <= nx; x++)
        run[x] = last[x];
      for (int x = 1; x <= nx; x++)
        if (matx[x] == -1)
          total += dfs(x);
    }
    return total;
  }
} match1;
struct vertex_cover {
  int n, m;
  vector<int> a[maxn];
  bool choosed[maxn];
  void init(int _n, int _m) {
    n = _n;
    m = _m;
    for (int i = 1; i <= n; i++) {
      a[i].clear();
      choosed[i] = 0;
    }
    match1.init(n, m);
  }
  void add(int u, int v) {
    a[u].push_back(v);
    match1.add(u, v);
  }
  void konig(int u) {
    choosed[u] = 1;
    for (int i = 0, v, x; i < a[u].size(); i++) {
      v = a[u][i];
      x = match1.maty[v];
      if (x != -1 && x != oo) {
        match1.maty[v] = oo;
        if (!choosed[x])
          konig(x);
      }
    }
  }
  void solve() {
    cout << match1.maxMatching();
    for (int i = 1; i <= n; i++)
      if (match1.matx[i] == -1)
        konig(i);
    for (int i = 1; i <= n; i++)
      if (match1.matx[i] != -1 && !choosed[i])
        cout << " r" << i;
    for (int j = 1; j <= m; j++)
      if (match1.maty[j] == oo)
        cout << " c" << j;
    cout << "\n";
  }
} cover;

// https://en.wikipedia.org/wiki/K%C5%91nig%27s_theorem_(graph_theory)
