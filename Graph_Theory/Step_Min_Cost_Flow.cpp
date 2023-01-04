#include <bits/extc++.h>
#include <bits/stdc++.h>
#define int long long
using namespace std;

using ll = long long;
const ll INF = numeric_limits<ll>::max() / 4;
typedef vector<ll> VL;

struct StepMCMF {
  int N;
  vector<vector<int>> ed, red;
  vector<VL> cap, flow, cost;
  vector<int> seen;
  VL dist, pi;
  vector<pair<int, int>> par;

  StepMCMF(int _N)
      : N(_N), ed(N), red(N), cap(N, VL(N)), flow(cap), cost(cap), seen(N),
        dist(N), pi(N), par(N) {}

  void addEdge(int from, int to, int capa, int costEd) {
    cap[from][to] = capa;
    cost[from][to] = costEd;
    ed[from].push_back(to);
    red[to].push_back(from);
  }

  void path(int s) {
    fill(seen.begin(), seen.end(), 0);
    fill(dist.begin(), dist.end(), INF);
    dist[s] = 0;
    int di;
    __gnu_pbds::priority_queue<pair<int, int>> q;
    vector<decltype(q)::point_iterator> its(N);
    q.push({0, s});

    auto relax = [&](int i, int capa, int costEd, int dir) {
      int val = di - pi[i] + costEd;
      if (capa and val < dist[i]) {
        dist[i] = val;
        par[i] = {s, dir};
        if (its[i] == q.end())
          its[i] = q.push({-dist[i], i});
        else
          q.modify(its[i], {-dist[i], i});
      }
    };

    while (!q.empty()) {
      s = q.top().second;
      q.pop();
      seen[s] = 1;
      di = dist[s] + pi[s];
      for (int i : ed[s])
        if (!seen[i])
          relax(i, cap[s][i] - flow[s][i], cost[s][i], 1);
      for (int i : red[s])
        if (!seen[i])
          relax(i, flow[i][s], -cost[i][s], 0);
    }
    for (int i = 0; i < N; ++i)
      pi[i] = min(pi[i] + dist[i], INF);
  }

  // returns the min cost for flow=1, flow=2, etc...
  vector<int> maxflow(int s, int t) {
    vector<int> ret;
    int totflow = 0, totcost = 0;
    while (path(s), seen[t]) {
      int fl = INF;
      for (int p, r, x = t; tie(p, r) = par[x], x != s; x = p)
        fl = min(fl, r ? cap[p][x] - flow[p][x] : flow[x][p]);
      assert(fl == 1);
      totflow += fl;
      for (int p, r, x = t; tie(p, r) = par[x], x != s; x = p)
        if (r) {
          flow[p][x] += fl;
          totcost += fl * cost[p][x];
        } else {
          flow[x][p] -= fl;
          totcost -= fl * cost[x][p];
        }
      ret.push_back(totcost);
    }
    return ret;
  }
};

signed main(void) {
  ios_base::sync_with_stdio(false);
  cin.tie(0);

  int N;
  cin >> N;
  vector<int> A(N), B(N), C(N);
  const int M = 1e9;
  const int X = 150;
  vector<vector<int>> best(X, vector<int>(X, 0));
  for (int i = 0; i < N; ++i) {
    cin >> A[i] >> B[i] >> C[i];
    --A[i], --B[i];
    best[A[i]][B[i]] = max(best[A[i]][B[i]], C[i]);
  }

  StepMCMF mcmf(2 * X + 2);

  int s = 2 * X, t = 2 * X + 1;

  for (int i = 0; i < X; ++i)
    mcmf.addEdge(s, i, 1, 0);
  for (int i = 0; i < X; ++i)
    mcmf.addEdge(X + i, t, 1, 0);
  for (int i = 0; i < X; ++i)
    for (int j = 0; j < X; ++j)
      if (best[i][j])
        mcmf.addEdge(i, X + j, 1, M - best[i][j]);
  auto sol = mcmf.maxflow(s, t);
  cout << sol.size() << endl;
  // sol[i] is negative so we make the trick to reduce M*flow
  for (int i = 0; i < (int)sol.size(); ++i)
    cout << (i + 1) * M - sol[i] << '\n';
}