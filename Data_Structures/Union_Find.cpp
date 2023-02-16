#include <algorithm>
#include <iostream>
#include <set>
#include <utility>
#include <vector>
using namespace std;

#define repeat(n) for (int repc = (n); repc > 0; --repc)
typedef long long int64;

const int MAXN = 100005;

struct UnionFind {
  int n;
  vector<int> dad, rank;

  UnionFind(int n) : n(n) {
    dad.resize(n);
    rank.resize(n);
    reset();
  }

  void reset() {
    for (int i=0; i<n; ++i) {
      dad[i] = i;
    }
    fill(rank.begin(), rank.end(), 0);
  }

  int find(int a) {
    int top;
    for (top=a; top != dad[top]; top=dad[top]) ;
    while (a != top) { int x = dad[a]; dad[a] = top; a = x; }
    return top;
  }

  int union_find(int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) {
      if (rank[a] > rank[b]) {
        dad[b] = a;
      } else {
        dad[a] = b;
        if (rank[a] == rank[b]) {
          ++rank[b];
        }
      }
      return 1;
    }
    return 0;
  }
};

int main() {
  cin.sync_with_stdio(0);

  int n, m, p, q;
  cin >> n >> m >> p >> q;

  UnionFind union_find(n);
  struct Edge { int a, b, w; };
  vector<Edge> edges(m);
  for (Edge &e : edges) {
    cin >> e.a >> e.b >> e.w;
    --e.a; --e.b;
    union_find.union_find(e.a, e.b);
  }

  vector<int64> sizes(n, 0);
  for (const Edge &e : edges) {
    sizes[union_find.find(e.a)] += e.w;
  }

  set<pair<int64, int> > pq;
  int ncomp = 0;
  for (int i=0; i<n; ++i) {
    if (union_find.find(i) == i) {
      ++ncomp;
      pq.insert(make_pair(sizes[i], i));
    }
  }

  if (q > ncomp) {
    cout << "NO\n";
    return 0;
  }

  vector<pair<int, int> > output;
  repeat (p) {
    if (pq.size() > 1 && ncomp > q) {
      auto first = pq.begin();
      auto next = first; ++next;
      int64 S = first->first + next->first;
      int64 newlen = min<int64>(1000000000, S + 1);
      Edge edge = { first->second, next->second, newlen };
      edges.push_back(edge);
      output.push_back(make_pair(edge.a, edge.b));

      pair<int64, int> newelem(S + newlen, first->second);
      pq.erase(first);
      pq.erase(next);
      pq.insert(newelem);
      --ncomp;
    } else {
      if (edges.empty()) {
        cout << "NO\n";
        return 0;
      }

      const Edge &e = edges.front();
      output.push_back(make_pair(e.a, e.b));
    }
  }

  if (ncomp != q) {
    cout << "NO\n";
    return 0;
  }

  cout << "YES\n";
  for (const auto &p : output) {
    cout << p.first+1 << ' ' << p.second+1 << '\n';
  }

  return 0;
}

// https://codeforces.com/contest/362/problem/D
