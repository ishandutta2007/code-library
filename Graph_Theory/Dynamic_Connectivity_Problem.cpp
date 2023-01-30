/**
 * Author: Sergey Kopeliovich (Burunduk30@gmail.com)
 * Date: 2013.09.08
 */

#include <cstdio>
#include <cassert>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;

#define forn(i, n) for (int i = 0; i < (int)(n); i++)
#define forit(i, a) for (__typeof((a).begin()) i = (a).begin(); i != (a).end(); i++)
#define pb push_back
#define mp make_pair
#define sz(a) (int)(a).size()

typedef vector <int> vi;
typedef pair <int, int> pii;

const int N = (int)3e5;

struct edge
{
  int a, b, l, r;
};
typedef vector <edge> List;

int cnt[N + 1], ans[N], u[N], color[N], deg[N];
vi g[N];

void add( int a, int b )
{
  g[a].pb(b), g[b].pb(a);
}

void dfs( int v, int value )
{
  u[v] = 1, color[v] = value;
  forn(i, sz(g[v]))
    if (!u[g[v][i]])
      dfs(g[v][i], value);
}

void go( int l, int r, const List &v, int vn, int add_vn ) // [l, r)
{
  if (cnt[l] == cnt[r])
    return; // no queries, only changes
  if (!sz(v))
  {
    while (l < r)
      ans[l++] = vn + add_vn;
    return;
  }  

  List v1;
  forn(i, vn)
    g[i].clear();
  forn(i, sz(v))
    if (v[i].a != v[i].b)
    {
      if (v[i].l <= l && r <= v[i].r)
        add(v[i].a, v[i].b);
      else if (l < v[i].r && v[i].l < r)
        v1.pb(v[i]);
    }

  int vn1 = 0;
  forn(i, vn)
    u[i] = 0;
  forn(i, vn)
    if (!u[i])
      deg[vn1] = 0, dfs(i, vn1++);
  
  forn(i, sz(v1))
  {
    v1[i].a = color[v1[i].a];
    v1[i].b = color[v1[i].b];
    if (v1[i].a != v1[i].b)
      deg[v1[i].a]++, deg[v1[i].b]++;
  }

  vn = vn1, vn1 = 0;
  forn(i, vn)
    u[i] = vn1, vn1 += (deg[i] > 0), add_vn += !deg[i];
  forn(i, sz(v1))
  {
    v1[i].a = u[v1[i].a];
    v1[i].b = u[v1[i].b];
  }

  int m = (l + r) / 2; // [l, m) [m, r)
  go(l, m, v1, vn1, add_vn);
  go(m, r, v1, vn1, add_vn);
}

int main()
{
  #define NAME "connect"
  assert(freopen(NAME ".in", "r", stdin));
  assert(freopen(NAME ".out", "w", stdout));

  map <pii, int> m;
  List v;
  int n, k;
  scanf("%d%d", &n, &k);
  forn(i, k)
  {
    char type;
    scanf(" %c", &type);
    if (type == '+' || type == '-')
    {
      int a, b;
      scanf("%d%d", &a, &b), a--, b--;
      if (a == b) // loop
      {
        k--, i--;
        continue;
      }
      if (a > b)
        swap(a, b);

      if (type == '+')
        m[mp(a, b)] = i;
      else
      {
        int &j = m[mp(a, b)];
        v.pb({a, b, j, i});
        j = -1;
      }
    }
    else
      cnt[i + 1]++;
    cnt[i + 1] += cnt[i];
  }
  forit(it, m)
    if (it->second != -1)
      v.pb({it->first.first, it->first.second, it->second, k});

  go(0, k, v, n, 0);
  
  forn(i, k)
    if (cnt[i + 1] != cnt[i])
      printf("%d\n", ans[i]);
  return 0;
}
//https://codeforces.com/blog/entry/15296
