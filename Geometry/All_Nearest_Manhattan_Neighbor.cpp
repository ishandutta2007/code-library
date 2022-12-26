// all-nearest-neighbor-manhattan
// all-nearest-neighbor problem is solved using KD Tree in O(nlog n) when
// distances are eucledian
// but if the distances are manhattan the problem can be solved in O(n)
#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;
inline void getmin(int &x, int y) {
  if (y < x)
    x = y;
}
#define debug() printf("%d\n", __LINE__)

struct p {
  int x, y, ix;
};

int lxly(const p &a, const p &b) {
  if (a.x == b.x && a.y == b.y)
    return a.ix < b.ix;
  if (a.x == b.x)
    return b.y > a.y;
  return b.x > a.x;
}
int lxgy(const p &a, const p &b) {
  if (a.x == b.x)
    return b.y < a.y;
  return b.x > a.x;
}

p ps1[200002], ps2[200002];
int pps[200002];
int ans[200002];

int ty[200002], y[200002];
int maxy;

int A[200002], B[200002], C[200002], D[200002];

void add(int i, int x, int *d) {
  for (; i <= maxy; i += i & -i)
    if (x > d[i])
      d[i] = x;
}

int query(int i, int *d) {
  int ret = d[0];
  for (; i > 0; i -= i & -i)
    if (ret < d[i])
      ret = d[i];
  return ret;
}

int bsearch(int v) {
  int l = 0, r = maxy, mid;
  while (l < r) {
    mid = (l + r) / 2;
    if (y[mid] == v)
      return mid;
    if (y[mid] < v)
      l = mid + 1;
    else
      r = mid;
  }
  return -1;
}

int main() {
  int i, j;
  int n;

  scanf("%d", &n);
  for (i = 0; i < n; i++) {
    scanf("%d %d", &ps1[i].x, &ty[i]);
    ps1[i].y = ty[i];
    ps1[i].ix = i;
    ps2[i] = ps1[i];
    ans[i] = 100000000;
  }

  sort(ty, ty + n);

  maxy = -1;
  for (i = 0; i < n; i++)
    if (maxy == -1 || ty[i] != y[maxy])
      y[++maxy] = ty[i];
  maxy++;

  for (i = 0; i < n; i++)
    pps[ps1[i].ix] = bsearch(ps1[i].y) + 1;
  maxy++;

  sort(ps1, ps1 + n, lxly);
  sort(ps2, ps2 + n, lxgy);

  for (i = 0; i <= maxy; i++)
    A[i] = B[i] = C[i] = D[i] = -100000000;

  for (i = 0, j = n - 1; i < n; i++, j--) {
    getmin(ans[ps1[i].ix], ps1[i].x + ps1[i].y - query(pps[ps1[i].ix] + 1, A));
    add(pps[ps1[i].ix] + 1, ps1[i].x + ps1[i].y, A);

    getmin(ans[ps2[i].ix],
           ps2[i].x - ps2[i].y - query(maxy - pps[ps2[i].ix], B));
    add(maxy - pps[ps2[i].ix], ps2[i].x - ps2[i].y, B);

    getmin(ans[ps2[j].ix], -ps2[j].x + ps2[j].y - query(pps[ps2[j].ix] + 1, C));
    add(pps[ps2[j].ix] + 1, -ps2[j].x + ps2[j].y, C);

    getmin(ans[ps1[j].ix],
           -ps1[j].x - ps1[j].y - query(maxy - pps[ps1[j].ix], D));
    add(maxy - pps[ps1[j].ix], -ps1[j].x - ps1[j].y, D);
  }

  for (i = 0; i < n; i++)
    printf("%d\n", ans[i]);

  return 0;
}