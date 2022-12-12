#include <bits/stdc++.h>
using namespace std;

int n, m;
vector<vector<long long>> dp;

void calc(int x = 0, int y = 0, int mask = 0, int next_mask = 0) {
  if (x == n)
    return;
  if (y >= m)
    dp[x + 1][next_mask] += dp[x][mask];
  else {
    int my_mask = 1 << y;
    if (mask & my_mask)
      calc(x, y + 1, mask, next_mask);
    else {
      calc(x, y + 1, mask, next_mask | my_mask);
      if (y + 1 < m && !(mask & my_mask) && !(mask & (my_mask << 1)))
        calc(x, y + 2, mask, next_mask);
    }
  }
}

int main() {
  cin >> n >> m;

  dp.resize(n + 1, vector<long long>(1 << m));
  dp[0][0] = 1;
  for (int x = 0; x < n; ++x)
    for (int mask = 0; mask < (1 << m); ++mask)
      calc(x, 0, mask, 0);

  cout << dp[n][0];
}
// https://cp-algorithms.com/dynamic_programming/profile-dynamics.html#implementation
// tiling problems
// https://projecteuler.net/problem=189
// 7255 - Land of Farms
// dp tiling
// http://fileadmin.cs.lth.se/contest/nwerc/Problemset_NWERC2004.pdf
// poj 1038
// sgu 131
// sgu 132
// sgu 223
// sgu 225
// zoj 1346
// poj 3254
// poj 1185
// poj 3311
// hdu 3001
// poj 2288
// zoj 4257
// hdu 3681
// poj 2430
// poj 2436
// poj 2541
// poj 2836
// poj 1699
// poj 2288
// poj 2688
// poj 3411
// poj 2686
// poj 1482
// poj 2690
// poj 3719
// poj 1795
// poj 1739
// poj 3593
// poj 2088
// UVA 10359 - Tiling
// UVA 10918 - Tri Tiling
// SPOJ GNY07H (Four Tiling)
// SPOJ M5TILE (Five Tiling)
// SPOJ MNTILE (MxN Tiling)
// SPOJ DOJ1
// SPOJ DOJ2
// SPOJ BTCODE_J
// SPOJ PBOARD
// ACM HDU 4285 - Circuits
// LiveArchive 4608 - Mosaic
// Timus 1519 - Formula 1
// Codeforces Parquet
