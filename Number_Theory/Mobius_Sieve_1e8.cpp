#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

// Sieves till 1e8 in 2 sec

using namespace std;
typedef long long ll;
const ll m = 100000001;
bool flag[m];
ll p[m], u[m], i, j, k, n, tot;
void mobius_sieve() {
  ll i, j, k;
  // memset(flag,0,sizeof(flag));
  tot = 0;
  for (i = 2, u[1] = 1; i < m; i++) {
    if (!flag[i])
      p[++tot] = i, u[i] = -1;
    for (j = 1; j <= tot; j++) {
      k = i * p[j];
      if (k >= m)
        break;
      flag[k] = 1;
      if (i % p[j] == 0) {
        u[k] = 0;
        break;
      }
      u[k] = -u[i];
    }
  }
}

int main() {
  auto start_time = clock();
  mobius_sieve();
  cout << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << " " << u[4]
       << " " << u[5] << endl;
  cout << "Time till prime and mobius combined sieve: "
       << (1.0 * (clock() - start_time) / CLOCKS_PER_SEC) << "s" << endl;
  return 0;
}
