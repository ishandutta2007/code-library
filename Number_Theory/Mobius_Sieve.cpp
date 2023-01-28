// computes till 10^6 in 1 sec (1000x slower than rust version)
#include <bits/stdc++.h>
using namespace std;
// http://oeis.org/A008683
using int64 = long long;
const int N = 100000 + 10;
std::vector<int> divs[N];
std::vector<std::pair<int, int>> edges[N];
int p[N], m;
int64 f[N], mu[N];
void mobius_sieve() {
  for (int i = 1; i < N; ++i) {
    for (int j = i; j < N; j += i) {
      divs[j].push_back(i);
    }
  }
  mu[1] = 1;
  m = 0;
  for (int i = 2; i < N; ++i) {
    if (!p[i])
      p[m++] = i, mu[i] = -1;
    for (int j = 0; i * p[j] < N && j < m; ++j) {
      p[i * p[j]] = 1;
      if (i % p[j])
        mu[i * p[j]] = -mu[i];
      else {
        mu[i * p[j]] = 0;
        break;
      }
    }
  }
}
int main() {
  mobius_sieve();
  for (int i = 0; i < 10; i++)
    cout << mu[i] << endl;
  return 0;
}
