// computes till 10^6 in 1 sec (1000x slower than rust version)
#include <bits/stdc++.h>
using namespace std;
// http://oeis.org/A008683
using int64 = long long;
const int N = 1e6 + 10;
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
  auto start_time = clock();
  mobius_sieve();
  cout << "Time till prime and mobius combined sieve: " << (1.0 * (clock() - start_time) / CLOCKS_PER_SEC) << "s" << endl;
  for (int i = 0; i < 10; i++)
    cout << mu[i] << " ";
  return 0;
}
