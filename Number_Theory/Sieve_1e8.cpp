#include <bits/stdc++.h>
using namespace std;

const int N = 1e8 + 9;

// takes 0.7s for n = 1e8
int spf[N];
vector<int> primes;
void sieve() {
  for (int i = 2; i < N; i++) {
    if (spf[i] == 0)
      spf[i] = i, primes.push_back(i);
    int sz = primes.size();
    for (int j = 0; j < sz && i * primes[j] < N && primes[j] <= spf[i]; j++) {
      spf[i * primes[j]] = primes[j];
    }
  }
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  sieve();
  cout << primes.size() << endl;
  cout << primes[primes.size() - 1] << endl;
  return 0;
}