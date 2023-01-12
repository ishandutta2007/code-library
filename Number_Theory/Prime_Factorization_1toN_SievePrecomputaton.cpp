// If you need to factorise all numbers from 1 to N .
// Precomputation via Sieve is faster than querrying Pollard Rho N times.
// For this N should be small enough such that it is possible to do
// precomputation.

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

#define MAXN 10000001
const int MOD = 1000000087;

// stores smallest prime factor for every number
int spf[MAXN];

// Calculating SPF (Smallest Prime Factor) for every number till MAXN.
// Time Complexity : O(nloglogn)
void sieve() {
  spf[1] = 1;
  for (int i = 2; i < MAXN; i++)
    spf[i] = i;
  for (int i = 4; i < MAXN; i += 2)
    spf[i] = 2;

  for (int i = 3; i * i < MAXN; i++) {
    if (spf[i] == i) {
      for (int j = i * i; j < MAXN; j += i)
        if (spf[j] == j)
          spf[j] = i;
    }
  }
}

// Time Complexity : O(logn) per query
vector<ll> getFactorization(int x) {
  vector<ll> ret;
  while (x != 1) {
    ret.push_back(spf[x]);
    x = x / spf[x];
  }
  return ret;
}

auto flat_to_map_format(vector<ll> v) {
  map<ll, int> mp;
  int n = v.size();
  for (int i = 0; i < n; i++) {
    int count = 1;
    while (i < n - 1 && v[i] == v[i + 1]) {
      count++;
      i++;
    }
    // cout << v[i] << ":" << count << ", ";
    mp[v[i]] = count;
  }
  return mp;
}

int main(int argc, char const *argv[]) {
  auto start_time = clock();
  sieve();
  cout << "Sieve time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC
       << endl;
  start_time = clock();
  int n = 20;
  for (int x = 1; x <= n; x++) {
    cout << "Factors of " << x << " = {";
    vector<ll> p = getFactorization(x);
    map<ll, int> facs = flat_to_map_format(p);
    for (const auto & [ key, val ] : facs) {
      cout << key << ":" << val << ", ";
    }
    cout << "}" << endl;
  }
  cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
  return 0;
}
