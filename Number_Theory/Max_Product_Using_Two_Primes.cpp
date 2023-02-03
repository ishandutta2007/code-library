#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// generate maximum product using two primes(or any number) and their powers
// under a limit
// O(log limit /log min(p1,p2))

ll gen_max_prod(ll p1, ll p2, ll limit) {
  ll product = p1 * p2;
  ll maxProduct = 1;
  do {
    ll current = product;
    while (current * p2 <= limit)
      current *= p2;

    if (maxProduct < current)
      maxProduct = current;

    product *= p1;
  } while (product <= limit);

  return maxProduct;
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto start = clock();
  cout << "gen_max_prod(" << 3 << "," << 7 << "," << 10000000
       << ")= " << gen_max_prod(3, 7, 10000000) << endl;
  cout << "time: " << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
