#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// generate maximum product using two primes(or any number) and their powers under a limit
// O(log limit /log min(p1,p2))

int gen_max_prod(int p1, int p2, int limit){
  auto product = p1 * p2;
  int maxProduct = 1;
  // for p1^1, p1^2, p1^3, ... find the maximum exponent for p2
  do
  {
  // increase exponent of p2 as much as possible
  auto current = product;
  while (current * p2 <= limit)
    current *= p2;

  // better than before ?
  if (maxProduct < current)
    maxProduct = current;

  // increment p1's exponent by one
  product *= p1;
  } while (product <= limit);

  return maxProduct;
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto start = clock();
  cout << "gen_max_prod(" << 3 << "," << 7 << "," << 10000000 << ")= " << gen_max_prod(3, 7, 10000000) << endl;
  cout << "time: " << (double)(clock() - start) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
