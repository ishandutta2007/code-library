#include <bits/stdc++.h>

using namespace std;
using ll = long long;

#define DEBUG

/**
 * An algorithm to calculate the number of factor p in factorial(n)
 *
 * Pre-condition: p is a factor of factorial(n)
 * @return the number of factor p in n!
 */
inline static ll countOfAPrimeFactorInFactorial(ll n, ll p) {

  ll result = 0;
  while (n != 0) {
    n /= p;
    result += n;
  }
  return result;
}

int main() {

#ifdef DEBUG

  assert(countOfAPrimeFactorInFactorial(20, 2) == 18);

#endif

  return 0;
}
