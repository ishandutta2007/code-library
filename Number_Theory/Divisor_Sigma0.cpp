#include <bits/stdc++.h>
using namespace std;

int number_of_divisors(long long x) {
  int t = 1;
  for (long long i = 2; i * i <= x; ++i) {
    int cnt = 1;
    for (; x % i == 0; x /= i)
      ++cnt;
    t *= cnt;
  }

  if (x != 1)
    t *= 2;
  return t;
}

int main() {
  int ans = 0;
  for (int n = 1; n <= 39; ++n)
    cout<< "number_of_divisors(" << n << ") = " << number_of_divisors(n) << endl;
}
