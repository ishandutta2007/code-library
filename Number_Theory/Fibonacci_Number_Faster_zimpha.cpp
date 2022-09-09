#include <bits/stdc++.h>
using namespace std;

using int64 = long long;

inline int64 mul_mod(int64 a, int64 b, int64 mod) {
  if (mod < int(1e9))
    return a * b % mod;
  int64 k = (int64)((long double)a * b / mod);
  int64 res = a * b - k * mod;
  res %= mod;
  if (res < 0)
    res += mod;
  return res;
}

void fib(int64 n, int64 &x, int64 &y, int64 mod) { // store in x, n-th
  if (n == 0) {
    x = 0, y = 1;
    return;
  } else if (n & 1) {
    fib(n - 1, y, x, mod);
    y += x;
    if (y >= mod)
      y -= mod;
  } else {
    int64 a, b;
    fib(n >> 1, a, b, mod);
    y = mul_mod(a, a, mod) + mul_mod(b, b, mod);
    x = mul_mod(a, b, mod) + mul_mod(a, b - a + mod, mod);
    if (y >= mod)
      y -= mod;
    if (x >= mod)
      x -= mod;
  }
}

int64 fib(int64 n, int64 mod) {
  int64 x, y;
  fib(n, x, y, mod);
  return x;
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  cout << fib(10, 100) << '\n';
  return 0;
}