#include <bits/stdc++.h>
using namespace std;

using namespace std;

const int MAX = 1000005;
const int MODULO = 1000000007;

int power(int a, long long b, int modulo) {
  if (b == 0)
    return 1 % modulo;
  if (b & 1)
    return 1LL * a * power(a, b - 1, modulo) % modulo;
  int half = power(a, b >> 1, modulo);
  return 1LL * half * half % modulo;
}

struct matrix {
  int a[2][2];
  matrix() {}
};

matrix operator*(matrix a, matrix b) {
  matrix c;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      c.a[i][j] = 0;
      for (int k = 0; k < 2; k++) {
        c.a[i][j] =
            (c.a[i][j] + (1LL * a.a[i][k] * b.a[k][j] % MODULO)) % MODULO;
      }
    }
  }
  return c;
}

matrix power(matrix a, int p) {
  if (p == 1) {
    return a;
  }
  if (p & 1) {
    return power(a, p - 1) * a;
  }
  matrix half = power(a, p >> 1);
  return half * half;
}
void print_matrix(matrix mat) {
  cout << mat.a[0][0] << " " << mat.a[0][1] << endl;
  cout << mat.a[1][0] << " " << mat.a[1][1] << endl;
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  int N;
  cin >> N;
  matrix fib;
  fib.a[0][0] = 1;
  fib.a[0][1] = 1;
  fib.a[1][0] = 1;
  fib.a[1][1] = 0;
  matrix result = power(fib, N);
  print_matrix(result);
  return 0;
}
