#include <cstdio>
#include <cctype>
#include <vector>

#define maxn 22

template <class T>

inline T read() {
  T r = 0, f = 0;
  char c;
  while (!isdigit(c = getchar()))
    f |= (c == '-');
  while (isdigit(c))
    r = (r << 1) + (r << 3) + (c ^ 48), c = getchar();
  return f ? -r : r;
}

const int mod = 998244353;

unsigned long long A[maxn][maxn];

int n;

struct Mat {

  int a[maxn][maxn];

  Mat() {
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        a[i][j] = 0;
    for (int i = 1; i <= n; i++)
      a[i][i] = 1;
  }

  Mat operator+(const Mat &b) const {
    Mat c;
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        c.a[i][j] = (a[i][j] + b.a[i][j]) % mod;
    return c;
  }

  Mat operator*(const Mat &b) const {
    Mat c;
    for (int i = 1; i <= n; i++)
      for (int k = 1; k <= n; k++)
        for (int j = 1; j <= n; j++)
          A[i][j] += 1ll * a[i][k] * b.a[k][j];
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        c.a[i][j] = A[i][j] % mod, A[i][j] = 0;
    return c;
  }
};

struct D {
  Mat a, b, s;
  D() {
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        s.a[i][j] = 0;
  }
  D operator+(const D &B) const {
    D C;
    C.a = a * B.a, C.b = b * B.b;
    C.s = s + (a * B.s * b);
    return C;
  }
};

template <class T>

inline T qmul(T a, int b) {
  T ans;
  for (; b; b >>= 1) {
    if (b & 1)
      ans = ans + a;
    a = a + a;
  }
  return ans;
}

#define ll long long

inline ll Div(ll p, ll x, ll r, ll q) { return ((__int128)x * p + r) / q; }

D solve(ll p, ll q, ll l, ll r, D U, D R) { //
  if (!l)
    return D();
  if (p >= q)
    return solve(p % q, q, l, r, U, qmul(U, p / q) + R);
  ll Maxy = Div(p, l, r, q);
  if (!Maxy)
    return qmul(R, l);
  ll cntr = l - Div(q, Maxy, -r - 1, p);
  return qmul(R, (q - r - 1) / p) + U +
         solve(q, p, Maxy - 1, (q - r - 1) % p, R, U) + qmul(R, cntr);
}

ll p, q, l, r;

int main() {
#ifndef ONLINE_JUDGE
  freopen("1.in", "r", stdin);
  freopen("1.out", "w", stdout);
#endif
  p = read<ll>();
  q = read<ll>();
  r = read<ll>();
  l = read<ll>();
  n = read<int>();
  D U, R;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      R.a.a[i][j] = read<int>();
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      U.b.a[i][j] = read<int>();
  R.s = R.a;
  D ans = qmul(U, r / q) + solve(p, q, l, r % q, U, R);
  for (int i = 1; i <= n; i++, puts(""))
    for (int j = 1; j <= n; j++)
      printf("%d ", ans.s.a[i][j]);
  return 0;
}
// https://loj.ac/p/6440
// https://loj.ac/p/138

