#include <bits/stdc++.h>
using namespace std;

const int mod = 998244353;

struct Mat {
  int n, m;
  vector<vector<int>> a;
  Mat() {}
  Mat(int _n, int _m) {
    n = _n;
    m = _m;
    a.assign(n, vector<int>(m, 0));
  }
  Mat(vector<vector<int>> v) {
    n = v.size();
    m = n ? v[0].size() : 0;
    a = v;
  }
  inline void make_unit() {
    assert(n == m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        a[i][j] = i == j;
    }
  }
  inline Mat operator+(const Mat &b) {
    assert(n == b.n && m == b.m);
    Mat ans = Mat(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        ans.a[i][j] = (a[i][j] + b.a[i][j]) % mod;
      }
    }
    return ans;
  }
  inline Mat operator-(const Mat &b) {
    assert(n == b.n && m == b.m);
    Mat ans = Mat(n, m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        ans.a[i][j] = (a[i][j] - b.a[i][j] + mod) % mod;
      }
    }
    return ans;
  }
  inline Mat operator*(const Mat &b) {
    assert(m == b.n);
    Mat ans = Mat(n, b.m);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < b.m; j++) {
        for (int k = 0; k < m; k++) {
          ans.a[i][j] = (ans.a[i][j] + 1LL * a[i][k] * b.a[k][j] % mod) % mod;
        }
      }
    }
    return ans;
  }

  inline void modinvert(const Mat &b) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        assert(mod - a[i][j] >= 0) if (a[i][j] < 0) a[i][j] = mod - a[i][j];
        else if (a[i][j] >= mod) a[i][j] = a[i][j] % mod;
      }
    }
  }

  inline Mat pow(long long k) {
    assert(n == m);
    Mat ans(n, n), t = a;
    ans.make_unit();
    while (k) {
      if (k & 1)
        ans = ans * t;
      t = t * t;
      k >>= 1;
    }
    return ans;
  }
  inline Mat &operator+=(const Mat &b) { return *this = (*this) + b; }
  inline Mat &operator-=(const Mat &b) { return *this = (*this) - b; }
  inline Mat &operator*=(const Mat &b) { return *this = (*this) * b; }
  inline bool operator==(const Mat &b) { return a == b.a; }
  inline bool operator!=(const Mat &b) { return a != b.a; }
};

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  Mat b(3, 3);

  b.a[0][0] = 1;
  b.a[0][1] = 1;
  b.a[0][2] = 1;

  b.a[1][0] = 1;
  b.a[1][1] = 0;
  b.a[1][2] = 0;

  b.a[2][0] = 0;
  b.a[2][1] = 1;
  b.a[2][2] = 0;

  Mat ans(3, 3);
  ans = b.pow(2);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      cout << ans.a[i][j] << ' ';
    }
    cout << '\n';
  }
  return 0;
}
