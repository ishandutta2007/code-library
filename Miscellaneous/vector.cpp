// This Code was made by Chinese_zjc_.
#include <bits/stdc++.h>
// #define debug
const long double PI = std::acos(-1);
const long double eps = 1e-10;
template <class T> long double sqr(const T &x) { return x * x; }
struct vector {
  long double x, y;
  vector(long double x_, long double y_) : x(x_), y(y_) {}
  vector operator+() const { return *this; }
  vector operator-() const { return vector(-this->x, -this->y); }
  friend vector operator+(const vector &A, const vector &B) {
    return vector(A.x + B.x, A.y + B.y);
  }
  friend vector operator-(const vector &A, const vector &B) {
    return vector(A.x - B.x, A.y - B.y);
  }
  friend vector operator*(const long double &A, const vector &B) {
    return vector(A * B.x, A * B.y);
  }
  friend vector operator*(const vector &A, const long double &B) {
    return vector(A.x * B, A.y * B);
  }
  friend vector operator/(const vector &A, const long double &B) {
    return vector(A.x / B, A.y / B);
  }
  friend long double operator*(const vector &A, const vector &B) {
    return A.x * B.x + A.y * B.y;
  }
  friend long double operator^(const vector &A, const vector &B) {
    return A.x * B.y - A.y * B.x;
  } long double length() const {
    return std::sqrt(*this * *this);
  }
  friend long double dis(const vector &A, const vector &B) {
    return (A - B).length();
  }
  friend long double outlen(const vector &A, const vector &B, const vector &C,
                            const long double &Al, const long double &Bl) {
    long double X = A ^ B;
    if (std::abs(X) < eps)
      return INFINITY;
    long double a = (-C.x * B.y + C.y * B.x) / X,
                b = (+C.x * A.y - C.y * A.x) / X;
    if (a < 0 || b < 0)
      return INFINITY;
    return Al * a + Bl * b;
  }
  friend long double outlen(const vector &A, const vector &B, const vector &C) {
    long double X = A ^ B;
    if (std::abs(X) < eps)
      return INFINITY;
    long double a = (-C.x * B.y + C.y * B.x) / X,
                b = (+C.x * A.y - C.y * A.x) / X;
    if (a < 0 || b < 0)
      return INFINITY;
    return A.length() * a + B.length() * b;
  }
  friend vector gettoge(const vector &A, const vector &B, const long double &l,
                        const long double &Al) {
    return (sqr(l) - sqr(B)) / (2 * (A * B + l * Al)) * A;
  }
  friend vector gettoge(const vector &A, const vector &B,
                        const long double &l) {
    return (sqr(l) - sqr(B)) / (2 * (A * B + l * A.length())) * A;
  }
  friend vector trans(const vector &A, const long double B11,
                      const long double B12, const long double B21,
                      const long double B22) {
    return vector(A.x * B11 + A.y * B21, A.x * B12 + A.y * B22);
  }
  friend long double getarea(const vector &C1, const vector &C2,
                             const vector &A, const vector &B,
                             const long double &l) {
    long double Cl = sqr(C1 - C2) / 4;
    long double l_ = l / 2 / std::sqrt(Cl);
    vector p(std::sqrt(sqr(l_) - 1), l_);
    vector O = (C1 + C2) / 2;
    vector C = C2 - O;
    long double T11 = C.y / Cl / p.x, T12 = C.x / Cl / p.y,
                T21 = -C.x / Cl / p.x, T22 = C.y / Cl / p.y;
    vector AA = trans(A - O, T11, T12, T21, T22),
           BB = trans(B - O, T11, T12, T21, T22);
    long double angle = std::acos(AA * BB);
    return (angle - (AA ^ BB)) / 2 * p.x * p.y * Cl;
  }
  friend std::istream &operator>>(std::istream &A, vector &B) {
    return A >> B.x >> B.y;
  }
  friend std::ostream &operator<<(std::ostream &A, const vector &B) {
    return A << '(' << B.x << ',' << B.y << ')';
  }
};
std::vector<vector> a, b;
std::vector<int> Ls, Rs;
int n;
long double len, Length[20005], pre[20005], suf[20005], ans;
inline long double totlen(int x, int y) {
  return pre[x] + suf[y] + outlen(a[x] - a[x - 1], a[y + 1] - a[y], a[x] - a[y],
                                  Length[x], Length[y + 1]);
}
signed main() {
  std::ios::sync_with_stdio(false);
  std::cout << std::fixed << std::setprecision(15);
  std::cin >> n >> len;
  a.reserve(2 * n + 1);
  b.reserve(2 * n + 1);
  Ls.reserve(2 * n + 1);
  Rs.reserve(2 * n + 1);
  for (int i = 0, x, y; i != n; ++i)
    std::cin >> x >> y, a.emplace_back(x, y);
  for (int i = 0; i != n; ++i)
    a.emplace_back(a[i]);
  a.emplace_back(a[0]);
  for (int i = 1; i <= 2 * n; ++i)
    Length[i] = dis(a[i], a[i - 1]);
  for (int i = 1; i <= 2 * n; ++i)
    pre[i] = pre[i - 1] + Length[i];
  for (int i = n - 1; i >= 0; --i)
    suf[i] = suf[i + 1] + Length[i + 1];
  for (int i = 1; i <= n; ++i)
    suf[i + n] = -pre[i];
  for (int l = 1, r = 1; l <= n; ++l) {
    while (r - l < n - 2 && (l == r || totlen(l, r) <= len))
      ++r;
    b.emplace_back(a[l] + gettoge(a[l] - a[l - 1], a[l] - a[r],
                                  len - pre[l] - suf[r], Length[l]));
    Ls.emplace_back(l);
    Rs.emplace_back(r);
    while (totlen(l + 1, r) <= len) {
      b.emplace_back(a[r] + gettoge(a[r] - a[r + 1], a[r] - a[l],
                                    len - pre[l] - suf[r], Length[r + 1]));
      Ls.emplace_back(l);
      Rs.emplace_back(++r);
    }
  }
  b.emplace_back(b.front());
  for (std::size_t i = 0; i != Ls.size(); ++i)
    ans += (b[i] ^ b[i + 1]) / 2 + getarea(a[Ls[i]], a[Rs[i]], b[i], b[i + 1],
                                           len - pre[Ls[i]] - suf[Rs[i]]);
  std::cout << ans << std::endl;
  return 0;
}
// https://codeforces.com/contest/1510/problem/F