#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
using namespace std;

template <unsigned Mod> // DIVISION ASSUMES PRIME MOD
class Mint {
public:
  unsigned v; // PRIVATE IF REPLACE .v WITH ()

public:
  using LL = int64_t;
  using ULL = uint64_t;
  using M = Mint;

  Mint(LL a = 0) : v(a >= 0 ? (a < Mod ? a : a % Mod) : a % Mod + Mod) {}
  int getInt() const { return int(v); }
  int operator()() const { return int(v); }
  M operator-() const { return M() -= *this; }
  M &operator+=(M r) {
    if ((v += r.v) >= Mod)
      v -= Mod;
    return *this;
  }
  M &operator-=(M r) {
    if ((v += Mod - r.v) >= Mod)
      v -= Mod;
    return *this;
  }
  M &operator*=(M r) {
    v = (LL)v * r.v % Mod;
    return *this;
  }
  M &operator/=(M r) { return *this *= power(r, Mod - 2); }
  M &operator++() { return *this += M(1); }
  M &operator--() { return *this -= M(1); }

  M pow(ULL k) const { return power(*this, k); }
  M inv() const { return pow(Mod - 2); }
  friend M operator+(M l, M r) { return l += r; }
  friend M operator-(M l, M r) { return l -= r; }
  friend M operator*(M l, M r) { return l *= r; }
  friend M operator/(M l, M r) { return l /= r; }
  friend bool operator==(M l, M r) { return l.v == r.v; }
  friend bool operator!=(M l, M r) { return l.v != r.v; }
  bool operator<(const M &rhs) const { return v < rhs.v; }
  bool operator>(const M &rhs) const { return v > rhs.v; }
  friend M power(M b, LL e) {
    M r = 1;
    for (; e > 0; e >>= 1LL) {
      if (e & 1LL)
        r *= b;
      b *= b;
    }
    return r;
  }
  int legendre() { return pow((Mod - 1) >> 1).v; }
  int tonelli() { // RETURNS -1 IF NO INV, ELSE CAST TO Mint
    if (v == 0)
      return 0;
    if (legendre() != 1)
      return -1;
    int s = 0, q = Mod - 1;
    while (q % 2 == 0) {
      q >>= 1;
      ++s;
    }
    if (s == 1)
      return pow((Mod + 1) >> 2).v;
    M z = 2, val = Mod - 1;
    for (; z != Mint(0); ++z)
      if (val == z.legendre())
        break;
    M c = power(z, q);
    M r = pow((q + 1) >> 1);
    M t = pow(q);
    int m = s;
    while (t.v != 1) {
      z = t * t;
      int i = 1;
      for (; i < s; ++i) {
        if (z.v == 1)
          break;
        z *= z;
      }
      auto b = power(c, 1 << (m - i - 1));
      r *= b;
      c = b * b;
      t *= c;
      m = i;
    }
    return min(r.v, Mod - r.v);
  }
  friend ostream &operator<<(ostream &os, M a) { return os << a.v; }
  friend istream &operator>>(istream &is, M &a) {
    ULL w;
    is >> w;
    a = M(w);
    return is;
  }
};

template <unsigned MOD, typename T> class Matrix {
private:
  using VT = vector<T>;
  using VVT = vector<vector<T>>;
  using ULL = uint64_t;

  size_t hgt, wth;
  VVT A;

public:
  Matrix(size_t n) : A(n, VT(n, 0)), hgt(n), wth(n) {}
  Matrix(size_t m, size_t n) : A(m, VT(n, 0)), hgt(m), wth(n) {}

  size_t height() const { return hgt; }
  size_t width() const { return wth; }

  inline const VT &operator[](int k) const { return A[k]; }
  inline VT &operator[](int k) { return A[k]; }

  Matrix &operator+=(const Matrix &B) {
    assert(hgt == B.height() && wth == B.width());
    for (size_t i = 0; i < hgt; ++i)
      for (size_t j = 0; j < wth; ++j)
        (*this)[i][j] += B[i][j];
    return *this;
  }
  Matrix &operator-=(const Matrix &B) {
    assert(hgt == B.height() && wth == B.width());
    for (size_t i = 0; i < hgt; ++i)
      for (size_t j = 0; j < wth; ++j)
        (*this)[i][j] -= B[i][j];
    return *this;
  }
  Matrix &operator*=(const Matrix &B) {
    assert(wth == B.height());
    size_t p = B.width();
    Matrix ret(hgt, p);
    for (size_t i = 0; i < hgt; ++i)
      for (size_t k = 0; k < wth; ++k)
        for (size_t j = 0; j < p; ++j)
          ret[i][j] += A[i][k] * B[k][j];
    swap(A, ret.A);
    return *this;
  }
  Matrix &operator*=(const T &v) {
    for (size_t i = 0; i < hgt; ++i)
      for (size_t j = 0; j < wth; ++j)
        A[i][j] *= v;
    return *this;
  }
  Matrix operator^=(uint64_t k) {
    assert(hgt == wth);
    Matrix work = I(hgt);
    for (; k > 0ULL; k >>= 1ULL, (*this) *= (*this))
      if (k & 1ULL)
        work *= (*this);
    return (*this) = work;
  }
  Matrix &operator/=(const T &v) { return (*this) *= T(1) / v; }

  Matrix operator+(const Matrix &B) const { return Matrix(*this) += B; }
  Matrix operator-(const Matrix &B) const { return Matrix(*this) -= B; }
  Matrix operator*(const Matrix &B) const { return Matrix(*this) *= B; }
  Matrix operator*(const T &v) const { return Matrix(*this) *= v; }
  Matrix operator/(const T &v) const { return Matrix(*this) /= v; }
  Matrix operator^(uint64_t k) const { return Matrix(*this) ^= k; }
  friend vector<T> operator*(const Matrix &B, const vector<T> &v) {
    size_t n = B.height(), m = B.width();
    assert(m == v.size());
    vector<T> ret(n, T(0));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < m; ++j)
        ret[i] += B[i][j] * v[j];
    return ret;
  }

private:
  pair<size_t, T> row_reduction(vector<T> &b) {
    bool neg = false;
    size_t check = 0, rank = 0;
    T det = 1;
    assert(b.size() == hgt);
    for (size_t j = 0; j < wth; ++j) {
      size_t pivot = check;
      for (size_t i = check; i < hgt; ++i)
        if (A[i][j] != T(0))
          pivot = i;
      if (check != pivot) {
        neg = !neg;
        swap(A[check], A[pivot]), swap(b[check], b[pivot]);
      }
      if (A[check][j] == T(0)) {
        det = T(0);
        continue;
      }

      ++rank;
      det *= A[check][j];
      T r = T(1) / A[check][j];
      for (size_t k = j + 1; k < wth; ++k)
        A[check][k] *= r;
      b[check] *= r;
      A[check][j] = T(1);
      for (size_t i = 0; i < hgt; ++i) {
        if (i == check)
          continue;
        if (A[i][j] != T(0)) {
          for (size_t k = j + 1; k < wth; ++k)
            A[i][k] -= A[i][j] * A[check][k];
          b[i] -= A[i][j] * b[check];
        }
        A[i][j] = T(0);
      }
      if (++check == hgt)
        break;
    }
    return make_pair(rank, det);
  }

public:
  static Matrix I(size_t n) {
    Matrix ret(n);
    for (size_t i = 0; i < n; ++i)
      ret[i][i] = 1;
    return ret;
  }

  Matrix inverse() { // RETURNS EMPTY MATRIX IF NOT FULL RANK
    assert(wth == hgt);
    if (hgt != wth)
      return Matrix(0);
    Matrix ret = I(hgt);

    for (size_t j = 0; j < hgt; ++j) {
      if (A[j][j] == T(0)) {
        size_t pivot = j;
        for (size_t i = j; i < hgt; ++i)
          if (A[i][j] != T(0)) {
            pivot = i;
            break;
          }
        swap(A[j], A[pivot]), swap(ret[j], ret[pivot]);
      }
      if (A[j][j] == T(0))
        return Matrix(0, 0);
      T r = T(1) / A[j][j];
      for (size_t k = j + 1; k < hgt; ++k)
        A[j][k] *= r;
      for (size_t k = 0; k < hgt; ++k)
        ret[j][k] *= r;
      A[j][j] = T(1);
      for (size_t i = 0; i < hgt; ++i) {
        if (i == j)
          continue;
        if (A[i][j] != T(0)) {
          for (size_t k = j + 1; k < hgt; ++k)
            A[i][k] -= A[i][j] * A[j][k];
          for (size_t k = 0; k < hgt; ++k)
            ret[i][k] -= A[i][j] * ret[j][k];
        }
        A[i][j] = T(0);
      }
    }
    return ret;
  }

  VVT gauss_elimination(VT b) {
    row_reduction(b);
    VVT ret;
    vector<size_t> p(hgt, wth);
    vector<bool> is_zero(wth, true);

    for (size_t i = 0; i < hgt; ++i) {
      for (size_t j = 0; j < wth; ++j) {
        if (A[i][j] != T(0)) {
          p[i] = j;
          break;
        }
      }
      if (p[i] < wth)
        is_zero[p[i]] = false;
      else if (b[i] != T(0))
        return {};
    }
    VT x(wth, T(0));

    for (size_t i = 0; i < hgt; ++i)
      if (p[i] < wth)
        x[p[i]] = b[i];

    ret.push_back(x);
    for (size_t j = 0; j < wth; ++j) {
      if (!is_zero[j])
        continue;
      x[j] = T(1);
      for (size_t i = 0; i < hgt; ++i)
        if (p[i] < wth)
          x[p[i]] = -A[i][j];
      ret.push_back(x), x[j] = T(0);
    }
    return ret;
  }

  T determinant() {
    T det = 1;
    bool neg = false;
    assert(hgt == wth);

    for (size_t i = 0; i < wth; ++i) {
      if (A[i][i] == T(0)) {
        size_t pivot = i;
        for (size_t j = i + 1; j < hgt; ++j)
          if (A[pivot][i] < A[j][i])
            pivot = j;
        if (i != pivot) {
          neg = !neg;
          swap(A[i], A[pivot]);
        }
      }
      if (A[i][i] == T(0))
        return T(0);

      det *= A[i][i];
      T r = T(-1) / A[i][i];

      for (size_t k = i + 1; k < hgt; ++k) {
        if (A[k][i] != T(0)) {
          T val = r * A[k][i];
          for (size_t j = i + 1; j < wth; ++j)
            A[k][j] += val * A[i][j];
          A[k][i] = T(0);
        }
      }
    }
    return neg ? T(-1) * det : det;
  }

  vector<T> charPoly() {
    assert(wth == hgt);
    size_t n = hgt;

    // REDUCE MATRIX TO UPPER HESSENBURG FORM
    for (size_t j = 0; n > 2 && j < n - 2; ++j) {
      for (size_t i = j + 2; i < n; ++i)
        if (A[i][j] != T(0)) {
          swap(A[j + 1], A[i]);
          for (int k = 0; k < n; ++k)
            swap(A[k][j + 1], A[k][i]);
          break;
        }
      if (A[j + 1][j] == 0)
        continue;

      T inv = T(1) / A[j + 1][j];
      for (size_t i = j + 2; i < n; ++i) {
        T coef1 = A[i][j] * inv;
        ULL coef2 = coef1.v;
        ULL coef3 = (-coef1).v;
        for (size_t k = j; k < n; ++k)
          A[i][k] = (coef3 * A[j + 1][k].v + A[i][k].v) % MOD;
        for (size_t k = 0; k < n; ++k)
          A[k][j + 1] = (coef2 * A[k][i].v + A[k][j + 1].v) % MOD;
      }
    }

    // COMPUTE CHAR POLY FROM UPPER HESSENBURG MATRIX
    VVT p(n + 1);
    p[0] = {T(1)};
    for (size_t i = 0; i < n; ++i) {
      p[i + 1].resize(i + 2);
      for (size_t j = 0; j <= i; ++j) {
        p[i + 1][j + 1] += p[i][j];
        p[i + 1][j] -= p[i][j] * A[i][i];
      }
      T bta = 1;
      for (size_t j = i; j > 0; --j) {
        bta *= A[j][j - 1];
        T coef = -bta * A[j - 1][i];
        for (size_t k = 0; k < j; ++k)
          p[i + 1][k] = (ULL(coef.v) * p[j - 1][k].v + p[i + 1][k].v) % MOD;
      }
    }
    return p[n];
  }

private:
  using poly = vector<T>;
  using polyMat = vector<vector<poly>>;

  polyMat hafnianMatrix() {
    polyMat h(hgt);
    for (size_t i = 0; i < hgt; ++i) {
      h[i].resize(i);
      //            assert( A[i][i] == T( 0 ) );       // ZERO DIAGONAL
      for (size_t j = 0; j < i; ++j) {
        h[i][j].resize(hgt / 2 + 1);
        h[i][j][0] = A[i][j];
        //                assert( A[i][j] == A[j][i] );  // SYMMETRIC
      }
    }
    return h;
  }

  poly hafnianVector(const polyMat &b) {
    if (b.size() == 0) {
      poly r(hgt / 2 + 1);
      r[0] = T(1);
      return r;
    }

    size_t m = b.size();
    polyMat c = b;
    c.resize(m - 2);

    poly r = hafnianVector(c);
    for (size_t i = 0; i + 2 < m; ++i) {
      for (size_t j = 0; j < i; ++j) {
        for (size_t k = 0; k < hgt / 2; ++k) {
          T v = 0;
          for (size_t q = 0; q <= k; ++q)
            v = (v.v + 1ULL * b[m - 2][i][q].v * b[m - 1][j][k - q].v +
                 1ULL * b[m - 1][i][q].v * b[m - 2][j][k - q].v) %
                MOD;
          c[i][j][k + 1] += v;
        }
      }
    }

    poly t = hafnianVector(c);
    for (size_t i = 0; i <= hgt / 2; ++i) {
      T v = t[i];
      for (size_t j = 0; j < i; ++j)
        v = (1ULL * t[j].v * b[m - 1][m - 2][i - j - 1].v + v.v) % MOD;
      r[i] = v - r[i];
    }
    return r;
  }

public:
  T numMatchings() { return hafnianVector(hafnianMatrix())[hgt / 2]; }

  friend std::istream &operator>>(std::istream &is, Matrix &x) {
    for (auto &v : x.A)
      for (auto &w : v)
        is >> w;
    return is;
  }
  friend std::ostream &operator<<(std::ostream &out, const Matrix &x) {
    for (auto &v : x.A) {
      for (auto w : v)
        cout << w << " ";
      cout << '\n';
    }
    return out;
  }
};

int main() {
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);

  const unsigned MOD = 998244353;

  int N;
  cin >> N;

  Matrix<MOD, Mint<MOD>> A(N);
  cin >> A;

  vector<Mint<MOD>> ans = A.charPoly();
  for (auto &v : ans)
    cout << v << " ";
  cout << '\n';

  return 0;
}
